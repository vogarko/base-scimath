/// @file
///
/// A linear solver for parameters from the normal equations
///
/// @copyright (c) 2007 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///
/// @author Tim Cornwell <tim.cornwell@csiro.au>
/// @author Vitaliy Ogarko <vogarko@gmail.com>
///
#include <askap/scimath/fitting/LinearSolver.h>
#include <askap/scimath/fitting/LinearSolverUtils.h>
#include <askap/scimath/fitting/LinearSolverLsqrUtils.h>
#include <askap/scimath/fitting/GenericNormalEquations.h>
#include <askap/scimath/fitting/CalParamNameHelper.h>

#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>
#include <boost/config.hpp>
#include <algorithm>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/Timer.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include <askap/scimath/lsqr_solver/LSQRSolver.h>
#include <askap/scimath/lsqr_solver/ModelDamping.h>
#include <askap/scimath/lsqr_solver/ParallelTools.h>

#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".linearsolver");

#include <iostream>

#include <string>
#include <map>

#include <cmath>
using std::abs;
using std::map;
using std::string;

namespace askap
{
  namespace scimath
  {
    BOOST_CONSTEXPR_OR_CONST double LinearSolver::KeepAllSingularValues;

/// @brief Constructor
/// @details Optionally, it is possible to limit the condition number of
/// normal equation matrix to a given number.
/// @param maxCondNumber maximum allowed condition number of the range
/// of the normal equation matrix for the SVD algorithm. Effectively this
/// puts the limit on the singular values, which are considered to be
/// non-zero (all greater than the largest singular value divided by this
/// condition number threshold). Default is 1e3. Put a negative number
/// if you don't want to drop any singular values (may be a not very wise
/// thing to do!). A very large threshold has the same effect. Zero
/// threshold is not allowed and will cause an exception.
LinearSolver::LinearSolver(double maxCondNumber) :
       itsMaxCondNumber(maxCondNumber),
       itsMajorLoopIterationNumber(0)
#ifdef HAVE_MPI
       ,itsWorkersComm(MPI_COMM_NULL)
#endif
{
    ASKAPASSERT(itsMaxCondNumber != 0);
};

LinearSolver::~LinearSolver()
{
#ifdef HAVE_MPI
    if (itsWorkersComm != MPI_COMM_NULL) {
        MPI_Comm_free(&itsWorkersComm);
    }
#endif
}

void LinearSolver::init()
{
    resetNormalEquations();
}

/// @brief test that all matrix elements are below tolerance by absolute value
/// @details This is a helper method to test all matrix elements
/// @param[in] matr matrix to test
/// @param[in] tolerance tolerance on the element absolute values
/// @return true if all elements are zero within the tolerance
static bool allMatrixElementsAreZeros(const casa::Matrix<double> &matr, const double tolerance)
{
  for (casacore::uInt row = 0; row < matr.nrow(); ++row) {
       for (casacore::uInt col = 0; col < matr.ncolumn(); ++col) {
            if (abs(matr(row,col)) > tolerance) {
                return false;
            }
       }
  }
  return true;
}

/// @brief extract an independent subset of parameters
/// @details This method analyses the normal equations and forms a subset of
/// parameters which can be solved for independently. Although the SVD is more than
/// capable of dealing with degeneracies, it is often too slow if the number of parameters is large.
/// This method essentially gives the solver a hint based on the structure of the equations
/// @param[in] names names for parameters to choose from
/// @param[in] tolerance tolerance on the matrix elements to decide whether they can be considered independent
/// @return names of parameters in this subset
std::vector<std::string> LinearSolver::getIndependentSubset(std::vector<std::string> &names, const double tolerance) const
{
    ASKAPTRACE("LinearSolver::getIndependentSubset");
    ASKAPDEBUGASSERT(names.size() > 0);
    std::vector<std::string> resultNames;
    resultNames.reserve(names.size());
    resultNames.push_back(names[0]);

    // this has been added to the subset, so delete it
    names.erase(names.begin());

    // for each name in subset (which grows as matches are found), check all remaining names for associates.
    for (std::vector<std::string>::const_iterator ciRes = resultNames.begin(); ciRes != resultNames.end(); ++ciRes) {

        // keep a temporary list of added names to erase from the main list
        std::vector<std::string> toErase;
        toErase.reserve(names.size());

        for (std::vector<std::string>::const_iterator ci = names.begin(); ci != names.end(); ++ci) {
            const casacore::Matrix<double>& nm1 = normalEquations().normalMatrix(*ci, *ciRes);
            const casacore::Matrix<double>& nm2 = normalEquations().normalMatrix(*ciRes, *ci);

            if (!allMatrixElementsAreZeros(nm1,tolerance) || !allMatrixElementsAreZeros(nm2,tolerance)) {
                // this parameter (iterated in the inner loop) belongs to the subset
                resultNames.push_back(*ci);
                // also add it to the temporary list of added names to erased
                toErase.push_back(*ci);
            }
        }

        // erase all of the names that have just been added from the main list
        for (std::vector<std::string>::iterator it0 = toErase.begin(); it0 != toErase.end(); ++it0) {
            std::vector<std::string>::iterator it = std::find(names.begin(), names.end(), *it0);
            if (it != names.end()) names.erase(it);
        }

    }
    return resultNames;
}

size_t LinearSolver::calculateGainNameIndices(const std::vector<std::string> &names,
                                              const Params &params,
                                              std::vector<std::pair<std::string, int> > &indices) const
{
    ASKAPCHECK(indices.size() == names.size(), "Wrong vector size in calculateGainNameIndices!");

    size_t nParameters = 0;
    std::vector<std::pair<std::string, int> >::iterator it = indices.begin();
    for (std::vector<std::string>::const_iterator cit = names.begin();
            cit != names.end(); ++cit, ++it) {
        ASKAPDEBUGASSERT(it != indices.end());
        it->second = nParameters;
        it->first = *cit;
        ASKAPLOG_DEBUG_STR(logger, "Processing " << *cit << " " << nParameters);

        const casa::uInt newParameters = normalEquations().dataVector(*cit).nelements();
        nParameters += newParameters;
        ASKAPDEBUGASSERT((params.isFree(*cit) ? params.value(*cit).nelements() : newParameters) == newParameters);
    }
    ASKAPLOG_DEBUG_STR(logger, "Done");

    return nParameters;
}

/// @brief solve for a subset of parameters
/// @details This method is used in solveNormalEquations
/// @param[in] params parameters to be updated
/// @param[in] quality Quality of the solution
/// @param[in] names names of the parameters to solve for
std::pair<double,double> LinearSolver::solveSubsetOfNormalEquations(Params &params, Quality& quality,
                   const std::vector<std::string> &names) const
{
    ASKAPTRACE("LinearSolver::solveSubsetOfNormalEquations");
#ifdef HAVE_MPI
    ASKAPLOG_INFO_STR(logger, "Started LinearSolver::solveSubsetOfNormalEquations, with MPI.");
#else
    ASKAPLOG_INFO_STR(logger, "Started LinearSolver::solveSubsetOfNormalEquations, without MPI.");
#endif

    std::pair<double,double> result(0.,0.);

    // Solving A^T Q^-1 V = (A^T Q^-1 A) P

    std::vector<std::pair<std::string, int> > indices(names.size());
    int nParameters = calculateGainNameIndices(names, params, indices);
    ASKAPCHECK(nParameters > 0, "No free parameters in a subset of normal equations!");

    // Convert the normal equations to gsl format.
    gsl_matrix * A = gsl_matrix_alloc(nParameters, nParameters);
    gsl_vector * B = gsl_vector_alloc(nParameters);
    gsl_vector * X = gsl_vector_alloc(nParameters);

    for (std::vector<std::pair<std::string, int> >::const_iterator indit2=indices.begin();indit2!=indices.end(); ++indit2)  {
        for (std::vector<std::pair<std::string, int> >::const_iterator indit1=indices.begin();indit1!=indices.end(); ++indit1)  {

             // Axes are dof, dof for each parameter.
             // Take a deep breath for const-safe indexing into the double layered map.
             const casa::Matrix<double>& nm = normalEquations().normalMatrix(indit1->first, indit2->first);

             if (&nm == &emptyMatrix) {
                 continue;
             }

             for (size_t row=0; row<nm.nrow(); ++row)  {
                  for (size_t col=0; col<nm.ncolumn(); ++col) {
                       const double elem = nm(row,col);
                       ASKAPCHECK(!std::isnan(elem), "Normal matrix seems to have NaN for row = "<<
                           row<<" and col = "<<col<<", this shouldn't happen!");
                       gsl_matrix_set(A, row+(indit1->second), col+(indit2->second), elem);
                       //std::cout<<"A "<<row+(indit1->second)<<" "<<col+(indit2->second)<<" "<<nm(row,col)<<std::endl;
                  }
             }
         }
    }

    // Populate the right-hand side vector B.
    solverutils::populate_B(normalEquations(), indices, B, solverutils::assign_to_gsl_vector);

    if (algorithm() == "SVD") {
        ASKAPLOG_INFO_STR(logger, "Solving normal equations using the SVD solver");

        gsl_matrix * V = gsl_matrix_alloc (nParameters, nParameters);
        ASKAPDEBUGASSERT(V!=NULL);
        gsl_vector * S = gsl_vector_alloc (nParameters);
        ASKAPDEBUGASSERT(S!=NULL);
        gsl_vector * work = gsl_vector_alloc (nParameters);
        ASKAPDEBUGASSERT(work!=NULL);

        gsl_error_handler_t *oldhandler=gsl_set_error_handler_off();
        ASKAPLOG_DEBUG_STR(logger, "Running SV decomp");
        
        int status = 1;
        int failures = 0;
        int failure_limit = 10;
        while (status != 0 && failures < failure_limit ) {
            status = gsl_linalg_SV_decomp (A, V, S, work);
            ++failures;
        }
        ASKAPCHECK(status == 0, "gsl_linalg_SV_decomp failed, status = "<<status << ". Tries = " << failures );
        gsl_set_error_handler(oldhandler);

        // A hack for now. For some reason, for some matrices gsl_linalg_SV_decomp may return NaN as singular value, perhaps some
        // numerical precision issue inside SVD. Although it needs to be investigated further  (see ASKAPSDP-2270), for now trying
        // to replace those singular values with zeros to exclude them from processing. Note, singular vectors may also contain NaNs.
        for (int i=0; i<nParameters; ++i) {
          if (std::isnan(gsl_vector_get(S,i))) {
              gsl_vector_set(S,i,0.);
          }
          for (int k=0; k < nParameters; ++k) {
               if(std::isnan(gsl_matrix_get(V,i,k))) {
                   gsl_matrix_set(V,i,k,0.);
               }
          }
        }
        // end of the hack

         //SVDecomp (A, V, S);

        // Code to put a limit on the condition number of the system.
        const double singularValueLimit = nParameters>1 ?
                    gsl_vector_get(S,0)/itsMaxCondNumber : -1.;
        for (int i=1; i<nParameters; ++i) {
             if (gsl_vector_get(S,i)<singularValueLimit) {
                 gsl_vector_set(S,i,0.);
             }
        }

        gsl_vector * X = gsl_vector_alloc(nParameters);
        ASKAPDEBUGASSERT(X!=NULL);

        const int solveStatus = gsl_linalg_SV_solve (A, V, S, B, X);
        ASKAPCHECK(solveStatus == 0, "gsl_linalg_SV_solve failed");

        // Now find the statistics for the decomposition.
        int rank=0;
        double smin = 1e50;
        double smax = 0.0;
        for (int i=0;i<nParameters; ++i) {
             const double sValue = std::abs(gsl_vector_get(S, i));
             ASKAPCHECK(!std::isnan(sValue), "Got NaN as a singular value for normal matrix, this shouldn't happen S[i]="<<gsl_vector_get(S,i)<<" parameter "<<i<<" singularValueLimit="<<singularValueLimit);
             if(sValue>0.0) {
                ++rank;
                if ((sValue>smax) || (i == 0)) {
                    smax=sValue;
                }
                if ((sValue<smin) || (i == 0)) {
                    smin=sValue;
                }
            }
        }
        result.first = smin;
        result.second = smax;
        quality.setDOF(nParameters);
        if (status != 0) {
            ASKAPLOG_WARN_STR(logger, "Solution is considered invalid due to gsl_linalg_SV_decomp failure, main matrix is effectively rank zero");
            quality.setRank(0);
            quality.setCond(0.);
            quality.setInfo("SVD decomposition rank deficient");
        } else {
            quality.setRank(rank);
            quality.setCond(smax/smin);
            if (rank==nParameters) {
                quality.setInfo("SVD decomposition rank complete");
            } else {
                quality.setInfo("SVD decomposition rank deficient");
            }
        }

        // Update the parameters for the calculated changes.
        solverutils::update_solution(indices, params, X, solverutils::retrieve_from_gsl_vector);

        gsl_vector_free(S);
        gsl_vector_free(work);
        gsl_matrix_free(V);
    }
    else if (algorithm() == "Chol") {
        // TODO: It seems this branch is never actually used. Need to remove it?
        ASKAPLOG_INFO_STR(logger, "Solving normal equations using the Cholesky decomposition solver");

        quality.setInfo("Cholesky decomposition");
        gsl_linalg_cholesky_decomp(A);
        gsl_linalg_cholesky_solve(A, B, X);

        // Update the parameters for the calculated changes.
        solverutils::update_solution(indices, params, X, solverutils::retrieve_from_gsl_vector);
    }
    else {
        ASKAPTHROW(AskapError, "Wrong calibration solver type: " << algorithm());
    }

    // Free up gsl storage.
    gsl_matrix_free(A);
    gsl_vector_free(B);
    gsl_vector_free(X);

    return result;
}

std::pair<double,double> LinearSolver::solveSubsetOfNormalEquationsLSQR(Params &params, Quality& quality,
                   const std::vector<std::string> &__names) const
{
    ASKAPTRACE("LinearSolver::solveSubsetOfNormalEquationsLSQR");
#ifdef HAVE_MPI
    ASKAPLOG_INFO_STR(logger, "Started LinearSolver::solveSubsetOfNormalEquationsLSQR, with MPI.");
#else
    ASKAPLOG_INFO_STR(logger, "Started LinearSolver::solveSubsetOfNormalEquationsLSQR, without MPI.");
#endif

    std::pair<double, double> result(0.,0.);

    std::vector<std::string> names = __names;
    std::vector<std::pair<std::string, int> > indices;
    size_t nParameters;

    const GenericNormalEquations& gne = dynamic_cast<const GenericNormalEquations&>(normalEquations());

    if (gne.indexedNormalMatrixInitialized()) {
    // new matrix format
        ASKAPCHECK(gne.indexedDataVectorInitialized(), "Indexed data vector is not initialized!");
        nParameters = 2 * names.size();
    } else {
    // old matrix format
        std::sort(names.begin(), names.end(), lsqrutils::compareGainNames);
        indices.resize(names.size());
        nParameters = calculateGainNameIndices(names, params, indices);
    }
    ASKAPCHECK(nParameters > 0, "No free parameters in a subset of normal equations!");

    // Solving A^T Q^-1 V = (A^T Q^-1 A) P

    //------------------------------------------------------------------------------
    // Define MPI partitioning.
    //------------------------------------------------------------------------------
    int myrank = 0;
    int nbproc = 1;

    bool matrixIsParallel = solverutils::getParameter("parallelMatrix", parameters(), false);

#ifdef HAVE_MPI
    if (matrixIsParallel) {
    // The parallel matrix case - need to define the parallel partitioning.
        ASKAPCHECK(itsWorkersComm != MPI_COMM_NULL, "Workers communicator is not defined!");
        MPI_Comm_rank(itsWorkersComm, &myrank);
        MPI_Comm_size(itsWorkersComm, &nbproc);
    }
#endif
    if (myrank == 0) {
        ASKAPLOG_DEBUG_STR(logger, "it, matrixIsParallel, nbproc = " << itsMajorLoopIterationNumber
                                    << ", " << matrixIsParallel << ", " << nbproc);
    }

    //------------------------------------------------------------------------------
    // Define LSQR solver sparse matrix.
    //------------------------------------------------------------------------------
#ifdef HAVE_MPI
    size_t nParametersTotal = lsqr::ParallelTools::get_total_number_elements(nParameters, nbproc, itsWorkersComm);
#else
    size_t nParametersTotal = nParameters;
#endif
    if (myrank == 0) ASKAPLOG_DEBUG_STR(logger, "nParameters = " << nParameters);
    if (myrank == 0) ASKAPLOG_DEBUG_STR(logger, "nParametersTotal = " << nParametersTotal);

    lsqr::SparseMatrix matrix;
#ifdef HAVE_MPI
    if (matrixIsParallel) {
        matrix = lsqr::SparseMatrix(nParametersTotal, itsWorkersComm);
    }
    else
#endif
    {
        matrix = lsqr::SparseMatrix(nParametersTotal);
    }

    // Copy matrix elements from normal matrix (map of map of matrixes) to the solver sparse matrix (in CSR format).
    lsqrutils::buildLSQRSparseMatrix(gne, indices, matrix, nParameters, matrixIsParallel);

    size_t nonzeros = matrix.GetNumberElements();
    double sparsity = (double)(nonzeros) / (double)(nParameters) / (double)(nParameters);
    ASKAPLOG_DEBUG_STR(logger, "Jacobian nonzeros, sparsity = " << nonzeros << ", " << sparsity << " on rank " << myrank);

    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Solving normal equations using the LSQR solver");

    //------------------------------------------------------------------
    // Define the right-hand side (the data misfit part).
    //------------------------------------------------------------------
    lsqr::Vector b_RHS(nParametersTotal, 0.);

    // Populate the right-hand side vector B.
    size_t nDataAdded;
    if (gne.indexedDataVectorInitialized()) {
    // new matrix format
        ASKAPCHECK(gne.indexedNormalMatrixInitialized(), "Indexed normal matrix is not initialized!");
        nDataAdded = gne.unrollIndexedDataVector(b_RHS);
    } else {
    // old matrix format
        nDataAdded = solverutils::populate_B(normalEquations(), indices, b_RHS, solverutils::assign_to_lsqr_vector);
    }
    ASKAPCHECK(nDataAdded == (size_t)(nParameters), "Wrong number of data added on rank " << myrank);

    if (matrixIsParallel) {
#ifdef HAVE_MPI
        lsqr::ParallelTools::get_full_array_in_place(nParameters, b_RHS, true, myrank, nbproc, itsWorkersComm);
#endif
    }

    double misfit_cost = lsqrutils::calculateCost(b_RHS);
    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Misfit cost = " << misfit_cost);

    //------------------------------------------------------------------
    // Adding smoothness constraints.
    //------------------------------------------------------------------
    bool addSmoothing = solverutils::getParameter("smoothing", parameters(), false);

    if (addSmoothing) {
    // Adding smoothing constraints into the system of equations.
        ASKAPCHECK(matrixIsParallel, "Smoothing constraints should be used in the parallel matrix mode!");
        // TODO: Remove non indexed normal matrix branches below, after all tests are passing.
        ASKAPCHECK(gne.indexedNormalMatrixInitialized(), "Smoothing constraints should be used with indexed normal matrix format!");

        //--------------------------------------------------------------
        // Extract the solution at the current major iteration (before the update).
        std::vector<double> x0(nParametersTotal);

        // Retrieve the local solution (at the current worker).
        if (gne.indexedNormalMatrixInitialized()) {
            lsqrutils::getCurrentSolutionVector(gne, params, x0);
        } else {
            lsqrutils::getCurrentSolutionVector(indices, params, x0);
        }

#ifdef HAVE_MPI
        // Retrieve the global solution (at all workers).
        lsqr::ParallelTools::get_full_array_in_place(nParameters, x0, true, myrank, nbproc, itsWorkersComm);
#endif

        // Reading the number of channels.
        size_t nChannels = solverutils::getParameter("nChan", parameters(), 0);

        if (!gne.indexedNormalMatrixInitialized()) {
        // Sanity check (for channel order consistency) for old matrix format.
            // Assume the same number of channels at every CPU.
            size_t nChannelsLocal = nChannels / (nParametersTotal / nParameters);
            // NOTE: Assume channels are ordered with the MPI rank order, i.e., the higher the rank the higher the channel number.
            ASKAPCHECK(lsqrutils::testMPIRankOrderWithChannels(myrank, nChannelsLocal, indices), "Channels are not ordered with MPI ranks!");
        }

        //--------------------------------------------------------------
        double smoothingWeight = lsqrutils::getSmoothingWeight(parameters(), itsMajorLoopIterationNumber);
        double smoothingLevel = solverutils::getParameter("smoothingLevel", parameters(), 0.);
        int smoothingType = solverutils::getParameter("smoothingType", parameters(), 2);
        bool addSpectralDiscont = solverutils::getParameter("spectralDiscont", parameters(), false);
        size_t spectralDiscontStep = (size_t)solverutils::getParameter("spectralDiscontStep", parameters(), 40);

        bool indexedNormalMatrixFormat = gne.indexedNormalMatrixInitialized();
        lsqrutils::addSmoothnessConstraints(matrix, b_RHS, x0, nParameters, nChannels,
                                            smoothingWeight, smoothingLevel, smoothingType,
                                            addSpectralDiscont, spectralDiscontStep,
                                            indexedNormalMatrixFormat);
    }
    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Matrix nelements = " << matrix.GetNumberElements());

    // A simple approximation for the upper bound of the rank of the  A'A matrix.
    size_t rank_approx = matrix.GetNumberNonemptyRows();

    //------------------------------------------------------------------
    // Adding damping.
    //------------------------------------------------------------------
    // Setting damping parameters.
    double alpha = solverutils::getParameter("alpha", parameters(), 0.01);
    double norm = solverutils::getParameter("norm", parameters(), 2.0);

    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Adding model damping, with alpha = " << alpha);

    lsqr::ModelDamping damping(nParameters);
    damping.Add(alpha, norm, matrix, b_RHS, NULL, NULL, NULL, myrank, nbproc);

    //------------------------------------------------------------------
    // Column normalisation.
    bool normalizeColumns = solverutils::getParameter("normalizeColumns", parameters(), false);

    std::vector<double> columnNorms;
    if (normalizeColumns) {
        if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Normalizing matrix columns");
        columnNorms.resize(nParameters);
        matrix.NormalizeColumns(columnNorms);
    }

    //------------------------------------------------------------------
    // Calculating the total cost.
    double total_cost = lsqrutils::calculateCost(b_RHS);
    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Total cost = " << total_cost);

    //------------------------------------------------------------------
    // Setting solver parameters.
    int niter = solverutils::getParameter("niter", parameters(), 100);
    double rmin = solverutils::getParameter("rmin", parameters(), 1.e-13);
    bool suppress_output = !(solverutils::getParameter("verbose", parameters(), false));

    //------------------------------------------------------------------
    // Solving the matrix system.
    //------------------------------------------------------------------
    casa::Timer timer;
    timer.mark();

    lsqr::Vector x(nParameters, 0.);
    lsqr::LSQRSolver solver(matrix.GetCurrentNumberRows(), nParameters);

    solver.Solve(niter, rmin, matrix, b_RHS, x, suppress_output);

    ASKAPLOG_INFO_STR(logger, "Completed LSQR in " << timer.real() << " seconds on rank " << myrank);

    //------------------------------------------------------------------
    // Solution re-scaling.
    if (normalizeColumns) {
        for (size_t i = 0; i < x.size(); i++) {
            assert(columnNorms[i] != 0.);
            x[i] /= columnNorms[i];
        }
    }

    //------------------------------------------------------------------
    // Update the parameters with the calculated changes.
    if (gne.indexedNormalMatrixInitialized()) {
    // TODO: Move to a function.
    // new matrix format
        size_t nChannelsLocal = gne.getNumberLocalChannels();
        size_t nBaseParameters = gne.getNumberBaseParameters();

        for (size_t i = 0; i < nBaseParameters; i++) {
            for (size_t chan = 0; chan < nChannelsLocal; chan++) {
                std::string paramName = gne.getFullParameterName(i, chan);

                auto *data = params.value(paramName).data();
                size_t index = i + nBaseParameters * chan;

                data[0] += x[2 * index];
                data[1] += x[2 * index + 1];
            }
        }

    } else {
    // old matrix format
        solverutils::update_solution(indices, params, x, solverutils::retrieve_from_lsqr_vector);
    }

    //------------------------------------------------------------------
    // Set approximate solution quality.
    quality.setDOF(nParameters);
    quality.setRank(rank_approx);

    return result;
}

/// @brief solve for parameters
/// The solution is constructed from the normal equations and given
/// parameters are updated. If there are no free parameters in the
/// given Params class, all unknowns in the normal
/// equatons will be solved for.
/// @param[in] params parameters to be updated
/// @param[in] quality Quality of solution
/// @note This is fully general solver for the normal equations for any shape
/// parameters.
bool LinearSolver::solveNormalEquations(Params &params, Quality& quality)
{
    ASKAPTRACE("LinearSolver::solveNormalEquations");

    // Solving A^T Q^-1 V = (A^T Q^-1 A) P
    // Find all the free parameters.
    std::vector<std::string> names;

    const GenericNormalEquations& gne = dynamic_cast<const GenericNormalEquations&>(normalEquations());

    if (gne.indexedNormalMatrixInitialized()) {
    // Solve for all unknowns in the indexed normal matrix case.
    // Currently, the indexed normal matrix format is activated when the parallel matrix is used.
    // Note that for the parallel matrix case, params.freeNames() should also contain all unknowns,
    // as we do not fix model parameters for the local models, and thus solve for all parameters (including flagged data).
        names = normalEquations().unknowns();
    } else {
        names = params.freeNames();
    }

    if (names.size() == 0) {
        // List of parameters is empty, will solve for all unknowns in the equation.
        names = normalEquations().unknowns();
    }
    ASKAPCHECK(names.size() > 0, "No free parameters in Linear Solver");

    ASKAPLOG_INFO_STR(logger, "LinearSolver free names: " << params.freeNames().size());
    ASKAPLOG_INFO_STR(logger, "LinearSolver fixed names: " << params.fixedNames().size());
    ASKAPLOG_INFO_STR(logger, "LinearSolver names: " << names.size());

    if (algorithm() == "LSQR") {
    // LSQR solver.
        solveSubsetOfNormalEquationsLSQR(params, quality, names);
    } else {
        while (names.size() > 0) {
            ASKAPLOG_INFO_STR(logger, "Solving independent subset of parameters");
            const std::vector<std::string> subsetNames = getIndependentSubset(names,1e-6);
            solveSubsetOfNormalEquations(params, quality, subsetNames);
        }
    }

    return true;
};

Solver::ShPtr LinearSolver::clone() const
{
    return Solver::ShPtr(new LinearSolver(*this));
}

#ifdef HAVE_MPI
void LinearSolver::setWorkersCommunicator(const MPI_Comm &comm)
{
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_dup(comm, &itsWorkersComm);
}
#endif

void LinearSolver::setMajorLoopIterationNumber(size_t it)
{
    itsMajorLoopIterationNumber = it;
}

}
}
