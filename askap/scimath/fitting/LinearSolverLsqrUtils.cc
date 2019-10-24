/// @file
///
/// Utils for the LSQR solver with smoothing constraints.
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
/// @author Vitaliy Ogarko <vogarko@gmail.com>
///
#include <fitting/LinearSolverLsqrUtils.h>
#include <fitting/LinearSolverUtils.h>
#include <fitting/GenericNormalEquations.h>
#include <lsqr_solver/ParallelTools.h>

#include <askap/AskapError.h>
#include <profile/AskapProfiler.h>

#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".lsqrutils");

#include <iostream>
#include <string>
#include <map>

namespace askap { namespace scimath { namespace lsqrutils {

// NOTE: Copied from "calibaccess/CalParamNameHelper.h", as currently accessors depends of scimath.
/// @brief extract coded channel and parameter name
/// @details This is a reverse operation to codeInChannel. Note, no checks are done that the name passed
/// has coded channel present.
/// @param[in] name full name of the parameter
/// @return a pair with extracted channel and the base parameter name
static std::pair<casa::uInt, std::string> extractChannelInfo(const std::string &name)
{
    size_t pos = name.rfind(".");
    ASKAPCHECK(pos != std::string::npos, "Expect dot in the parameter name passed to extractChannelInfo, name=" << name);
    ASKAPCHECK(pos + 1 != name.size(), "Parameter name=" << name << " ends with a dot");
    return std::pair<casa::uInt, std::string>(utility::fromString<casa::uInt>(name.substr(pos + 1)), name.substr(0, pos));
}

bool compareGainNames(const std::string& gainA, const std::string& gainB) {
    try {
        std::pair<casa::uInt, std::string> paramInfoA = extractChannelInfo(gainA);
        std::pair<casa::uInt, std::string> paramInfoB = extractChannelInfo(gainB);

        // Parameter name excluding channel number.
        std::string parNameA = paramInfoA.second;
        std::string parNameB = paramInfoB.second;

        casa::uInt chanA = paramInfoA.first;
        casa::uInt chanB = paramInfoB.first;

        int res = parNameA.compare(parNameB);

        if (res == 0) {
        // Same names, sort by channel number.
            return (chanA <= chanB);
        } else {
            if (res < 0) {
                return true;
            } else {
                return false;
            }
        }
    }
    catch (AskapError &e) {
        return (gainA.compare(gainB) < 0);
    }
}

void buildLSQRSparseMatrix(const INormalEquations &ne,
                           const std::vector<std::pair<string, int> > &indices,
                           lsqr::SparseMatrix &matrix,
                           size_t nParameters,
                           bool matrixIsParallel)
{
#ifdef HAVE_MPI
    MPI_Comm workersComm = matrix.GetComm();
    ASKAPCHECK(workersComm != MPI_COMM_NULL, "Workers communicator is not defined!");

    int myrank, nbproc;
    MPI_Comm_rank(workersComm, &myrank);
    MPI_Comm_size(workersComm, &nbproc);

    size_t nParametersTotal = lsqr::ParallelTools::get_total_number_elements(nParameters, nbproc, workersComm);
    size_t nParametersSmaller = lsqr::ParallelTools::get_nsmaller(nParameters, myrank, nbproc, workersComm);
#else
    size_t nParametersTotal = nParameters;
    size_t nParametersSmaller = 0;
#endif

    std::map<string, size_t> indicesMap;
    for (std::vector<std::pair<string, int> >::const_iterator it = indices.begin();
         it != indices.end(); ++it) {
        indicesMap[it->first] = (size_t)(it->second);
    }

    //------------------------------------------------------------------------------------------------------------------------
    const GenericNormalEquations& gne = dynamic_cast<const GenericNormalEquations&>(ne);

    if (matrixIsParallel) {
        // Adding starting matrix empty rows, i.e., the rows in a big block-diagonal matrix above the current block.
        for (size_t i = 0; i < nParametersSmaller; ++i) {
            matrix.NewRow();
        }
    }

    // Loop over matrix rows.
    for (std::vector<std::pair<string, int> >::const_iterator indit1 = indices.begin();
            indit1 != indices.end(); ++indit1) {

        const std::map<string, casa::Matrix<double> >::const_iterator colItBeg = gne.getNormalMatrixRowBegin(indit1->first);
        const std::map<string, casa::Matrix<double> >::const_iterator colItEnd = gne.getNormalMatrixRowEnd(indit1->first);

        if (colItBeg != colItEnd) {
            const size_t nrow = colItBeg->second.nrow();
            for (size_t row = 0; row < nrow; ++row) {
                matrix.NewRow();
                // Loop over column elements.
                for (std::map<string, casa::Matrix<double> >::const_iterator colIt = colItBeg;
                        colIt != colItEnd; ++colIt) {

                    const std::map<string, size_t>::const_iterator indicesMapIt = indicesMap.find(colIt->first);
                    if (indicesMapIt != indicesMap.end()) {
                    // It is a parameter to solve for, adding it to the matrix.

                        const size_t colIndex = indicesMapIt->second;
                        const casa::Matrix<double>& nm = colIt->second;

                        ASKAPCHECK(nrow == nm.nrow(), "Not consistent normal matrix element element dimension!");
                        const size_t ncolumn = nm.ncolumn();
                        for (size_t col = 0; col < ncolumn; ++col) {
                             const double elem = nm(row, col);
                             ASKAPCHECK(!std::isnan(elem), "Normal matrix seems to have NaN for row = "<< row << " and col = " << col << ", this shouldn't happen!");
                             matrix.Add(elem, col + colIndex);
                        }
                    }
                }
            }
        } else {
        // Adding empty matrix rows.
            // Need to add corresponding empty rows in the sparse matrix.
            const size_t nrow = gne.dataVector(indit1->first).nelements();
            for (size_t i = 0; i < nrow; ++i) {
                matrix.NewRow();
            }
        }
    }

    ASKAPCHECK(matrix.GetCurrentNumberRows() == (nParametersSmaller + nParameters), "Wrong number of matrix rows!");
    if (matrixIsParallel) {
        // Adding ending matrix empty rows, i.e., the rows in a big block-diagonal matrix below the current block.
        size_t nEndRows = nParametersTotal - nParametersSmaller - nParameters;
        for (size_t i = 0; i < nEndRows; ++i) {
            matrix.NewRow();
        }
    }
    ASKAPCHECK(matrix.GetCurrentNumberRows() == nParametersTotal, "Wrong number of matrix rows!");
    matrix.Finalize(nParameters);
}

void getCurrentSolutionVector(const std::vector<std::pair<std::string, int> >& indices,
                              const Params& params,
                              std::vector<double>& solution)
{
    ASKAPCHECK(solution.size() >= 2 * indices.size(), "Wrong size of the solution vector in getCurrentSolutionVector!");

    for (std::vector<std::pair<string, int> >::const_iterator indit = indices.begin();
            indit != indices.end(); ++indit) {
        casa::IPosition vecShape(1, params.value(indit->first).nelements());
        casa::Vector<double> value(params.value(indit->first).reform(vecShape));
        for (size_t i = 0; i < value.nelements(); ++i) {
            solution[indit->second + i] = value(i);
        }
    }
}

double getSmoothingWeight(const std::map<std::string, std::string>& parameters,
                          size_t majorLoopIterationNumber) {

    double smoothingWeight = 0.;

    double smoothingMinWeight = solverutils::getParameter("smoothingMinWeight", parameters, 0.);
    double smoothingMaxWeight = solverutils::getParameter("smoothingMaxWeight", parameters, 3.e+6);
    size_t smoothingNsteps = solverutils::getParameter("smoothingNsteps", parameters, 10);

    if (majorLoopIterationNumber < smoothingNsteps) {
        if (smoothingMinWeight == smoothingMaxWeight) {
            smoothingWeight = smoothingMaxWeight;
        } else {
            double span = smoothingMaxWeight - smoothingMinWeight;
            ASKAPCHECK(span > 0, "Wrong smoothing weight!");

            // Logarithmic sweep (between the minimum and maximum weights).
            smoothingWeight = smoothingMinWeight + std::pow(10., log10(span) / (double)(smoothingNsteps) * (double)(majorLoopIterationNumber));
        }
    } else {
        // Relaxation with constant weight.
        smoothingWeight = smoothingMaxWeight;
    }
    return smoothingWeight;
}

/// @brief Adding smoothness constraints to the system of equations.
void addSmoothnessConstraints(lsqr::SparseMatrix& matrix,
                              lsqr::Vector& b_RHS,
                              const std::vector<std::pair<std::string, int> >& indices,
                              const std::vector<double>& x0,
                              size_t nParameters,
                              size_t nChannels,
                              double smoothingWeight,
                              int gradientType)
{
    ASKAPCHECK(nChannels > 1, "Wrong number of channels for smoothness constraints!");

#ifdef HAVE_MPI
    MPI_Comm workersComm = matrix.GetComm();
    ASKAPCHECK(workersComm != MPI_COMM_NULL, "Workers communicator is not defined!");

    int myrank, nbproc;
    MPI_Comm_rank(workersComm, &myrank);
    MPI_Comm_size(workersComm, &nbproc);

    size_t nParametersTotal = lsqr::ParallelTools::get_total_number_elements(nParameters, nbproc, workersComm);
    size_t nParametersSmaller = lsqr::ParallelTools::get_nsmaller(nParameters, myrank, nbproc, workersComm);
#else
    int myrank = 0;
    size_t nParametersTotal = nParameters;
    size_t nParametersSmaller = 0;
#endif

    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Adding smoothness constraints, with weight = " << smoothingWeight);

    matrix.Extend(nParametersTotal);
    b_RHS.resize(b_RHS.size() + nParametersTotal);

    //-----------------------------------------------------------------------------
    // Assume the same number of channels at every CPU.
    size_t nChannelsLocal = nChannels / (nParametersTotal / nParameters);
    size_t nextChannelIndexShift = nParameters - (nChannelsLocal - 1) * 2;

    if (myrank == 0) ASKAPLOG_DEBUG_STR(logger, "nChannelsLocal = " << nChannelsLocal);

    std::vector<int> leftIndexGlobal(nParametersTotal);
    std::vector<int> rightIndexGlobal(nParametersTotal);

    // NOTE: Assume channels are ordered with the MPI rank order, i.e., the higher the rank the higher the channel number.
    // E.g.: for 40 channels and 4 workers, rank 0 has channels 0-9, rank 1: 10-19, rank 2: 20-29, and rank 3: 30-39.

    if (gradientType == 0) {
    // Forawrd difference scheme.

        size_t localChannelNumber = 0;
        for (size_t i = 0; i < nParametersTotal; i += 2) {

            bool lastLocalChannel = (localChannelNumber == nChannelsLocal - 1);

            if (lastLocalChannel) {
                size_t shiftedIndex = i + nextChannelIndexShift;

                if (shiftedIndex < nParametersTotal) {
                // Reached last local channel - shift the 'next' index.
                    // Real part.
                    leftIndexGlobal[i] = i;
                    rightIndexGlobal[i] = shiftedIndex;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
                    rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;
                } else {
                // Reached last global channel - do not add constraints - it is already coupled with previous one.
                    // Real part.
                    leftIndexGlobal[i] = -1;
                    rightIndexGlobal[i] = -1;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = -1;
                    rightIndexGlobal[i + 1] = -1;
                }

            } else {
                // Real part.
                leftIndexGlobal[i] = i;
                rightIndexGlobal[i] = i + 2;
                // Imaginary part.
                leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
                rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;
            }

            if (lastLocalChannel) {
                // Reset local channel counter.
                localChannelNumber = 0;
            } else {
                localChannelNumber++;
            }
        }
    } else {
    // Central difference scheme.

        size_t localChannelNumber = 0;
        for (size_t i = 0; i < nParametersTotal; i += 2) {

            bool firstLocalChannel = (localChannelNumber == 0);
            bool lastLocalChannel = (localChannelNumber == nChannelsLocal - 1);

            int shiftedLeftIndex = i - nextChannelIndexShift;
            size_t shiftedRightIndex = i + nextChannelIndexShift;

            if (firstLocalChannel && lastLocalChannel) {
            // One local channel.
                if ((shiftedLeftIndex >= 0) && (shiftedRightIndex < nParametersTotal)) {
                    // Real part.
                    leftIndexGlobal[i] = shiftedLeftIndex;
                    rightIndexGlobal[i] = shiftedRightIndex;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
                    rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;

                } else {
                // First/last global channel - do not add constraints - it is already coupled with next/previous one.
                    // Real part.
                    leftIndexGlobal[i] = -1;
                    rightIndexGlobal[i] = -1;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = -1;
                    rightIndexGlobal[i + 1] = -1;
                }

            } else if (firstLocalChannel) {

                if (shiftedLeftIndex >= 0) {
                // First local channel - shift the 'left' index.
                    // Real part.
                    leftIndexGlobal[i] = shiftedLeftIndex;
                    rightIndexGlobal[i] = i + 2;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
                    rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;

                } else {
                // First/last global channel - do not add constraints - it is already coupled with next/previous one.
                    // Real part.
                    leftIndexGlobal[i] = -1;
                    rightIndexGlobal[i] = -1;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = -1;
                    rightIndexGlobal[i + 1] = -1;
                }
            } else if (lastLocalChannel) {

                if (shiftedRightIndex < nParametersTotal) {
                // Last local channel - shift the 'right' index.
                    // Real part.
                    leftIndexGlobal[i] = i - 2;
                    rightIndexGlobal[i] = shiftedRightIndex;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
                    rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;
                } else {
                // Reached last global channel - do not add constraints - it is already coupled with previous one.
                    // Real part.
                    leftIndexGlobal[i] = -1;
                    rightIndexGlobal[i] = -1;
                    // Imaginary part.
                    leftIndexGlobal[i + 1] = -1;
                    rightIndexGlobal[i + 1] = -1;
                }

            } else {
                // Real part.
                leftIndexGlobal[i] = i - 2;
                rightIndexGlobal[i] = i + 2;
                // Imaginary part.
                leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
                rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;
            }

            if (lastLocalChannel) {
                // Reset local channel counter.
                localChannelNumber = 0;
            } else {
                localChannelNumber++;
            }
        }
    }

    //-----------------------------------------------------------------------------
    // Adding Jacobian of the gradient to the matrix.
    //-----------------------------------------------------------------------------
    double cost = 0.;

    double matrix_val[2];
    matrix_val[0] = - smoothingWeight;
    matrix_val[1] = + smoothingWeight;

    for (size_t i = 0; i < nParametersTotal; i++) {

        //----------------------------------------------------
        // Adding matrix values.
        //----------------------------------------------------
        matrix.NewRow();

        // Global matrix column indexes.
        size_t globIndex[2];
        globIndex[0] = leftIndexGlobal[i];  // Left index: "minus" term in finite difference (e.g. -x[i] for x' =  x[i+1] - x[i])
        globIndex[1] = rightIndexGlobal[i]; // Right index: "plus" term in finite difference (e.g. x[i+1] for x' =  x[i+1] - x[i]).

        for (size_t k = 0; k < 2; k++) {
            if (globIndex[k] >= 0
                && globIndex[k] >= nParametersSmaller
                && globIndex[k] < nParametersSmaller + nParameters) {

                // Local matrix column index (at the current CPU).
                size_t localColumnIndex = globIndex[k] - nParametersSmaller;
                matrix.Add(matrix_val[k], localColumnIndex);
            }
        }

        //----------------------------------------------------
        // Adding the Right-Hand Side.
        //----------------------------------------------------
        double b_RHS_value = 0.;
        if (leftIndexGlobal[i] >= 0 && rightIndexGlobal[i] >= 0) {
            b_RHS_value = - smoothingWeight * (x0[rightIndexGlobal[i]] - x0[leftIndexGlobal[i]]);
        } else {
            ASKAPCHECK(leftIndexGlobal[i] == -1 && rightIndexGlobal[i] == -1,
                    "Wrong finite difference indexes: " << i << ", " << leftIndexGlobal[i] << ", " << rightIndexGlobal[i]);
        }
        size_t b_index = matrix.GetCurrentNumberRows() - 1;
        b_RHS[b_index] = b_RHS_value;

        cost += b_RHS_value * b_RHS_value;
    }

    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Smoothness constraints cost = " << cost / (smoothingWeight * smoothingWeight));

    matrix.Finalize(nParameters);
}

}}}
