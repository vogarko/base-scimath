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
#include <askap/scimath/fitting/LinearSolverLsqrUtils.h>
#include <askap/scimath/fitting/LinearSolverUtils.h>
#include <askap/scimath/fitting/CalParamNameHelper.h>
#include <askap/scimath/lsqr_solver/ParallelTools.h>

#include <askap/AskapError.h>
#include <profile/AskapProfiler.h>

#include <askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".lsqrutils");

#include <iostream>
#include <string>
#include <map>

namespace askap { namespace scimath { namespace lsqrutils {

bool compareGainNames(const std::string& gainA, const std::string& gainB) {
    try {
        std::pair<casa::uInt, std::string> paramInfoA = CalParamNameHelper::extractChannelInfo(gainA);
        std::pair<casa::uInt, std::string> paramInfoB = CalParamNameHelper::extractChannelInfo(gainB);

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

void buildLSQRSparseMatrix(const GenericNormalEquations& gne,
                           const std::vector<std::pair<std::string, int> > &indices,
                           lsqr::SparseMatrix &matrix,
                           size_t nParameters,
                           bool matrixIsParallel)
{
    size_t nParametersTotal;
    size_t nParametersSmaller;

#ifdef HAVE_MPI
    MPI_Comm workersComm = matrix.GetComm();
    if (workersComm != MPI_COMM_NULL) {
        int myrank, nbproc;
        MPI_Comm_rank(workersComm, &myrank);
        MPI_Comm_size(workersComm, &nbproc);
        nParametersTotal = lsqr::ParallelTools::get_total_number_elements(nParameters, nbproc, workersComm);
        nParametersSmaller = lsqr::ParallelTools::get_nsmaller(nParameters, myrank, nbproc, workersComm);
    }
    else
#endif
    {
        nParametersTotal = nParameters;
        nParametersSmaller = 0;
    }

    //------------------------------------------------------------------------------------------------------------------------

    if (matrixIsParallel) {
        // Adding starting matrix empty rows, i.e., the rows in a big block-diagonal matrix above the current block.
        for (size_t i = 0; i < nParametersSmaller; ++i) {
            matrix.NewRow();
        }
    }

    if (gne.indexedNormalMatrixInitialized()) {
    // new matrix format
        ASKAPCHECK(gne.indexedDataVectorInitialized(), "Indexed data vector is not initialized!");

        size_t nChannelsLocal = gne.getNumberLocalChannels();
        size_t nBaseParameters = gne.getNumberBaseParameters();

        ASKAPCHECK(2 * nBaseParameters * nChannelsLocal == nParameters, "Wrong number of parameters in buildLSQRSparseMatrix!");

        for (size_t chan = 0; chan < nChannelsLocal; chan++) {
            for (size_t row = 0; row < nBaseParameters; row++) {
                // Adding 2x2 complex elements (indexed with i, j).
                for (size_t i = 0; i < 2; i++) {
                    matrix.NewRow();
                    // Loop over matrix columns.
                    for (size_t col = 0; col < nBaseParameters; col++) {
                        const IndexedNormalMatrix::elem_type& elem = gne.indexedNormalMatrix(col, row, chan);
                        for (size_t j = 0; j < 2; j++) {
                            size_t colIndex = 2 * col + j + chan * (2 * nBaseParameters);
                            matrix.Add(elem.data[i][j], colIndex);
                        }
                    }
                }
            }
        }

    } else {
    // old matrix format

        std::map<std::string, size_t> indicesMap;
        for (std::vector<std::pair<std::string, int> >::const_iterator it = indices.begin();
             it != indices.end(); ++it) {
            indicesMap[it->first] = (size_t)(it->second);
        }

        // Loop over matrix rows.
        for (std::vector<std::pair<std::string, int> >::const_iterator indit1 = indices.begin();
                indit1 != indices.end(); ++indit1) {

            const auto colItBeg = gne.getNormalMatrixRowBegin(indit1->first);
            const auto colItEnd = gne.getNormalMatrixRowEnd(indit1->first);

            if (colItBeg != colItEnd) {
                const size_t nrow = colItBeg->second.nrow();
                for (size_t row = 0; row < nrow; ++row) {
                    matrix.NewRow();
                    // Loop over column elements.
                    for (auto colIt = colItBeg;
                            colIt != colItEnd; ++colIt) {

                        const std::map<std::string, size_t>::const_iterator indicesMapIt = indicesMap.find(colIt->first);
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

    for (std::vector<std::pair<std::string, int> >::const_iterator indit = indices.begin();
            indit != indices.end(); ++indit) {
        casa::IPosition vecShape(1, params.value(indit->first).nelements());
        casa::Vector<double> value(params.value(indit->first).reform(vecShape));
        for (size_t i = 0; i < value.nelements(); ++i) {
            solution[indit->second + i] = value(i);
        }
    }
}

void getCurrentSolutionVector(const GenericNormalEquations& gne,
                              const Params& params,
                              std::vector<double>& solution)
{
    size_t nChannelsLocal = gne.getNumberLocalChannels();
    size_t nBaseParameters = gne.getNumberBaseParameters();

    for (size_t i = 0; i < nBaseParameters; i++) {
        for (size_t chan = 0; chan < nChannelsLocal; chan++) {
            std::string paramName = gne.getFullParameterName(i, chan);

            auto *data = params.value(paramName).data();
            size_t index = i + nBaseParameters * chan;

            solution[2 * index] = data[0];
            solution[2 * index + 1] = data[1];
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

bool testMPIRankOrderWithChannels(int workerRank, size_t nChannelsLocal,
                                  const std::vector<std::pair<std::string, int> >& indices)
{
    for (std::vector<std::pair<std::string, int> >::const_iterator indit = indices.begin();
         indit != indices.end(); ++indit) {
        // Extracting channel number from gain name.
        std::string gainName = indit->first;
        std::pair<casa::uInt, std::string> paramInfo = CalParamNameHelper::extractChannelInfo(gainName);
        casa::uInt chan = paramInfo.first;

        // Testing that the channel number is aligned with MPI partitioning:
        // E.g.: for 40 channels and 4 workers, rank 0 has channels 0-9, rank 1: 10-19, rank 2: 20-29, and rank 3: 30-39.
        if (chan < workerRank * nChannelsLocal
            || chan >= (workerRank + 1) * nChannelsLocal) {
            return false;
        }
    }
    return true;
}

// Calculates the index shift for the next channel located on the next MPI rank.
// Assumes channels are ordered with MPI ranks.
size_t getNextChannelIndexShift(size_t nParametersLocal, size_t nChannelsLocal)
{
    return nParametersLocal - (nChannelsLocal - 1) * 2;
}

// Calculates matrix indexes for the Forward Difference (FD) gradient operator.
// For old normal matrix format.
void calculateIndexesFWD_old(size_t nParametersTotal,
                             size_t nParametersLocal,
                             size_t nChannelsLocal,
                             std::vector<int>& leftIndexGlobal,
                             std::vector<int>& rightIndexGlobal)
{
    size_t nextChannelIndexShift = getNextChannelIndexShift(nParametersLocal, nChannelsLocal);

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
}

// Calculates matrix indexes for the Forward Difference (FD) gradient operator.
// For new indexed normal matrix format (it has different column order).
void calculateIndexesFWD_new(size_t nParametersTotal,
                             size_t nParametersLocal,
                             size_t nChannelsLocal,
                             std::vector<int>& leftIndexGlobal,
                             std::vector<int>& rightIndexGlobal)
{
    size_t nextChannelIndexShift = nParametersLocal / nChannelsLocal;

    for (size_t i = 0; i < nParametersTotal; i += 2) {
        size_t shiftedIndex = i + nextChannelIndexShift;

        if (shiftedIndex < nParametersTotal) {
            // Real part.
            leftIndexGlobal[i] = i;
            rightIndexGlobal[i] = shiftedIndex;
            // Imaginary part.
            leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
            rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;

        } else {
        // Last global channel - do not add constraints - it is already coupled with previous one.
            // Real part.
            leftIndexGlobal[i] = -1;
            rightIndexGlobal[i] = -1;
            // Imaginary part.
            leftIndexGlobal[i + 1] = -1;
            rightIndexGlobal[i + 1] = -1;
        }
    }
}

// Calculates matrix indexes for the Forward Difference (FD) gradient operator.
void calculateIndexesFWD(size_t nParametersTotal,
                         size_t nParametersLocal,
                         size_t nChannelsLocal,
                         std::vector<int>& leftIndexGlobal,
                         std::vector<int>& rightIndexGlobal,
                         bool indexedNormalMatrixFormat)
{
    ASKAPCHECK(nParametersTotal == leftIndexGlobal.size()
               && nParametersTotal == rightIndexGlobal.size(), "Wrong vector size in calculateIndexesFWD!");

    if (indexedNormalMatrixFormat) {
        calculateIndexesFWD_new(nParametersTotal, nParametersLocal, nChannelsLocal, leftIndexGlobal, rightIndexGlobal);
    } else {
        calculateIndexesFWD_old(nParametersTotal, nParametersLocal, nChannelsLocal, leftIndexGlobal, rightIndexGlobal);
    }
}

// Calculates matrix indexes for the Central Difference (CD) gradient operator.
void calculateIndexesCD_old(size_t nParametersTotal,
                            size_t nParametersLocal,
                            size_t nChannelsLocal,
                            std::vector<int>& leftIndexGlobal,
                            std::vector<int>& rightIndexGlobal)
{
    size_t nextChannelIndexShift = getNextChannelIndexShift(nParametersLocal, nChannelsLocal);

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

// Calculates matrix indexes for the Central Difference (CD) gradient operator.
void calculateIndexesCD_new(size_t nParametersTotal,
                            size_t nParametersLocal,
                            size_t nChannelsLocal,
                            std::vector<int>& leftIndexGlobal,
                            std::vector<int>& rightIndexGlobal)
{
    size_t nextChannelIndexShift = nParametersLocal / nChannelsLocal;

    for (size_t i = 0; i < nParametersTotal; i += 2) {

        int shiftedLeftIndex = i - nextChannelIndexShift;
        size_t shiftedRightIndex = i + nextChannelIndexShift;

        if (shiftedLeftIndex >= 0 && shiftedRightIndex < nParametersTotal) {
            // Real part.
            leftIndexGlobal[i] = shiftedLeftIndex;
            rightIndexGlobal[i] = shiftedRightIndex;
            // Imaginary part.
            leftIndexGlobal[i + 1] = leftIndexGlobal[i] + 1;
            rightIndexGlobal[i + 1] = rightIndexGlobal[i] + 1;

        } else {
        // First or last global channel - do not add constraints - it is already coupled with the next one.
            // Real part.
            leftIndexGlobal[i] = -1;
            rightIndexGlobal[i] = -1;
            // Imaginary part.
            leftIndexGlobal[i + 1] = -1;
            rightIndexGlobal[i + 1] = -1;
        }
    }
}

// Calculates matrix indexes for the Central Difference (CD) gradient operator.
void calculateIndexesCD(size_t nParametersTotal,
                         size_t nParametersLocal,
                         size_t nChannelsLocal,
                         std::vector<int>& leftIndexGlobal,
                         std::vector<int>& rightIndexGlobal,
                         bool indexedNormalMatrixFormat)
{
    ASKAPCHECK(nParametersTotal == leftIndexGlobal.size()
               && nParametersTotal == rightIndexGlobal.size(), "Wrong vector size in calculateIndexesCD!");

    if (indexedNormalMatrixFormat) {
        calculateIndexesCD_new(nParametersTotal, nParametersLocal, nChannelsLocal, leftIndexGlobal, rightIndexGlobal);
    } else {
        calculateIndexesCD_old(nParametersTotal, nParametersLocal, nChannelsLocal, leftIndexGlobal, rightIndexGlobal);
    }
}

/// @brief Adding smoothness constraints to the system of equations.
void addSmoothnessConstraints(lsqr::SparseMatrix& matrix,
                              lsqr::Vector& b_RHS,
                              const std::vector<double>& x0,
                              size_t nParameters,
                              size_t nChannels,
                              double smoothingWeight,
                              int smoothingType,
                              bool indexedNormalMatrixFormat)
{
    ASKAPCHECK(nChannels > 1, "Wrong number of channels for smoothness constraints!");

    int myrank;
    size_t nParametersTotal;

#ifdef HAVE_MPI
    MPI_Comm workersComm = matrix.GetComm();
    if (workersComm != MPI_COMM_NULL) {
        int nbproc;
        MPI_Comm_rank(workersComm, &myrank);
        MPI_Comm_size(workersComm, &nbproc);
        nParametersTotal = lsqr::ParallelTools::get_total_number_elements(nParameters, nbproc, workersComm);
    }
    else
#endif
    {
        myrank = 0;
        nParametersTotal = nParameters;
    }

    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Adding smoothness constraints, with weight = " << smoothingWeight << ", smoothingType = " << smoothingType);

    //-----------------------------------------------------------------------------
    // Assume the same number of channels at every CPU.
    size_t nChannelsLocal = nChannels / (nParametersTotal / nParameters);

    if (myrank == 0) ASKAPLOG_DEBUG_STR(logger, "nChannelsLocal = " << nChannelsLocal);

    //-----------------------------------------------------------------------------
    size_t nDiag;
    if (smoothingType == 0 || smoothingType == 1) {
        nDiag = 2;
    } else if (smoothingType == 2) {
        nDiag = 3;
    } else {
        throw std::invalid_argument("Unknown smoothing type!");
    }

    std::vector<std::vector<int> > columnIndexGlobal(nDiag, std::vector<int>(nParametersTotal));
    std::vector<double> matrixValue(nDiag);
    {
        std::vector<int> leftIndexGlobal(nParametersTotal);
        std::vector<int> rightIndexGlobal(nParametersTotal);

        if (smoothingType == 0 || smoothingType == 1) {
            if (smoothingType == 0) {
            // Forawrd difference scheme.
                calculateIndexesFWD(nParametersTotal, nParameters, nChannelsLocal, leftIndexGlobal, rightIndexGlobal, indexedNormalMatrixFormat);
            }
            else if (smoothingType == 1) {
            // Central difference scheme.
                calculateIndexesCD(nParametersTotal, nParameters, nChannelsLocal, leftIndexGlobal, rightIndexGlobal, indexedNormalMatrixFormat);
            }

            columnIndexGlobal[0] = leftIndexGlobal;
            columnIndexGlobal[1] = rightIndexGlobal;

            matrixValue[0] = - smoothingWeight;
            matrixValue[1] = + smoothingWeight;
        }
        else if (smoothingType == 2) {
        // Laplacian.
            // Utilize that the left and right indexes in Laplacian are the same as in the Central Difference (CD) scheme.
            calculateIndexesCD(nParametersTotal, nParameters, nChannelsLocal, leftIndexGlobal, rightIndexGlobal, indexedNormalMatrixFormat);

            std::vector<int> middleIndexGlobal(nParametersTotal);
            for (size_t i = 0; i < nParametersTotal; i++) {
                if (leftIndexGlobal[i] >= 0) {
                    assert(rightIndexGlobal[i] >= 0);
                    middleIndexGlobal[i] = i;
                } else {
                    assert(leftIndexGlobal[i] == -1);
                    assert(rightIndexGlobal[i] == -1);
                    middleIndexGlobal[i] = -1;
                }
            }

            columnIndexGlobal[0] = leftIndexGlobal;
            columnIndexGlobal[1] = middleIndexGlobal;
            columnIndexGlobal[2] = rightIndexGlobal;

            // 1D Laplacian kernel = [1 -2 1].
            matrixValue[0] = smoothingWeight;
            matrixValue[1] = - 2. * smoothingWeight;
            matrixValue[2] = smoothingWeight;
        }
    }

    //-----------------------------------------------------------------------------
    // Adding Jacobian of the gradient/Laplacian to the matrix.
    //-----------------------------------------------------------------------------
    matrix.addParallelSparseOperator(nDiag, nParameters, columnIndexGlobal, matrixValue);

    //-----------------------------------------------------------------------------
    // Adding the Right-Hand Side.
    //-----------------------------------------------------------------------------
    size_t b_size0 = b_RHS.size();
    b_RHS.resize(b_RHS.size() + nParametersTotal);
    assert(matrix.GetCurrentNumberRows() == b_RHS.size());

    double cost = 0.;
    for (size_t i = 0; i < nParametersTotal; i++) {
        double Ax0 = 0.;
        for (size_t k = 0; k < nDiag; k++) {
            Ax0 += matrixValue[k] * x0[columnIndexGlobal[k][i]];
        }

        // b = - F(x0) = - A.x0.
        double b_RHS_value = - Ax0;

        size_t b_index = b_size0 + i;
        b_RHS[b_index] = b_RHS_value;

        cost += b_RHS_value * b_RHS_value;
    }
    if (myrank == 0) ASKAPLOG_INFO_STR(logger, "Smoothness constraints cost = " << cost / (smoothingWeight * smoothingWeight));
}

double calculateCost(const std::vector<double> &b_RHS)
{
    double cost = 0.;
    for (size_t i = 0; i < b_RHS.size(); ++i) {
        cost += b_RHS[i] * b_RHS[i];
    }
    return cost;
}

}}}
