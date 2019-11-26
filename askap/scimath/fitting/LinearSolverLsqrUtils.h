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
#ifndef LinearSolverLsqrUtils_h_
#define LinearSolverLsqrUtils_h_

#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/lsqr_solver/SparseMatrix.h>

#include <vector>
#include <string>

namespace askap { namespace scimath { namespace lsqrutils {

    /// @brief Used for sorting gain names to have continuous channel number.
    /// @details Continuous channel number is needed to apply smoothing constraints.
    bool compareGainNames(const std::string& gainA, const std::string& gainB);

    /// @brief Builds a sparse matrix for the LSQR solver.
    /// @details Copy matrix elements from normal matrix (map of map of matrixes) to the LSQR solver sparse matrix (in CSR format).
    /// @param[in] ne Normal equation.
    /// @param[in] indices List of gain name/index pairs.
    /// @param[in] matrix The output sparse matrix.
    /// @param[in] nParameters Local number of parameters (at the current worker).
    /// @param[in] matrixIsParallel Flag for whether the matrix is parallel (columns distributed among workers).
    void buildLSQRSparseMatrix(const INormalEquations &ne,
                               const std::vector<std::pair<string, int> > &indices,
                               lsqr::SparseMatrix &matrix,
                               size_t nParameters,
                               bool matrixIsParallel);

    /// @brief Returns a current solution vector of doubles.
    /// @param[in] indices List of gain name/index pairs (note two solution parameters per gain - real & imaginary part).
    /// @param[in] params Normal equation parameters.
    /// @param[in] solution A container where the solution will be returned.
    void getCurrentSolutionVector(const std::vector<std::pair<std::string, int> >& indices,
                                  const Params& params,
                                  std::vector<double>& solution);

    /// @brief Calculates the smoothing weight for given major loop iteration.
    /// @details This smoothing weight is used to weight the smoothing constraints in the cost function.
    /// @param[in] parameters Configuration parameters.
    /// @param[in] majorLoopIterationNumber Major loop iteration number.
    /// @return The smoothing weight.
    double getSmoothingWeight(const std::map<std::string, std::string>& parameters,
                              size_t majorLoopIterationNumber);

    /// @brief Adding smoothness constraints to the system of equations.
    /// @details Extends the matrix and right-hand side with smoothness constraints,
    /// which are defined for the least squares minimization framework.
    /// @param[in] matrix The matrix where constraints will be added.
    /// @param[in] b_RHS The right-hand side where constraints will be added.
    /// @param[in] indices List of gain name/index pairs (note two parameters in x0 per gain: real & imaginary parts).
    /// @param[in] x0 The current global solution (at all workers).
    /// @param[in] nParameters Local number of parameters (at the current worker).
    /// @param[in] nChannels The total number of channels.
    /// @param[in] smoothingWeight The smoothing weight.
    /// @param[in] smoothingType The type of gradient approximation (0 - forward difference, 1- central difference, 2 - Laplacian).
    void addSmoothnessConstraints(lsqr::SparseMatrix& matrix,
                                  lsqr::Vector& b_RHS,
                                  const std::vector<std::pair<std::string, int> >& indices,
                                  const std::vector<double>& x0,
                                  size_t nParameters,
                                  size_t nChannels,
                                  double smoothingWeight,
                                  int smoothingType);
}}}
#endif
