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
#ifndef LinearSolverUtils_h_
#define LinearSolverUtils_h_

#include <fitting/INormalEquations.h>
#include <fitting/Params.h>

#include <vector>
#include <string>

#include <gsl/gsl_vector.h>

namespace askap { namespace scimath { namespace solverutils {

    typedef std::vector<double> Vector;

    /// @brief Assignment routines for populate_B function.
    void assign_to_lsqr_vector(Vector &B, std::size_t index, double elem);
    void assign_to_gsl_vector(gsl_vector *B, std::size_t index, double elem);

    /// @brief Retrieval routines for update_solution function.
    double retrieve_from_lsqr_vector(const Vector &X, std::size_t index);
    double retrieve_from_gsl_vector(const gsl_vector *X, std::size_t index);

    /// @brief Populate the right-hand side vector b (in the system of equations Ax = b).
    /// @param[in] ne Normal equations.
    /// @param[in] indices List of gain name/index pairs.
    /// @param[in] B Data holder.
    /// @param[in] assignment Assignment function.
    template <typename DataHolder, typename AssignmentFunc>
    size_t populate_B(const INormalEquations &ne,
                      const std::vector<std::pair<string, int> > &indices,
                      DataHolder &B,
                      AssignmentFunc assignment);

    /// @brief Update the solution with given perturbation.
    /// @param[in] indices List of gain name/index pairs.
    /// @param[in] params Solution accessor.
    /// @param[in] delta_X Solution perturbation (obtained at the current major iteration).
    /// @param[in] retrieval Retrieval function.
    template <typename DataHolder, typename RetrievalFunc>
    void update_solution(const std::vector<std::pair<string, int> > &indices,
                         Params& params,
                         const DataHolder &delta_X,
                         RetrievalFunc retrieval);

}}}
#endif
