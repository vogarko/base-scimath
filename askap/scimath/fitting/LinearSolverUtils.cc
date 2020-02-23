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
#include <askap/scimath/fitting/LinearSolverUtils.h>

#include <askap/askap/AskapError.h>

namespace askap { namespace scimath { namespace solverutils {

double getParameter(std::string parname, const SolverParameters& parameters, double defaultValue)
{
    double value = defaultValue;
    if (parameters.count(parname) > 0) {
        value = std::atof(parameters.at(parname).c_str());
    }
    return value;
}

int getParameter(std::string parname, const SolverParameters& parameters, int defaultValue)
{
    int value = defaultValue;
    if (parameters.count(parname) > 0) {
        value = std::atoi(parameters.at(parname).c_str());
    }
    return value;
}

bool getParameter(std::string parname, const SolverParameters& parameters, bool defaultValue)
{
    bool value = defaultValue;
    if (parameters.count(parname) > 0) {
        if (parameters.at(parname) == "true") {
            value = true;
        } else if (parameters.at(parname) == "false") {
            value = false;
        }
    }
    return value;
}

void assign_to_lsqr_vector(Vector &B, std::size_t index, double elem)
{
    B[index] = elem;
}

void assign_to_gsl_vector(gsl_vector *B, std::size_t index, double elem)
{
    gsl_vector_set(B, index, elem);
}

double retrieve_from_lsqr_vector(const Vector &X, std::size_t index)
{
    return X[index];
}

double retrieve_from_gsl_vector(const gsl_vector *X, std::size_t index)
{
    return gsl_vector_get(X, index);
}

template <typename DataHolder, typename AssignmentFunc>
size_t populate_B(const INormalEquations &ne,
                  const std::vector<std::pair<std::string, int> > &indices,
                  DataHolder &B,
                  AssignmentFunc assignment)
{
    size_t counter = 0;
    for (std::vector<std::pair<std::string, int> >::const_iterator indit = indices.begin();
            indit != indices.end(); ++indit) {
        const casa::Vector<double> &dv = ne.dataVector(indit->first);
        for (size_t row = 0; row < dv.nelements(); ++row) {
             const double elem = dv(row);
             ASKAPCHECK(!std::isnan(elem), "Data vector seems to have NaN for row = " << row << ", this shouldn't happen!");
             assignment(B, row + indit->second, elem);
             counter++;
        }
    }
    return counter;
}

template
size_t populate_B<Vector, typeof(&assign_to_lsqr_vector)>(const INormalEquations &,
                                                          const std::vector<std::pair<std::string, int> > &,
                                                          Vector &,
                                                          typeof(&assign_to_lsqr_vector));

template
size_t populate_B<gsl_vector *, typeof(&assign_to_gsl_vector)>(const INormalEquations &,
                                                               const std::vector<std::pair<std::string, int> > &,
                                                               gsl_vector *&,
                                                               typeof(&assign_to_gsl_vector));

template <typename DataHolder, typename RetrievalFunc>
void update_solution(const std::vector<std::pair<std::string, int> > &indices,
                     Params &params,
                     const DataHolder &delta_X,
                     RetrievalFunc retrieval)
{
    // Exploit reference semantics of casa::Array.
    std::vector<std::pair<std::string, int> >::const_iterator indit;
    for (indit = indices.begin(); indit != indices.end(); ++indit) {
        casa::IPosition vecShape(1, params.value(indit->first).nelements());
        casa::Vector<double> value(params.value(indit->first).reform(vecShape));
        for (size_t i = 0; i < value.nelements(); ++i) {
            double adjustment = retrieval(delta_X, indit->second + i);
            //ASKAPCHECK(!std::isnan(adjustment), "Solution resulted in NaN as an update for parameter " << (indit->second + i));
            if (!std::isnan(adjustment)) value(i) += adjustment;
        }
    }
}

template
void update_solution<Vector, typeof(&retrieve_from_lsqr_vector)>(const std::vector<std::pair<std::string, int> > &,
                                                                 Params &,
                                                                 const Vector &,
                                                                 typeof(&retrieve_from_lsqr_vector));

template
void update_solution<gsl_vector *, typeof(&retrieve_from_gsl_vector)>(const std::vector<std::pair<std::string, int> > &,
                                                                      Params &,
                                                                      gsl_vector * const &,
                                                                      typeof(&retrieve_from_gsl_vector));


}}}
