/// @file
///
/// Expresses a polynomial equation
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
///

#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/PolynomialEquation.h>
#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/DesignMatrix.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/scimath/Mathematics/AutoDiff.h>
#include <casacore/scimath/Mathematics/AutoDiffMath.h>

#include <cmath>

using askap::scimath::DesignMatrix;

namespace askap
{
  namespace scimath
  {

    PolynomialEquation::PolynomialEquation(const Params& ip,
      casacore::Vector<double>& data,
      casacore::Vector<double>& weights,
      casacore::Vector<double>& arguments,
      casacore::Vector<double>& model) : Equation(ip), GenericEquation(ip), itsData(data),
      itsWeights(weights), itsArguments(arguments), itsModel(model)
    {
    };

    PolynomialEquation::PolynomialEquation(casacore::Vector<double>& data,
      casacore::Vector<double>& weights,
      casacore::Vector<double>& arguments,
      casacore::Vector<double>& model) :  itsData(data),
      itsWeights(weights), itsArguments(arguments), itsModel(model)
    {
    };

    PolynomialEquation::PolynomialEquation(const PolynomialEquation& other) :
       Equation(other), GenericEquation(other)
    {
      operator=(other);
    }

    PolynomialEquation& PolynomialEquation::operator=(const PolynomialEquation& other)
    {
      if(this!=&other)
      {
        static_cast<Equation*>(this)->operator=(other);
        itsData=other.itsData;
        itsWeights=other.itsWeights;
        itsArguments=other.itsArguments;
        itsModel=other.itsModel;
      }
      return *this;
    }

    PolynomialEquation::~PolynomialEquation()
    {
    }
    
    Params PolynomialEquation::defaultParameters()
    {
      Params ip;
      ip.add("poly");
      return ip;
    }
     
    void PolynomialEquation::init()
    {
    }

    void PolynomialEquation::predict() const
    {
    	itsModel.set(0.0);

      std::vector<std::string> completions(parameters().completions("poly"));
      if(completions.size()>0) {
// Loop over all polynomials adding to the values. 
        for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
        {
          string polyName("poly"+(*it));
          const casacore::Vector<double> par(parameters().value(polyName));
          this->calcPoly(itsArguments, par, itsModel);
        }
      }
    }

    void PolynomialEquation::calcGenericEquations(GenericNormalEquations& ne) const
    {
    	itsModel.set(0.0);

    	std::vector<std::string> completions(parameters().completions("poly"));
      if(completions.size()>0) {
        DesignMatrix designmatrix; // old parameters: parameters();

        casacore::Vector<double> values(itsData.size());
        values=0.0;

// Loop over all polynomials adding to the values. 
        for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
        {
          string polyName("poly"+(*it));
          const casacore::Vector<double> par(parameters().value(polyName));
          casacore::Matrix<double> valueDerivs(itsData.size(), par.size());
          this->calcPoly(itsArguments, par, itsModel);
          this->calcPolyDeriv(itsArguments, par, valueDerivs);
          designmatrix.addDerivative(polyName, valueDerivs);
        }
        casacore::Vector<double> residual(itsData.copy());
  
        residual-=itsModel;
        designmatrix.addResidual(residual, itsWeights);
        ne.add(designmatrix);
      }
    };

    void PolynomialEquation::calcPoly(const casacore::Vector<double>& x,
      const casacore::Vector<double>& parameters,
      casacore::Vector<double>& values)
    {
      for (size_t ix=0; ix<x.size(); ++ix)
      {
        for (size_t ipar=0; ipar<parameters.size(); ++ipar)
        {
          /// @todo Optimize calculation of values
          values[ix]+=parameters[ipar]*std::pow(x[ix], int(ipar));
        }
      }
    }
    void PolynomialEquation::calcPolyDeriv(const casacore::Vector<double>& x,
      const casacore::Vector<double>& parameters,
      casacore::Matrix<double>& valueDerivs)
    {
      const uint nPoly=parameters.size();
      const uint n=x.size();
      for (uint ix=0; ix<n; ++ix)
      {
        for (uint ipar=0; ipar<nPoly; ++ipar)
        {
          /// @todo Optimize calculation of derivatives
          valueDerivs(ix,ipar)=std::pow(x[ix], int(ipar));
        }
      }
    }

    Equation::ShPtr PolynomialEquation::clone() const
    {
      return Equation::ShPtr(new PolynomialEquation(*this));
    }

  }

}
