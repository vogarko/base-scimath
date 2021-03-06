/// @file
/// @brief Naming convention for calibratable parameters
/// @details It is handy to use the same names of the calibratable parameters
/// in different parts of the code, e.g. when they're written to a parset file or
/// added as a model parameter. This class holds methods forming the name out of
/// antenna/beam/polarisation indices and parsing the string name to get these
/// indices back. 
///
/// @copyright (c) 2011 CSIRO
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>

#include <askap/scimath/fitting/CalParamNameHelper.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>

namespace askap {

namespace scimath {

std::string CalParamNameHelper::paramName(const JonesIndex &index, casacore::Stokes::StokesTypes par, bool isBP)
{
   std::string res;   
   if ((par == casacore::Stokes::XX) || (par == casacore::Stokes::YY)) {
       res = "gain."; 
       res += (par == casacore::Stokes::XX ? "g11." : "g22.");
   } else if ((par == casacore::Stokes::XY) || (par == casacore::Stokes::YX)) {
       res = "leakage."; 
       res += (par == casacore::Stokes::XY ? "d12." : "d21.");
   } else {
       ASKAPTHROW(AskapError, 
           "Unsupported polarisation descriptor passed to ParsetCalSolutionAccessor::paramName, only XX,XY,YX and YY are allowed");
   }
   
   return (isBP ? bpPrefix() + res : res) + utility::toString<casacore::Short>(index.antenna()) +"."+ utility::toString<casacore::Short>(index.beam());
}

std::string CalParamNameHelper::ionoParamName(const JonesIndex &index, const int paramType)
{
   std::string res;   
   if (paramType == 0) {
       res = "ionosphere.linear."; 
   } else {
       ASKAPTHROW(AskapError,  "Unsupported parameterisation type, only type 0 is allowed (linear)");
   }
   
   return res+utility::toString<casacore::Short>(index.antenna())+"."+utility::toString<casacore::Short>(index.beam());
}

std::pair<JonesIndex, casacore::Stokes::StokesTypes> CalParamNameHelper::parseParam(const std::string &name)
{

  size_t startPos = 0;
  if (bpParam(name)) {
      // to ignore bandpass prefix
      startPos = bpPrefix().size();
  }
  const size_t pos = name.find(".",startPos);
  ASKAPCHECK((pos != std::string::npos) && (pos + 1 != name.size()), 
             "Parameter name should be in the form something.something.ant.beam; you have "<<name);
  const std::string what = name.substr(startPos,pos - startPos);
  ASKAPCHECK((what == "gain") || (what == "leakage") || (what == "ionosphere"),
             "Only gain, leakage and ionospheric parameters are supported, you have "<<name);
  const size_t pos2 = name.find(".", pos + 1);
  ASKAPCHECK((pos2 != std::string::npos) && (pos2 + 1 != name.size()) && (pos + 1 != pos2), 
             "Parameter name should be in the form something.something.ant.beam; you have "<<name);
  const std::string type = name.substr(pos + 1, pos2 - pos - 1);
  casacore::Stokes::StokesTypes polDescriptor;
  if (what  == "gain") {
      ASKAPCHECK((type == "g11") || (type == "g22"), "Unrecognised polarisation product "<<type<<" in "<<name);
      polDescriptor = (type == "g11" ? casacore::Stokes::XX : casacore::Stokes::YY);
  } else if (what == "leakage") {
      ASKAPCHECK((type == "d12") || (type == "d21"), "Unrecognised polarisation product "<<type<<" in "<<name);
      polDescriptor = (type == "d12" ? casacore::Stokes::XY : casacore::Stokes::YX);
  } else if (what == "ionosphere") {
      ASKAPCHECK(type == "linear", "Unrecognised coefficient type "<<type<<" in "<<name);
      polDescriptor = casacore::Stokes::Undefined;
  } else {
     ASKAPTHROW(AskapError, "This line should never be executed!");
  }
  const size_t pos3 = name.find(".", pos2 + 1);
  ASKAPCHECK((pos3 != std::string::npos) && (pos3 + 1 != name.size()) && (pos2 + 1 != pos3), 
             "Parameter name should be in the form something.something.ant.beam; you have "<<name);
  const casacore::Short ant = utility::fromString<casacore::Short>(name.substr(pos2 + 1, pos3 - pos2 - 1));
  const casacore::Short beam = utility::fromString<casacore::Short>(name.substr(pos3 + 1));
  return std::pair<JonesIndex, casacore::Stokes::StokesTypes>(JonesIndex(ant,beam),polDescriptor);
}

bool CalParamNameHelper::bpParam(const std::string &name)
{
  return name.find(bpPrefix()) == 0;
}

bool CalParamNameHelper::ionoParam(const std::string &name)
{
  return name.find(ionoPrefix()) == 0;
}

std::string CalParamNameHelper::addChannelInfo(const std::string &name, casacore::uInt chan)
{
  return name + "." + utility::toString<casacore::uInt>(chan);
}

std::pair<casacore::uInt, std::string> CalParamNameHelper::extractChannelInfo(const std::string &name)
{
  size_t pos = name.rfind(".");
  ASKAPCHECK(pos != std::string::npos, "Expect dot in the parameter name passed to extractChannelInfo, name="<<name);
  ASKAPCHECK(pos + 1 != name.size(), "Parameter name="<<name<<" ends with a dot");
  return std::pair<casacore::uInt, std::string>(utility::fromString<casacore::uInt>(name.substr(pos+1)),name.substr(0,pos));
}

std::string CalParamNameHelper::extractBaseParamName(const std::string &name)
{
    if (!CalParamNameHelper::bpParam(name)) {
        return name;
    } else {
        return extractChannelInfo(name).second;
    }
}

} // namespace scimath

} // namespace askap
