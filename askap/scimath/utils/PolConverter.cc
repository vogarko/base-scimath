/// @file
/// @brief Converter of polarisation frames
/// @details This is the class which handles polarisation frame conversion and contains some
/// helper methods related to it (i.e. converting strings into Stokes enums). It may eventually
/// replace or become derived from IPolSelector, which is not used at the moment.
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
///


#include <askap/scimath/utils/PolConverter.h>
#include <askap/askap/AskapError.h>

using namespace askap;
using namespace askap::scimath;

/// @brief constructor of the converter between two frames
/// @details
/// @param[in] polFrameIn input polarisation frame defined as a vector of Stokes enums
/// @param[in] polFrameOut output polarisation frame defined as a vector of Stokes enums
/// @param[in] checkUnspecifiedProducts if true (default), the code checks that all
///            polarisation products missing in the input frame are multiplied by 0 (and
///            therefore don't affect the result), see itsCheckUnspecifiedProducts for more info
PolConverter::PolConverter(const casacore::Vector<casacore::Stokes::StokesTypes> &polFrameIn,
               const casacore::Vector<casacore::Stokes::StokesTypes> &polFrameOut,
               bool checkUnspecifiedProducts) : itsVoid(false),
               itsTransform(polFrameOut.nelements(),polFrameIn.nelements(),casacore::Complex(0.,0.)),
               itsPolFrameIn(polFrameIn), itsPolFrameOut(polFrameOut),
               itsCheckUnspecifiedProducts(checkUnspecifiedProducts),
               itsSwapPols(false),itsInitialised(false)
{
  if (equal(polFrameIn, polFrameOut)) {
      itsVoid = true;
  } else {
    for (casacore::uInt pol=0; pol<polFrameIn.nelements(); ++pol) {
         ASKAPCHECK(isValid(polFrameIn[pol]), "Conversion is unsupported for polarisation product "<<
                    int(polFrameIn[pol])<<" ("<<casacore::Stokes::type(polFrameIn[pol])<<")");
    }
    for (casacore::uInt pol=0; pol<polFrameOut.nelements(); ++pol) {
         ASKAPCHECK(isValid(polFrameOut[pol]), "Conversion is unsupported for polarisation product "<<
                    int(polFrameOut[pol])<<" ("<<casacore::Stokes::type(polFrameOut[pol])<<")");
    }
    itsInputLinear = isLinear(itsPolFrameIn);
    itsOutputLinear = isLinear(itsPolFrameOut);
  }
}

/// @brief default constructor - no conversion
/// @details Constructed via this method the object passes all visibilities intact
PolConverter::PolConverter() : itsVoid(true), itsCheckUnspecifiedProducts(false)
{
}

/// @brief copy constructor
/// @details
/// @param[in] other, object to construct from
PolConverter::PolConverter(const PolConverter & other):
itsVoid(other.itsVoid),
itsTransform(other.itsTransform.copy()),
itsPARotation(other.itsPARotation.copy()),
itsPolFrameIn(other.itsPolFrameIn.copy()),
itsPolFrameOut(other.itsPolFrameOut.copy()),
itsCheckUnspecifiedProducts(other.itsCheckUnspecifiedProducts),
itsInputLinear(other.itsInputLinear),itsOutputLinear(other.itsOutputLinear),
itsSwapPols(other.itsSwapPols),itsInitialised(other.itsInitialised)
{
}

/// @brief assignment operator
/// @details
/// @param[in] other, object to assign from
PolConverter & PolConverter::operator=(const PolConverter & other)
{
    if (this == &other) return *this;
    itsVoid = other.itsVoid;
    itsTransform.assign(other.itsTransform);
    itsPARotation.assign(other.itsPARotation);
    itsPolFrameIn.assign(other.itsPolFrameIn);
    itsPolFrameOut.assign(other.itsPolFrameOut);
    itsCheckUnspecifiedProducts = other.itsCheckUnspecifiedProducts;
    itsInputLinear = other.itsInputLinear;
    itsOutputLinear = other.itsOutputLinear;
    itsSwapPols = other.itsSwapPols;
    itsInitialised = other.itsInitialised;
    return *this;
}


/// @brief compare two vectors of Stokes enums
/// @param[in] first first polarisation frame
/// @param[in] second second polarisation frame
/// @return true if two given frames are the same, false if not.
bool PolConverter::equal(const casacore::Vector<casacore::Stokes::StokesTypes> &first,
                    const casacore::Vector<casacore::Stokes::StokesTypes> &second)
{
  if (first.nelements() != second.nelements()) {
      return false;
  }
  for (casacore::uInt pol = 0; pol<first.nelements(); ++pol) {
       if (first[pol]!=second[pol]) {
           return false;
       }
  }
  return true;
}

/// @brief main method doing conversion
/// @details Convert the given visibility vector between two polarisation frames supplied
/// in the constructor.
/// @param[in] vis visibility vector
/// @return converted visibility vector
/// @note vis should have the same size (<=4) as both polFrames passed in the constructor,
/// the output vector will have the same size.
casacore::Vector<casacore::Complex> PolConverter::operator()(casacore::Vector<casacore::Complex> vis) const
{
  if (itsVoid) {
      return vis;
  }
  casacore::Vector<casacore::Complex> res(itsTransform.nrow(),0.);
  convert(res,vis);

  return res;
}

/// @brief alternative method doing conversion
/// @details Convert the given visibility vector between two polarisation frames supplied
/// in the constructor.
/// @param[in] vis visibility vector
/// @param[out] out converted visibility vector
/// @return converted visibility vector
/// @note vis should have the same size (<=4) as both polFrames passed in the constructor,
/// the output vector should have the correct size.
casacore::Vector<casacore::Complex>& PolConverter::convert(casacore::Vector<casacore::Complex>& out,
    const casacore::Vector<casacore::Complex>& vis) const
{
    if (itsVoid) {
        out = vis;
        return out;
    }
    if (!itsInitialised) {
        fillMatrix();
    }
    ASKAPDEBUGASSERT(vis.nelements() == itsTransform.ncolumn());
    ASKAPDEBUGASSERT(out.nelements() == itsTransform.nrow());

    for (casacore::uInt row = 0; row<out.nelements(); ++row) {
         out[row] = 0.;
         for (casacore::uInt col = 0; col<vis.nelements(); ++col) {
              out[row] += itsTransform(row,col)*vis[col];
         }
    }
    return out;
}


/// @brief propagate noise through conversion
/// @details Convert the visibility noise visibility vector between two
/// polarisation frames supplied in the constructor.
/// @param[in] visNoise visibility noise vector
/// @return converted noise vector
/// @note visNoise should have the same size (<=4) as both polFrames passed in the constructor,
/// the output vector will have the same size. Real and imaginary parts of the output vectors are the noise
/// levels of real and imaginary parts of the visibility.
casacore::Vector<casacore::Complex> PolConverter::noise(casacore::Vector<casacore::Complex> visNoise) const
{
  if (itsVoid) {
      return visNoise;
  }
  casacore::Vector<casacore::Complex> res(itsTransform.nrow(),0.);
  noise(res,visNoise);

  return res;
}

/// @brief propagate noise through conversion
/// @details Convert the visibility noise visibility vector between two
/// polarisation frames supplied in the constructor.
/// @param[in] visNoise visibility noise vector
/// @param[out] converted noise vector
/// @note visNoise should have the same size (<=4) as both polFrames passed in the constructor,
/// the output vector will have the same size. Real and imaginary parts of the output vectors are the noise
/// levels of real and imaginary parts of the visibility.
casacore::Vector<casacore::Complex>& PolConverter::noise(casacore::Vector<casacore::Complex>& outNoise, const casacore::Vector<casacore::Complex>& visNoise) const
{
  if (itsVoid) {
      outNoise = visNoise;
      return outNoise;
  }
  if (!itsInitialised) {
      fillMatrix();
  }
  ASKAPDEBUGASSERT(visNoise.nelements() == itsTransform.ncolumn());
  ASKAPDEBUGASSERT(outNoise.nelements() == itsTransform.nrow());

  for (casacore::uInt row = 0; row<outNoise.nelements(); ++row) {
       float reNoise = 0.;
       float imNoise = 0.;
       for (casacore::uInt col = 0; col<visNoise.nelements(); ++col) {
            const casacore::Complex coeff = itsTransform(row,col);
            const casacore::Complex val = visNoise[col];
            reNoise += casacore::square(casacore::real(coeff)*casacore::real(val)) +
                       casacore::square(casacore::imag(coeff)*casacore::imag(val));
            imNoise += casacore::square(casacore::imag(coeff)*casacore::real(val)) +
                       casacore::square(casacore::real(coeff)*casacore::imag(val));
       }
       ASKAPDEBUGASSERT(reNoise >= 0.);
       ASKAPDEBUGASSERT(imNoise >= 0.);
       outNoise[row] = casacore::Complex(sqrt(reNoise),sqrt(imNoise));
  }

  return outNoise;
}

/// @brief build transformation matrix
/// @details This is the core of the algorithm, this method builds the transformation matrix
/// given the two frames. It now also includes optional rotation.
void PolConverter::fillMatrix() const
{
  ASKAPDEBUGASSERT(itsTransform.nrow() == itsPolFrameOut.nelements());
  ASKAPDEBUGASSERT(itsTransform.ncolumn() == itsPolFrameIn.nelements());
  // See Hamaker, Bregman and Sault, 1996, A&ASS, 117, 137 for matrix formalism of
  // the polarisation conversion

  // todo, check whether we can do the same in a more elegant and general way.
  // todo 2, could cache the stokes conversion matrix separately for multiplication with rotation matrix

  casacore::Matrix<casacore::Complex> T(4,4,0.);
  if (isStokes(itsPolFrameOut)) {
      if (itsInputLinear) {
          // linear to stokes
          T(0,0)=1.; T(0,3)=1.;
          T(1,0)=1.; T(1,3)=-1.;
          T(2,1)=1.; T(2,2)=1.;
          T(3,1)=casacore::Complex(0.,-1.); T(3,2)=casacore::Complex(0.,1.);
      } else if (isCircular(itsPolFrameIn)) {
          // circular to stokes
          T(0,0)=1.; T(0,3)=1.;
          T(1,1)=casacore::Complex(0.,-1.); T(1,2)=casacore::Complex(0.,1.);
          T(2,0)=1.; T(2,3)=-1.;
          T(3,1)=1.; T(3,2)=1.;
      } else if (isStokes(itsPolFrameIn)) {
          T.diagonal() = 1.;
      } else {
          ASKAPTHROW(AskapError, "Conversion of the selected input polarisation frame into stokes parameters is not supported");
      }
  } else if (isStokes(itsPolFrameIn)) {
      if (itsOutputLinear) {
          // stokes to linear
          T(0,0)=0.5; T(0,1)=0.5;
          T(1,2)=0.5; T(1,3)=casacore::Complex(0.,0.5);
          T(2,2)=0.5; T(2,3)=casacore::Complex(0.,-0.5);
          T(3,0)=0.5; T(3,1)=-0.5;
      } else if (isCircular(itsPolFrameOut)) {
          // stokes to circular
          T(0,0)=0.5; T(0,2)=0.5;
          T(1,1)=casacore::Complex(0.,0.5); T(1,3)=0.5;
          T(2,1)=casacore::Complex(0.,-0.5); T(2,3)=0.5;
          T(3,0)=0.5; T(3,2)=-0.5;
      } else {
          ASKAPTHROW(AskapError, "Conversion of stokes parameters into the selected output polarisation frame is not supported");
      }
  } else if ((itsInputLinear && itsOutputLinear) ||
             (isCircular(itsPolFrameIn) && isCircular(itsPolFrameOut)))  {
             T.diagonal() = 1.;
  } else {
      ASKAPTHROW(AskapError, "Unsupported combination of input and output polarisation frames");
  }

  ASKAPDEBUGASSERT(itsPolFrameIn.nelements()>0);
  ASKAPDEBUGASSERT(itsPolFrameOut.nelements()>0);

  // if we are rotating inputs we require 4 linear pols to be present and
  // rearranged into standard order
  bool rotateInput = itsPARotation.nelements()>0 && itsInputLinear;
  bool rotateOutput = itsPARotation.nelements()>0 && itsOutputLinear;
  // have to copy, because the transformation may not preserve dimensionality
  for (casacore::uInt row = 0; row<itsTransform.nrow(); ++row) {
       const casacore::uInt rowIndex = (rotateOutput ? row : getIndex(itsPolFrameOut[row]));
       ASKAPDEBUGASSERT(rowIndex<4);
       // vector of flags, true if a particular product is present in the input for the given row
       // it is used to check that all required data are present
       std::vector<bool> presentPols(4,false);
       for (casacore::uInt col = 0; col<itsTransform.ncolumn(); ++col) {
            const casacore::uInt colIndex = (rotateInput ? col : getIndex(itsPolFrameIn[col]));
            ASKAPDEBUGASSERT(colIndex<4);
            presentPols[colIndex] = true;
            itsTransform(row,col) = T(rowIndex,colIndex);
       }
       // now check that nothing depends on all products that are absent in the input
       if (itsCheckUnspecifiedProducts) {
           for (casacore::uInt pol = 0; pol<presentPols.size(); ++pol) {
                if (!presentPols[pol]) {
                    ASKAPDEBUGASSERT(pol<T.ncolumn());
                    ASKAPCHECK(abs(T(rowIndex,pol))<1e-5, "Polarisation product "<<
                          casacore::Stokes::name(stokesFromIndex(pol,itsPolFrameIn[0]))<<
                          " is required to get "<< casacore::Stokes::name(itsPolFrameOut[row])<<
                          " polarisation");
                }
           }
       }
  }
  // Now multiply in the rotation matrix if we have one
  if (itsPARotation.nelements()>0) {
      if (itsInputLinear) {
          // form S * R (S=Stokes conversion, R=Rotation)
          // SR(row,col) += S(row,k) * R(k,col)
          casa::Matrix<casa::Complex> SR(itsTransform.shape());
          for (casacore::uInt row = 0; row < itsTransform.nrow(); ++row) {
               for (casacore::uInt col = 0; col < 4; ++col) {
                   for (casacore::uInt k = 0; k < 4; ++k) {
                       SR(row,col) += itsTransform(row,k) * itsPARotation(k,col);
                   }
               }
           }
           itsTransform = SR;
      } else if (itsOutputLinear) {
           // form R * S (S=Stokes conversion, R=Rotation)
           // RS(row,col) += R(row,k) * S(k,col)
           casa::Matrix<casa::Complex> RS(itsTransform.shape());
           for (casacore::uInt row = 0; row < 4; ++row) {
                for (casacore::uInt col = 0; col < itsTransform.ncolumn(); ++col) {
                    for (casacore::uInt k = 0; k < 4; ++k) {
                        RS(row,col) += itsPARotation(row,k) * itsTransform(k,col);
                    }
                }
            }
           itsTransform = RS;
      }
  }
  itsInitialised = true;
}

/// @brief fill matrix describing parallactic angle rotation for linears
/// @details (cos(pa1)  sin(pa1)) * (cos(pa2) sin(pa2)) (outerproduct)
///          (-sin(pa1) cos(pa1))   (-sin(pa2) cos(pa2))
/// @param[in] pa1 parallactic angle on the first antenna
/// @param[in] pa2 parallactic angle on the second antenna
void PolConverter::fillPARotationMatrix(double pa1, double pa2)
{
    // Linear to linear rotation not implemented yet
    ASKAPDEBUGASSERT(!(itsInputLinear && itsOutputLinear));
    if (itsInputLinear || itsOutputLinear ) {
      if (itsPARotation.nelements() == 0) {
          itsPARotation.resize(4,4,0.);
      } else {
          // reuse matrix if possible
          if (pa1 == itsLastpa1 && pa2 == itsLastpa2) return;
      }
      itsLastpa1 = pa1;
      itsLastpa2 = pa2;
      // need to recompute the product of the stokes matrix and the rotation matrix
      itsInitialised = false;

      const double cpa1 = cos(pa1);
      const double cpa2 = cos(pa2);
      double spa1 = sin(pa1);
      double spa2 = sin(pa2);
      // Use inverse matrix for conversion to linears in antenna frame
      // This assumes we're always degridding Stokes to linear or gridding linear to Stokes, not linear to linear
      // The latter would require explicit knowledge of gridding/degridding direction
      if (itsOutputLinear) {
          spa1 = -spa1;
          spa2 = -spa2;
      }
      // sort out the indexing (for non standard orders)
      casacore::Vector<casacore::uInt> rowidx(4),colidx(4);
      for (casacore::uInt col = 0; col<4; ++col) {
          casacore::uInt polIndex =
            (itsInputLinear ? getIndex(itsPolFrameIn[col]) : getIndex(itsPolFrameOut[col]));
          // swap X-Y polarisations if needed
          if (itsSwapPols) polIndex = 3 - polIndex;
          if (itsInputLinear) {
              rowidx(col) = col;
              colidx(polIndex) = col;
          } else {
              // output linear
              rowidx(polIndex) = col;
              colidx(col) = col;
          }
      }
      itsPARotation(rowidx(0),colidx(0)) = cpa1 * cpa2;
      itsPARotation(rowidx(0),colidx(1)) = cpa1 * spa2;
      itsPARotation(rowidx(0),colidx(2)) = spa1 * cpa2;
      itsPARotation(rowidx(0),colidx(3)) = spa1 * spa2;
      itsPARotation(rowidx(1),colidx(0)) = -cpa1 * spa2;
      itsPARotation(rowidx(1),colidx(1)) = cpa1 * cpa2;
      itsPARotation(rowidx(1),colidx(2)) = -spa1 * spa2;
      itsPARotation(rowidx(1),colidx(3)) = spa1 * cpa2;
      itsPARotation(rowidx(2),colidx(0)) = -spa1 * cpa2;
      itsPARotation(rowidx(2),colidx(1)) = -spa1 * spa2;
      itsPARotation(rowidx(2),colidx(2)) = cpa1 * cpa2;
      itsPARotation(rowidx(2),colidx(3)) = cpa1 * spa2;
      itsPARotation(rowidx(3),colidx(0)) = spa1 * spa2;
      itsPARotation(rowidx(3),colidx(1)) = -spa1 * cpa2;
      itsPARotation(rowidx(3),colidx(2)) = -cpa1 * spa2;
      itsPARotation(rowidx(3),colidx(3)) = cpa1 * cpa2;
    } // For circulars this turns into phase and needs to be handled differently
}

/// @brief reverse method for getIndex
/// @details convert index into stokes enum. Because the same index can correspond to a number
/// of polarisation products (meaning of index is frame-dependent), a second parameter is
/// required to unambiguate it. It can be any stokes enum of the frame, not necessarily the
/// first one.
/// @param[in] index an index to convert
/// @param[in] stokes any stokes enum from the working frame
/// @note This method is actually used only to provide a sensible message in the exception. No
/// other code depends on it.
casacore::Stokes::StokesTypes PolConverter::stokesFromIndex(casacore::uInt index, casacore::Stokes::StokesTypes stokes)
{
   const casacore::Vector<casacore::Stokes::StokesTypes> buf(1,stokes);
   casacore::uInt stokesAsuInt = index;
   if (isCircular(buf)) {
       stokesAsuInt += casacore::uInt(casacore::Stokes::RR);
   } else if (isLinear(buf)) {
       stokesAsuInt += casacore::uInt(casacore::Stokes::XX);
   } else if (isStokes(buf)) {
       stokesAsuInt += casacore::uInt(casacore::Stokes::I);
   }
   return casacore::Stokes::StokesTypes(stokesAsuInt);
}

/// @brief test if frame matches a given stokes enum
/// @param[in] polFrame polarisation frame defined as a vector of Stokes enums
/// @param[in] stokes a single stokes enum defining the frame (should be the first in the set)
/// @return true, if the given vector and one stokes enum belong to the same frame
bool PolConverter::sameFrame(const casacore::Vector<casacore::Stokes::StokesTypes> &polFrame,
                        casacore::Stokes::StokesTypes stokes)
{
  ASKAPASSERT(polFrame.nelements()!=0);
  for (casacore::uInt pol = 0; pol<polFrame.nelements(); ++pol) {
       const int index = (int(polFrame[pol]) - int(stokes));
       if ((index<0) || (index>=4)) {
           return false;
       }
  }
  return true;
}

/// @brief return index of a particular polarisation
/// @details To be able to fill matrices efficiently we want to convert, say IQUV into 0,1,2,3.
/// This method does it for all supported types of polarisation products
/// @param[in] stokes a single stokes enum of the polarisation product to convert
/// @return unsigned index
casacore::uInt PolConverter::getIndex(casacore::Stokes::StokesTypes stokes)
{
  casacore::Vector<casacore::Stokes::StokesTypes> buf(1,stokes);
  if (isCircular(buf)) {
     return casacore::uInt(stokes)-casacore::uInt(casacore::Stokes::RR);
  } else if (isLinear(buf)) {
     return casacore::uInt(stokes)-casacore::uInt(casacore::Stokes::XX);
  } else if (isStokes(buf)) {
     return casacore::uInt(stokes)-casacore::uInt(casacore::Stokes::I);
  }
  ASKAPTHROW(AskapError, "Unsupported type of polarisation product in PolConverter::getIndex "<<
             int(stokes));
}

/// @brief check whether stokes parameter correspond to cross-correlation
/// @details casacore allows to code single-dish polarisation and there are some reserved codes
/// as well. As we're doing lots of indexing, it is good to check that given parameter is
/// valid before doing any further work.
/// @note Technically, this and a few other helper methods should be part of casacore::Stokes
/// @param[in] pol polarisation type
/// @return true, if it is a normal cross-correlation or I,Q,U or V.
bool PolConverter::isValid(casacore::Stokes::StokesTypes pol)
{
  // casacore's order is checked by unit test
  if ((int(pol) >= int(casacore::Stokes::I)) || (int(pol) <= int(casacore::Stokes::V))) {
      return true;
  }
  if ((int(pol) >= int(casacore::Stokes::RR)) || (int(pol) <= int(casacore::Stokes::LL))) {
      return true;
  }
  if ((int(pol) >= int(casacore::Stokes::XX)) || (int(pol) <= int(casacore::Stokes::YY))) {
      return true;
  }
  if ((int(pol) >= int(casacore::Stokes::RX)) || (int(pol) <= int(casacore::Stokes::YL))) {
      return true;
  }
  return false;
}

/// @brief convert string representation into a vector of Stokes enums
/// @details It is convenient to define polarisation frames like "xx,xy,yx,yy" or "iquv".
/// This method does it and return a vector of Stokes enums. Comma symbol is ignored. i.e.
/// "iquv" and "i,q,u,v" are equivalent.
/// @param[in] frame a string representation of the frame
/// @return vector with Stokes enums
casacore::Vector<casacore::Stokes::StokesTypes> PolConverter::fromString(const std::string &frame)
{
  if (frame.size() == 0) {
      return casacore::Vector<casacore::Stokes::StokesTypes>();
  }
  // parse the string, it is certainly not empty at this stage
  std::vector<std::string> products;
  products.reserve(4);
  for(size_t pos=0; pos<frame.size(); ++pos) {
     if (frame[pos]!=',' && frame[pos]!=' ') {
         if (frame.find_first_of("iquvIQUV",pos) == pos) {
             products.push_back(frame.substr(pos,1));
             continue;
         }
         ASKAPCHECK(pos+1<frame.size(), "Unable to interpret polarisation product "<<frame[pos]);
         const std::string polProduct = frame.substr(pos,2);
         ASKAPCHECK(polProduct.find_first_not_of("xyrlXYRL") == std::string::npos,
                    "Unknown polarisation product "<<polProduct);
         products.push_back(polProduct);
         ++pos; // two-symbol descriptor has been extracted
     }
  }
  return fromString(products);
}

/// @brief convert string representation into a vector of Stokes enums
/// @details This version of the method accept string representations in a vector and doesn't
/// parse the concatenated string.
/// @param[in] products vector of strings representation of the frame
/// @return vector with Stokes enums
casacore::Vector<casacore::Stokes::StokesTypes> PolConverter::fromString(const std::vector<std::string> &products)
{
  casacore::Vector<casacore::Stokes::StokesTypes> res(products.size());
  for (size_t pol=0;pol<products.size();++pol) {
       res[pol] = casacore::Stokes::type(products[pol]);
  }
  return res;
}

/// @brief convert a vector of Stokes enums into a vector of strings
/// @details This method does a reverse job to fromString. It converts a vector of stokes enums
/// into a vector of strings (with one to one correspondence between elements)
/// @param[in] frame vector of stokes enums
/// @return vector with string represenation
std::vector<std::string> PolConverter::toString(const casacore::Vector<casacore::Stokes::StokesTypes> &frame)
{
  std::vector<std::string> res(frame.nelements());
  for (size_t pol=0; pol<res.size(); ++pol) {
       res[pol] = casacore::Stokes::name(frame[pol]);
  }
  return res;
}

/// @brief helper method to get the elements of transform matrix in the sparse form
/// @details For calibration code when the model is setup with the descrete components there
/// is potentially a heavy calculation for each polarisation product. Therefore, we want to be
/// able to skip it if the coefficient is zero. This method returns non-zero elements of a selected column
/// of the transform matrix.
/// @param[in] pol polarisation type (should be part of the input frame)
/// @return map with the polarisation product as the key and the complex coefficient as the value
std::map<casacore::Stokes::StokesTypes, casacore::Complex> PolConverter::getSparseTransform(const casacore::Stokes::StokesTypes pol) const
{
  if (!itsInitialised) fillMatrix();
  std::map<casacore::Stokes::StokesTypes, casacore::Complex> result;
  casacore::uInt product = 0;
  for (; product < itsPolFrameIn.nelements(); ++product) {
       if (itsPolFrameIn[product] == pol) {
           break;
       }
  }
  ASKAPCHECK(product != itsPolFrameIn.nelements(), "Requested polarisation product "<<casacore::Stokes::name(pol)<<
       " is not found in the input frame the converter is set up with");
  if (itsVoid) {
      result[pol] = casacore::Complex(1.,0.);
  } else {
      ASKAPDEBUGASSERT(product < itsTransform.ncolumn());
      ASKAPDEBUGASSERT(itsPolFrameOut.nelements() <= itsTransform.nrow());
      for (casacore::uInt index = 0; index < itsPolFrameOut.nelements(); ++index) {
           const casacore::Complex val = itsTransform(index,product);
           if (abs(val) > 1e-5) {
               result[itsPolFrameOut[index]] = val;
           }
      }
  }

  ASKAPCHECK(result.size(), "A non-zero output is expected from getSparseTransform, but all selected matrix elements seem to be zeros");
  return result;
}

/// @brief all polarisation products in the stokes frame in the canonical order
/// @details This is a helper method to get 4-element vector filled with polarisation
/// products in the canonical order. This particular version returns them in the Stokes frame.
/// It can be used to setup the converter when either input or output frame (or both) are to be
/// given explicitly
/// @return 4-element vector with polarisation products
casacore::Vector<casacore::Stokes::StokesTypes> PolConverter::canonicStokes()
{
  casacore::Vector<casacore::Stokes::StokesTypes> result(4);
  result[0] = casacore::Stokes::I;
  result[1] = casacore::Stokes::Q;
  result[2] = casacore::Stokes::U;
  result[3] = casacore::Stokes::V;
  return result;
}

/// @brief all polarisation products in the linear frame in the canonical order
/// @details This is a helper method to get 4-element vector filled with polarisation
/// products in the canonical order. This particular version returns them in the linear frame.
/// It can be used to setup the converter when either input or output frame (or both) are to be
/// given explicitly
/// @return 4-element vector with polarisation products
casacore::Vector<casacore::Stokes::StokesTypes> PolConverter::canonicLinear()
{
  casacore::Vector<casacore::Stokes::StokesTypes> result(4);
  result[0] = casacore::Stokes::XX;
  result[1] = casacore::Stokes::XY;
  result[2] = casacore::Stokes::YX;
  result[3] = casacore::Stokes::YY;
  return result;
}

/// @brief all polarisation products in the circular frame in the canonical order
/// @details This is a helper method to get 4-element vector filled with polarisation
/// products in the canonical order. This particular version returns them in the circular frame.
/// It can be used to setup the converter when either input or output frame (or both) are to be
/// given explicitly
/// @return 4-element vector with polarisation products
casacore::Vector<casacore::Stokes::StokesTypes> PolConverter::canonicCircular()
{
  casacore::Vector<casacore::Stokes::StokesTypes> result(4);
  result[0] = casacore::Stokes::RR;
  result[1] = casacore::Stokes::RL;
  result[2] = casacore::Stokes::LR;
  result[3] = casacore::Stokes::LL;
  return result;
}
