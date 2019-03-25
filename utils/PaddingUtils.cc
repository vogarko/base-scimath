/// @file
///
/// PaddingUtils a class containing utilities used for FFT padding in preconditioners. Code like this
/// can probably be moved to a higer level. At this stage we just need to make these methods available not
/// just to the WienerPreconditioner, but for other classes as well.
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

#include <utils/PaddingUtils.h>

#include <askap/AskapError.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LEL/LatticeExpr.h>
#include <casacore/lattices/Lattices/SubLattice.h>

#include <fft/FFTWrapper.h>
#include <casacore/casa/Arrays/ArrayIter.h>


using namespace askap;
using namespace askap::scimath;

/// @brief Inject source into the centre quarter of the target
/// @details 
/// @param[in] target target array to alter, the source will be converted to Complex and stored in the 
/// inner quarter of the target
/// @param[in] source input array
void PaddingUtils::inject(casacore::Lattice<casacore::Complex>& target, casacore::Lattice<float>& source)
{
  target.set(0.0);
  casacore::IPosition corner(target.shape().nelements(),0);
  ASKAPDEBUGASSERT(corner.nelements()>=2);
  ASKAPDEBUGASSERT(target.shape()(0) == source.shape()(0)*2);
  ASKAPDEBUGASSERT(target.shape()(1) == source.shape()(1)*2);
      
  corner(0) = target.shape()(0)/4;
  corner(1) = target.shape()(1)/4;
  casacore::Slicer slicer(corner, source.shape());
  casacore::SubLattice<casacore::Complex> inner(target, slicer, casacore::True);
  inner.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(source)));
}
      
/// @brief Extract target from the center quarter of the source 
/// @details
/// @param[in] target target array to save the reslt, a real part of the inner quarter of the the source array 
/// will be extracted
/// @param[in] source input array
void PaddingUtils::extract(casacore::Lattice<float>& target, casacore::Lattice<casacore::Complex>& source)
{
  target.set(0.0);
  casacore::IPosition corner(source.shape().nelements(),0);
  ASKAPDEBUGASSERT(corner.nelements()>=2);
  ASKAPDEBUGASSERT(source.shape()(0) == target.shape()(0)*2);
  ASKAPDEBUGASSERT(source.shape()(1) == target.shape()(1)*2);
  corner(0) = source.shape()(0)/4;
  corner(1) = source.shape()(1)/4;
  casacore::Slicer slicer(corner, target.shape());
  casacore::SubLattice<casacore::Complex> inner(source, slicer, casacore::True);
  target.copyData(casacore::LatticeExpr<float>(real(inner)));
}


/// @brief helper method to get padded shape
/// @details Most padding applications in the ASKAPsoft require operations on just two
/// axes. This method froms a shape of an array padded on first two axes with the given factor.
/// @param[in] shape shape of the original array
/// @param[in] padding padding factor
/// @return shape of the padded array
casacore::IPosition PaddingUtils::paddedShape(const casacore::IPosition &shape, const float padding)
{
  casacore::IPosition result(shape);
  ASKAPDEBUGASSERT(result.nelements()>=2);
  ASKAPDEBUGASSERT(padding>0);
  result(0) = int(result(0) * padding);
  result(1) = int(result(1) * padding);
  return result;
}

/// @brief padding with fft
/// @details Sometimes it is necessary to do padding in the other domain. This routine 
/// does the Fourier transform, pad the result to the size of the output and then transforms 
/// back to the original domain. It is done if the size of the output array along the first two
/// axes is larger than the size of the input array. If the output array size is smaller, just 
/// the inner subimage is copied and no fft is done. Equal size result in no operation.
/// @note Both input and output arrays should be at least 2-dimensional, otherwise an exception
/// is thrown (in debug mode)
/// @param[in] in input array 
/// @param[in] out output array (should already be resized to a desired size) 
void PaddingUtils::fftPad(const casacore::Array<double>& in, casacore::Array<double>& out)
{
  ASKAPDEBUGASSERT(in.shape().nelements()>=2);    
  const int inx=in.shape()(0);
  const int iny=in.shape()(1);

  ASKAPDEBUGASSERT(out.shape().nelements()>=2);    
      
  const int onx=out.shape()(0);
  const int ony=out.shape()(1);
            
  // Shortcut no-op
  if ((inx == onx) && (iny == ony)) {
       out = in.copy();
       return;
  }
      
  ASKAPCHECK((onx >= inx) == (ony >= iny), 
             "Attempting to pad to a rectangular array smaller on one axis");
  if (onx<inx) {
      // no fft padding required, the output array is smaller.
      casacore::Array<double> tempIn(in); // in is a conceptual const array here
      out = centeredSubArray(tempIn,out.shape()).copy();
      return;
  }
      
      
  /// Make an iterator that returns plane by plane
  casacore::ReadOnlyArrayIterator<double> inIt(in, 2);
  casacore::ArrayIterator<double> outIt(out, 2);
  while (!inIt.pastEnd()&&!outIt.pastEnd()) {
         casacore::Matrix<casacore::DComplex> inPlane(inx, iny);
         casacore::Matrix<casacore::DComplex> outPlane(onx, ony);
         casacore::convertArray(inPlane, inIt.array());
         outPlane.set(0.0);
         fft2d(inPlane, false);
         for (int iy=0; iy<iny; ++iy) {
              for (int ix=0; ix<inx; ++ix) {
                   outPlane(ix+(onx-inx)/2, iy+(ony-iny)/2) = inPlane(ix, iy);
              }
         }
         
         fft2d(outPlane, true);
         const casacore::Array<casacore::DComplex> constOutPlane(outPlane);
         casacore::Array<double> outArray(outIt.array());
	
         casacore::real(outArray, constOutPlane);
	
         inIt.next();
         outIt.next();
  }
}

/// @brief padding with fft
/// @details This variant of the method is intended for the case where internal padding is used.
/// The input array is padded with fft to a larger area (size is controlled by the given factor)
/// and then an inner sub-array is extracted. With factor==1, there is no difference from the
/// other fftPad methid
/// @note Both input and output arrays should be at least 2-dimensional, otherwise an exception
/// is thrown (in debug mode)
/// @param[in] in input array 
/// @param[in] out output array (should already be resized to a desired size) 
/// @param[in] factor additional padding of the output array
void PaddingUtils::fftPad(const casacore::Array<double>& in, casacore::Array<double>& out, float factor)
{
  ASKAPDEBUGASSERT(factor>0);
  const casacore::IPosition shape = paddedShape(out.shape(),factor);
  if (shape.isEqual(out.shape())) {
      // factor == 1 case, by comparing shapes we can avoid problems with rounding off errors
      fftPad(in,out);
  } else {
      casacore::Array<double> tempOut(shape);
      fftPad(in,tempOut);
      out = centeredSubArray(tempOut,out.shape()).copy();
  }
}

/// @brief helper method to get shape before padding
/// @details Most padding applications in the ASKAPsoft require operations on just two
/// axes. This method froms a shape of an array before padding from the padded shape
/// @param[in] shape shape of the padded array
/// @param[in] padding padding factor (should be a positive number)
/// @return shape before padding
casacore::IPosition PaddingUtils::unpadShape(const casacore::IPosition &shape, const float padding)
{
   ASKAPDEBUGASSERT(shape.nelements()>=2);
   ASKAPDEBUGASSERT(padding>0);
   casacore::IPosition outShape(shape);
   // form desired shape
   for (size_t dim=0; dim<2; ++dim) {
        outShape(dim) = int(outShape(dim) / padding);
        // rounding off operation does not commute with division/multiplication, hence an extra check is required
        if (int(padding*outShape(dim))<shape(dim)) {
            ++outShape(dim);
        }
   }
   // ensure that this method does reverse operation to paddedShape
   ASKAPDEBUGASSERT(paddedShape(outShape,padding) == shape);
   return outShape;
}


