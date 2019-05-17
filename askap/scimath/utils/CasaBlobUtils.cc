/// @file CasaBlobUtils.cc
///
/// @copyright (c) 2016 CSIRO
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
/// @author Daniel Mitchell <daniel.mitchell@csiro.au>
/// measures stuff - Max Voronkov
///

#include "askap/scimath/utils/CasaBlobUtils.h"


// casacore includes
#include "casacore/casa/aips.h"

// for coordinate system
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>

#include <Blob/BlobSTL.h>

/// @brief increment this if there is any change to the stuff written into blob
const int coordSysBlobVersion = 2;

using namespace askap;

namespace LOFAR
{

    /// @brief blob support for casacore::CoordinateSystem
    LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os,
                                   const casacore::CoordinateSystem& cSys)
    {
        os.putStart("CoordinateSystem",coordSysBlobVersion);

        int nCoordinates = cSys.nCoordinates();
        int coordinateCount = 0;
        int dcPos = cSys.findCoordinate(casacore::Coordinate::DIRECTION,-1);
        int fcPos = cSys.findCoordinate(casacore::Coordinate::SPECTRAL,-1);
        int pcPos = cSys.findCoordinate(casacore::Coordinate::STOKES,-1);

        os << nCoordinates;
        os << dcPos << fcPos << pcPos;
        if (dcPos >= 0) {
          const casacore::DirectionCoordinate& dc = cSys.directionCoordinate(dcPos);
          os << casacore::MDirection::showType(dc.directionType()) <<
                dc.projection().name() << dc.referenceValue() <<
                dc.increment() << dc.linearTransform() <<
                dc.referencePixel() << dc.worldAxisUnits();
          coordinateCount++;
        }
        if (fcPos >= 0) {
          const casacore::SpectralCoordinate& fc = cSys.spectralCoordinate(fcPos);
          os << casacore::MFrequency::showType(fc.frequencySystem()) <<
                fc.referenceValue() << fc.increment() << fc.referencePixel() <<
                fc.restFrequency() << fc.worldAxisUnits();
          coordinateCount++;
        }
        if (pcPos >= 0) {
          const casacore::StokesCoordinate& pc = cSys.stokesCoordinate(pcPos);
          os << pc.stokes();
          coordinateCount++;
        }

        ASKAPCHECK(coordinateCount == nCoordinates,
            "BlobOStream currently only supports DIRECTION, SPECTRAL and STOKES coordinates");

        os.putEnd();
        return os;
    }

    /// @brief blob support for casacore::CoordinateSystem
    LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is,
                                   casacore::CoordinateSystem& cSys)
    {
        const int version = is.getStart("CoordinateSystem");
        ASKAPCHECK(version == coordSysBlobVersion, "Attempting to read from" <<
            " a blob stream a CoordinateSystem object of the wrong version," <<
            " expect "<< coordSysBlobVersion<<" got "<<version);

        int nCoordinates;
        is >> nCoordinates;
        int dcPos, fcPos, pcPos;
        is >> dcPos >> fcPos >> pcPos;

        cSys = casacore::CoordinateSystem();

        if (dcPos >= 0) {
          casacore::String dirTypeStr;
          is >> dirTypeStr;
          casacore::MDirection::Types dirType;
          casacore::MDirection::getType(dirType, dirTypeStr);
          casacore::String projectionName;
          is >> projectionName;
          casacore::Vector<casacore::Double> increment, refPix, refVal;
          casacore::Matrix<casacore::Double> xform;
          is >> refVal >> increment >> xform >> refPix;
          ASKAPCHECK(increment.nelements() == 2,
              "Direction axis increment should be a vector of size 2");
          ASKAPCHECK(refPix.nelements() == 2,
              "Direction axis reference pixel should be a vector of size 2");
          ASKAPCHECK(refVal.nelements() == 2,
              "Direction axis reference value should be a vector of size 2");
          ASKAPCHECK(xform.shape() == casacore::IPosition(2,2,2),
              "Direction axis transform matrix should be 2x2");
          casacore::DirectionCoordinate dc(dirType,
              casacore::Projection::type(projectionName), refVal[0],refVal[1],
              increment[0],increment[1], xform, refPix[0],refPix[1]);
          casacore::Vector<casacore::String> worldAxisUnits;
          is >> worldAxisUnits;
          dc.setWorldAxisUnits(worldAxisUnits);
          cSys.addCoordinate(dc);
        }
        if (fcPos >= 0) {
          casacore::String freqTypeStr;
          is >> freqTypeStr;
          casacore::MFrequency::Types freqType;
          casacore::MFrequency::getType(freqType, freqTypeStr);
          casacore::Double restFreq;
          casacore::Vector<casacore::Double> increment, refPix, refVal;
          is >> refVal >> increment >> refPix >> restFreq;
          ASKAPCHECK(increment.nelements() == 1,
              "Spectral axis increment should be a vector of size 1");
          ASKAPCHECK(refPix.nelements() == 1,
              "Spectral axis reference pixel should be a vector of size 1");
          ASKAPCHECK(refVal.nelements() == 1,
              "Spectral axis reference value should be a vector of size 1");
          casacore::SpectralCoordinate fc(freqType, refVal[0], increment[0],
              refPix[0], restFreq);
          casacore::Vector<casacore::String> worldAxisUnits;
          is >> worldAxisUnits;
          fc.setWorldAxisUnits(worldAxisUnits);
          cSys.addCoordinate(fc);
        }
        if (pcPos >= 0) {
          casacore::Vector<casacore::Int> whichStokes;
          is >> whichStokes;
          casacore::StokesCoordinate pc(whichStokes);
          cSys.addCoordinate(pc);
        }

        is.getEnd();
        return is;
    }

    // MV: shift operators for measures-related casacore types
 
    /// @brief output operator for casacore::Quantity
    /// @param[in] os output stream
    /// @param[in] q Quantity to serialise
    /// @return output stream for chaining
    LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os, const casacore::Quantity& q)
    {
      os<<q.getFullUnit().getName()<<q.getValue();
      return os;
    }

    /// @brief input operator for casacore::Quantity
    /// @param[in] is input stream
    /// @param[in] q quantity object to populate
    /// @return input stream for chaining
    LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is, casacore::Quantity& q)
    {
      std::string unitName;
      casacore::Double val;
      is>>unitName>>val;
      q=casacore::Quantity(val, casacore::Unit(unitName));
      return is;
    }

    /// @brief output operator for casacore::MDirection::Ref
    /// @param[in] os output stream
    /// @param[in] ref object to serialise
    /// @return output stream for chaining
    LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os, const casacore::MDirection::Ref& ref)
    {
      const casacore::uInt refType = ref.getType();
      // for now ignore frame and offset - we're not using them in ingest anyway
      // but do check that the user didn't set them. If someone sees the exception later on, they can add 
      // the required logic in this and the following method
      ASKAPCHECK(const_cast<casacore::MDirection::Ref&>(ref).getFrame().empty(), 
                 "Serialisation of frame information attached to measures is not implemented");
      ASKAPCHECK(ref.offset() == NULL, "Serialisation of frame offset in measures is not implemented");
      os << refType;
      return os;
    }

    /// @brief input operator for casacore::MDirection::Ref
    /// @param[in] is input stream
    /// @param[in] ref object to populate
    /// @return input stream for chaining
    LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is, casacore::MDirection::Ref& ref)
    {
      casacore::uInt refType;
      is >> refType;
      // for now ignore frame and offset - we're not using them in ingest anyway
      // see output operator for cross checks.
      ref = casacore::MDirection::Ref(refType);
      return is;
    }

    /// @brief output operator for casacore::MVDirection
    /// @param[in] os output stream
    /// @param[in] dir object to serialise
    /// @return output stream for chaining
    LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os, const casacore::MVDirection& dir)
    {
      casacore::Vector<casacore::Double> angles = dir.get();
      ASKAPDEBUGASSERT(angles.nelements() == 2);
      os << angles;
      return os;
    }

    /// @brief input operator for casacore::MVDirection
    /// @param[in] is input stream
    /// @param[in] dir object to populate
    /// @return input stream for chaining
    LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is, casacore::MVDirection& dir)
    {
      casacore::Vector<casacore::Double> angles;
      is >> angles;
      ASKAPDEBUGASSERT(angles.nelements() == 2);
      dir = casacore::MVDirection(angles);
      return is;
    }

    /// @brief output operator for casacore::MDirection
    /// @param[in] os output stream
    /// @param[in] dir object to serialise
    /// @return output stream for chaining
    LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os, const casacore::MDirection& dir)
    {
      os<<dir.getValue()<<dir.getRef();
      return os;
    }

    /// @brief input operator for casacore::MDirection
    /// @param[in] is input stream
    /// @param[in] dir object to populate
    /// @return input stream for chaining
    LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is, casacore::MDirection& dir)
    {
      casacore::MVDirection val;
      casacore::MDirection::Ref ref;
      is >> val >> ref;
      dir = casacore::MDirection(val,ref);
      return is;
    }

    /// @brief output operator for casacore::Stokes::StokesTypes
    /// @param[in] os output stream
    /// @param[in] pol object to serialise
    /// @return output stream for chaining
    LOFAR::BlobOStream& operator<<(LOFAR::BlobOStream& os, const casacore::Stokes::StokesTypes& pol)
    {
      os<<(int)pol;
      return os;
    }

    /// @brief input operator for casacore::Stokes::StokesTypes
    /// @param[in] is input stream
    /// @param[in] pol object to populate
    /// @return input stream for chaining
    LOFAR::BlobIStream& operator>>(LOFAR::BlobIStream& is, casacore::Stokes::StokesTypes& pol)
    {
      int intPol;
      is >> intPol;
      pol = casacore::Stokes::StokesTypes(intPol);
      return is;
    }

} // End namespace LOFAR

