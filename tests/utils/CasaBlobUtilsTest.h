/// @file CasaBlobUtilsTest.cc
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

// CPPUnit includes
#include <cppunit/extensions/HelperMacros.h>

// Support classes
#include <askap/askap/AskapError.h>
#include "casacore/casa/Arrays/Vector.h"
#include "casacore/casa/Arrays/Cube.h"
#include "Blob/BlobIStream.h"
#include "Blob/BlobIBufVector.h"
#include "Blob/BlobOStream.h"

// for coordinate system
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include "Blob/BlobOBufVector.h"

// Classes to test
#include "askap/scimath/utils/CasaBlobUtils.h"

using namespace casa;
using namespace LOFAR;

namespace askap {

namespace scimath {

class CasaBlobUtilsTest : public CppUnit::TestFixture {
        CPPUNIT_TEST_SUITE(CasaBlobUtilsTest);
        CPPUNIT_TEST(testCoordinateSystem);
        CPPUNIT_TEST(testQuantity);
        CPPUNIT_TEST(testMDirectionRef);
        CPPUNIT_TEST(testMVDirection);
        CPPUNIT_TEST(testMDirection);
        CPPUNIT_TEST(testStokesTypes);
        CPPUNIT_TEST_SUITE_END();

    public:
        void setUp()
        {
        };

        void tearDown()
        {
        }

        void testCoordinateSystem()
        {

          casa::CoordinateSystem source;
          casa::CoordinateSystem target;

          // DirectionCoordinate
          {
            casa::Vector<casa::Double> refPix(2,512.);
            casa::Vector<casa::Double> refVal(2,0.0);
            casa::Vector<casa::Double> increment(2,0.1);
            casa::Matrix<casa::Double> xform(2,2,0.0);
            xform.row(0)(0) = 1.0;
            xform.row(1)(1) = 1.0;
            casa::DirectionCoordinate dc(casa::MDirection::B1950,
                casa::Projection(casa::Projection::TAN), refVal[0],refVal[1],
                increment[0],increment[1], xform, refPix[0],refPix[1]);
            source.addCoordinate(dc);
          }

          // SpectralCoordinate
          {
            casa::Vector<casa::Double> refPix(1,512.);
            casa::Vector<casa::Double> refVal(1,0.0);
            casa::Vector<casa::Double> increment(1,0.1);
            casa::Double restFreq = 0.0;
            casa::SpectralCoordinate sc(casa::MFrequency::GALACTO,
                refVal[0], increment[0], refPix[0], restFreq);
            source.addCoordinate(sc);
          }

          // StokesCoordinate
          {
            casa::Vector<casa::Int> whichStokes(5);
            whichStokes(0) = casa::Stokes::I;
            whichStokes(1) = casa::Stokes::XX;
            whichStokes(2) = casa::Stokes::RR;
            whichStokes(3) = casa::Stokes::RX;
            whichStokes(4) = casa::Stokes::XR;
            casa::StokesCoordinate sc(whichStokes);
            source.addCoordinate(sc);
          }

          // Encode
          std::vector<int8_t> buf;
          LOFAR::BlobOBufVector<int8_t> obv(buf);
          LOFAR::BlobOStream out(obv);
          out.putStart("CoordinateSystem", 1);
          out << source;
          out.putEnd();

          // Decode
          LOFAR::BlobIBufVector<int8_t> ibv(buf);
          LOFAR::BlobIStream in(ibv);
          int version = in.getStart("CoordinateSystem");
          ASKAPASSERT(version == 1);
          in >> target;
          in.getEnd();

          // Tolerance for double equality
          const double tol = 1.0E-8;

          // CoordinateSystem tests
            CPPUNIT_ASSERT_EQUAL(source.nCoordinates(),target.nCoordinates());

          // DirectionCoordinate tests
          {
            int dcPos;
            dcPos = source.findCoordinate(casa::Coordinate::DIRECTION,-1);
            const casa::DirectionCoordinate& sourceDC =
                source.directionCoordinate(dcPos);
            dcPos = target.findCoordinate(casa::Coordinate::DIRECTION,-1);
            const casa::DirectionCoordinate& targetDC =
                target.directionCoordinate(dcPos);
            CPPUNIT_ASSERT_EQUAL(sourceDC.directionType(),
                                 targetDC.directionType());
            CPPUNIT_ASSERT_EQUAL(sourceDC.projection().name(),
                                 targetDC.projection().name());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceDC.referencePixel()(0),
                                         targetDC.referencePixel()(0), tol);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceDC.referencePixel()(1),
                                         targetDC.referencePixel()(1), tol);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceDC.increment()(0),
                                         targetDC.increment()(0), tol);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceDC.increment()(1),
                                         targetDC.increment()(1), tol);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceDC.referenceValue()(0),
                                         targetDC.referenceValue()(0), tol);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceDC.referenceValue()(1),
                                         targetDC.referenceValue()(1), tol);
          }

          // SpectralCoordinate tests
          {
            int scPos;
            scPos = source.findCoordinate(casa::Coordinate::SPECTRAL,-1);
            const casa::SpectralCoordinate& sourceSC =
                source.spectralCoordinate(scPos);
            scPos = target.findCoordinate(casa::Coordinate::SPECTRAL,-1);
            const casa::SpectralCoordinate& targetSC =
                target.spectralCoordinate(scPos);
            CPPUNIT_ASSERT_EQUAL(sourceSC.frequencySystem(),
                                 targetSC.frequencySystem());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceSC.referencePixel()(0),
                                         targetSC.referencePixel()(0), tol);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceSC.increment()(0),
                                         targetSC.increment()(0), tol);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sourceSC.referenceValue()(0),
                                         targetSC.referenceValue()(0), tol);
          }
          // StokesCoordinate
          {
            int scPos;
            scPos = source.findCoordinate(casa::Coordinate::STOKES,-1);
            const casa::StokesCoordinate& sourceSC =
                source.stokesCoordinate(scPos);
            scPos = target.findCoordinate(casa::Coordinate::STOKES,-1);
            const casa::StokesCoordinate& targetSC =
                target.stokesCoordinate(scPos);
            CPPUNIT_ASSERT_EQUAL(sourceSC.stokes().shape(),
                                 targetSC.stokes().shape());
            for (casa::Int k = 0; k<sourceSC.stokes().shape(); ++k) {
              CPPUNIT_ASSERT_EQUAL(sourceSC.stokes()(k),
                                   targetSC.stokes()(k));
            }
          }

        }

        // MV: serialisation for the measures-related casacore types
        
        // helper method for serialisation and deserialisation of any type
        // it simply takes the first variable, serialises it, then deserialises and
        // writes the result into the second variable.
        template<typename T> 
        static void copyViaBlob(const T &source, T &target) {
          // Encode
          std::vector<int8_t> buf;
          const int datagramVersion = 10;
          {
             LOFAR::BlobOBufVector<int8_t> obv(buf);
             LOFAR::BlobOStream out(obv);
             out.putStart("copyViaBlobTest", datagramVersion);
             out << source;
             out.putEnd();
          }
          // Decode
          LOFAR::BlobIBufVector<int8_t> ibv(buf);
          LOFAR::BlobIStream in(ibv);
          const int version = in.getStart("copyViaBlobTest");
          CPPUNIT_ASSERT_EQUAL(datagramVersion,version);
          in >> target;
          in.getEnd();
        }

        void testQuantity() {
           const casa::Quantity q(3.1415, "km/s");
           casa::Quantity receivedQ(0.1, "MHz");
           copyViaBlob(q, receivedQ);
           CPPUNIT_ASSERT_DOUBLES_EQUAL(q.getValue(), receivedQ.getValue(),1e-6);
           CPPUNIT_ASSERT_EQUAL(q.getFullUnit().getName(), receivedQ.getFullUnit().getName());
        }

        void testMDirectionRef() {
          const casa::MDirection::Ref ref(casa::MDirection::AZEL);
          casa::MDirection::Ref receivedRef(casa::MDirection::J2000);
          copyViaBlob(ref,receivedRef);
          CPPUNIT_ASSERT_EQUAL(ref.getType(), receivedRef.getType());
        }

        void testMVDirection() {
          const casa::MVDirection dir(casa::Quantity(135.0, "deg"), casa::Quantity(-31.0, "deg"));
          casa::MVDirection receivedDir(casa::Quantity(1.0, "rad"), casa::Quantity(0, "arcsec"));
          copyViaBlob(dir, receivedDir);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(dir.getLong(), receivedDir.getLong(), 1e-6);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(dir.getLat(), receivedDir.getLat(), 1e-6);
        }

        void testMDirection() {
          const casa::MDirection dir(casa::MVDirection(casa::Quantity(135.0, "deg"), casa::Quantity(+31.0, "deg")), 
                                     casa::MDirection::AZEL);
          casa::MDirection receivedDir(casa::MVDirection(casa::Quantity(1.0, "rad"), casa::Quantity(0, "arcsec")), 
                                     casa::MDirection::J2000);
          copyViaBlob(dir, receivedDir);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(dir.getValue().getLong(), receivedDir.getValue().getLong(), 1e-6);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(dir.getValue().getLat(), receivedDir.getValue().getLat(), 1e-6);
          CPPUNIT_ASSERT_EQUAL(dir.getRef().getType(), receivedDir.getRef().getType());
        }

        void testStokesTypes() {
          const casa::Stokes::StokesTypes pol = casa::Stokes::XY;
          casa::Stokes::StokesTypes receivedPol(casa::Stokes::RR);
          copyViaBlob(pol, receivedPol);
          CPPUNIT_ASSERT_EQUAL(pol, receivedPol);
        }
};

}   // End namespace scimath
}   // End namespace askap
