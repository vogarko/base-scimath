/// @file
/// @brief Normal equations with an approximation for imaging
/// @details There are two kinds of normal equations currently supported. The
/// first one is a generic case, where the full normal matrix is retained. It
/// is used for calibration. The second one is intended for imaging, where we
/// can't afford to keep the whole normal matrix. In the latter approach, the
/// matrix is approximated by a sum of diagonal and shift invariant matrices.
/// This class represents the approximated case, and is used with imaging
/// algorithms.
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

#include <askap/scimath/fitting/DesignMatrix.h>
#include <askap/scimath/fitting/ImagingNormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/Vector.h>

#include <casacore/coordinates/Coordinates/CoordinateUtil.h>

#include <askap/scimath/utils/CasaBlobUtils.h>
#include <Blob/BlobArray.h>
#include <Blob/BlobSTL.h>

#include <askap/askap/AskapError.h>

#include <askap/scimath/utils/DeepCopyUtils.h>

#include <linmos/LinmosAccumulator.h>

#include <stdexcept>
#include <string>
#include <map>
#include <cmath>
#include <vector>

using namespace LOFAR;

using std::abs;
using std::map;
using std::string;
using std::vector;

ASKAP_LOGGER(nelogger, ".nelogger");
namespace askap
{
  namespace scimath
  {


    ImagingNormalEquations::ImagingNormalEquations() : itsWeightType(FROM_BP_MODEL), itsWeightState(INHERENT) {};

    ImagingNormalEquations::ImagingNormalEquations(const Params& ip) : itsWeightState(INHERENT),itsWeightType(FROM_BP_MODEL)
    {
      std::vector<std::string> names=ip.freeNames();
      std::vector<std::string>::iterator iterRow;
      std::vector<std::string>::iterator iterCol;
      for (iterRow=names.begin();iterRow!=names.end();++iterRow)
      {
        itsDataVector[*iterRow]=casacore::Vector<double>(0);
        itsShape[*iterRow]=casacore::IPosition();
        itsReference[*iterRow]=casacore::IPosition();
        itsCoordSys[*iterRow]=casacore::CoordinateSystem();
        itsNormalMatrixSlice[*iterRow]=casacore::Vector<double>(0);
        itsNormalMatrixDiagonal[*iterRow]=casacore::Vector<double>(0);
        itsPreconditionerSlice[*iterRow]=casacore::Vector<double>(0);
      }
    }


    /// @brief copy constructor
    /// @details Data members of this class are non-trivial types including
    /// std containers of casa containers. The letter are copied by reference by default. We,
    /// therefore, need this copy constructor to achieve proper copying.
    /// @param[in] src input measurement equations to copy from
    ImagingNormalEquations::ImagingNormalEquations(const ImagingNormalEquations &src) :
         INormalEquations(src),itsShape(src.itsShape), itsReference(src.itsReference),
         itsCoordSys(src.itsCoordSys), itsWeightType(src.itsWeightType), itsWeightState(src.itsWeightState)
    {
      deepCopyOfSTDMap(src.itsNormalMatrixSlice, itsNormalMatrixSlice);
      deepCopyOfSTDMap(src.itsNormalMatrixDiagonal, itsNormalMatrixDiagonal);
      deepCopyOfSTDMap(src.itsPreconditionerSlice, itsPreconditionerSlice);
      deepCopyOfSTDMap(src.itsDataVector, itsDataVector);
    }

    /// @brief assignment operator
    /// @details Data members of this class are non-trivial types including
    /// std containers of casa containers. The letter are copied by reference by default. We,
    /// therefore, need this copy constructor to achieve proper copying.
    /// @param[in] src input measurement equations to copy from
    /// @return reference to this object
    ImagingNormalEquations& ImagingNormalEquations::operator=(const ImagingNormalEquations &src)
    {
      if (&src != this) {
        itsShape = src.itsShape;
        itsReference = src.itsReference;
        itsCoordSys = src.itsCoordSys;
        itsWeightType = src.itsWeightType;
        itsWeightState = src.itsWeightState;

        deepCopyOfSTDMap(src.itsNormalMatrixSlice, itsNormalMatrixSlice);
        deepCopyOfSTDMap(src.itsNormalMatrixDiagonal, itsNormalMatrixDiagonal);
        deepCopyOfSTDMap(src.itsPreconditionerSlice, itsPreconditionerSlice);
        deepCopyOfSTDMap(src.itsDataVector, itsDataVector);
      }
      return *this;
    }
    /// @brief zero all the data elements but keep thier size and coordinate systems
    /// @details I have found I need a utility function that allows me to zero the contents
    /// of the ImagingNormalEquations but to keep the meta information for the params
    /// @param[in] ip the params
    /// @return void

    void ImagingNormalEquations::zero(const Params& ip) {


        std::vector<std::string> names=ip.freeNames();
        std::vector<std::string>::iterator iterRow;
        std::vector<std::string>::iterator iterCol;
        for (iterRow=names.begin();iterRow!=names.end();++iterRow)
        {
          std::fill(itsDataVector[*iterRow].begin(),itsDataVector[*iterRow].end(),0); //data
          // std::fill(itsNormalMatrixSlice[*iterRow].begin(),itsNormalMatrixSlice[*iterRow].end(),0); //PSF?
          std::fill(itsNormalMatrixDiagonal[*iterRow].begin(),itsNormalMatrixDiagonal[*iterRow].end(),0); // weights
          std::fill(itsPreconditionerSlice[*iterRow].begin(),itsPreconditionerSlice[*iterRow].end(),0); // preconditioner
        }
    }

    ImagingNormalEquations::~ImagingNormalEquations()
    {
      reset();
    }

    /// @brief Merge these normal equations with another
    /// @details Combining two normal equations depends on the actual class type
    /// (different work is required for a full matrix and for an approximation).
    /// This method must be overriden in the derived classes for correct
    /// implementation.
    /// This means that we just add
    /// @param[in] src an object to get the normal equations from
    void ImagingNormalEquations::merge(const INormalEquations& src)
    {
      ASKAPTRACE("ImagingNormalEquations::merge");
      // std::cout << "In merge" << std::endl;

      try {
        const ImagingNormalEquations &other =
                           dynamic_cast<const ImagingNormalEquations&>(src);

        // list of parameters covered by input normal equations
        const std::vector<std::string> otherParams = other.unknowns();
        if (!otherParams.size()) {
          // do nothing, src is empty
          return;
        }

        // initialise an image accumulator
        imagemath::LinmosAccumulator<double> accumulator;

        std::vector<std::string> names = unknowns();

        if (!names.size()) {
          // this object is empty, just do an assignment
          std::cout << "In merge - this object is empty just assign" << std::endl;
          operator=(other);
          return;
        }

        // concatenate unique parameter names
        std::vector<std::string>::const_iterator iterCol = otherParams.begin();
        for (; iterCol != otherParams.end(); ++iterCol) {
          if (std::find(names.begin(),names.end(),*iterCol) == names.end()) {
            names.push_back(*iterCol);
          }
        }

        // step through parameter names and add/merge new ones
        for (iterCol = names.begin(); iterCol != names.end(); ++iterCol)
        {

          // check dataVector
          if (other.itsDataVector.find(*iterCol) == other.itsDataVector.end()) {
            // no new data for this parameter
            continue;
          }

          // record how data are updated
          enum updateType_t{ overwrite, add, linmos };
          int updateType;

          // check coordinate systems and record how data are updated

          ASKAPDEBUGASSERT(other.itsDataVector.find(*iterCol) != other.itsDataVector.end());

          const casacore::Vector<double> &newDataVec = other.itsDataVector.find(*iterCol)->second;
          const casacore::CoordinateSystem &newCoordSys = other.itsCoordSys.find(*iterCol)->second;
          const casacore::IPosition &newShape = other.itsShape.find(*iterCol)->second;

          if(itsDataVector[*iterCol].size()==0)
          {
            // no old data for this parameter
            itsDataVector[*iterCol].assign(newDataVec);
            updateType = overwrite;
            std::cout << "In merge - no old data to merge" << std::endl;
          }
          else if(accumulator.coordinatesAreEqual(itsCoordSys[*iterCol],newCoordSys,
                                                  itsShape[*iterCol],newShape))
          {
            // new and old data can be added directly for this parameter
            itsDataVector[*iterCol] += newDataVec;
            updateType = add;
            ASKAPLOG_INFO_STR(nelogger,"In merge of " << *iterCol << "- old data to merge but coordinates match");


          }
          else
          {
            // new and old data cannot be added directly for this parameter
            if (itsCoordSys[*iterCol].nCoordinates() == 0) {
              // no coordinate information, so just use the new data
              itsDataVector[*iterCol].assign(newDataVec);
              updateType = overwrite;
              ASKAPLOG_INFO_STR(nelogger,"In merge of " << *iterCol <<  " old data to merge - no old coordinate info");
            } else if (itsCoordSys[*iterCol].nCoordinates() != newCoordSys.nCoordinates()) {
              // different dimensions, so just use the new data
              itsDataVector[*iterCol].assign(newDataVec);
              updateType = overwrite;
              ASKAPLOG_INFO_STR(nelogger,"In merge of " << *iterCol <<  " - old data to merge - dimension mismatch");
            } else {
              // regrid then add (using weights)
              ASKAPLOG_INFO_STR(nelogger,"In merge of " << *iterCol <<  " old data to merge coordinates present - linmos merge beginning");
              linmosMerge(other, *iterCol);
              updateType = linmos;
            }
          }

          // update shape and reference.
          if (updateType == overwrite)
          {
            itsShape[*iterCol].resize(0);
            ASKAPDEBUGASSERT(other.itsShape.find(*iterCol) != other.itsShape.end());
            itsShape[*iterCol] = other.itsShape.find(*iterCol)->second;

            itsReference[*iterCol].resize(0);
            ASKAPDEBUGASSERT(other.itsReference.find(*iterCol) != other.itsReference.end());
            itsReference[*iterCol] = other.itsReference.find(*iterCol)->second;

            ASKAPDEBUGASSERT(other.itsCoordSys.find(*iterCol) != other.itsCoordSys.end());
            // leave "keep" and "replace" axis vectors empty, so they are all removed.
            Vector<Int> worldAxes;
            Vector<Double> worldRep;
            CoordinateUtil::removeAxes(itsCoordSys[*iterCol], worldRep, worldAxes, False);
            itsCoordSys[*iterCol] = other.itsCoordSys.find(*iterCol)->second;
          }

          // linmos uses itsNormalMatrixDiagonal to store weights. Otherwise, leave this as it was.
          if (updateType != linmos)
          {
            // check NormalMatrixSlice
            ASKAPDEBUGASSERT(other.itsNormalMatrixSlice.find(*iterCol) != other.itsNormalMatrixSlice.end());
            if(itsNormalMatrixSlice[*iterCol].shape()==0)
            {
              itsNormalMatrixSlice[*iterCol].assign(other.itsNormalMatrixSlice.find(*iterCol)->second);
            }
            else if(itsNormalMatrixSlice[*iterCol].shape() !=
                    other.itsNormalMatrixSlice.find(*iterCol)->second.shape())
            {
              itsNormalMatrixSlice[*iterCol].assign(other.itsNormalMatrixSlice.find(*iterCol)->second);
            }
            else
            {
              itsNormalMatrixSlice[*iterCol] += other.itsNormalMatrixSlice.find(*iterCol)->second;
            }
            // check NormalMatrixDiagonal
            ASKAPDEBUGASSERT(other.itsNormalMatrixDiagonal.find(*iterCol) !=
                             other.itsNormalMatrixDiagonal.end());
            if(itsNormalMatrixDiagonal[*iterCol].shape()==0)
            {
                itsNormalMatrixDiagonal[*iterCol].assign(other.itsNormalMatrixDiagonal.find(*iterCol)->second);
            }
            else if(itsNormalMatrixDiagonal[*iterCol].shape() !=
                    other.itsNormalMatrixDiagonal.find(*iterCol)->second.shape())
            {
              itsNormalMatrixDiagonal[*iterCol].assign(other.itsNormalMatrixDiagonal.find(*iterCol)->second);
            }
            else
            {
              itsNormalMatrixDiagonal[*iterCol] += other.itsNormalMatrixDiagonal.find(*iterCol)->second;
            }
            // check PreconditionerSlice
            ASKAPDEBUGASSERT(other.itsPreconditionerSlice.find(*iterCol) != other.itsPreconditionerSlice.end());
            if(itsPreconditionerSlice[*iterCol].shape()==0)
            {
              itsPreconditionerSlice[*iterCol].assign(other.itsPreconditionerSlice.find(*iterCol)->second);
            }
            else if(itsPreconditionerSlice[*iterCol].shape() !=
                    other.itsPreconditionerSlice.find(*iterCol)->second.shape())
            {
              itsPreconditionerSlice[*iterCol].assign(other.itsPreconditionerSlice.find(*iterCol)->second);
            }
            else
            {
              itsPreconditionerSlice[*iterCol] += other.itsPreconditionerSlice.find(*iterCol)->second;
            }
          }
        }
      }
      catch (const std::bad_cast &bc) {
        ASKAPTHROW(AskapError, "An attempt to merge NormalEquations with an "
                   "equation of incompatible type");
      }
    }


    /// @brief Regrid and add new parameter
    /// @details Regrid new image parameter, which is assume to be and image,
    /// onto the current image grid, which is assumed to be of an appropriate
    /// extent.
    /// @param[in] other normal equations
    /// @param[in] name of parameter under consideration
    void ImagingNormalEquations::linmosMerge(const ImagingNormalEquations &other, const string col)
    {

      ASKAPASSERT(itsShape[col].nelements() >= 2);
      ASKAPASSERT(itsShape[col].nelements() == other.itsShape.find(col)->second.nelements());

      ASKAPASSERT(itsWeightState == other.itsWeightState);
      ASKAPASSERT(itsWeightType == other.itsWeightType);


      // initialise an image accumulator
      imagemath::LinmosAccumulator<double> accumulator;
      accumulator.weightType(itsWeightType);
      accumulator.weightState(itsWeightState);

      std::cout << "linmosMerge - Weight state/type = " << accumulator.weightState() << "/" << accumulator.weightType() << std::endl;

      accumulator.setDefaultPB(); // this will probably need to be smarter
      // these inputs should be set up to take the full mosaic.
      accumulator.setOutputParameters(itsShape[col], itsCoordSys[col]);

      // I think these are essentially copy by reference which means the
      // accumulation is happening in place. This needs to deal with the weightState-

      casacore::Array<double> outPix(itsDataVector[col].reform(accumulator.outShape()));
      casacore::Array<double> outWgtPix(itsNormalMatrixDiagonal[col].reform(accumulator.outShape()));
      casacore::Array<double> outSenPix(accumulator.outShape(),0.);

      accumulator.setInputParameters(other.itsShape.find(col)->second,
                                     other.itsCoordSys.find(col)->second);

      casacore::Array<double> inPix(other.itsDataVector.find(col)->second.reform(accumulator.inShape()));
      casacore::Array<double> inWgtPix(other.itsNormalMatrixDiagonal.find(col)->second.reform(accumulator.inShape()));
      casacore::Array<double> inSenPix(accumulator.inShape(),1.);

      if ( accumulator.outputBufferSetupRequired() ) {
        accumulator.initialiseRegridder();
      }
      accumulator.initialiseOutputBuffers();
      accumulator.initialiseInputBuffers();

      // loop over non-direction axes (e.g. spectral and/or polarisation)
      IPosition curpos(accumulator.inShape());
      for (uInt dim=0; dim<curpos.nelements(); ++dim) {
        curpos[dim] = 0;
      }
      scimath::MultiDimArrayPlaneIter planeIter(accumulator.inShape());
      for (; planeIter.hasMore(); planeIter.next()) {
        curpos = planeIter.position();

        if (itsWeightState == CORRECTED) {
          std::cout << "weighting stored image plane" << std::endl;
          accumulator.weightPlane(outPix, outWgtPix, outSenPix,curpos);
        }
        // load input buffer for the current plane
        accumulator.loadAndWeightInputBuffers(curpos, inPix, inWgtPix, inSenPix);
        // call regrid for any buffered images
        accumulator.regrid();
        // update the accululation arrays for this plane
        accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, curpos);
      }
      std::cout << "Deweighting accumulated images" << std::endl;
      scimath::MultiDimArrayPlaneIter deweightIter(accumulator.outShape());
      for (; deweightIter.hasMore(); deweightIter.next()) {
        curpos = deweightIter.position();
        accumulator.deweightPlane(outPix, outWgtPix, outSenPix, curpos);
      }

    }

    // normalMatrixDiagonal
    const casacore::Vector<double>& ImagingNormalEquations::normalMatrixDiagonal(const std::string &par) const
    {
      std::map<std::string, casacore::Vector<double> >::const_iterator cIt =
                                        itsNormalMatrixDiagonal.find(par);
      ASKAPASSERT(cIt != itsNormalMatrixDiagonal.end());
      return cIt->second;
    }

    const std::map<std::string, casacore::Vector<double> >& ImagingNormalEquations::normalMatrixDiagonal() const
    {
      return itsNormalMatrixDiagonal;
    }

    // normalMatrixSlice
    const casacore::Vector<double>& ImagingNormalEquations::normalMatrixSlice(const std::string &par) const
    {
      std::map<std::string, casacore::Vector<double> >::const_iterator cIt =
                                        itsNormalMatrixSlice.find(par);
      ASKAPASSERT(cIt != itsNormalMatrixSlice.end());
      return cIt->second;
    }

    const std::map<std::string, casacore::Vector<double> >& ImagingNormalEquations::normalMatrixSlice() const
    {
      return itsNormalMatrixSlice;
    }

    // preconditionerSlice
    const casacore::Vector<double>& ImagingNormalEquations::preconditionerSlice(const std::string &par) const
    {
      std::map<std::string, casacore::Vector<double> >::const_iterator cIt =
                                        itsPreconditionerSlice.find(par);
      ASKAPASSERT(cIt != itsPreconditionerSlice.end());
      return cIt->second;
    }

    const std::map<std::string, casacore::Vector<double> >& ImagingNormalEquations::preconditionerSlice() const
    {
      return itsPreconditionerSlice;
    }

/// @brief normal equations for given parameters
/// @details In the current framework, parameters are essentially
/// vectors, not scalars. Each element of such vector is treated
/// independently (but only vector as a whole can be fixed). As a
/// result any element of the normal matrix is another matrix for
/// all non-scalar parameters. For scalar parameters each such
/// matrix has a shape of [1,1].
/// @param[in] par1 the name of the first parameter
/// @param[in] par2 the name of the second parameter
/// @return one element of the sparse normal matrix (a dense matrix)
const casacore::Matrix<double>& ImagingNormalEquations::normalMatrix(const std::string &par1,
                        const std::string &par2) const
{
   ASKAPTHROW(AskapError,
               "ImagingNormalEquations::normalMatrix has not yet been implemented, attempted access to elements par1="<<
               par1<<" and par2="<<par2);
}

/// @brief data vector for a given parameter
/// @details In the current framework, parameters are essentially
/// vectors, not scalars. Each element of such vector is treated
/// independently (but only vector as a whole can be fixed). As a
/// result any element of the normal matrix as well as an element of the
/// data vector are, in general, matrices, not scalar. For the scalar
/// parameter each element of data vector is a vector of unit length.
/// @param[in] par the name of the parameter of interest
/// @return one element of the sparse data vector (a dense vector)
const casacore::Vector<double>& ImagingNormalEquations::dataVector(const std::string &par) const
{
   std::map<std::string, casacore::Vector<double> >::const_iterator cIt =
                                     itsDataVector.find(par);
   ASKAPASSERT(cIt != itsDataVector.end());
   return cIt->second;
}

const std::map<std::string, casacore::Vector<double> >& ImagingNormalEquations::dataVector() const
{
  return itsDataVector;
}

/// Return shape
    const std::map<std::string, casacore::IPosition >& ImagingNormalEquations::shape() const
    {
      return itsShape;
    }

/// Return reference
    const std::map<std::string, casacore::IPosition >& ImagingNormalEquations::reference() const
    {
      return itsReference;
    }

/// Return coordinate system
    const std::map<std::string, casacore::CoordinateSystem >& ImagingNormalEquations::coordSys() const
    {
      return itsCoordSys;
    }

    void ImagingNormalEquations::reset()
    {
      map<std::string, casacore::Vector<double> >::iterator iterRow;
      for (iterRow=itsDataVector.begin();iterRow!=itsDataVector.end();iterRow++)
      {
        itsDataVector[iterRow->first].resize();
        itsDataVector[iterRow->first]=casacore::Vector<double>(0);
        itsShape[iterRow->first].resize(0);
        itsShape[iterRow->first]=casacore::IPosition();
        itsReference[iterRow->first].resize(0);
        itsReference[iterRow->first]=casacore::IPosition();
        // leave "keep" and "replace" axis vectors empty, so they are all removed.
        Vector<Int> worldAxes;
        Vector<Double> worldRep;
        CoordinateUtil::removeAxes(itsCoordSys[iterRow->first], worldRep, worldAxes, False);
        itsCoordSys[iterRow->first] = casacore::CoordinateSystem();
        itsNormalMatrixSlice[iterRow->first].resize();
        itsNormalMatrixSlice[iterRow->first]=casacore::Vector<double>(0);
        itsNormalMatrixDiagonal[iterRow->first].resize();
        itsNormalMatrixDiagonal[iterRow->first]=casacore::Vector<double>(0);
        itsPreconditionerSlice[iterRow->first].resize();
        itsPreconditionerSlice[iterRow->first]=casacore::Vector<double>(0);
      }
    }

    void ImagingNormalEquations::addSlice(const string& name,
      const casacore::Vector<double>& normalmatrixslice,
      const casacore::Vector<double>& normalmatrixdiagonal,
      const casacore::Vector<double>& preconditionerslice,
      const casacore::Vector<double>& datavector,
      const casacore::IPosition& shape,
      const casacore::IPosition& reference,
      const casacore::CoordinateSystem& coordSys)
    {
      ASKAPTRACE("ImagingNormalEquations::addSlice");

      // if coordinate systems exist, make sure they're equal
      if(itsCoordSys[name].nCoordinates()>0)
      {
        // initialise an image accumulator
        imagemath::LinmosAccumulator<double> accumulator;
        ASKAPCHECK(accumulator.coordinatesAreEqual(itsCoordSys[name], coordSys,
                                                   itsShape[name], shape),
            "Cannot combine slices with different coord systems using addSlice. Use merge.");
      }

      if(datavector.size()!=itsDataVector[name].size())
      {
        ASKAPDEBUGASSERT(itsDataVector[name].size() == 0);
        itsDataVector[name]=datavector;
      }
      else
      {
        itsDataVector[name]+=datavector;
      }
      if(normalmatrixdiagonal.shape()!=itsNormalMatrixDiagonal[name].shape())
      {
        ASKAPDEBUGASSERT(itsNormalMatrixDiagonal[name].size() == 0);
        itsNormalMatrixDiagonal[name]=normalmatrixdiagonal;
      }
      else
      {
        itsNormalMatrixDiagonal[name]+=normalmatrixdiagonal;
      }
      if(normalmatrixslice.shape()!=itsNormalMatrixSlice[name].shape())
      {
        itsNormalMatrixSlice[name]=normalmatrixslice;
      }
      else
      {
        itsNormalMatrixSlice[name]+=normalmatrixslice;
      }
      if(preconditionerslice.shape()!=itsPreconditionerSlice[name].shape())
      {
        ASKAPDEBUGASSERT(itsPreconditionerSlice[name].size() == 0);
        itsPreconditionerSlice[name]=preconditionerslice;
      }
      else
      {
        itsPreconditionerSlice[name]+=preconditionerslice;
      }
      itsShape[name].resize(0);
      itsShape[name]=shape;
      itsReference[name].resize(0);
      itsReference[name]=reference;
      if(itsCoordSys[name].nCoordinates()==0)
      {
          itsCoordSys[name]=coordSys;
      }

    }

    /// @brief Store slice of the normal matrix for a given parameter.
    ///
    /// This means
    /// that the cross terms between parameters are excluded and only
    /// a slice of the normal matrix is retained.
    /// @param name Name of parameter
    /// @param normalmatrixslice Slice of normal matrix for this parameter
    /// @param normalmatrixdiagonal Diagonal of normal matrix for
    ///        this parameter
    /// @param datavector Data vector for this parameter
    /// @param reference Reference point for the slice
    void ImagingNormalEquations::addSlice(const string& name,
                    const casacore::Vector<double>& normalmatrixslice,
                    const casacore::Vector<double>& normalmatrixdiagonal,
                    const casacore::Vector<double>& preconditionerslice,
                    const casacore::Vector<double>& datavector,
                    const casacore::IPosition& reference)
    {
      addSlice(name,normalmatrixslice,normalmatrixdiagonal,preconditionerslice,
               datavector,casacore::IPosition(1,datavector.nelements()),reference);
    }

    void ImagingNormalEquations::addDiagonal(const string& name, const casacore::Vector<double>& normalmatrixdiagonal,
      const casacore::Vector<double>& datavector, const casacore::IPosition& shape)
    {
      ASKAPTRACE("ImagingNormalEquations::addDiagonal");

      if(datavector.size()!=itsDataVector[name].size())
      {
        ASKAPDEBUGASSERT(itsDataVector[name].size() == 0);
        itsDataVector[name]=datavector;
      }
      else
      {
        itsDataVector[name]+=datavector;
      }
      if(normalmatrixdiagonal.shape()!=itsNormalMatrixDiagonal[name].shape())
      {
        ASKAPDEBUGASSERT(itsNormalMatrixDiagonal[name].size() == 0);
        itsNormalMatrixDiagonal[name]=normalmatrixdiagonal;
      }
      else
      {
        itsNormalMatrixDiagonal[name]+=normalmatrixdiagonal;
      }
      itsShape[name].resize(0);
      itsShape[name]=shape;
    }

    void ImagingNormalEquations::addDiagonal(const string& name, const casacore::Vector<double>& normalmatrixdiagonal,
      const casacore::Vector<double>& datavector)
    {
      casacore::IPosition shape(1, datavector.nelements());
      addDiagonal(name, normalmatrixdiagonal, datavector, shape);
    }

    INormalEquations::ShPtr ImagingNormalEquations::clone() const
    {
      return INormalEquations::ShPtr(new ImagingNormalEquations(*this));
    }

    /// @brief write the object to a blob stream
    /// @param[in] os the output stream
    void ImagingNormalEquations::writeToBlob(LOFAR::BlobOStream& os) const
    {
      os << itsNormalMatrixSlice << itsNormalMatrixDiagonal
         << itsPreconditionerSlice << itsShape
         << itsReference << itsCoordSys << itsDataVector;
    }

    /// @brief read the object from a blob stream
    /// @param[in] is the input stream
    /// @note Not sure whether the parameter should be made const or not
    void ImagingNormalEquations::readFromBlob(LOFAR::BlobIStream& is)
    {
      is >> itsNormalMatrixSlice >> itsNormalMatrixDiagonal
         >> itsPreconditionerSlice >> itsShape
         >> itsReference >> itsCoordSys >> itsDataVector;
    }

    /// @brief obtain all parameters dealt with by these normal equations
    /// @details Normal equations provide constraints for a number of
    /// parameters (i.e. unknowns of these equations). This method returns
    /// a vector with the string names of all parameters mentioned in the
    /// normal equations represented by the given object.
    /// @return a vector listing the names of all parameters (unknowns of these equations)
    /// @note if ASKAP_DEBUG is set some extra checks on consistency of these
    /// equations are done
    std::vector<std::string> ImagingNormalEquations::unknowns() const {
      std::vector<std::string> result;
      result.reserve(itsNormalMatrixSlice.size());
      for (std::map<std::string, casacore::Vector<double> >::const_iterator ci = itsNormalMatrixSlice.begin();
           ci!=itsNormalMatrixSlice.end(); ++ci) {
           result.push_back(ci->first);

// extra checks in debug mode
#ifdef ASKAP_DEBUG
           ASKAPCHECK(itsNormalMatrixDiagonal.find(ci->first) != itsNormalMatrixDiagonal.end(),
                      "Parameter "<<ci->first<<" is present in the matrix slice but is missing in the diagonal");
           ASKAPCHECK(itsShape.find(ci->first) != itsShape.end(),
                      "Parameter "<<ci->first<<" is present in the matrix slice but is missing in the shape map");
           ASKAPCHECK(itsReference.find(ci->first) != itsReference.end(),
                      "Parameter "<<ci->first<<" is present in the matrix slice but is missing in the reference map");
           ASKAPCHECK(itsDataVector.find(ci->first) != itsDataVector.end(),
                      "Parameter "<<ci->first<<" is present in the matrix slice but is missing in the data vector");
#endif // #ifdef ASKAP_DEBUG
      }
      return result;
    } // unknowns method
    /// get the weightstate
    int ImagingNormalEquations::weightState() {
      return itsWeightState;
    } // get the weightstate

    /// set the weightState
    void ImagingNormalEquations::weightState(int theState) {
      this->itsWeightState = theState;
    }

    /// get the weightType
    int ImagingNormalEquations::weightType() {
      return itsWeightType;
    }

    /// set the weightType.
    void ImagingNormalEquations::weightType(int theType){
      this->itsWeightType = theType;
    }
  } // namespace scimath
} // namespace askap
