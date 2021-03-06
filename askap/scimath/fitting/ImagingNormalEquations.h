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
#ifndef IMAGING_NORMAL_EQUATIONS_H
#define IMAGING_NORMAL_EQUATIONS_H

#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/INormalEquations.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Slice.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>

#include <Blob/BlobOStream.h>
#include <Blob/BlobIStream.h>

#include <boost/shared_ptr.hpp>

namespace askap
{
  namespace scimath
  {

    class DesignMatrix;

/// @brief Normal equations with an approximation for imaging
/// @details There are two kinds of normal equations currently supported. The
/// first one is a generic case, where the full normal matrix is retained. It
/// is used for calibration. The second one is intended for imaging, where we
/// can't afford to keep the whole normal matrix. In the latter approach, the
/// matrix is approximated by a sum of diagonal and shift invariant matrices.
/// This class represents the approximated case, and is used with imaging
/// algorithms.
/// @ingroup fitting
    class ImagingNormalEquations : public INormalEquations
    {
    public:

      ImagingNormalEquations();

      /// @brief Construct for the specified parameters
      ///
      /// Initialisation does not allocate much memory.
      /// @param ip Parameters
      ImagingNormalEquations(const Params& ip);

      /// @brief copy constructor
      /// @details Data members of this class are non-trivial types including
      /// std containers of casa containers. The letter are copied by reference by default. We,
      /// therefore, need this copy constructor to achieve proper copying.
      /// @param[in] src input measurement equations to copy from
      ImagingNormalEquations(const ImagingNormalEquations &src);

      /// @brief assignment operator
      /// @details Data members of this class are non-trivial types including
      /// std containers of casa containers. The letter are copied by reference by default. We,
      /// therefore, need this copy constructor to achieve proper copying.
      /// @param[in] src input measurement equations to copy from
      /// @return reference to this object
      ImagingNormalEquations& operator=(const ImagingNormalEquations &src);

      virtual ~ImagingNormalEquations();

      /// @brief Store slice of the normal matrix for a given parameter.
      /// This uses the shape parametr to determince how much of the datavector you
      /// take.
      ///
      /// @param[in] name parameter name
      /// @param[in] normalmatrixslice Slice of normal matrix for this parameter
      /// @param[in] normalmatrixdiagonal Diagonal of normal matrix for
      ///        this parameter
      /// @param[in] datavector Data vector for this parameter
      /// @param[in] shape Shape of this parameter
      /// @param[in] reference Reference point for the slice
      /// @param[in] the coordinate system of the slice

      void addSlice(const std::string& name,
                    const casacore::Vector<double>& normalmatrixslice,
                    const casacore::Vector<double>& normalmatrixdiagonal,
                    const casacore::Vector<double>& preconditionerslice,
                    const casacore::Vector<double>& datavector,
                    const casacore::IPosition& shape,
                    const casacore::IPosition& reference,
                    const casacore::CoordinateSystem& coordSys = casacore::CoordinateSystem());

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
      void addSlice(const std::string& name,
                    const casacore::Vector<double>& normalmatrixslice,
                    const casacore::Vector<double>& normalmatrixdiagonal,
                    const casacore::Vector<double>& preconditionerslice,
                    const casacore::Vector<double>& datavector,
                    const casacore::IPosition& reference);

      /// @brief Store diagonal of the normal matrix for a given parameter.
      ///
      /// This means
      /// that the cross terms between parameters are excluded and only
      /// the diagonal inside a parameter is kept.
      /// @param name Name of parameter
      /// @param normalmatrixdiagonal Normal Matrix for this parameter
      /// @param datavector Data vector for this parameter
      /// @param shape Shape of this parameter
      void addDiagonal(const std::string& name,
                       const casacore::Vector<double>& normalmatrixdiagonal,
                       const casacore::Vector<double>& datavector,
                       const casacore::IPosition& shape);

      /// @brief Store diagonal of the normal matrix for a given parameter.
      ///
      /// This means
      /// that the cross terms between parameters are excluded and only
      /// the diagonal inside a parameter is kept.
      /// @param name Name of parameter
      /// @param normalmatrixdiagonal Normal Matrix for this parameter
      /// @param datavector Data vector for this parameter
      void addDiagonal(const std::string& name,
                       const casacore::Vector<double>& normalmatrixdiagonal,
                       const casacore::Vector<double>& datavector);



      /// @brief zero all the data elements but keep thier size and coordinate systems
      /// @details I have ound I need a utility function that allows me to zero the contents
      /// of the ImagingNormalEquations but to keep the meta information for the params
      /// @param[in] ip the params
      /// @return void

      void zero(const Params& ip);
      /// @brief Merge these normal equations with another
      /// @details Combining two normal equations depends on the actual class type
      /// (different work is required for a full matrix and for an approximation).
      /// This method must be overriden in the derived classes for correct
      /// implementation.
      /// This means that we just add
      /// @param[in] src an object to get the normal equations from
      virtual void merge(const INormalEquations& src);

      virtual void linmosMerge(const ImagingNormalEquations &other, const string col);

      /// @brief normal equations for given parameters
      /// @details In the current framework, parameters are essentially
      /// vectors, not scalars. Each element of such vector is treated
      /// independently (but only vector as a whole can be fixed). As a
      /// result the element of the normal matrix is another matrix for
      /// all non-scalar parameters. For scalar parameters each such
      /// matrix has a shape of [1,1].
      /// @param[in] par1 the name of the first parameter
      /// @param[in] par2 the name of the second parameter
      /// @return one element of the sparse normal matrix (a dense matrix)
      virtual const casacore::Matrix<double>& normalMatrix(const std::string &par1,
                       const std::string &par2) const;

      /// @brief Return normal equations slice
      virtual const casacore::Vector<double>& normalMatrixSlice(const std::string &par) const;

      /// @brief Return normal equations diagonal
      virtual const casacore::Vector<double>& normalMatrixDiagonal(const std::string &par) const;

      /// @brief Return preconditioner slice (otherwise normalMatrixSlice is used for preconditioning)
      virtual const casacore::Vector<double>& preconditionerSlice(const std::string &par) const;

      /// @brief data vector for a given parameter
      /// @details In the current framework, parameters are essentially
      /// vectors, not scalars. Each element of such vector is treated
      /// independently (but only vector as a whole can be fixed). As a
      /// result any element of the normal matrix as well as an element of the
      /// data vector are, in general, matrices, not scalar. For the scalar
      /// parameter each element of data vector is a vector of unit length.
      /// @param[in] par the name of the parameter of interest
      /// @return one element of the sparse data vector (a dense vector)
      virtual const casacore::Vector<double>& dataVector(const std::string &par) const;

      /// @brief obtain all parameters dealt with by these normal equations
      /// @details Normal equations provide constraints for a number of
      /// parameters (i.e. unknowns of these equations). This method returns
      /// a vector with the string names of all parameters mentioned in the
      /// normal equations represented by the given object.
      /// @return a vector listing the names of all parameters (unknowns of these equations)
      /// @note if ASKAP_DEBUG is set some extra checks on consistency of these
      /// equations are done
      virtual std::vector<std::string> unknowns() const;

      /// Return normal equations slice map
      const std::map<std::string, casacore::Vector<double> >& normalMatrixSlice() const;

      /// Return normal equations diagonal map
      const std::map<std::string, casacore::Vector<double> >& normalMatrixDiagonal() const;

      /// Return preconditioner slice map
      const std::map<std::string, casacore::Vector<double> >& preconditionerSlice() const;

      /// Return data vector map
      const std::map<std::string, casacore::Vector<double> >& dataVector() const;

      /// Return shape map
      const std::map<std::string, casacore::IPosition>& shape() const;

      /// Return references map
      const std::map<std::string, casacore::IPosition >& reference() const;

      /// Return coordinate system map
      const std::map<std::string, casacore::CoordinateSystem >& coordSys() const;

      /// Reset to empty
      virtual void reset();

      /// Shared pointer definition
      typedef boost::shared_ptr<ImagingNormalEquations> ShPtr;

      /// Clone this into a shared pointer
      virtual INormalEquations::ShPtr clone() const;

      /// @brief write the object to a blob stream
      /// @param[in] os the output stream
      virtual void writeToBlob(LOFAR::BlobOStream& os) const;

      /// @brief read the object from a blob stream
      /// @param[in] is the input stream
      /// @note Not sure whether the parameter should be made const or not
      virtual void readFromBlob(LOFAR::BlobIStream& is);

      /// get the weightstate
      int weightState();

      /// set the weightState
      void weightState(int theState);

      /// get the weightType
      int weightType();

      /// set the weightType.
      void weightType(int theType);

    private:
      /// A slice through a specified plane
      std::map<std::string, casacore::Vector<double> > itsNormalMatrixSlice;
      /// The diagonal
      std::map<std::string, casacore::Vector<double> > itsNormalMatrixDiagonal;
      /// A slice through a specified plane
      std::map<std::string, casacore::Vector<double> > itsPreconditionerSlice;
      /// The shape
      std::map<std::string, casacore::IPosition> itsShape;
      /// The Reference point for the slice
      std::map<std::string, casacore::IPosition> itsReference;
      /// The coordinate system for the slice
      std::map<std::string, casacore::CoordinateSystem> itsCoordSys;
      /// The data vectors
      std::map<std::string, casacore::Vector<double> > itsDataVector;

      /// The weighting type and state of the params stored in the NE
      int itsWeightType;
      int itsWeightState;

    };

  }  // namespace scimath
} // namespace askap
#endif // #ifndef IMAGING_NORMAL_EQUATIONS_H
