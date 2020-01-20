/// @file
/// @brief Normal equations without any approximation
/// @details There are two kinds of normal equations currently supported. The
/// first one is a generic case, where the full normal matrix is retained. It
/// is used for calibration. The second one is intended for imaging, where we
/// can't afford to keep the whole normal matrix. In this approach, the matrix
/// is approximated by a sum of diagonal and shift invariant matrices. This
/// class represents the generic case, where no approximation to the normal
/// matrix is done. Implementation of this class is largely taken from the
/// old NormalEquation (revisions up to 4637) written by Tim Cornwell.
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
/// @author Vitaliy Ogarko <vogarko@gmail.com>
///

// own includes
#include <askap/scimath/fitting/GenericNormalEquations.h>
#include <askap/scimath/fitting/DesignMatrix.h>
#include <askap/scimath/fitting/CalParamNameHelper.h>
#include <askap/askap/AskapError.h>
#include <askap/scimath/utils/DeepCopyUtils.h>

#include <Blob/BlobArray.h>
#include <Blob/BlobSTL.h>

// std includes
#include <utility>
#include <set>
#include <stdexcept>

// casa includes
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/OS/Timer.h>

// logging stuff
#include <askap_scimath.h>
#include <askap/askap/AskapLogging.h>

using namespace LOFAR;

using casacore::product;
using casacore::transpose;

namespace askap { namespace scimath {

casacore::Matrix<double> emptyMatrix;

ASKAP_LOGGER(logger, ".genericne");

/// @brief a default constructor
/// @details It creates an empty normal equations class
GenericNormalEquations::GenericNormalEquations() {}
  
/// @brief constructor from a design matrix
/// @details This version of the constructor is equivalent to an
/// empty constructor plus a call to add method with the given
/// design matrix
/// @param[in] dm Design matrix to use
GenericNormalEquations::GenericNormalEquations(const DesignMatrix& dm)
{
  add(dm);
}

/// @brief copy constructor
/// @details It is required because this class has non-trivial types (std containers
/// of casa containers)
/// @param[in] src other class
GenericNormalEquations::GenericNormalEquations(const GenericNormalEquations &src) :
        INormalEquations(src),
        itsParameterNameToIndex(src.itsParameterNameToIndex),
        itsParameterIndexToBaseName(src.itsParameterIndexToBaseName),
        itsMetadata(src.itsMetadata)
{
  deepCopyOfSTDMap(src.itsDataVector, itsDataVector);
  for (std::map<std::string, MapOfMatrices>::const_iterator ci = src.itsNormalMatrix.begin();
       ci!=src.itsNormalMatrix.end(); ++ci) {
       deepCopyOfSTDMap(ci->second, itsNormalMatrix[ci->first]);
  }
}

/// @brief assignment operator
/// @details It is required because this class has non-trivial types (std containers
/// of casa containers)
/// @param[in] src other class
/// @return reference to this object
GenericNormalEquations& GenericNormalEquations::operator=(const GenericNormalEquations &src)
{
  if (&src != this) {
      itsDataVector.clear();
      deepCopyOfSTDMap(src.itsDataVector, itsDataVector);
      itsNormalMatrix.clear();
      for (std::map<std::string, MapOfMatrices>::const_iterator ci = src.itsNormalMatrix.begin();
           ci!=src.itsNormalMatrix.end(); ++ci) {
           deepCopyOfSTDMap(ci->second, itsNormalMatrix[ci->first]);
      }
      itsParameterNameToIndex = src.itsParameterNameToIndex;
      itsParameterIndexToBaseName = src.itsParameterIndexToBaseName;
      itsMetadata = src.itsMetadata;
  }
  return *this;
}

/// @brief reset the normal equation object
/// @details After a call to this method the object has the same pristine
/// state as immediately after creation with the default constructor
void GenericNormalEquations::reset()
{
  itsDataVector.clear();
  itsNormalMatrix.clear();
  itsParameterNameToIndex.clear();
  itsParameterIndexToBaseName.clear();
  itsMetadata.reset();
}

/// @brief Clone this into a shared pointer
/// @details "Virtual constructor" - creates a copy of this object. Derived
/// classes must override this method to instantiate the object of a proper 
/// type.
/// @return shared pointer on INormalEquation class
GenericNormalEquations::ShPtr GenericNormalEquations::clone() const
{
  return ShPtr(new GenericNormalEquations(*this));
}

/// @brief Merge these normal equations with another
/// @details Combining two normal equations depends on the actual class type
/// (different work is required for a full matrix and for an approximation).
/// This method must be overriden in the derived classes for correct 
/// implementation. 
/// This means that we just add
/// @param[in] src an object to get the normal equations from
void GenericNormalEquations::merge(const INormalEquations& src) 
{
   try {
      casacore::Timer timer;
      timer.mark();
      ASKAPLOG_INFO_STR(logger, "Merging normal equations");

      const GenericNormalEquations &gne = 
                dynamic_cast<const GenericNormalEquations&>(src);

      // loop over all parameters, add them one by one.
      // We could have passed iterator directly to mergeParameter and it
      // would work faster (no extra search accross the map). But current
      // code is more readable.
      for (MapOfVectors::const_iterator ci = gne.itsDataVector.begin(); 
           ci != gne.itsDataVector.end(); ++ci) {
           mergeParameter(ci->first, gne);
      }
      itsMetadata.merge(gne.metadata());

      // TODO: We also need to merge itsParameterNameToIndex & itsParameterIndexToName.
      //       Let's keep it in mind, and add this once we really need it.

      ASKAPLOG_INFO_STR(logger, "Merged normal equations in "<< timer.real() << " seconds");
   }
   catch (const AskapError &) {
      throw;
   }
   catch (const std::bad_cast &) {
      throw AskapError("Attempt to use GenericNormalEquation::merge with "
                        "incompatible type of the normal equation class");
   }
}

/// @brief Add one parameter from another normal equations class
/// @details This helper method is used in merging two normal equations.
/// It processes just one parameter.
/// @note This helper method works with instances of this class only (as
/// only then it knows how the actual normal matrix is handled). One could
/// have a general code which would work for every possible normal equation,
/// but in some cases it would be very inefficient. Therefore, the decision
/// has been made to throw an exception if incompatible operation is requested
/// and add the code to handle this situation later, if it appears to be 
/// necessary.
/// @param[in] par name of the parameter to copy
/// @param[in] src an object to get the normal equations from
void GenericNormalEquations::mergeParameter(const std::string &par, 
                              const GenericNormalEquations& src)
{
   // srcItRow is an iterator over rows in source matrix;
   // by analogy, destItCol is an iterator over columns in the destination matrix
   const std::map<std::string, MapOfMatrices>::const_iterator srcItRow = 
                         src.itsNormalMatrix.find(par);
   ASKAPDEBUGASSERT(srcItRow != src.itsNormalMatrix.end());

   const MapOfVectors::const_iterator srcItData = src.itsDataVector.find(par);
   ASKAPDEBUGASSERT(srcItData != src.itsDataVector.end());

   addParameterSparsely(par, srcItRow->second, srcItData->second);
}

/// @brief Add/update one parameter using given matrix and data vector
/// @details This helper method is the main workhorse used in merging two
/// normal equations, adding an independent parameter or a design matrix.
/// The normal matrix to be integrated with this class is given in the form
/// of map of matrices (effectively a sparse matrix). Each element of the map
/// corresponds to a cross- or parallel term in the normal equations. Data
/// vector is given simply as a casacore::Vector, rather than the map of vectors,
/// because only one parameter is concerned here. If a parameter with the given
/// name doesn't exist, the method adds it to both normal matrix and data vector,
/// populating correctly all required cross-terms with 0-matrices of an 
/// appropriate shape.
/// @param[in] par name of the parameter to work with
/// @param[in] inNM input normal matrix
/// @param[in] inDV input data vector 
void GenericNormalEquations::addParameter(const std::string &par, 
           const MapOfMatrices &inNM, const casacore::Vector<double>& inDV)
{
  // before processing this row, check that the columns already exist as rows
  for (MapOfMatrices::const_iterator nmColIt = inNM.begin();
                      nmColIt != inNM.end(); ++nmColIt) {
      if (itsNormalMatrix.find(nmColIt->first) == itsNormalMatrix.end()) {
          // this is a new parameter. Add it to the NM
          std::map<std::string, MapOfMatrices>::iterator
                  newRowIt = itsNormalMatrix.insert(std::make_pair(
                                 nmColIt->first,MapOfMatrices())).first;
          const casacore::uInt rowParDim = nmColIt->second.ncolumn();
          // set all of the new parameter columns to zero
          std::map<std::string, MapOfMatrices>::iterator newColIt;
          for (newColIt = itsNormalMatrix.begin();
                     newColIt != itsNormalMatrix.end(); ++newColIt) {
              casacore::uInt colParDim;
              if (newColIt->first == newRowIt->first) {
                  // this is the new parameter.
                  colParDim = rowParDim;
              } else {
                  // this is an old parameter. Get size from itsNormalMatrix
                  colParDim = parameterDimension(newColIt->second);
              }
              newRowIt->second.insert(std::make_pair(newColIt->first,
                  casacore::Matrix<double>(rowParDim, colParDim, 0.)));
              // fill in the symmetric term
              if (newRowIt->first != newColIt->first) {
                  newColIt->second.insert(std::make_pair(newRowIt->first,
                      casacore::Matrix<double>(colParDim, rowParDim, 0.)));
              }
          }
      }
  }

  // first, process normal matrix 
  // nmRowIt is an iterator over rows (the outer map) of the normal matrix
  // stored in this class 
  std::map<std::string, MapOfMatrices>::iterator nmRowIt =
                         itsNormalMatrix.find(par);
  ASKAPDEBUGASSERT(nmRowIt != itsNormalMatrix.end());
  ASKAPDEBUGASSERT(nmRowIt->second.find(par) != nmRowIt->second.end());

  for (MapOfMatrices::const_iterator inNMIt = inNM.begin();
          inNMIt != inNM.end(); ++inNMIt) {

      // search for an appropriate parameter in the normal matrix
      MapOfMatrices::iterator nmColIt = nmRowIt->second.find(inNMIt->first);
      // work with cross-terms if the input matrix have them
      if (nmColIt != nmRowIt->second.end()) {
          ASKAPCHECK(inNMIt->second.shape() == nmColIt->second.shape(),
                   "shape mismatch for normal matrix, parameters ("<<
                   nmColIt->first<<" , "<<nmRowIt->first<<"). "<<
                   nmColIt->second.shape()<<" != "<<inNMIt->second.shape());
          nmColIt->second += inNMIt->second; // add up a matrix
      }
  }

  // Processing the data vector.
  addDataVector(par, inDV);
}

void GenericNormalEquations::addParameterSparselyToMatrix(const std::string &par, const MapOfMatrices &inNM)
{
    // First, process normal matrix.
    // nmRowIt is an iterator over rows (the outer map) of the normal matrix stored in this class
    std::map<std::string, MapOfMatrices>::iterator nmRowIt = itsNormalMatrix.find(par);
    if (nmRowIt == itsNormalMatrix.end()) {
    // Parameter not found in the matrix, adding the corresponding row.
        nmRowIt = itsNormalMatrix.insert(std::make_pair(par, MapOfMatrices())).first;
    }

    for (MapOfMatrices::const_iterator inNMIt = inNM.begin();
        inNMIt != inNM.end(); ++inNMIt) {

        // Search for an appropriate parameter in the normal matrix.
        MapOfMatrices::iterator nmColIt = nmRowIt->second.find(inNMIt->first);
        if (nmColIt == nmRowIt->second.end()) {

            // Extract parameter dimensions.
            const casacore::uInt rowParDim = inNMIt->second.nrow();
            const casacore::uInt colParDim = inNMIt->second.ncolumn();

            // Initialize a column with an empty matrix.
            // Note: only the columns used in the calculation get initialized.
            nmColIt = nmRowIt->second.insert(
                            std::make_pair(inNMIt->first,
                                        casacore::Matrix<double>(rowParDim, colParDim, 0.))).first;
        }

        // Work with cross-terms if the input matrix have them.
        ASKAPCHECK(inNMIt->second.shape() == nmColIt->second.shape(),
                "shape mismatch for normal matrix, parameters ("<<
                nmColIt->first<<" , "<<nmRowIt->first<<"). "<<
                nmColIt->second.shape()<<" != "<<inNMIt->second.shape());
        nmColIt->second += inNMIt->second; // add up a matrix
    }
}

void GenericNormalEquations::addParameterSparsely(const std::string &par,
           const MapOfMatrices &inNM, const casacore::Vector<double>& inDV)
{
    addParameterSparselyToMatrix(par, inNM);

    // Processing the data vector.
    addDataVector(par, inDV);
}

size_t GenericNormalEquations::getNumberElements() const {
    size_t nElements = 0;
    for (std::map<std::string, MapOfMatrices>::const_iterator it = itsNormalMatrix.begin();
         it != itsNormalMatrix.end(); ++it) {

        if (it->second.begin() != it->second.end()) {
            size_t parRowDim = it->second.begin()->second.nrow();
            size_t parColDim = it->second.begin()->second.ncolumn();

            // All row elements should have the same dimension, hence multiply by the number of elements.
            nElements += (parRowDim * parColDim) * it->second.size();
        }
    }
    return nElements;
}

void GenericNormalEquations::addDataVector(const std::string &par, const casacore::Vector<double>& inDV)
{
    MapOfVectors::const_iterator dvIt = itsDataVector.find(par);
    if (dvIt != itsDataVector.end()) {
        ASKAPCHECK(inDV.shape() == dvIt->second.shape(),
                 "shape mismatch for data vector, parameter: "<<dvIt->first<<
                 ". "<<inDV.shape()<<" != "<<dvIt->second.shape());
        // we have to instantiate explicitly a casacore::Vector object because
        // otherwise, for some reason, the compiler can't figure out the type
        // properly at the += operator. Exploit reference semantics - no copying!
        casacore::Vector<double> destVec = dvIt->second;
        destVec += inDV; // add up a vector
    } else {
        itsDataVector.insert(std::make_pair(par, inDV));
    }
}

/// @brief extract dimension of a parameter from the given row
/// @details This helper method analyses the matrices stored in the supplied
/// map (effectively a row of a sparse matrix) and extracts the dimension of
/// the parameter this row corresponds to. If compiled with ASKAP_DEBUG, 
/// this method does an additional consistency check that all elements of
/// the sparse matrix give the same dimension (number of rows is the same for
/// all elements).
/// @param[in] nmRow a row of the sparse normal matrix to work with
/// @return dimension of the corresponding parameter
casacore::uInt GenericNormalEquations::parameterDimension(const MapOfMatrices &nmRow)
{
  const MapOfMatrices::const_iterator it = nmRow.begin();
  ASKAPDEBUGASSERT(it != nmRow.end());
  const casacore::uInt dim = it->second.nrow();
#ifdef ASKAP_DEBUG
  for (MapOfMatrices::const_iterator cur = it; cur != nmRow.end(); ++cur) {
       ASKAPASSERT(cur->second.nrow() == dim);
  }     
#endif // ASKAP_DEBUG
  return dim;
}

void GenericNormalEquations::add(const ComplexDiffMatrix &cdm, const PolXProducts &pxp, size_t columnOffset)
{
    if (pxp.nPol() == 0) {
        return; // nothing to process
    }
    ASKAPDEBUGASSERT(pxp.nPol() == cdm.nRow());
    const casacore::uInt nDataPoints = pxp.nPol();

    // Pre-calculate the array frequently used in the embedded loops.
    std::vector<casacore::DComplex> modelProductMatrix(nDataPoints * nDataPoints);
    for (casacore::uInt p1 = 0; p1 < nDataPoints; ++p1) {
        for (casacore::uInt p2 = 0; p2 < nDataPoints; ++p2) {
            modelProductMatrix[p2 + nDataPoints * p1] = pxp.getModelProduct(p1, p2);
        }
    }

    const std::complex<double> czero(0, 0);

    //-----------------------------------------------------------------------------
    // Form the data vector.
    //-----------------------------------------------------------------------------

    // The first loop is over polarisations, essentially summing over
    // data points in the calculation of normal matrix
    for (casacore::uInt p = 0; p < nDataPoints; ++p) {

        // Two inner loops are from matrix multiplication of cdm to a vector
        // for Y and A^T parts in the product forming the element of the data vector (projected).
        for (casacore::uInt p1 = 0; p1 < nDataPoints; ++p1) {

            const casacore::DComplex measProduct = pxp.getModelMeasProduct(p1, p);

            const ComplexDiff& cd1 = cdm(p, p1 + columnOffset);

            auto itRe1 = cd1.derivReBegin();
            auto itIm1 = cd1.derivImBegin();

            // Loop over the list of derivatives corresponding to the ComplexDiff cd1.
            for (; itRe1 != cd1.derivReEnd(); itRe1++, itIm1++) {

                ASKAPDEBUGASSERT(itIm1 != cd1.derivImEnd());

                const casacore::DComplex rowParDerivRe1 = itRe1->second;
                const casacore::DComplex rowParDerivIm1 = itIm1->second;

                // TODO: Perhaps after we inverted the loops to iterate only over derivatives of the current complex diff,
                // we dont need this check anymore, as there wont be zero derivatives after the first major iteration?
                if (rowParDerivRe1 == czero && rowParDerivIm1 == czero) {
                    continue;
                }

                // Partial projected data vector buffer for this row (size of 2 - all parameters are complex).
                // Elements correspond to real and imaginary part derivatives of projected residual.
                casacore::Vector<double> dataVector(2, 0.);

                // Parameter name corresponding to the row of the normal matrix.
                std::string rowName = itRe1->first;

                if (measProduct != czero) {
                    dataVector[0] += real(conj(rowParDerivRe1) * measProduct);
                    dataVector[1] += real(conj(rowParDerivIm1) * measProduct);
                }

                for (casacore::uInt p2 = 0; p2 < nDataPoints; ++p2) {

                    const casacore::DComplex modelProduct = modelProductMatrix[p2 + nDataPoints * p1];
                    if (modelProduct == czero) {
                        continue;
                    }

                    const ComplexDiff &cd2 = cdm(p, p2 + columnOffset);
                    const casacore::DComplex val2 = cd2.value();
                    const casacore::DComplex val2_modelProduct = val2 * modelProduct;

                    dataVector[0] -= real(conj(rowParDerivRe1) * val2_modelProduct);
                    dataVector[1] -= real(conj(rowParDerivIm1) * val2_modelProduct);
                }

                // Adding partial data vector.
                addDataVector(rowName, dataVector);
            }
        }
    }

    //-----------------------------------------------------------------------------
    // Now form the row of normal matrix.
    //-----------------------------------------------------------------------------

    // Buffer for the element of normal matrix
    // treat all parameters as complex here to simplify the logic
    // (and they're complex anyway) -> 2x2 matrix
    casacore::Matrix<casacore::Double> nmElementBuf(2, 2, 0.);

    // The first loop is over polarisations, essentially summing over
    // data points in the calculation of normal matrix.
    for (casacore::uInt p = 0; p < nDataPoints; ++p) {

        // Two inner loops are from matrix multiplication of cdm to a vector
        // for A and A^T parts in the product forming the element of the normal
        // matrix. We can optimise this for speed later, if proved to be a problem.
        for (casacore::uInt p1 = 0; p1 < nDataPoints; ++p1) {

            const ComplexDiff& cd1 = cdm(p, p1 + columnOffset);

            auto itRe1 = cd1.derivReBegin();
            auto itIm1 = cd1.derivImBegin();

            // Loop over the list of derivatives corresponding to the ComplexDiff cd1.
            for (; itRe1 != cd1.derivReEnd(); itRe1++, itIm1++) {

                ASKAPDEBUGASSERT(itIm1 != cd1.derivImEnd());

                const casacore::DComplex rowParDerivRe1 = itRe1->second;
                const casacore::DComplex rowParDerivIm1 = itIm1->second;

                // TODO: Perhaps after we inverted the loops to iterate only over derivatives of the current complex diff,
                // we dont need this check anymore, as there wont be zero derivatives after the first major iteration?
                if (rowParDerivRe1 == czero && rowParDerivIm1 == czero) {
                    continue;
                }

                // Partial normal matrix buffer for this row.
                MapOfMatrices normalMatrix;

                // Parameter name corresponding to the row of the normal matrix.
                std::string rowName = itRe1->first;

                for (casacore::uInt p2 = 0; p2 < nDataPoints; ++p2) {

                    const casacore::DComplex modelProduct = modelProductMatrix[p2 + nDataPoints * p1];
                    if (modelProduct == czero) {
                        continue;
                    }

                    const ComplexDiff& cd2 = cdm(p, p2 + columnOffset);

                    auto itRe2 = cd2.derivReBegin();
                    auto itIm2 = cd2.derivImBegin();

                    // Loop over the list of derivatives corresponding to the ComplexDiff cd2.
                    for (; itRe2 != cd2.derivReEnd(); itRe2++, itIm2++) {

                        ASKAPDEBUGASSERT(itIm2 != cd2.derivImEnd());

                        const casacore::DComplex colParDerivRe2 = itRe2->second;
                        const casacore::DComplex colParDerivIm2 = itIm2->second;

                        if (colParDerivRe2 == czero && colParDerivIm2 == czero) {
                            continue;
                        }

                        const casacore::DComplex colParDerivRe2_modelProduct = colParDerivRe2 * modelProduct;
                        const casacore::DComplex colParDerivIm2_modelProduct = colParDerivIm2 * modelProduct;

                        nmElementBuf(0, 0) = real(conj(rowParDerivRe1) * colParDerivRe2_modelProduct);
                        nmElementBuf(0, 1) = real(conj(rowParDerivRe1) * colParDerivIm2_modelProduct);
                        nmElementBuf(1, 0) = real(conj(rowParDerivIm1) * colParDerivRe2_modelProduct);
                        nmElementBuf(1, 1) = real(conj(rowParDerivIm1) * colParDerivIm2_modelProduct);

                        std::string colName = itRe2->first;

                        if (normalMatrix.find(colName) == normalMatrix.end()) {
                            normalMatrix[colName] = nmElementBuf;
                        } else {
                            normalMatrix[colName] += nmElementBuf;
                        }
                    }
                }

                // Now add this partial row to the normal matrix.
                addParameterSparselyToMatrix(rowName, normalMatrix);

                //-----------------------------------------------------
                // TODO: Added temporarily for tests.
                size_t parIndex = getParameterIndexByName(rowName);
                std::string parBaseNameMapped = getParameterBaseNameByIndex(parIndex);
                ASKAPCHECK(parBaseNameMapped == CalParamNameHelper::extractBaseParamName(rowName), "Wrong parameter mapping for: " << rowName);
            }
        }
    }
}

/// @brief Add a design matrix to the normal equations
/// @details This method computes the contribution to the normal matrix 
/// using a given design matrix and adds it.
/// @param[in] dm Design matrix to use
void GenericNormalEquations::add(const DesignMatrix& dm)
{
  std::set<std::string> names=dm.parameterNames();
  const casacore::uInt nDataSet=dm.residual().size();
  if (!nDataSet) {
      return; // nothing to process
  }

  // Loop over all parameters defined by the design matrix.
  // It may be better to write an iterator over parameters defined in
  // the design matrix instead of building a set or list. 
  for (std::set<std::string>::const_iterator iterRow = names.begin(); 
       iterRow != names.end(); ++iterRow) {
       const DMAMatrix &derivMatrices = dm.derivative(*iterRow);
       DMAMatrix::const_iterator derivMatricesIt = derivMatrices.begin();
       ASKAPDEBUGASSERT(derivMatricesIt != derivMatrices.end());
       DMBVector::const_iterator residualIt = dm.residual().begin();
       ASKAPDEBUGASSERT(residualIt != dm.residual().end());
       ASKAPDEBUGASSERT(derivMatricesIt->ncolumn());
       
       casacore::Vector<double> dataVector; // data vector buffer for this row 
       
       // it looks unnecessary from the first glance to fill the map
       // of matrices for the whole row. However, the design matrix can
       // be defined for a subset of parameters used by this normal equation
       // class. Therefore, one must resize appropriate elements of 
       // itsNormalMatrix to have there zero matrix of appropriate shape.
       // it requires access to the size of the result anyway, therefore
       // it is not too bad to calculate all elements in the row before
       // merging them with itsNormalMatrix
       MapOfMatrices normalMatrix; // normal matrix buffer for this row

       // the first contribution
       dataVector = dvElement(*derivMatricesIt, *residualIt);

       for (std::set<std::string>::const_iterator iterCol = names.begin();
                iterCol != names.end(); ++iterCol) {

           const casacore::Matrix<double> elem = nmElement(*derivMatricesIt, extractDerivatives(dm, *iterCol, 0));
           normalMatrix.insert(std::make_pair(*iterCol, elem));
       }

       // now add up all other data points
       ++derivMatricesIt;
       for(casacore::uInt dataPoint = 1; derivMatricesIt != derivMatrices.end() ;
                               ++dataPoint,++derivMatricesIt) {
           dataVector += dvElement(*derivMatricesIt, *residualIt);

           for (MapOfMatrices::iterator iterCol = normalMatrix.begin();
                               iterCol != normalMatrix.end(); ++iterCol) {

                iterCol->second += nmElement(*derivMatricesIt,
                           extractDerivatives(dm,iterCol->first,dataPoint));
           }
       }

       // Erase zero elements. We don't want to store them in the sparse normal matrix.
       for (MapOfMatrices::iterator iterCol = normalMatrix.begin();
            iterCol != normalMatrix.end();) {
           if (allMatrixElementsAreZeros(iterCol->second)) {
               normalMatrix.erase(iterCol++);
           } else {
               iterCol++;
           }
       }
       addParameterSparsely(*iterRow, normalMatrix, dataVector);
  }
}

bool GenericNormalEquations::allMatrixElementsAreZeros(const casacore::Matrix<double>& matrix)
{
    for (size_t row = 0; row < matrix.nrow(); ++row) {
         for (size_t col = 0; col < matrix.ncolumn(); ++col) {
              if (matrix(row, col) != 0.) {
                  return false;
              }
         }
    }
    return true;
}

/// @brief Extract derivatives from design matrix
/// @details This method extracts an appropriate derivative matrix
/// from the given design matrix. Effectively, it implements
/// dm.derivative(par)[dataPoint] with some additional validity checks
/// @param[in] dm Design matrix to work with
/// @param[in] par parameter name of interest
/// @param[in] dataPoint a sequence number of the data point, for which 
/// the derivatives are returned
/// @return matrix of derivatives
const casacore::Matrix<double>& 
     GenericNormalEquations::extractDerivatives(const DesignMatrix &dm,
             const std::string &par, casacore::uInt dataPoint)
{
  const DMAMatrix &derivMatrices = dm.derivative(par);
  // there is no benefit here from introducing an iterator as
  // only one specific offset is always taken
  ASKAPDEBUGASSERT(dataPoint < derivMatrices.size());
  return derivMatrices[dataPoint];
}  
  
/// @brief Calculate an element of A^tA
/// @details Each element of a sparse normal matrix is also a matrix
/// in general. However, due to some limitations of CASA operators, a
/// separate treatment is required for degenerate cases. This method
/// calculates an element of the normal matrix (effectively an element of
/// a product of A transposed and A, where A is the whole design matrix)
/// @param[in] matrix1 the first element of a sparse normal matrix
/// @param[in] matrix2 the second element of a sparse normal matrix
/// @return a product of matrix1 transposed to matrix2
casacore::Matrix<double> GenericNormalEquations::nmElement(const casacore::Matrix<double> &matrix1,
               const casacore::Matrix<double> &matrix2)
{
  ASKAPDEBUGASSERT(matrix1.ncolumn() && matrix2.ncolumn());
  ASKAPDEBUGASSERT(matrix1.nrow() == matrix2.nrow());
  if (matrix1.ncolumn() == 1 && matrix2.ncolumn() == 1) {
      const casacore::Vector<double> &m1ColVec = matrix1.column(0);
      const casacore::Vector<double> &m2ColVec = matrix2.column(0);
      return casacore::Matrix<double>(1,1,sum(m1ColVec*m2ColVec));
  }
  
  // at least one of the matrices is non-degenerate
  return product(transpose(matrix1),matrix2);
}               
  
/// @brief Calculate an element of A^tB
/// @details Each element of a sparse normal matrix is also a matrix
/// in general. However, due to some limitations of CASA operators, a
/// separate treatment is required for degenerate cases. This method
/// calculates an element of the right-hand side of the normal equation
/// (effectively an element of a product of A transposed and the data
/// vector, where A is the whole design matrix)
/// @param[in] dm an element of the design matrix
/// @param[in] dv an element of the data vector
casacore::Vector<double> GenericNormalEquations::dvElement(const casacore::Matrix<double> &dm,
              const casacore::Vector<double> &dv)
{
  ASKAPDEBUGASSERT(dm.ncolumn() && dv.nelements());
  ASKAPDEBUGASSERT(dm.nrow() == dv.nelements());
  if (dm.ncolumn() == 1) {
      const casacore::Vector<double> &dmColVec = dm.column(0);
      return casacore::Vector<double>(1, sum(dmColVec*dv));
  }
  // dm is non-degenerate
  return product(transpose(dm), dv);
}               
  
/// @brief add normal matrix for a given parameter
/// @details This means that the cross terms between parameters 
/// are excluded. However the terms inside a parameter are retained.
/// @param[in] name Name of the parameter
/// @param[in] normalmatrix Normal Matrix for this parameter
/// @param[in] datavector Data vector for this parameter
void GenericNormalEquations::add(const string& name, 
                               const casacore::Matrix<double>& normalmatrix,
                               const casacore::Vector<double>& datavector)
{
  MapOfMatrices tempSparseMatrix;
  tempSparseMatrix[name] = normalmatrix;
  
  addParameter(name, tempSparseMatrix, datavector);
}  
  
/// @brief normal equations for given parameters
/// @details In the current framework, parameters are essentially 
/// vectors, not scalars. Each element of such vector is treated
/// independently (but only vector as a whole can be fixed). As a 
/// result the element of the normal matrix is another matrix for
/// all non-scalar parameters. For scalar parameters each such
/// matrix has a shape of [1,1].
/// @param[in] par1 the name of the first parameter
/// @param[in] par2 the name of the second parameter
const casacore::Matrix<double>& GenericNormalEquations::normalMatrix(const std::string &par1, 
                          const std::string &par2) const
{
  auto cIt1 = itsNormalMatrix.find(par1);
  ASKAPCHECK(cIt1 != itsNormalMatrix.end(), "Missing first parameter "<<par1<<" is requested from the normal matrix");
  auto cIt2 = cIt1->second.find(par2);
  //ASKAPCHECK(cIt2 != cIt1->second.end(), "Missing second parameter "<<par2<<" is requested from the normal matrix");
  if (cIt2 == cIt1->second.end()) {
  // Added this to allow for sparse matrix.
      return emptyMatrix;
  }
  return cIt2->second;
}

std::map<std::string, casacore::Matrix<double> >::const_iterator GenericNormalEquations::getNormalMatrixRowBegin(const std::string &par) const
{
    auto cIt = itsNormalMatrix.find(par);
    ASKAPCHECK(cIt != itsNormalMatrix.end(), "Missing parameter " << par << " is requested from the normal matrix");

    return cIt->second.begin();
}

std::map<std::string, casacore::Matrix<double> >::const_iterator GenericNormalEquations::getNormalMatrixRowEnd(const std::string &par) const
{
    auto cIt = itsNormalMatrix.find(par);
    ASKAPCHECK(cIt != itsNormalMatrix.end(), "Missing parameter " << par << " is requested from the normal matrix");

    return cIt->second.end();
}

/// @brief data vector for a given parameter
/// @details In the current framework, parameters are essentially 
/// vectors, not scalars. Each element of such vector is treated
/// independently (but only vector as a whole can be fixed). As a 
/// result any element of the normal matrix as well as an element of the
/// data vector are, in general, matrices, not scalar. For the scalar 
/// parameter each element of data vector is a vector of unit length.
/// @param[in] par the name of the parameter of interest
const casacore::Vector<double>& GenericNormalEquations::dataVector(const std::string &par) const
{
  auto cIt = itsDataVector.find(par);
  ASKAPCHECK(cIt != itsDataVector.end(),"Parameter "<<par<<" is not found in the normal equations");
  return cIt->second;                                  
}                          
  
/// @brief write the object to a blob stream
/// @param[in] os the output stream
void GenericNormalEquations::writeToBlob(LOFAR::BlobOStream& os) const
{ 
  std::map<std::string, casacore::Vector<double> > dataVectorTemp;
  dataVectorTemp.insert(itsDataVector.begin(), itsDataVector.end());

  // increment version number on the next line and in the next method
  // if any new data members are added  
  os.putStart("GenericNormalEquations",2);
  os << itsNormalMatrix << dataVectorTemp << itsMetadata;
  os.putEnd();
}

/// @brief read the object from a blob stream
/// @param[in] is the input stream
/// @note Not sure whether the parameter should be made const or not 
void GenericNormalEquations::readFromBlob(LOFAR::BlobIStream& is)
{ 
  const int version = is.getStart("GenericNormalEquations");
  ASKAPCHECK(version == 2, 
              "Attempting to read from a blob stream an object of the wrong "
              "version: expect version 2, found version "<<version);

  std::map<std::string, casacore::Vector<double> > dataVectorTemp;

  is >> itsNormalMatrix >> dataVectorTemp >> itsMetadata;
  is.getEnd();

  itsDataVector.insert(dataVectorTemp.begin(), dataVectorTemp.end());
}

/// @brief obtain all parameters dealt with by these normal equations
/// @details Normal equations provide constraints for a number of 
/// parameters (i.e. unknowns of these equations). This method returns
/// a vector with the string names of all parameters mentioned in the
/// normal equations represented by the given object.
/// @return a vector listing the names of all parameters (unknowns of these equations)
/// @note if ASKAP_DEBUG is set some extra checks on consistency of these 
/// equations are done
std::vector<std::string> GenericNormalEquations::unknowns() const
{
  std::vector<std::string> result;
  result.reserve(itsNormalMatrix.size());
  for (auto ci=itsNormalMatrix.begin();
       ci!=itsNormalMatrix.end(); ++ci) {
       result.push_back(ci->first);
// extra consistency checks in the debug mode
#ifdef ASKAP_DEBUG
       ASKAPCHECK(itsDataVector.find(ci->first) != itsDataVector.end(), 
                  "The parameter "<<ci->first<<" is present in the normal matrix but missing in the data vector");
#endif // #ifdef ASKAP_DEBUG
  }
  return result;
} // unknowns method

bool GenericNormalEquations::addParameterNameToIndexMap(const std::string &parName)
{
    if (itsParameterNameToIndex.find(parName) != itsParameterNameToIndex.end()) {
        // Parameter already exists.
        return false;
    } else {
        size_t parIndex = itsParameterIndexToBaseName.size();

        itsParameterNameToIndex[parName] = parIndex;

        std::string baseParName = CalParamNameHelper::extractBaseParamName(parName);
        if (baseParName != parName) {
            itsParameterNameToIndex[baseParName] = parIndex;
        }
        itsParameterIndexToBaseName.push_back(baseParName);

        return true;
    }
}

ssize_t GenericNormalEquations::getParameterIndexByName(const std::string &parName) const
{
    std::map<std::string, size_t>::const_iterator it = itsParameterNameToIndex.find(parName);
    if (it != itsParameterNameToIndex.end()) {
        return it->second;
    } else {
        // Parameter does not exist.
        return -1;
    }
}

std::string GenericNormalEquations::getParameterBaseNameByIndex(size_t parIndex) const
{
    if (parIndex < itsParameterIndexToBaseName.size()) {
        return itsParameterIndexToBaseName[parIndex];
    } else {
        ASKAPCHECK(false, "The parameter index is out of range: " << parIndex);
    }
}

}}

