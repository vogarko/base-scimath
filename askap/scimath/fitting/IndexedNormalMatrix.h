/// @file
/// @brief Normal equations without any approximation
/// @details There are two kinds of normal equations currently supported. The
/// first one is a generic case, where the full normal matrix is retained. It
/// is used for calibration. The second one is intended for imaging, where we
/// can't afford to keep the whole normal matrix. In the latter approach, the 
/// matrix is approximated by a sum of diagonal and shift invariant matrices. 
/// This class represents the generic case, where no approximation to the normal
/// matrix is done.
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

#ifndef INDEXED_NORMAL_MATRIX_H
#define INDEXED_NORMAL_MATRIX_H

// casa includes
#include <casacore/casa/Arrays/Matrix.h>

// own includes
#include <askap/askap/AskapError.h>

// std includes
#include <vector>
#include <string>

namespace askap { namespace scimath {

struct IndexedNormalMatrix
{
public:
    /// @brief default constructor
    IndexedNormalMatrix();

    /// @brief assignment operator
    IndexedNormalMatrix& operator=(const IndexedNormalMatrix &src);

    // Allocate memory for matrix elements, and set default value to zero.
    void initialize(size_t nChannelsLocal_, size_t nBaseParameters_, size_t chanOffset_);

    // Resets the object to its initial state.
    void reset();

    const casacore::Matrix<double>& getValue(size_t col, size_t row, size_t chan) const;

    void addValue(size_t col, size_t row, size_t chan, const casacore::Matrix<double>& value);

    inline size_t get1Dindex(size_t col, size_t row, size_t chan) const
    {
        ASKAPDEBUGASSERT(col < nBaseParameters);
        ASKAPDEBUGASSERT(row < nBaseParameters);
        ASKAPDEBUGASSERT(chan < nChannelsLocal);
        return col + nBaseParameters * row + nBaseParameters * nBaseParameters * chan;
    }

    inline bool initialized() const
    {
        return isInitialized;
    }

    inline size_t getChanOffset() const
    {
        return chanOffset;
    }

private:
    // Flag for whether the matrix is initialized.
    bool isInitialized;
    // Number of channels at the current worker.
    size_t nChannelsLocal;
    // Number of parameters at one channel number.
    size_t nBaseParameters;
    // Channel offset (store it here for convenience).
    size_t chanOffset;
    // Matrix elements.
    std::vector<casacore::Matrix<casacore::Double>> elements;
};

}}

#endif // #ifndef INDEXED_NORMAL_MATRIX_H
