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

// own includes
#include <askap/askap/AskapError.h>

// std includes
#include <complex>
#include <string>
#include <vector>

namespace askap { namespace scimath {

/// @brief Data holder for the indexed normal matrix element.
struct IndexedMatrixElelment {
    IndexedMatrixElelment(double a);
    IndexedMatrixElelment& operator+=(const IndexedMatrixElelment& rhs);

    // Complex value stored as 2x2 matrix.
    double data[2][2];
};

struct IndexedNormalMatrix
{
public:
    /// @brief default constructor
    IndexedNormalMatrix();

    /// @brief assignment operator
    IndexedNormalMatrix& operator=(const IndexedNormalMatrix &src);

    /// @brief Allocate memory for matrix elements, and set default values to zero.
    void initialize(size_t nBaseParameters_, size_t nChannelsLocal_, size_t chanOffset_);

    /// @brief Resets the object to its initial state.
    void reset();

    /// @brief Returns the matrix element by its 3D index.
    const IndexedMatrixElelment& getValue(size_t col, size_t row, size_t chan) const;

    /// @brief Increments (+=) the matrix element by a given value
    void addValue(size_t col, size_t row, size_t chan, const IndexedMatrixElelment& value);

    /// @brief Returns whether the matrix is initialized.
    inline bool initialized() const
    {
        return isInitialized;
    }

    /// @brief Returns the channel offset.
    /// @details This is needed to retrieve the true channel numbers on the workers in the parallel case.
    inline size_t getChanOffset() const
    {
        return chanOffset;
    }

private:
    /// @brief Converts the 3D index to 1D index.
    inline size_t get1Dindex(size_t col, size_t row, size_t chan) const
    {
        ASKAPDEBUGASSERT(col < nBaseParameters);
        ASKAPDEBUGASSERT(row < nBaseParameters);
        ASKAPDEBUGASSERT(chan < nChannelsLocal);
        return col + nBaseParameters * row + nBaseParameters * nBaseParameters * chan;
    }

    // Flag for whether the matrix is initialized.
    bool isInitialized;
    // Number of parameters at one channel number.
    size_t nBaseParameters;
    // Number of channels at the current worker.
    size_t nChannelsLocal;
    // Channel offset (store it here for convenience).
    size_t chanOffset;
    // Matrix elements.
    std::vector<IndexedMatrixElelment> elements;
};

struct IndexedDataVector
{
    using element_type = std::complex<double>;
public:
    /// @brief default constructor
    IndexedDataVector();

    /// @brief assignment operator
    IndexedDataVector& operator=(const IndexedDataVector &src);

    /// @brief Allocate memory for vector elements, and set default values to zero.
    void initialize(size_t nBaseParameters_, size_t nChannelsLocal_);

    /// @brief Resets the object to its initial state.
    void reset();

    /// @brief Returns the data vector element by its 2D index.
    const element_type& getValue(size_t row, size_t chan) const;

    /// @brief Increments (+=) the data vector element by a given value
    void addValue(size_t row, size_t chan, const element_type& value);

    /// @brief Returns whether the vector is initialized.
    inline bool initialized() const
    {
        return isInitialized;
    }

private:
    /// @brief Converts the 2D index to 1D index.
    inline size_t get1Dindex(size_t row, size_t chan) const
    {
        ASKAPDEBUGASSERT(row < nBaseParameters);
        ASKAPDEBUGASSERT(chan < nChannelsLocal);
        return row + nBaseParameters * chan;
    }

    // Flag for whether the matrix is initialized.
    bool isInitialized;
    // Number of parameters at one channel number.
    size_t nBaseParameters;
    // Number of channels at the current worker.
    size_t nChannelsLocal;
    // Vector elements.
    std::vector<element_type> elements;
};

}}

#endif // #ifndef INDEXED_NORMAL_MATRIX_H
