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

// own includes
#include <askap/scimath/fitting/IndexedNormalMatrix.h>

namespace askap { namespace scimath {

IndexedMatrixElelment::IndexedMatrixElelment(double a)
{
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            data[i][j] = a;
        }
    }
}

IndexedMatrixElelment& IndexedMatrixElelment::operator+=(const IndexedMatrixElelment& rhs)
{
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            this->data[i][j] += rhs.data[i][j];
        }
    }
    return *this;
}

IndexedNormalMatrix::IndexedNormalMatrix() :
    isInitialized(false),
    nChannelsLocal(0),
    nBaseParameters(0),
    chanOffset(0)
{}

IndexedNormalMatrix& IndexedNormalMatrix::operator=(const IndexedNormalMatrix &src)
{
    if (&src != this) {
        reset();
        initialize(src.nChannelsLocal, src.nBaseParameters, src.chanOffset);
        elements = src.elements;
    }
    return *this;
}

void IndexedNormalMatrix::initialize(size_t nChannelsLocal_, size_t nBaseParameters_, size_t chanOffset_)
{
    if (!initialized()) {
        nChannelsLocal = nChannelsLocal_;
        nBaseParameters = nBaseParameters_;
        chanOffset = chanOffset_;

        size_t nElements = nChannelsLocal * nBaseParameters * nBaseParameters;

        IndexedMatrixElelment zero = IndexedMatrixElelment(0.);
        elements.resize(nElements, zero);

        isInitialized = true;
    } else {
        throw AskapError("Attempt initialize an already initialized normal matrix!");
    }
}

void IndexedNormalMatrix::reset() {
    isInitialized = false;
    nChannelsLocal = 0;
    nBaseParameters = 0;
    chanOffset = 0;
    elements.clear();
}

const IndexedMatrixElelment& IndexedNormalMatrix::getValue(size_t col, size_t row, size_t chan) const
{
    if (initialized()) {
        size_t index = get1Dindex(col, row, chan);
        return elements[index];
    } else {
        throw AskapError("Attempt to get an element of non-initialized normal matrix!");
    }
}

void IndexedNormalMatrix::addValue(size_t col, size_t row, size_t chan, const IndexedMatrixElelment& value)
{
    if (initialized()) {
        size_t index = get1Dindex(col, row, chan);
        elements[index] += value;
    } else {
        throw AskapError("Attempt to set an element of non-initialized normal matrix!");
    }
}

}}

