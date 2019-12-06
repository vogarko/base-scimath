/*
 * SparseMatrix.cpp
 *
 * @author Vitaliy Ogarko <vogarko@gmail.com>
 */

#include <stdexcept>
#include <cassert>
#include <cmath>

#include <askap/scimath/lsqr_solver/SparseMatrix.h>
#include <askap/scimath/lsqr_solver/ParallelTools.h>

namespace askap { namespace lsqr {

SparseMatrix::SparseMatrix(size_t nl) :
    finalized(false),
    nel(0),
    nl(nl),
    nl_current(0),
    sa(),
    ija(),
    ijl(nl + 1)
#ifdef HAVE_MPI
    ,itsComm(MPI_COMM_NULL)
#endif
{
}

#ifdef HAVE_MPI
SparseMatrix::SparseMatrix(size_t nl, const MPI_Comm &comm) :
    finalized(false),
    nel(0),
    nl(nl),
    nl_current(0),
    sa(),
    ija(),
    ijl(nl + 1)
{
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_dup(comm, &itsComm);
}
#endif

SparseMatrix::~SparseMatrix()
{
#ifdef HAVE_MPI
    if (itsComm != MPI_COMM_NULL) {
        MPI_Comm_free(&itsComm);
    }
#endif
}

void SparseMatrix::Add(double value, size_t column)
{
    // Sanity check.
    if (finalized) {
        throw std::runtime_error("Matrix has already been finalized in SparseMatrix::Add!");
    }

    // Do not add zero values to a sparse matrix.
    if (value == 0.) return;

    // Sanity check for the number of lines.
    if (nl_current == 0) {
        throw std::runtime_error("Error in the number of lines in SparseMatrix::Add!");
    }

    nel += 1;
    sa.push_back(value);
    ija.push_back(column);
}

void SparseMatrix::NewRow()
{
    // Sanity check.
    if (finalized) {
        throw std::runtime_error("Matrix has already been finalized in SparseMatrix::NewRow!");
    }

    // Sanity check.
    if (nl_current >= nl) {
        throw std::runtime_error("Error in number of rows in SparseMatrix::NewRow!");
    }

    nl_current += 1;
    ijl[nl_current - 1] = nel;
}

size_t SparseMatrix::GetNumberElements() const
{
    return nel;
}

size_t SparseMatrix::GetCurrentNumberRows() const
{
    return nl_current;
}

size_t SparseMatrix::GetTotalNumberRows() const
{
    return nl;
}

size_t SparseMatrix::GetNumberNonemptyRows() const
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::GetNumberNonemptyRows!");
    }

    size_t number_empty_rows = 0;
    for (size_t i = 0; i < nl; i++) {
        if (ijl[i] == ijl[i + 1]) number_empty_rows++;
    }
    return (nl - number_empty_rows);
}

void SparseMatrix::GetColumnNorms(Vector& columnNorms) const
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::GetColumnNorms!");
    }

    // Set all elements to zero.
    std::fill(columnNorms.begin(), columnNorms.end(), 0.);

    for (size_t i = 0; i < nl; i++) {
        for (size_t k = ijl[i]; k < ijl[i + 1]; k++) {
            size_t j = ija[k];
            columnNorms.at(j) += sa[k] * sa[k];
        }
    }

    for (size_t j = 0; j < columnNorms.size(); j++) {
        columnNorms[j] = sqrt(columnNorms[j]);
    }
}

void SparseMatrix::ScaleColumns(Vector& columnWeight)
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::ScaleColumns!");
    }

    for (size_t i = 0; i < nl; i++) {
        for (size_t k = ijl[i]; k < ijl[i + 1]; k++) {
            sa[k] *= columnWeight.at(ija[k]);
        }
    }
}

void SparseMatrix::NormalizeColumns(Vector& columnNorms)
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::NormalizeColumns!");
    }

    GetColumnNorms(columnNorms);

    size_t nParameters = columnNorms.size();
    std::vector<double> columnWeight(nParameters);

    for (size_t i = 0; i < columnWeight.size(); i++) {
        assert(columnNorms[i] != 0.);
        columnWeight[i] = 1. / columnNorms[i];
    }

    ScaleColumns(columnWeight);
}

bool SparseMatrix::Finalized() const
{
    return finalized;
}

double SparseMatrix::GetValue(size_t i, size_t j) const
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::GetValue!");
    }

    double res = 0.0;
    for (size_t k = ijl[j]; k < ijl[j + 1]; k++) {
        if (ija[k] == i) {
        // Found a non-zero element at the column i.
            res = sa[k];
            break;
        }
    }
    return res;
}

void SparseMatrix::Reset()
{
    finalized = false;

    nl_current = 0;
    nel = 0;

    std::fill(sa.begin(), sa.end(), 0.);
    std::fill(ija.begin(), ija.end(), 0);
    std::fill(ijl.begin(), ijl.end(), 0);
}

void SparseMatrix::MultVector(const Vector& x, Vector& b) const
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::MultVector!");
    }

    // Set all elements to zero.
    std::fill(b.begin(), b.end(), 0.);

    for (size_t i = 0; i < nl; i++) {
        for (size_t k = ijl[i]; k < ijl[i + 1]; k++) {
            b[i] += sa[k] * x[ija[k]];
        }
    }
}

void SparseMatrix::TransMultVector(const Vector& x, Vector& b) const
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::TransMultVector!");
    }

    // Set all elements to zero.
    std::fill(b.begin(), b.end(), 0.);

    for (size_t i = 0; i < nl; i++) {
// Compiler directive to vectorize the following loop.
//#pragma ivdep       // For Intel compiler.
#pragma GCC ivdep   // For GCC compiler.

        for (size_t k = ijl[i]; k < ijl[i + 1]; k++) {
            size_t j = ija[k];
            b[j] += sa[k] * x[i];
        }
    }
}

void SparseMatrix::addParallelSparseOperator(size_t nDiag,
                                             size_t nParametersLocal,
                                             const std::vector<std::vector<int> >& columnIndexGlobal,
                                             const std::vector<double>& matrixValue)
{
#ifdef HAVE_MPI
    assert(itsComm != MPI_COMM_NULL);

    int myrank, nbproc;
    MPI_Comm_rank(itsComm, &myrank);
    MPI_Comm_size(itsComm, &nbproc);

    size_t nParametersTotal = ParallelTools::get_total_number_elements(nParametersLocal, nbproc, itsComm);
    size_t nParametersSmaller = ParallelTools::get_nsmaller(nParametersLocal, myrank, nbproc, itsComm);
#else
    int myrank = 0;
    size_t nParametersTotal = nParametersLocal;
    size_t nParametersSmaller = 0;
#endif

    Extend(nParametersTotal);

    for (size_t i = 0; i < nParametersTotal; i++) {
        NewRow();

        // For sanity check.
        bool allNegative = true;
        bool allNonNegative = true;

        // Loop over diagonals.
        for (size_t k = 0; k < nDiag; k++) {
            if (columnIndexGlobal[k][i] >= 0
                && columnIndexGlobal[k][i] >= nParametersSmaller
                && columnIndexGlobal[k][i] < nParametersSmaller + nParametersLocal) {

                // Local matrix column index (at the current MPI rank).
                size_t localColumnIndex = columnIndexGlobal[k][i] - nParametersSmaller;

                Add(matrixValue[k], localColumnIndex);
            }
            // For sanity check.
            if (columnIndexGlobal[k][i] >= 0) {
                allNegative = false;
            } else {
                allNonNegative = false;
            }
        }
        assert(allNegative || allNonNegative);
    }
    Finalize(nParametersLocal);
}

void SparseMatrix::Extend(size_t extra_nl)
{
    // Sanity check.
    if (!finalized) {
        throw std::runtime_error("Matrix has not been finalized yet in SparseMatrix::Extend!");
    }

    finalized = false;

    // Reset the last element index.
    ijl[nl] = 0;

    nl += extra_nl;
    ijl.resize(nl + 1);
}

bool SparseMatrix::Finalize(size_t ncolumns)
{
    // Sanity check.
    if (nl_current != nl) {
        throw std::runtime_error("Wrong total number of rows in SparseMatrix::Finalize!");
    }

    // Store index of the last element.
    ijl[nl] = nel;

    if (!ValidateIndexBoundaries(ncolumns)) {
        throw std::runtime_error("Sparse matrix validation failed!");
    }

    finalized = true;

    return true;
}

bool SparseMatrix::ValidateIndexBoundaries(size_t ncolumns)
{
    // Use the same loop as in A'x multiplication.
    for (size_t i = 0; i < nl; i++) {
        for (size_t k = ijl[i]; k < ijl[i + 1]; k++) {
            if (k > nel) {
                throw std::runtime_error("Sparse matrix validation failed for k-index!");
            }

            size_t j = ija[k];

            if (j > ncolumns - 1) {
                throw std::runtime_error("Sparse matrix validation failed for j-index!");
            }
        }
    }
    return true;
}

#ifdef HAVE_MPI
const MPI_Comm& SparseMatrix::GetComm() const
{
    return itsComm;
}
#endif

}} // namespace askap.lsqr

