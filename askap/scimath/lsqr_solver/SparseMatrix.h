/*
 * SparseMatrix.h
 *
 * @author Vitaliy Ogarko <vogarko@gmail.com>
 */

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <vector>

namespace askap { namespace lsqr {

typedef std::vector<double> Vector;

/*
 * A class to work with sparse matrices that are stored using Compressed Sparse Row (CSR) format.
 */
class SparseMatrix {
public:

    /**
     * Constructs an empty sparse matrix
     */
    SparseMatrix();

    /*
     * Constructor: initializes the sparse matrix.
     * nl - number of matrix lines.
     * comm - MPI communicator (used when the matrix is split among CPUs).
     */
    SparseMatrix(size_t nl);
#ifdef HAVE_MPI
    SparseMatrix(size_t nl, const MPI_Comm &comm);
#endif

    virtual ~SparseMatrix();

    /*
     * Adds one element at specified position (column index).
     * It does not check if the element has already been added at this position.
     */
    void Add(double value, size_t column);

    /*
     * Adds a new row.
     */
    void NewRow();

    /*
     * Returns the number of (non-zero) elements added into the sparse matrix.
     */
    size_t GetNumberElements() const;

    /*
     * Returns the current number of rows in the matrix.
     */
    size_t GetCurrentNumberRows() const;

    /*
     * Returns the total number of rows in the matrix.
     */
    size_t GetTotalNumberRows() const;

    /*
     * Returns the number of nonempty rows.
     */
    size_t GetNumberNonemptyRows() const;

    /*
     * Returns matrix column norms.
     */
    void GetColumnNorms(Vector& columnNorms) const;

    /*
     * Scales matrix columns.
     */
    void ScaleColumns(Vector& columnWeight);

    /*
     * Normalizes matrix columns.
     */
    void NormalizeColumns(Vector& columnNorms);

    /*
     * Returns a flag whether the matrix has been finalized.
     */
    bool Finalized() const;

    /*
     * Returns the (i, j)-element's value of the matrix,
     * where i - column number, j - row number.
     * This function is mainly used for testing this class,
     * and should not be used for general purpose as it has low performance.
     */
    double GetValue(size_t i, size_t j) const;

    /*
     * Resets the matrix to the initial state, i.e., removes all elements.
     * Does not release the reserved memory, so a matrix can be reused.
     */
    void Reset();

    /*
     * Computes the product between a sparse matrix and vector x, and stores the result in b.
     */
    void MultVector(const Vector &x, Vector &b) const;

    /*
     * Computes the product between the transpose of sparse matrix and vector x.
     */
    void TransMultVector(const Vector &x, Vector &b) const;

    /*
     * Adds sparse operator to the parallel matrix (e.g. gradient, or Laplacian).
     * Note: Assumes the full matrix is split between MPI ranks by columns.
     * nDiag - number of diagonals of the sparse operator (e.g. two diagonals for forward difference, three for Laplacian).
     * nParametersLocal - number of parameters at the current MPI rank.
     * columnIndexGlobal - column index (of the nonzero value) for each diagonal.
     * matrixValue - matrix value for each diagonal (constant for all rows).
     */
    void addParallelSparseOperator(size_t nDiag,
                                   size_t nParametersLocal,
                                   const std::vector<std::vector<int> >& columnIndexGlobal,
                                   const std::vector<double>& matrixValue);

    /*
     * Extends (finalized) matrix for adding more elements.
     * Makes matrix non-finalized.
     */
    void Extend(size_t extra_nl);

    /*
     * Should be called when all elements of the matrix have been added.
     * It stores the index of last element, and validates the matrix indexes.
     */
    bool Finalize(size_t ncolumns);

    /// Copy assignment operator
    SparseMatrix &operator=(const SparseMatrix &other)
    {
#ifdef HAVE_MPI
        if (itsComm != MPI_COMM_NULL) {
            MPI_Comm_free(&itsComm);
        }
        if (other.itsComm != MPI_COMM_NULL) {
            MPI_Comm_dup(other.itsComm, &itsComm);
        }
#endif // HAVE_MPI
        finalized = other.finalized;
        nel = other.nel;
        nl = other.nl;
        nl_current = other.nl_current;
        sa = other.sa;
        ija = other.ija;
        ijl = other.ijl;
        return *this;
    }

#ifdef HAVE_MPI
    /*
     * Returns the MPI communicator. If the communicator is MPI_COMM_NULL it
     * means that this object is not using MPI, and hence it's not split across
     * many ranks.
     */
    const MPI_Comm& GetComm() const;
#endif

private:
    // Flag for whether the matrix has been finalized.
    bool finalized;

    // Actual number of nonzero elements.
    size_t nel;
    // Total number of rows in the matrix.
    size_t nl;
    // Current number of the added lines (rows) in the matrix.
    size_t nl_current;

    // An array of the (left-to-right, then top-to-bottom) non-zero values of the matrix.
    Vector sa;
    // The column indexes corresponding to the values.
    std::vector<size_t> ija;
    // The list of 'sa' indexes where each row starts.
    std::vector<size_t> ijl;

#ifdef HAVE_MPI
    // MPI communicator.
    MPI_Comm itsComm;
#endif

    /*
     * Validates the boundaries of column indexes.
     */
    bool ValidateIndexBoundaries(size_t ncolumns) const;
};

}} // namespace askap.lsqr

#endif /* SPARSEMATRIX_H_ */
