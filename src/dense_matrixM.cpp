#include "dense_matrixM.hpp"

#include "Typedefs_conf.hpp"

// C++ Includes
#include <cmath> // for sqrt

// Local Includes

#include "dense_vectorM.hpp"


#ifdef MATRIX_HAVE_PETSC    //?
#include "petsc_macroM.h"

EXTERN_C_FOR_PETSC_BEGIN
#include <petscblaslapack.h>
EXTERN_C_FOR_PETSC_END
#endif

// ------------------------------------------------------------
// Dense Matrix member functions




void DenseMatrixM::left_multiply (const DenseMatrixBaseM& M2)
{
    if (this->use_blas)
        this->_multiply_blas(M2, LEFT_MULTIPLY);
    else
    {  // (*this) <- M2 * (*this) Where: (*this) = (m x n), M2      = (m x p), M3      = (p x n)
        // M3 is a copy of *this before it gets resize()d
        DenseMatrixM M3(*this);
        // Resize *this so that the result can fit
        this->resize (M2.m(), M3.n());
        // Call the multiply function in the base class
        this->multiply(*this, M2, M3);
    }
}

// =================================================
void DenseMatrixM::left_multiply_transpose(const DenseMatrixM& A)
{
    if (this->use_blas)  this->_multiply_blas(A, LEFT_MULTIPLY_TRANSPOSE);
    else    {
        //Check to see if we are doing (A^T)*A
        if (this == &A)  {
            //libmesh_here();
            DenseMatrixM B(*this);
            // Simple but inefficient way
            // return this->left_multiply_transpose(B);

            // More efficient, but more code way
            // If A is mxn, the result will be a square matrix of Size n x n.
            const unsigned int m = A.m();
            const unsigned int n = A.n();

            // resize() *this and also zero out all entries.
            this->resize(n,n);

            // Compute the lower-triangular part
            for (unsigned int i=0; i<n; ++i)
                for (unsigned int j=0; j<=i; ++j)
                    for (unsigned int k=0; k<m; ++k) // inner products are over m
                        (*this)(i,j) += B(k,i)*B(k,j);

            // Copy lower-triangular part into upper-triangular part
            for (unsigned int i=0; i<n; ++i)
                for (unsigned int j=i+1; j<n; ++j)
                    (*this)(i,j) = (*this)(j,i);
        }

        else
        {
            DenseMatrixM B(*this);
            this->resize (A.n(), B.n());

            assert (A.m() == B.m());
            assert (this->m() == A.n());
            assert (this->n() == B.n());

            const unsigned int m_s = A.n();
            const unsigned int p_s = A.m();
            const unsigned int n_s = this->n();

            // Do it this way because there is a
            // decent chance (at least for constraint matrices)
            // that A.transpose(i,k) = 0.
            for (unsigned int i=0; i<m_s; i++)
                for (unsigned int k=0; k<p_s; k++)
                    if (A.transpose(i,k) != 0.)
                        for (unsigned int j=0; j<n_s; j++)
                            (*this)(i,j) += A.transpose(i,k)*B(k,j);
        }
    }

}






void DenseMatrixM::right_multiply (const DenseMatrixBaseM& M3)
{
    if (this->use_blas)
        this->_multiply_blas(M3, RIGHT_MULTIPLY);
    else
    {
        // (*this) <- M3 * (*this)
        // Where:
        // (*this) = (m x n),
        // M2      = (m x p),
        // M3      = (p x n)

        // M2 is a copy of *this before it gets resize()d
        DenseMatrixM M2(*this);

        // Resize *this so that the result can fit
        this->resize (M2.m(), M3.n());

        this->multiply(*this, M2, M3);
    }
}





void DenseMatrixM::right_multiply_transpose (const DenseMatrixM& B)
{
    if (this->use_blas)
        this->_multiply_blas(B, RIGHT_MULTIPLY_TRANSPOSE);
    else
    {
        //Check to see if we are doing B*(B^T)
        if (this == &B)
        {
            //libmesh_here();
            DenseMatrixM A(*this);

            // Simple but inefficient way
            // return this->right_multiply_transpose(A);

            // More efficient, more code
            // If B is mxn, the result will be a square matrix of Size m x m.
            const unsigned int m = B.m();
            const unsigned int n = B.n();

            // resize() *this and also zero out all entries.
            this->resize(m,m);

            // Compute the lower-triangular part
            for (unsigned int i=0; i<m; ++i)
                for (unsigned int j=0; j<=i; ++j)
                    for (unsigned int k=0; k<n; ++k) // inner products are over n
                        (*this)(i,j) += A(i,k)*A(j,k);

            // Copy lower-triangular part into upper-triangular part
            for (unsigned int i=0; i<m; ++i)
                for (unsigned int j=i+1; j<m; ++j)
                    (*this)(i,j) = (*this)(j,i);
        }

        else
        {
            DenseMatrixM A(*this);

            this->resize (A.m(), B.m());

            assert (A.n() == B.n());
            assert (this->m() == A.m());
            assert (this->n() == B.m());

            const unsigned int m_s = A.m();
            const unsigned int p_s = A.n();
            const unsigned int n_s = this->n();

            // Do it this way because there is a
            // decent chance (at least for constraint matrices)
            // that B.transpose(k,j) = 0.
            for (unsigned int j=0; j<n_s; j++)
                for (unsigned int k=0; k<p_s; k++)
                    if (B.transpose(k,j) != 0.)
                        for (unsigned int i=0; i<m_s; i++)
                            (*this)(i,j) += A(i,k)*B.transpose(k,j);
        }
    }
}


void DenseMatrixM::vector_mult (DenseVectorM& dest,
                                const DenseVectorM& arg) const
{
    const unsigned int n_rows = this->m();
    const unsigned int n_cols = this->n();

    // Make sure the sizes are compatible
    assert(n_cols == arg.size());
    assert(n_rows == dest.size());

    dest.zero();
    DenseMatrixM A(*this);

    for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
            dest(i) += A(i,j)*arg(j);
}


void DenseMatrixM::vector_mult_add (DenseVectorM& dest,
                                    const Real factor,
                                    const DenseVectorM& arg) const
{
    DenseVectorM temp(arg.size());
    this->vector_mult(temp, arg);
    dest.add(factor, temp);
}

// =============================================
void DenseMatrixM::get_principal_submatrix (unsigned int sub_m,
        unsigned int sub_n,
        DenseMatrixM& dest) const
{
    assert( (sub_m <= this->m()) && (sub_n <= this->n()) );
    dest.resize(sub_m, sub_n);
    for (unsigned int i=0; i<sub_m; i++) for (unsigned int j=0; j<sub_n; j++)  dest(i,j) = (*this)(i,j);
}

// ================================================
void DenseMatrixM::get_principal_submatrix (unsigned int sub_m, DenseMatrixM& dest) const
{
    get_principal_submatrix(sub_m, sub_m, dest);
}


void DenseMatrixM::get_transpose (DenseMatrixM& dest) const
{
    dest.resize(this->n(), this->m());
    for (unsigned int i=0; i<dest.m(); i++) for (unsigned int j=0; j<dest.n(); j++) dest(i,j) = (*this)(j,i);
}


// ===================================================
void DenseMatrixM::lu_solve (DenseVectorM& b,
                             DenseVectorM& x,
                             const bool partial_pivot)
{
    // Check for a previous decomposition
    switch (this->_decomposition_type)    {
    case NONE:
        this->_lu_decompose (partial_pivot);
        break;
    case LU:
        // Already factored, just need to call back_substitute.
        break;
    default:
        std::cerr << "Error! This matrix already has a "
                  << "different decomposition..."
                  << std::endl;
        exit(0);
    }
    // Perform back substitution
    this->_lu_back_substitute (b, x, partial_pivot);
}




void DenseMatrixM::_lu_back_substitute (DenseVectorM& b,
                                        DenseVectorM& x,
                                        const bool ) const
{
    const unsigned int
    n = this->n();
    assert (this->m() == n);
    assert (this->m() == b.size());

    x.resize (n);
    // A convenient reference to *this
    const DenseMatrixM& A = *this;

    // Transform the RHS vector
    for (unsigned int i=0; i<(n-1); i++) {
        // Get the diagonal entry and take its inverse
        const Real diag = A(i,i);
        assert (fabs(diag) < 1.e-20);
        const Real diag_inv = 1./diag;
        // Get the entry b(i) and store it
        const Real bi = b(i);
        for (unsigned int j=(i+1); j<n; j++) b(j) -= bi*A(j,i)*diag_inv;
    }


    // Perform back-substitution
    {
        x(n-1) = b(n-1)/A(n-1,n-1);
        for (unsigned int i=0; i<=(n-1); i++)  {
            const unsigned int ib = (n-1)-i;
            // Get the diagonal and take its inverse
            const Real diag = A(ib,ib);
            assert (fabs(diag) < 1.e-20);
            const Real diag_inv = 1./diag;

            for (unsigned int j=(ib+1); j<n; j++)  {
                b(ib) -= A(ib,j)*x(j);
		x(ib)  = b(ib)*diag_inv;
            }
            
        }
    }
    return;
}


// ===========================================================
void DenseMatrixM::_lu_decompose (const bool partial_pivot)
{    // If this function was called, there better not be any
    // previous decomposition of the matrix.
    assert(this->_decomposition_type == NONE);
    // Get the matrix size and make sure it is square
    const unsigned int
    m = this->m();   assert (m == this->n());
    // A convenient reference to *this
    DenseMatrixM& A = *this;

    // Straight, vanilla LU factorization without pivoting
    if (!partial_pivot){
        // For each row in the matrix
        for (unsigned int i=0; i<m; i++)  {
            // Get the diagonal entry and take its inverse
            const Real diag = A(i,i);
            assert (fabs(diag) < 1.e-20);
            const Real diag_inv = 1./diag;

            // For each row in the submatrix
            for (unsigned int j=i+1; j<m; j++)      {
                // Get the scale factor for this row
                const Real fact = A(j,i)*diag_inv;
                // For each column in the subrow scale it by the factor
                for (unsigned int k=i+1; k<m; k++)  A(j,k) -= fact*A(i,k);
            }
        }
    }
    else
    {   // Do partial pivoting.
        exit(0);  
    }
    // Set the flag for LU decomposition
    this->_decomposition_type = LU;
}

// ========================================
Real DenseMatrixM::det ()
{
    // First LU decompose the matrix (without partial pivoting).
    // Note that the lu_decompose routine will check to see if the
    // matrix is square so we don't worry about it.
    if (this->_decomposition_type == NONE)
        this->_lu_decompose(false);
    else if (this->_decomposition_type != LU)
    {
        std::cerr << "Error! Can't compute the determinant under "
                  << "the current decomposition."
                  << std::endl;
        exit(0);
    }

    // A variable to keep track of the running product of diagonal terms.
    Real determinant = 1.;
    // Loop over diagonal terms, computing the product. 
    for (unsigned int i=0; i<this->m(); i++) determinant *= (*this)(i,i);
    // Return the determinant
    return determinant;
}



// The cholesky solve function first decomposes the matrix
// with cholesky_decompose and then uses the cholesky_back_substitute
// routine to find the solution x.

void DenseMatrixM::cholesky_solve (DenseVectorM& b,
                                   DenseVectorM& x)
{
    // Check for a previous decomposition
    switch (this->_decomposition_type)
    {
    case NONE:
    {
        this->_cholesky_decompose ();
        break;
    }

    case CHOLESKY:
    {
        // Already factored, just need to call back_substitute.
        break;
    }

    default:
    {
        std::cerr << "Error! This matrix already has a "
                  << "different decomposition..."
                  << std::endl;
        exit(0);
    }
    }

    // Perform back substitution
    this->_cholesky_back_substitute (b, x);
}




// This algorithm is based on the Cholesky decomposition in
// the Numerical Recipes in C book.

void DenseMatrixM::_cholesky_decompose ()
{
    // If we called this function, there better not be any
    // previous decomposition of the matrix.
    assert(this->_decomposition_type == NONE);

    // Shorthand notation for number of rows and columns.
    const unsigned int
    m = this->m(),
        n = this->n();

    // Just to be really sure...
    assert(m==n);

    // A convenient reference to *this
    DenseMatrixM& A = *this;

    for (unsigned int i=0; i<m; ++i)
    {
        for (unsigned int j=i; j<n; ++j)
        {
            for (unsigned int k=0; k<i; ++k)
                A(i,j) -= A(i,k) * A(j,k);

            if (i == j)
            {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
                if (A(i,j) <= 0.0)
                {
                    std::cerr << "Error! Can only use Cholesky decomposition "
                              << "with symmetric positive definite matrices."
                              << std::endl;
                    exit(0);
                }
#endif

                A(i,i) = std::sqrt(A(i,j));
            }
            else
                A(j,i) = A(i,j) / A(i,i);
        }
    }

    // Set the flag for CHOLESKY decomposition
    this->_decomposition_type = CHOLESKY;
}



void DenseMatrixM::_cholesky_back_substitute (DenseVectorM& b,
        DenseVectorM& x) const
{
    // Shorthand notation for number of rows and columns.
    const unsigned int
    m = this->m(),
        n = this->n();

    // Just to be really sure...
    assert(m==n);

    // A convenient reference to *this
    const DenseMatrixM& A = *this;

    // Now compute the solution to Ax =b using the factorization.
    x.resize(m);

    // Solve for Ly=b
    for (unsigned int i=0; i<n; ++i)
    {
        for (unsigned int k=0; k<i; ++k)
            b(i) -= A(i,k)*x(k);

        x(i) = b(i) / A(i,i);
    }

    // Solve for L^T x = y
    for (unsigned int i=0; i<n; ++i)
    {
        const unsigned int ib = (n-1)-i;

        for (unsigned int k=(ib+1); k<n; ++k)
            x(ib) -= A(k,ib) * x(k);

        x(ib) /= A(ib,ib);
    }
}





#ifdef MATRIX_HAVE_PETSC


void DenseMatrixM::_multiply_blas(const DenseMatrixBaseM& other,
                                  _BLAS_Multiply_Flag flag)
{
    int result_size = 0;

    // For each case, determine the size of the final result make sure
    // that the inner dimensions match
    switch (flag)
    {
    case LEFT_MULTIPLY:
    {
        result_size = other.m() * this->n();
        if (other.n() == this->m())
            break;
    }
    case RIGHT_MULTIPLY:
    {
        result_size = other.n() * this->m();
        if (other.m() == this->n())
            break;
    }
    case LEFT_MULTIPLY_TRANSPOSE:
    {
        result_size = other.n() * this->n();
        if (other.m() == this->m())
            break;
    }
    case RIGHT_MULTIPLY_TRANSPOSE:
    {
        result_size = other.m() * this->m();
        if (other.n() == this->n())
            break;
    }
    default:
    {
        std::cout << "Unknown flag selected or matrices are ";
        std::cout << "incompatible for multiplication." << std::endl;
        exit(0);
    }
    }

    // For this to work, the passed arg. must actually be a DenseMatrixM
    const DenseMatrixM* const_that = dynamic_cast< const DenseMatrixM* >(&other);
    if (!const_that)
    {
        std::cerr << "Unable to cast input matrix to usable type." << std::endl;
        exit(0);
    }

    // Also, although 'that' is logically const in this BLAS routine,
    // the PETSc BLAS interface does not specify that any of the inputs are
    // const.  To use it, I must cast away const-ness.
    DenseMatrixM* that = const_cast< DenseMatrixM* > (const_that);

    // Initialize A, B pointers for LEFT_MULTIPLY* cases
    DenseMatrixM
    *A = this,
         *B = that;

    // For RIGHT_MULTIPLY* cases, swap the meaning of A and B.
    // Here is a full table of combinations we can pass to BLASgemm, and what the answer is when finished:
    // pass A B   -> (Fortran) -> A^T B^T -> (C++) -> (A^T B^T)^T -> (identity) -> B A   "lt multiply"
    // pass B A   -> (Fortran) -> B^T A^T -> (C++) -> (B^T A^T)^T -> (identity) -> A B   "rt multiply"
    // pass A B^T -> (Fortran) -> A^T B   -> (C++) -> (A^T B)^T   -> (identity) -> B^T A "lt multiply t"
    // pass B^T A -> (Fortran) -> B A^T   -> (C++) -> (B A^T)^T   -> (identity) -> A B^T "rt multiply t"
    if (flag==RIGHT_MULTIPLY || flag==RIGHT_MULTIPLY_TRANSPOSE)
        std::swap(A,B);

    // transa, transb values to pass to blas
    char
    transa[] = "n",
               transb[] = "n";

    // Integer values to pass to BLAS:
    //
    // M
    // In Fortran, the number of rows of op(A),
    // In the BLAS documentation, typically known as 'M'.
    //
    // In C/C++, we set:
    // M = n_cols(A) if (transa='n')
    //     n_rows(A) if (transa='t')
    int M = static_cast<int>( A->n() );

    // N
    // In Fortran, the number of cols of op(B), and also the number of cols of C.
    // In the BLAS documentation, typically known as 'N'.
    //
    // In C/C++, we set: N = n_rows(B) if (transb='n')     n_cols(B) if (transb='t')
    int N = static_cast<int>( B->m() );

    // K
    // In Fortran, the number of cols of op(A), and also
    // the number of rows of op(B). In the BLAS documentation,
    // typically known as 'K'.
    //
    // In C/C++, we set: K = n_rows(A) if (transa='n')     n_cols(A) if (transa='t')
    int K = static_cast<int>( A->m() );
    // LDA (leading dimension of A). In our cases LDA is always the number of columns of A.
    int LDA = static_cast<int>( A->n() );
    // LDB (leading dimension of B).  In our cases, LDB is always the number of columns of B.
    int LDB = static_cast<int>( B->n() );

    if (flag == LEFT_MULTIPLY_TRANSPOSE)    {
        transb[0] = 't';  N = static_cast<int>( B->n() );
    }

    else if (flag == RIGHT_MULTIPLY_TRANSPOSE)
    {
        transa[0] = 't';  std::swap(M,K);
    }

    // LDC (leading dimension of C).  LDC is the number of columns in the solution matrix.
    int LDC = M;

    // Scalar values to pass to BLAS    //
    // scalar multiplying the whole product AB
    Real alpha = 1.;
    // scalar multiplying C, which is the original matrix.
    Real beta  = 0.;
    // Storage for the result
    std::vector  <Real> result (result_size);
    // Finally ready to call the BLAS
    BLASgemm_(transa, transb, &M, &N, &K, &alpha, &(A->_val[0]), &LDA, &(B->_val[0]), &LDB, &beta, &result[0], &LDC);

    // Update the relevant dimension for this matrix.
    switch (flag)    {
    case LEFT_MULTIPLY:            
        this->_m = other.m();
        break;
    case RIGHT_MULTIPLY:        
        this->_n = other.n();
        break;
    case LEFT_MULTIPLY_TRANSPOSE:  
        this->_m = other.n();
        break;
    case RIGHT_MULTIPLY_TRANSPOSE: 
        this->_n = other.m();
        break;
    default:
        std::cout << "Unknown flag selected." << std::endl;
        exit(0);
    }

    // Swap my data vector with the result
    this->_val.swap(result);
}

#else
// ==========================================================
void DenseMatrixM::_multiply_blas(const DenseMatrixBaseM& ,
                                  _BLAS_Multiply_Flag ){
    std::cerr << "No PETSc-provided BLAS available!" << std::endl;
    exit(0);
}

#endif


// This routine is commented out since it is not really a memory
// efficient implementation.  Also, you don't *need* the inverse
// for anything, instead just use lu_solve to solve Ax=b.
//
// void DenseMatrixM::inverse ()
// {
//   // First LU decompose the matrix (without partial pivoting).
//   // Note that the lu_decompose routine will check to see if the
//   // matrix is square so we don't worry about it.
//   if (!this->_lu_decomposed)
//     this->_lu_decompose(false);

//   // A unit vector which will be used as a rhs
//   // to pick off a single value each time.
//   DenseVector  <Real> e;
//   e.resize(this->m());

//   // An empty vector which will be used to hold the solution
//   // to the back substitutions.
//   DenseVector  <Real> x;
//   x.resize(this->m());

//   // An empty dense matrix to store the resulting inverse
//   // temporarily until we can overwrite A.
//   DenseMatrixM inv;
//   inv.resize(this->m(), this->n());

//   // Resize the passed in matrix to hold the inverse
//   inv.resize(this->m(), this->n());

//   for (unsigned int j=0; j<this->n(); ++j)
//     {
//       e.zero();
//       e(j) = 1.;
//       this->_lu_back_substitute(e, x, false);
//       for (unsigned int i=0; i<this->n(); ++i)
// 	inv(i,j) = x(i);
//     }

//   // Now overwrite all the entries
//   *this = inv;
// }

