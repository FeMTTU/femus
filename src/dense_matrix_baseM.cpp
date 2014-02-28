#include "dense_matrix_baseM.hpp"

#include "Typedefs_conf.hpp"
#include "DenseVectorBase.hpp"

#include <iomanip> // for std::setw()

// ===========================================================
void DenseMatrixBaseM::multiply (DenseMatrixBaseM& M1,
                                 const DenseMatrixBaseM& M2,
                                 const DenseMatrixBaseM& M3){
  
    // Assertions to make sure we have been
    // passed matrices of the correct dimension.
    assert (M1.m() == M2.m());
    assert (M1.n() == M3.n());
    assert (M2.n() == M3.m());

    const unsigned int m_s = M2.m();
    const unsigned int p_s = M2.n();
    const unsigned int n_s = M1.n();

    // Do it this way because there is a
    // decent chance (at least for constraint matrices)
    // that M3(k,j) = 0. when right-multiplying.
    for (unsigned int k=0; k<p_s; k++)
        for (unsigned int j=0; j<n_s; j++)
            if (M3.el(k,j) != 0.)
                for (unsigned int i=0; i<m_s; i++)
                    M1.el(i,j) += M2.el(i,k) * M3.el(k,j);
}

// ============================================================
void DenseMatrixBaseM::condense(const unsigned int iv,
                                const unsigned int jv,
                                const Real val,
                                DenseVectorBase& rhs)
{
    assert (this->_m == rhs.size());assert (iv == jv);


    // move the known value into the RHS
    // and zero the column
    for (unsigned int i=0; i<this->m(); i++) {
        rhs.el(i) -= this->el(i,jv)*val;
        this->el(i,jv) = 0.;
    }

    // zero the row
    for (unsigned int j=0; j<this->n(); j++) this->el(iv,j) = 0.;
    this->el(iv,jv) = 1.;
    rhs.el(iv) = val;

}

// ============================================================
void DenseMatrixBaseM::print_scientific (std::ostream& os) const
{
    // save the initial format flags
    std::ios_base::fmtflags os_flags = os.flags();

    // Print the matrix entries.
    for (unsigned int i=0; i<this->m(); i++)    {
        for (unsigned int j=0; j<this->n(); j++)
            os << std::setw(15) << std::scientific
            << std::setprecision(8)  << this->el(i,j) << " ";
        os << std::endl;
    }

    // reset the original format flags
    os.flags(os_flags);
}

// ====================================================
void DenseMatrixBaseM::print (std::ostream& os) const
{
    for (unsigned int i=0; i<this->m(); i++)  {
        for (unsigned int j=0; j<this->n(); j++)
            os << std::setw(8)  << this->el(i,j) << " ";
        os << std::endl;
    }
    return;
}



