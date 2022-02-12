/*=========================================================================

 Program: FEMUS
 Module: SparseMatrix
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <fstream>
#include <sstream>
#include "SparseMatrix.hpp"
#include "FemusConfig.hpp"
#include "NumericVector.hpp"
#include "PetscMatrix.hpp"


#ifdef HAVE_HDF5
#include "hdf5.h"
#endif

namespace femus {



// ====================================================================
//            SparseMatrix Methods: Constructor/Destructor/init/build
// ===================================================================

// =====================================================================================
/// This function builds a  SparseMatrix using the linear solver
/// package specified by  solver_package
  std::unique_ptr<SparseMatrix > SparseMatrix::build (// -----
    const SolverPackage solver_package //  solver_package
  ) { // =================================================================================
    // Build the appropriate vector
    switch (solver_package) {
#ifdef HAVE_PETSC // ------------------------------
      case PETSC_SOLVERS: {
        std::unique_ptr<SparseMatrix > ap (new PetscMatrix);
        return ap;
      }
#endif
#ifdef HAVE_TRILINOS // ----------------------------
      case TRILINOS_SOLVERSM: {
        std::unique_ptr<SparseMatrix > ap (new EpetraMatrix<double>);
        return ap;
      }
#endif
      default:
        std::cerr << "SolverPackage solver_package:  Unrecognized: " << solver_package;
        abort();
    }
  }

// =================================================
//            SparseMatrix Methods: Add/mult
// =================================================

// =========================================================
  void SparseMatrix::vector_mult (NumericVector& dest,
                                  const NumericVector& arg) const {
    dest.zero();
    this->vector_mult_add (dest, arg);
  }

// ===========================================================
  void SparseMatrix::vector_mult_add (NumericVector& dest,
                                      const NumericVector& arg) const {
    /* This functionality is actually implemented in the \p NumericVector class.  */
    dest.add_vector (arg, *this);
  }

// =================================================
//            SparseMatrix Methods: setting/return
// =================================================

// ==========================================================
/// This function sets all row entries to 0 then puts diag_value in the diagonal entry
  void SparseMatrix::zero_rows (std::vector<int>&, double) {
    std::cerr << " SparseMatrix::zero_rows: Not implememnted ";
    abort();
  }


// ===============================================
//            SparseMatrix Methods: Print/Read
// ==============================================
// // ----------------------------------------------
// void SparseMatrix::print(const std::string& name) const{
// std::cerr << "SparseMatrix::print:not implemented";
// //    std::ofstream outfile(name.c_str());
// // //    outfile<<" Matrix Level "<<std::endl <<0<<"  "<<m()<< std::endl;
// //    // data
// //    for (uint i=1; i<= m(); i++) {
// //      std::vector<uint> row;
// //      for(uint j=0; j < n(); j++)
// //        if((*this)(i-1,j) > 0.)
// //    row.push_back(j);
// //      outfile << row.size() << " ";
// //      for (uint j=0; j< row.size(); j++) outfile  << row[j] << "  ";
// //      outfile << std::endl;
// //    }
// //    outfile.close();
// //   //  std::cout  << Level << " level matrix written in  " <<name.c_str()<< std::endl;  std::cout  << " dim  " << Q_GetDim(&matrix)<< std::endl;
// //   return;
// }


// void SparseMatrix::print(std::ostream& os) const
//{
//   parallel_only();
//
//   libmesh_assert (this->initialized());
//
//   // We'll print the matrix from processor 0 to make sure
//   // it's serialized properly
//   if (libMesh::processor_id() == 0){
// //       libmesh_assert(this->_dof_map->first_dof() == 0);
//       for (unsigned int i=this->_dof_map->first_dof();
//            i!=this->_dof_map->end_dof(); ++i){
//           for (unsigned int j=0; j<this->n(); j++)   os << (*this)(i,j) << " ";
//           os << std::endl;
//         }
//
//       std::vector<unsigned int> ibuf, jbuf;
//       std::vector<double> cbuf;
//       unsigned int currenti = this->_dof_map->end_dof();
//       for (unsigned int p=1; p < libMesh::n_processors(); ++p)  {
//           Parallel::receive(p, ibuf);
//           Parallel::receive(p, jbuf);
//           Parallel::receive(p, cbuf);
//           libmesh_assert(ibuf.size() == jbuf.size());
//           libmesh_assert(ibuf.size() == cbuf.size());
//
//           if (ibuf.empty())
//             continue;
//           libmesh_assert(ibuf.front() >= currenti);
//           libmesh_assert(ibuf.back() >= ibuf.front());
//
//           unsigned int currentb = 0;
//           for (;currenti <= ibuf.back(); ++currenti)            {
//               for (unsigned int j=0; j<this->n(); j++)                {
//                   if (currentb < ibuf.size() && ibuf[currentb] == currenti && jbuf[currentb] == j) {
//                os << cbuf[currentb] << " ";
//                currentb++;
//                     }
//                   else os << static_cast<double>(0.0) << " ";
//                 }
//               os << std::endl;
//             }
//         }
//       for (; currenti != this->m(); ++currenti)   {
//           for (unsigned int j=0; j<this->n(); j++)  os << static_cast<double>(0.0) << " ";
//           os << std::endl;
//         }
//     }
//   else
//     {
//       std::vector<unsigned int> ibuf, jbuf;
//       std::vector<double> cbuf;
//
//       // We'll assume each processor has access to entire
//       // matrix rows, so (*this)(i,j) is valid if i is a local index.
//       for (unsigned int i=this->_dof_map->first_dof();
//            i!=this->_dof_map->end_dof(); ++i)        {
//           for (unsigned int j=0; j<this->n(); j++)     {
//               double c = (*this)(i,j);
//               if (c != static_cast<double>(0.0))                {
//                   ibuf.push_back(i); jbuf.push_back(j); cbuf.push_back(c);
//                 }
//      }
//         }
//       Parallel::send(0,ibuf); Parallel::send(0,jbuf); Parallel::send(0,cbuf);
//     }
//}
// =====================================================
  void SparseMatrix::print (const std::string& name) const {
    std::ofstream out (name.c_str());
    print (out);
  }

// =====================================================
/// This function reads len sparse matrix structures
  void SparseMatrix::read_len_hdf5 (const std::string namefile, // file name
                                    const int mode,              // linear or quadratic
                                    int len_row[],               // row lengths
                                    int len_off_row[]            // row off entry lengths
                                   ) {

#ifdef HAVE_HDF5

    hid_t  file_id = H5Fopen (namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    // row lengths -------------------------------------------------------------
    std::ostringstream name_dst;
    name_dst.str ("");
    name_dst << "LEN" << mode;
    hid_t dataset = H5Dopen (file_id, name_dst.str().c_str(), H5P_DEFAULT);
    hid_t status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, len_row);
    // matrix row off entry lengths --------------------------------------------
    name_dst.str ("");
    name_dst << "OFFLEN" << mode;
    dataset = H5Dopen (file_id, name_dst.str().c_str(), H5P_DEFAULT);
    status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, len_off_row);

    H5Fclose (file_id);
    return;
#else
    std::cout << "FEMuS has been installed with no HDF5" << std::endl;
    abort();
#endif
  }

// =====================================================
/// This function reads pos sparse matrix structure
  void SparseMatrix::read_pos_hdf5 (const std::string namefile, // file name
                                    const int mode,              // linear or quadratic
                                    int pos_row[]                // compressed row positions
                                   ) {
#ifdef HAVE_HDF5
    hid_t  file_id = H5Fopen (namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    std::ostringstream name_dst;
    name_dst.str ("");
    name_dst << "POS" << mode;
    hid_t dataset = H5Dopen (file_id, name_dst.str().c_str(), H5P_DEFAULT);
    H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_row);
    H5Fclose (file_id);
    return;
#else
    std::cout << "FEMuS has been installed with no HDF5" << std::endl;
    abort();
#endif
  }
// =====================================================
/// This function reads quad-linear matrix dimensions
  void SparseMatrix::read_dim_hdf5 (const std::string namefile, // file name
                                    int dim[]                   // dimensions
                                   ) {
#ifdef HAVE_HDF5
    hid_t  file_id = H5Fopen (namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    int ldim[2];
    // quadratic-quadratic -> dim[0]   ------------------------------------------------
    hid_t dataset = H5Dopen (file_id, "DIM0", H5P_DEFAULT);
    hid_t status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ldim);
    dim[0] = ldim[0];
    // linear-linear  -> dim[3]  ------------------------------------------------
    dataset = H5Dopen (file_id, "DIM3", H5P_DEFAULT);
    status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ldim);
    dim[3] = ldim[0];
    // quad-linear  -> dim[1]  ------------------------------------------------
    dataset = H5Dopen (file_id, "DIM1", H5P_DEFAULT);
    status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ldim);
    dim[1] = ldim[0];
    // linear-quad   -> dim[2] ------------------------------------------------
    dataset = H5Dopen (file_id, "DIM2", H5P_DEFAULT);
    status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ldim);
    dim[2] = ldim[0];

    H5Fclose (file_id);
    return;

#else
    std::cout << "FEMuS has been installed with no HDF5" << std::endl;
    abort();
#endif
  }


} //end namespace femus


