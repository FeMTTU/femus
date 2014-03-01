#include "sparse_matrixM.hpp"

#include "FemusExtLib_conf.hpp"

#include "Typedefs_conf.hpp"

// Local Includes
// #include "dof_map.h"
#include "Parallel.hpp"
#include "petsc_matrixM.hpp"
#include "numeric_vectorM.hpp"
 
 // C++ includes
#include <cstring>
#include <fstream>
 
 
//------------------------------------------------------------------
// SparseMatrix Methods

// Full specialization for Real datatypes 

std::auto_ptr<SparseMatrixM >
SparseMatrixM::build(const SolverPackage solver_package)
{
  // Build the appropriate vector
  switch (solver_package){
#ifdef FEMUS_HAVE_LASPACK
    case LASPACK_SOLVERSM:
      {	std::auto_ptr<SparseMatrixM> ap(new LaspackMatrixM);
	return ap;
      }
#endif
#ifdef FEMUS_HAVE_PETSC
    case PETSC_SOLVERS:
      {	std::auto_ptr<SparseMatrixM > ap(new PetscMatrixM);
	return ap;
      }
#endif
#ifdef HAVE_TRILINOSM
    case TRILINOS_SOLVERSM:
      {	std::auto_ptr<SparseMatrixM > ap(new EpetraMatrix<Real>);
	return ap;
      }
#endif
    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      abort();
    }
  std::auto_ptr<SparseMatrixM > ap(NULL);
  return ap;    
}


// =========================================================
void SparseMatrixM::vector_mult (NumericVectorM& dest,
				   const NumericVectorM& arg) const{
  dest.zero();
  this->vector_mult_add(dest,arg);
}

// ===========================================================
void SparseMatrixM::vector_mult_add (NumericVectorM& dest,
				       const NumericVectorM& arg) const{
  /* This functionality is actually implemented in the \p
     NumericVector class.  */
  dest.add_vector(arg,*this);
}

// ==========================================================
void SparseMatrixM::zero_rows (std::vector<int> &, Real){
  std::cerr << " This functionality isn't implemented or stubbed in every subclass yet " ;
  abort();
}

// ==============================================
void SparseMatrixM::read(const std::string& name) {
   int iproc=0;
   #ifdef HAVE_MPI
             MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
//             uint iproc= static_cast<unsigned int>(i);
 #endif
  std::ifstream infile(name.c_str());
  if (!infile) {std::cerr << " SparseMatrixM::read filename " 
    << name << " not found "<<std::endl;abort();}
  
  char  buf[256]="";  uint nrow,len,ipos,n_sdom,nodes_proc_in;
  
  infile >> buf >> buf >> buf;
  infile >> buf >> buf >> buf  >> nrow;  Graph graph;  graph.resize(nrow);
  infile >> n_sdom >> buf;
  uint *n_loc=new uint[n_sdom+1];
  for(uint ip=0;ip<=n_sdom;ip++)infile >>n_loc[ip];
  graph._m=nrow; graph._n=nrow;graph._ml=n_loc[iproc+1]-n_loc[iproc]; 
  graph._nl=n_loc[iproc+1]-n_loc[iproc];
  graph._ml_start=n_loc[iproc];
  
  // data
  for (uint i=0; i<(uint)nrow; i++) {
    
    infile >> len;  graph[i].resize(len+1);
    infile >> nodes_proc_in; 
    for (uint j=0; j< len; j++) {infile  >>ipos; graph[i][j] = ipos-1;}
    // last stored value is the number of in-matrix nonzero values
    graph[i][len] = nodes_proc_in;
//     printf(" grPH %d %d %d \n",libMesh::processor_id(), i,nodes_proc_in);
  }  
  infile.close();
  update_sparsity_pattern(graph); graph.clear();
  // std::cout << "   read matrix  " << std::endl;
}
// // ----------------------------------------------
// void SparseMatrixM::print(const std::string& name) const{
// std::cerr << "SparseMatrixM::print:not implemented";
// //    std::ofstream outfile(name.c_str());
// // //    outfile<<" Matrix Level "<<std::endl <<0<<"  "<<m()<< std::endl;
// //    // data
// //    for (uint i=1; i<= m(); i++) {
// //      std::vector<uint> row; 
// //      for(uint j=0; j < n(); j++) 
// //        if((*this)(i-1,j) > 0.) 
// // 	 row.push_back(j);
// //      outfile << row.size() << " ";
// //      for (uint j=0; j< row.size(); j++) outfile  << row[j] << "  ";
// //      outfile << std::endl;
// //    }
// //    outfile.close();
// //   //  std::cout  << Level << " level matrix written in  " <<name.c_str()<< std::endl;  std::cout  << " dim  " << Q_GetDim(&matrix)<< std::endl;
// //   return;
// }  


// void SparseMatrixM::print(std::ostream& os) const
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
//       std::vector<Real> cbuf;
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
// 	              os << cbuf[currentb] << " ";
// 	              currentb++;
//                     }
//                   else	os << static_cast<Real>(0.0) << " ";
//                 }
//               os << std::endl;
//             }
//         }
//       for (; currenti != this->m(); ++currenti)   {
//           for (unsigned int j=0; j<this->n(); j++)  os << static_cast<Real>(0.0) << " ";
//           os << std::endl;
//         }
//     }
//   else
//     {
//       std::vector<unsigned int> ibuf, jbuf;
//       std::vector<Real> cbuf;
// 
//       // We'll assume each processor has access to entire
//       // matrix rows, so (*this)(i,j) is valid if i is a local index.
//       for (unsigned int i=this->_dof_map->first_dof();
//            i!=this->_dof_map->end_dof(); ++i)        {
//           for (unsigned int j=0; j<this->n(); j++)	    {
//               Real c = (*this)(i,j); 
//               if (c != static_cast<Real>(0.0))                {
//                   ibuf.push_back(i); jbuf.push_back(j); cbuf.push_back(c);
//                 }
// 	    }
//         }
//       Parallel::send(0,ibuf); Parallel::send(0,jbuf); Parallel::send(0,cbuf);
//     }
//}



//------------------------------------------------------------------
// Explicit instantiations
// template class SparseMatrix<Number>;
