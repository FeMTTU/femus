#ifndef __numeric_vectorM_h__
#define __numeric_vectorM_h__

#include "FemusExtLib_conf.hpp"
#include "FemusDefault.hpp"

#include "Typedefs_conf.hpp" 
#include "SolverPackageEnum.hpp"
#include "Paralleltype_enum.hpp"


// C++ includes
#include <vector>
#include <set>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <iostream>



// forward declarations
class DenseVector;
class DenseSubVector;
class SparseMatrixM;



// ==========================================
// Numeric vector. Provides a uniform interface
// to vector storage schemes for different linear
// algebra libraries.
// ===============================================

class NumericVectorM
#ifdef LM_REFCOUNT
: public ReferenceCountedObject<NumericVectorM>
#endif
{
  // =====================================
  // DATA
  // =====================================
  protected:
  
  /// Flag to see if the Numeric assemble routines have been called yet
  bool _is_closed;
  /// Flag to tell if init  has been called yet
  bool _is_initialized;
  /// Type of vector
  ParallelTypeM _type;
  
  // =====================================
  // Constructor /Destructor
  // =====================================
public:
  /// Dummy-Constructor. Dimension=0
  explicit
  NumericVectorM (const ParallelTypeM = AUTOMATICM);
  
  /// Constructor. Set dimension to \p n and initialize all elements with zero.
  explicit
  NumericVectorM (const unsigned int n,
                 const ParallelTypeM = AUTOMATICM);
    
  /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
  NumericVectorM (const unsigned int n,
		 const unsigned int n_local,
                 const ParallelTypeM = AUTOMATICM);
  
  /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
  NumericVectorM (const unsigned int N,
		 const unsigned int n_local,
		 const std::vector<unsigned int>& ghost,
                 const ParallelTypeM = AUTOMATICM);
    
 /// Builds a \p NumericVectorM using the linear solver package 
 /// specified by \p solver_package
  static std::auto_ptr<NumericVectorM>
  build(const SolverPackage solver_package = LSOLVER); 
  /// Creates a copy of this vector and returns it in an \p AutoPtr.
  virtual std::auto_ptr<NumericVectorM > clone () const = 0;  
  
  /// Destructor, deallocates memory. 
  virtual ~NumericVectorM (){  clear ();}
  /// @returns the \p NumericalVectorM to a pristine state.
   virtual void clear (){ _is_closed= false;_is_initialized = false;}
   
 /// Call the assemble functions
   virtual void close () = 0; 
  /// Change the dimension of the vector to \p N. The reserved memory for
  /// this vector remains unchanged if possible, to make things faster, but
  /// this may waste some memory, so take this in the back of your head.
  /// However, if \p N==0 all memory is freed, i.e. if you want to resize
  /// the vector and release the memory not needed, you have to first call
  /// \p init(0) and then \p init(N). 
  virtual void init (const unsigned int,
		     const unsigned int,
		     const bool = false,
                     const ParallelTypeM = AUTOMATICM) = 0;
  
  /// call init with n_local = N,
  virtual void init (const unsigned int,
		     const bool = false,
                     const ParallelTypeM = AUTOMATICM) = 0;
    
  /// Create a vector that holds tha local indices plus those specified
  /// in the \p ghost argument.
  virtual void init (const unsigned int /*N*/,
		     const unsigned int /*n_local*/,
		     const std::vector<unsigned int>& /*ghost*/,
		     const bool /*fast*/ = false,
                     const ParallelTypeM = AUTOMATICM) = 0;

  /// Creates a vector that has the same dimension and storage type as
  /// \p other, including ghost dofs.
  virtual void init (const NumericVectorM& other,
                     const bool fast = false) = 0;
		     
		     
 ///  Creates the subvector "subvector" from the indices in the
 /// "rows" array.  Similar to the create_submatrix routine for
 /// the SparseMatrix class, it is currently only implemented for
 /// PetscVectors.
  virtual void create_subvector(NumericVectorM&,const std::vector<unsigned int>&) const
  {
    std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
    exit(0);
  }		     
    
  // =====================================
  // SETTING FUNCTIONS
  // =====================================  
  /// v(i) = value
  virtual void set (const unsigned int i, const Real value) = 0;
  /// v(i) += value
  virtual void add (const unsigned int i, const Real value) = 0;
  
 /// Set all entries to zero. Equivalent to \p v = 0
  virtual void zero () = 0; 
  /// \f$U(0-N) = s\f$: fill all components.
  virtual NumericVectorM & operator= (const Real s) = 0;
  ///  \f$U = V\f$: copy all components.
  virtual NumericVectorM & operator= (const NumericVectorM &V) = 0;
  ///  \f$U = V\f$: copy all components.
  virtual NumericVectorM & operator= (const std::vector<Real> &v) = 0;
  
    /// \f$ U=v \f$ where v is a std::vector<Real> andspecify WHERE to insert it  
  virtual void insert (const std::vector<Real>& v,
		       const std::vector<unsigned int>& dof_indices) = 0;
  /// \f$U=V\f$, and specify WHERE to insert
  virtual void insert (const NumericVectorM& V,
		       const std::vector<unsigned int>& dof_indices) = 0;
  /// \f$ U=V \f$ insert
  virtual void insert (const DenseVector& V,
		       const std::vector<unsigned int>& dof_indices) = 0;
  /// \f$ U=V \f$ and specify WHERE to insert it
  virtual void insert (const DenseSubVector& V,
		       const std::vector<unsigned int>& dof_indices) = 0;
  // =====================================
  // RETURN FUNCTIONS
  // =====================================
  /// @returns true if the vector has been initialized,
  virtual bool initialized() const { return _is_initialized; }
  /// @returns true if the vector is closed and ready 
  virtual bool closed() const { return _is_closed; }
  
  /// @returns the type (SERIAL, PARALLEL, GHOSTED) of the vector.
  ParallelTypeM type() const { return _type; }
  /// @returns the type (SERIAL, PARALLEL, GHOSTED) of the vector.
  ParallelTypeM & type() { return _type; }


  /// @returns the minimum element in the vector.
  virtual Real min () const = 0;  
  /// @returns the maximum element in the vector.
  virtual Real max () const = 0;
  /// returns the sum of the elements in a vector
  virtual Real sum() const = 0;

  /// @returns the \f$l_1\f$-norm of the vector, i.e.
  virtual Real l1_norm () const = 0;
  /// @returns the \f$l_2\f$-norm of the vector, i.e.
  virtual Real l2_norm () const = 0;
  /// @returns the maximum absolute value of the
  virtual Real linfty_norm () const = 0;

  /// @returns the \f$l_1\f$-norm of the vector, i.e.
  virtual Real subset_l1_norm (const std::set<unsigned int> & indices);
  /// @returns the \f$l_2\f$-norm of the vector, i.e.
  virtual Real subset_l2_norm (const std::set<unsigned int> & indices);
  /// @returns the maximum absolute value of the
  virtual Real subset_linfty_norm (const std::set<unsigned int> & indices);

  
  /// @returns dimension of the vector. 
  virtual unsigned int size () const = 0;
  /// @returns the local size of the vector (index_stop-index_start).
  virtual unsigned int local_size() const = 0;
  /// @returns the index of the first vector element
  virtual unsigned int first_local_index() const = 0;
  /// @returns the index+1 of the last vector element
  virtual unsigned int last_local_index() const = 0;
    
  ///Access components, returns \p U(i).
  virtual Real operator() (const unsigned int i) const = 0;
  /// @returns the element \p U(i)
  virtual Real el(const unsigned int i) const { return (*this)(i); }

  /**
   * Access multiple components at once.  \p values will be resized,
   * if necessary, and filled.  The default implementation calls \p
   * operator() for each index, but some implementations may supply
   * faster methods here.
   */
  virtual void get(const std::vector<unsigned int>& index, std::vector<Real>& values) const;

  // =====================================
  // algebra FUNCTIONS
  // =====================================   
  
  
  /// Addition operator. Fast equivalent to \p U.add(1, V).
  virtual NumericVectorM & operator += (const NumericVectorM &V) = 0;
  /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
  virtual NumericVectorM & operator -= (const NumericVectorM &V) = 0;  
  /// Multiplication operator. Equivalent to \p U.scale(a)
  NumericVectorM & operator *= (const Real a) { this->scale(a); return *this; }
  /// Division operator. Equivalent to \p U.scale(1./a)
  NumericVectorM & operator /= (const Real a) { this->scale(1./a); return *this; }


    
  /// \f$U(0-LIBMESH_DIM)+=s\f$. Addition of \p s to all components. Note
  virtual void add (const Real s) = 0;
  /// \f$U+=V\f$: Simple vector addition, equal to the
  virtual void add (const NumericVectorM& V) = 0;
  /// \f$U+=a*V\f$. Simple vector addition, equal to the
  virtual void add (const Real a, const NumericVectorM& v) = 0;
  
  /// \f$ U+=v \f$ where v is a DenseVector 
  virtual void add_vector (const std::vector<Real>& v,
			   const std::vector<unsigned int>& dof_indices) = 0;
  /// \f$U+=V\f$, where U and V are type 
  virtual void add_vector (const NumericVectorM& V,
			   const std::vector<unsigned int>& dof_indices) = 0;
  /// \f$U+=A*V\f$, add the product A*v
   virtual void add_vector (const NumericVectorM& /*A*/,const SparseMatrixM& /*V*/) = 0;

   virtual void resid (const NumericVectorM &/*rhs_in*/,const NumericVectorM& /*A*/,const SparseMatrixM& /*V*/) = 0;
    virtual     void matrix_mult (const NumericVectorM &vec_in,const SparseMatrixM &mat_in) = 0;
  /// \f$U+=A*V\f$, add the product of a \p ShellMatrix \p
//   void add_vector (const NumericVectorM& v,
// 		   const ShellMatrix<Real>& a); 
//        
  /// \f$ U+=V \f$ where U and V are type 
  virtual void add_vector (const DenseVector& V,
			   const std::vector<unsigned int>& dof_indices) = 0;

  /// Scale each element 
  virtual void scale (const Real factor) = 0;
  /// v = abs(v)... that is, each entry in v is replaced
  virtual void abs() = 0;
  /// Computes the dot product, p = U.V
  virtual Real dot(const NumericVectorM&) const = 0;
      
  /// Exchanges the values/sizes of two vectors.  
  virtual void swap (NumericVectorM &v);
  
  // =====================================
  // PARALLEL OPERATIONS
  // =====================================
  
  /// Creates a copy of the global vector in the local vector \p v_local.
  virtual void localize (std::vector<Real>& v_local) const = 0;
  /// Same, but fills a \p NumericVectorM instead
  virtual void localize (NumericVectorM& v_local) const = 0;
  /// Creates a local vector \p v_local 
  virtual void localize (NumericVectorM& v_local,
			 const std::vector<unsigned int>& send_list) const = 0;
  /// Updates a local vector with selected values from neighboring
  virtual void localize (const unsigned int first_local_idx,
			 const unsigned int last_local_idx,
			 const std::vector<unsigned int>& send_list) = 0;
  /// Creates a local copy of the global vector 
  virtual void localize_to_one (std::vector<Real>& v_local,
				const unsigned int proc_id=0) const = 0; 
  /// @returns \p -1 when \p this is equivalent to \p other_vector,
  virtual int compare (const NumericVectorM &other_vector,
		       const Real threshold = DEFAULT_NUMVEC_TOLERANCE) const;
  /// Computes the pointwise (i.e. component-wise) product of \p vec1
  virtual void pointwise_mult (const NumericVectorM& vec1,
			       const NumericVectorM& vec2) = 0;

			       
  // =====================================
  // PRINTING FUNCTIONS
  // =====================================  			       
		 /// Prints the local contents of the vector to the screen.
  virtual void print(std::ostream& os=std::cout) const;	       
  /// Prints the local contents of the vector to the screen.
  virtual void print_personal(std::ostream& os=std::cout) const=0;
  /// Prints the global contents of the vector to the screen.
  virtual void print_global(std::ostream& os=std::cout) const;
  /// Same as above but allows you to use stream syntax.
  friend std::ostream& operator << (std::ostream& os, const NumericVectorM& v){
    v.print_global(os); return os; }
  
  /// Print the contents of the matrix in hdf5 sparse matrix format. 
  virtual void print_hdf5(const std::string /*name="NULL"*/) const=0;


};


/*----------------------- Inline functions ----------------------------------*/



// ==============================================
inline NumericVectorM::NumericVectorM (const ParallelTypeM type) :
  _is_closed(false),  _is_initialized(false),  _type(type){}
// ==============================================
inline NumericVectorM::NumericVectorM (const unsigned int /*n*/,
                                 const ParallelTypeM type) :
  _is_closed(false),_is_initialized(false), _type(type){
   std::cout<< "Abstract base class! ";exit(0); // Abstract base class!
  // init(n, n, false, type);
}

// ==============================================
inline NumericVectorM::NumericVectorM (const uint /*n*/,const uint /*n_local*/,
                                 const ParallelTypeM type) :
  _is_closed(false),  _is_initialized(false),  _type(type){
 std::cout<< "Abstract base class! "; exit(0); // Abstract base class!
  // init(n, n_local, false, type);
}

// ==============================================
inline NumericVectorM::NumericVectorM (const uint /*n*/,const uint /*n_local*/,
				 const std::vector<uint>& /*ghost*/,
                                 const ParallelTypeM type) :
  _is_closed(false),  _is_initialized(false),  _type(type){
 std::cout<< "Abstract base class! "; exit(0); // Abstract base class!
  // init(n, n_local, ghost, false, type);
}


// ==============================================
inline void NumericVectorM::get(const std::vector<uint>& index, 
				std::vector<Real>& values) const{
  const unsigned int num = index.size();
  values.resize(num);
  for(unsigned int i=0; i<num; i++) values[i] = (*this)(index[i]);
}

// ==============================================
inline void  NumericVectorM::swap (NumericVectorM &v){
  std::swap(_is_closed, v._is_closed);
  std::swap(_is_initialized, v._is_initialized);
  std::swap(_type, v._type);
}

// ==============================================
inline void NumericVectorM::print(std::ostream& os) const {
  assert (this->initialized());
  os << "Size\tglobal =  " << this->size()
     << "\t\tlocal =  " << this->local_size() << std::endl;
  os << "#\tValue" << std::endl;
  for (unsigned int i=this->first_local_index(); i<this->last_local_index(); i++)
    os << i << "\t" << (*this)(i) << std::endl;
}
// ==============================================
inline void NumericVectorM::print_global(std::ostream& os) const{
  assert (this->initialized());
  std::vector<Real> v(this->size());
  this->localize(v);
  // Right now we only want one copy of the output
#ifdef HAVE_MPIM   //TODO
  if(static_cast<uint>(libMeshPrivateData::_processor_id)) return ;
#else
  if(0)return;
#endif
  os << "Size\tglobal =  " << this->size() << std::endl;
  os << "#\tValue" << std::endl;
  for (uint i=0; i!=v.size(); i++) os << i << "\t" << v[i] << std::endl;
//   printf("Size\tglobal =  %d \n", this->size());
//   printf("#\tValue \n");
//   for (uint i=0; i!=v.size(); i++) printf(" %d \t %g \n", i, v[i]);
}




#endif  // #ifdef __numeric_vector_h__
