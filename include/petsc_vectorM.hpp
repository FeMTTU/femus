#ifndef __petsc_vectorMa_h__
#define __petsc_vectorMa_h__

#include "FemusExtLib_conf.hpp"

// #ifdef LIBMESH_HAVE_PETSC
#ifdef FEMUS_HAVE_PETSC

#include "Typedefs_conf.hpp"
#include "ParalleltypeEnum.hpp"


//#include "libmesh.h" //TODO
// C++ includes
#include <map>
#include <vector>
#include <cstdio> 
// Local includes
#include "numeric_vectorM.hpp"
#include "petsc_macroM.hpp"

/// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
# include "petscvec.h"
EXTERN_C_FOR_PETSC_END


#include "Casts.hpp"







// forward declarations
class SparseMatrixM;

// ====================================================
// Petsc vector. Provides a nice interface to the
// Petsc C-based data structures for parallel vectors.
// ===================================================


class PetscVectorM : public NumericVectorM
{
    // =====================================
    // DATA
    // =====================================
private:
  
    /// Actual Petsc vector datatype
    Vec _vec;
  
    /// If \p true, the actual Petsc array of the values of the vector is currently accessible.
    /// Realhat means that the members \p _local_form and \p _values are valid.
    mutable bool _array_is_present;

    /// Petsc vector datatype to hold the local form of a ghosted vector.
    /// Realhe contents of this field are only valid if the vector is
    /// ghosted and \p _array_is_present is \p true.
    mutable Vec _local_form;

    /// Pointer to the actual Petsc array of the values of the vector.
    /// Realhis pointer is only valid if \p _array_is_present is \p true.
    mutable PetscScalar* _values;
    
      /// Realype for map that maps global to local ghost cells.
    typedef std::map<unsigned int,unsigned int> GlobalToLocalMap;
 
    /// Map that maps global to local ghost cells (will be empty if not in ghost cell mode
    GlobalToLocalMap _global_to_local_map;
    
       /// Realhis boolean value should only be set to false
    /// for the constructor which takes a PETSc Vec object.
    bool _destroy_vec_on_exit;
    
#ifndef NDEBUG
    ///Size of the local form, for being used in assertations.  Realhe
    /// contents of this field are only valid if the vector is ghosted
    /// and \p _array_is_present is \p true.
    mutable unsigned int _local_size;
#endif

 

  

 

public:
    // =====================================
    // Constructor /Destructor
    // =====================================
    ///  Dummy-Constructor. Dimension=0
    explicit
    PetscVectorM (const ParallelType type = AUTOMATIC);
    /// Constructor. Set dimension to \p n and initialize all elements with zero.
    explicit
    PetscVectorM (const unsigned int n,
                  const ParallelType type = AUTOMATIC);
    /// Constructor. Set local dimension to \p n_local, the global dimension
    /// to \p n, and initialize all elements with zero.
    PetscVectorM (const unsigned int n,
                  const unsigned int n_local,
                  const ParallelType type = AUTOMATIC);
    /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
    PetscVectorM (const unsigned int N,
                  const unsigned int n_local,
                  const std::vector<unsigned int>& ghost,
                  const ParallelType type = AUTOMATIC);
    /// Constructor.  Creates a PetscVectorM assuming a valid PETSc Vec object
    PetscVectorM(Vec v);

    /// Destructor, deallocates memory
    ~PetscVectorM ();


    /// Creates a copy of this vector and returns it in an \p AutoPtr.
    std::auto_ptr<NumericVectorM > clone () const;


    /// Call the assemble functions
    void close ();
    /// @returns the \p PetscVectorM to a pristine state.
    void clear ();


    /// Change the dimension of the vector to \p N. Realhe reserved memory for
    /// this vector remains unchanged if possible, to make things faster, but
    /// this may waste some memory, so take this in the back of your head.
    /// However, if \p N==0 all memory is freed, i.e. if you want to resize
    /// the vector and release the memory not needed, you have to first call
    /// \p init(0) and then \p init(N).
    void init (const unsigned int N,
               const unsigned int n_local,
               const bool         fast=false,
               const ParallelType type=AUTOMATIC);

    /// call init with n_local = N,
    void init (const unsigned int N,
               const bool         fast=false,
               const ParallelType type=AUTOMATIC);

    /// Create a vector that holds tha local indices plus those specified in the \p ghost argument.
    void init (const unsigned int /*N*/,
                       const unsigned int /*n_local*/,
                       const std::vector<unsigned int>& /*ghost*/,
                       const bool /*fast*/ = false,
                       const ParallelType = AUTOMATIC);
    /// Creates a vector that has the same dimension
    void init (const NumericVectorM& other,
                       const bool fast = false);

    //   /**
    //    * Change the dimension to that of the
    //    * vector \p V. Realhe same applies as for
    //    * the other \p init function.
    //    *
    //    * Realhe elements of \p V are not copied, i.e.
    //    * this function is the same as calling
    //    * \p init(V.size(),fast).
    //    */
    //   void init (const NumericVectorM& V,
    // 	     const bool fast=false);

    // ===========================
    // SETTING FUNCTIONS
    // ===========================

    /// v(i) = value
    void set (const unsigned int i, const Real value);
    /// v(i) += value
    void add (const unsigned int i, const Real value);

    /// Set all entries to zero. Equivalent to \p v = 0
    void zero ();

    /// \f$U(0-N) = s\f$: fill all components.
    NumericVectorM & operator= (const Real s);
    ///  \f$U = V\f$: copy all components.
    NumericVectorM & operator= (const NumericVectorM &V);
    /// \f$U = V\f$: copy all components.
    PetscVectorM & operator= (const PetscVectorM &V);
    ///  \f$U = V\f$: copy all components.
    NumericVectorM & operator= (const std::vector<Real> &v);

    /// \f$ U=v \f$ where v is a std::vector<Real>
    void insert (const std::vector<Real>& v,const std::vector<unsigned int>& dof_indices);
    /// \f$U=V\f$, where U and V are type NumericVectorM and you
    void insert(const NumericVectorM& V,const std::vector<unsigned int>& dof_indices);
    /// \f$ U=V \f$ where V is type DenseVector
    void insert (const DenseVector& V, const std::vector<unsigned int>& dof_indices);
    /// \f$ U=V \f$ where V is type DenseSubVector
     void insert (const DenseSubVector& V, const std::vector<unsigned int>& dof_indices);

    // ===========================
    // RETURN FUNCTIONS
    // ===========================

    /// @returns the minimum element in the vector.
    Real min () const;
    /// @returns the maximum element in the vector.
    Real max () const;
    /// @returns the sum of values in a vector
    Real sum () const;

    /// @returns the \f$l_1\f$-norm of the vector
    Real l1_norm () const;
    /// @returns the \f$l_2\f$-norm of the vector
    Real l2_norm () const;
    /// @returns the maximum absolute value of the elements of this vector
    Real linfty_norm () const;

    /// @returns dimension of the vector
    unsigned int size () const;

    /// @returns the local size of the vector (index_stop-index_start)
    unsigned int local_size() const;
    /// @returns the index of the first vector element
    unsigned int first_local_index() const;
    /// @returns the index of the last vector element
    unsigned int last_local_index() const;

    /// Maps the global index \p i to the corresponding global index. If
    /// the index is not a ghost cell, this is done by subtraction the
    /// number of the first local index.  If it is a ghost cell, it has
    /// to be looked up in the map.
    unsigned int map_global_to_local_index(const unsigned int i) const;

    /// Access components, returns \p U(i).
    Real operator() (const unsigned int i) const;

    /// Access multiple components at once.  Overloaded method that
    /// should be faster (probably much faster) than calling \p
    /// operator() individually for each index.
    void get(const std::vector<unsigned int>& index, std::vector<Real>& values) const;

    // ===========================
    // ALGEBRA FUNCTIONS
    // ===========================

    /// Addition operator. Fast equivalent to \p U.add(1, V).
    NumericVectorM & operator += (const NumericVectorM &V);
    /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
    NumericVectorM & operator -= (const NumericVectorM &V);

    /// \f$U(0-LIBMESH_DIM)+=s\f$. Addition of \p s to all components.
    void add (const Real s);
    /// \f$ U+=V \f$. Simple vector addition, equal to the
    void add (const NumericVectorM& V);
    /// \f$ U+=a*V \f$. Simple vector addition, equal to the
    void add (const Real a, const NumericVectorM&  v);
    
    
  
    
    /// \f$ U+=v \f$ where \p v is a std::vector
    void add_vector (const std::vector<Real>& v,
                     const std::vector<unsigned int>& dof_indices);

    /// \f$ U+=V \f$ where U and V are type  \p NumericVectorM and you
    void add_vector (const NumericVectorM& V,
                     const std::vector<unsigned int>& dof_indices);
    ///\f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
    void add_vector (const NumericVectorM &V,const SparseMatrixM &A);
    ///\f$U+=V \f$ where U and V are type DenseVector
    void add_vector (const DenseVector& V,const std::vector<unsigned int>& dof_indices);
  void resid (const NumericVectorM & rhs_in,const NumericVectorM& x,const SparseMatrixM& A);
     /// \f$U+=A*V\f$, add the product A*v
     void matrix_mult (const NumericVectorM &vec_in,const SparseMatrixM &mat_in);
   
    /// Scale each element of the vector by the given factor.
    void scale (const Real factor);

    /// v = abs(v)... that is, each entry in v is replaced
    void abs();

    /// Computes the dot product, p = U.V
     Real dot(const NumericVectorM& V) const;
    /// Computes the pointwise (i.e. component-wise) product of \p vec1
    /// and \p vec2 and stores the result in \p *this.
     void pointwise_mult (const NumericVectorM& vec1,
                                 const NumericVectorM& vec2);


    // ===========================
    // PARALLEL OPERATIONS
    // ===========================
     
    /// Creates a copy of the global vector in the local vector \p v_local.
    void localize (std::vector<Real>& v_local) const;
    /// Same, but fills a \p NumericVectorM instead of a \p std::vector.
    void localize (NumericVectorM& v_local) const;
    /// Creates a local vector \p v_local containing
    void localize (NumericVectorM& v_local,const std::vector<unsigned int>& send_list) const;
    /// Updates a local vector with selected values from neighboring processors
    void localize (const unsigned int first_local_idx,
                   const unsigned int last_local_idx,
                   const std::vector<unsigned int>& send_list);
    /// Creates a local copy of the global vector in
    /// \p v_local only on processor \p proc_id.  By
    /// default the data is sent to processor 0.
    void localize_to_one (std::vector<Real>& v_local,
                          const unsigned int proc_id=0) const;



    /// Creates a "subvector" from this vector using the rows indices of the "rows" array.
    void create_subvector(NumericVectorM& subvector,
                                  const std::vector<unsigned int>& rows) const;

    /// Swaps the raw PETSc vector context pointers.
    void swap (NumericVectorM &v);

    /// Returns the raw PETSc vector context pointer.  Note this is generally
    /// not required in user-level code. Just don't do anything crazy like calling VecDestroy()!
    Vec vec () { assert (_vec != NULL); return _vec;}

    // ===========================
    // PRINT
    // ===========================

    /**
     * Print the contents of the vector in Matlab
     * format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void print_personal(std::ostream& os=std::cout) const;
    void print_hdf5 (const std::string name="NULL") const;

    // ===========================
    // PROTECTED
    // ===========================
protected:
    /// Queries the array (and the local form if the vector is ghosted) from Petsc.
    void _get_array(void) const;
    ///  Restores the array (and the local form if the vector is ghosted) to Petsc.
    void _restore_array(void) const;
};


/*----------------------- Inline functions ----------------------------------*/

// ===========================================================
inline PetscVectorM::PetscVectorM (const ParallelType type)
        : _array_is_present(false),
        _local_form(NULL),
        _values(NULL),
        _global_to_local_map(),
        _destroy_vec_on_exit(true){  
	  this->_type = type;}

// =========================================================
inline PetscVectorM::PetscVectorM (const unsigned int n,const ParallelType type)
        : _array_is_present(false),
        _local_form(NULL),
        _values(NULL),
        _global_to_local_map(),
        _destroy_vec_on_exit(true)
{    this->init(n, n, false, type);}

// ===========================================================
inline PetscVectorM::PetscVectorM (const unsigned int n,
                                   const unsigned int n_local,
                                   const ParallelType type)
        : _array_is_present(false),
        _local_form(NULL),
        _values(NULL),
        _global_to_local_map(),
        _destroy_vec_on_exit(true){
    this->init(n, n_local, false, type);
}

// ===========================================================
inline PetscVectorM::PetscVectorM (const unsigned int n,
                                   const unsigned int n_local,
                                   const std::vector<unsigned int>& ghost,
                                   const ParallelType type)
        : _array_is_present(false),
        _local_form(NULL),
        _values(NULL),
        _global_to_local_map(),
        _destroy_vec_on_exit(true){
    this->init(n, n_local, ghost, false, type);
}

// ===========================================
inline PetscVectorM::PetscVectorM (Vec v)
        : _array_is_present(false),
        _local_form(NULL),
        _values(NULL),
        _global_to_local_map(),
        _destroy_vec_on_exit(false)
{
    this->_vec = v;
    this->_is_closed = true;
    this->_is_initialized = true;

    /* We need to ask PETSc about the (local to global) ghost value
       mapping and create the inverse mapping out of it.  */
    int ierr=0;
    int petsc_local_size=0;
// // //     int petsc_size=0;
// // //     ierr = VecGetSize(_vec, &petsc_size);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGetLocalSize(_vec, &petsc_local_size);
     CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Get the vector type from PETSc.
  // As of Petsc 3.0.0, the VecType #define lost its const-ness, so we
  // need to have it in the code
#if PETSC_VERSION_LESS_THAN(3,0,0)  || !PETSC_VERSION_LESS_THAN(3,4,0)
  VecType type;
#else
  const VecType type;
#endif    
  ierr = VecGetType(_vec, &type);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
    
  if((strcmp(type,VECSHARED) == 0) || (strcmp(type,VECMPI) == 0))
  {
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
    ISLocalToGlobalMapping mapping = _vec->mapping;
#else
    ISLocalToGlobalMapping mapping;
    ierr = VecGetLocalToGlobalMapping(_vec, &mapping);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif    
    
    // If is a sparsely stored vector, set up our new mapping
    if (mapping)
    {
      const unsigned int local_size = static_cast<unsigned int>(petsc_local_size);
      const unsigned int ghost_begin = static_cast<unsigned int>(petsc_local_size);
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,4,0)
      const unsigned int ghost_end = static_cast<unsigned int>(mapping->n);
//       const int *indices = mapping->indices;
#else
      PetscInt n;
//       const int *indices = mapping->indices;
      ierr = ISLocalToGlobalMappingGetSize(mapping,&n);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      const unsigned int ghost_end = static_cast<unsigned int>(n);
#endif
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
      const PetscInt *indices = mapping->indices;
#else
      const PetscInt *indices;
      ierr = ISLocalToGlobalMappingGetIndices(mapping,&indices);
      CHKERRABORT(MPI_COMM_WORLD,ierr); 
#endif
      
      for(unsigned int i=ghost_begin; i<ghost_end; i++)
        _global_to_local_map[indices[i]] = i-local_size;
      this->_type = GHOSTED;
    }
    else
      this->_type = PARALLEL;
  }
  else
    this->_type = SERIAL;

  this->close();
  
// // // /*  
// // //     /* \p petsc_local_size is the number of non-ghost values.  If it
// // //        equals the global size, then we are a serial vector, and there
// // //        are no ghost values.  */
// // //     if (petsc_size!=petsc_local_size)  {
// // //         ISLocalToGlobalMapping mapping = _vec->mapping;
// // // 
// // //         // If is a sparsely stored vector, set up our new mapping
// // //         if (mapping)        {
// // //             const unsigned int local_size = static_cast<unsigned int>(petsc_local_size);
// // //             const unsigned int ghost_begin = static_cast<unsigned int>(petsc_local_size);
// // //             const unsigned int ghost_end = static_cast<unsigned int>(mapping->n);
// // //             for (unsigned int i=ghost_begin; i<ghost_end; i++)
// // //                 _global_to_local_map[mapping->indices[i]] = i-local_size;
// // //             this->_type = GHOSTED;
// // //         }
// // //         else this->_type = PARALLEL;
// // //     }
// // //     else    this->_type = SERIAL;*/
      
}

// ======================================================
inline PetscVectorM::~PetscVectorM () { this->clear ();}

// =-=======================================================
inline void PetscVectorM::init (const unsigned int n,
                                const unsigned int n_local,
                                const bool fast,
                                const ParallelType type){
  
  
  int ierr=0;
    int petsc_n=static_cast<int>(n);
    int petsc_n_local=static_cast<int>(n_local);
    // Clear initialized vectors
    if (this->initialized())  this->clear();
    if (type == AUTOMATIC)    {
        if (n == n_local)   this->_type = SERIAL;
        else    this->_type = PARALLEL;
    }
    else    this->_type = type;
    
//  int size0; MPI_Comm_size(MPI_COMM_WORLD,&size0);
//  if(size0>1) this->_type = PARALLEL;
    assert ((this->_type==SERIAL && n==n_local) || this->_type==PARALLEL);

    // create a sequential vector if on only 1 processor
    if (this->_type == SERIAL)    {
        ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &_vec); CHKERRABORT(PETSC_COMM_SELF,ierr);
        ierr = VecSetFromOptions (_vec);        CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
    // otherwise create an MPI-enabled vector
    else if (this->_type == PARALLEL) {
        assert (n_local <= n);
       ierr = VecCreateMPI (MPI_COMM_WORLD, petsc_n_local, petsc_n, &_vec);
//         ierr = VecCreateMPI (MPI_COMM_WORLD,PETSC_DECIDE, petsc_n, &_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr); ierr = VecSetFromOptions (_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    else   {std::cout << "Not good" <<std::endl;abort();}
//   int rank,size,Istart,Iend;
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);//CHKERRQ(ierr);
//   ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);//CHKERRQ(ierr);
//   ierr = VecGetOwnershipRange(_vec,&Istart,&Iend);//CHKERRQ(ierr);
//   printf(" Vec Rostartend %d %d %d \n ",Istart,Iend,rank);
    this->_is_initialized = true;
//  this->_is_closed = true;
    if (fast == false)  this->zero ();
}

// ==================================================
inline void PetscVectorM::init (const unsigned int n,
                                const bool fast,
                                const ParallelType type){
    this->init(n,n,fast,type);
}

// ==================================================
inline void PetscVectorM::init (const unsigned int n,
                                const unsigned int n_local,
                                const std::vector<unsigned int>& ghost,
                                const bool fast,
                                const ParallelType type){
    int ierr=0;
    int petsc_n=static_cast<int>(n);
    int petsc_n_local=static_cast<int>(n_local);
    int petsc_n_ghost=static_cast<int>(ghost.size());

    // If the mesh is disjoint, the following assertion will fail.
    // If the mesh is not disjoint, every processor will either have
    // all the dofs, none of the dofs, or some non-zero dofs at the
    // boundary between processors.
    assert(n_local == 0 || n_local == n || !ghost.empty());

    int* petsc_ghost = ghost.empty() ? PETSC_NULL :
                       const_cast<int*>(reinterpret_cast<const int*>(&ghost[0]));

    // Clear initialized vectors
    if (this->initialized())   this->clear();

    assert(type == AUTOMATIC || type == GHOSTED);
    this->_type = GHOSTED;

    /* Make the global-to-local ghost cell map.  */
    for (unsigned int i=0; i<ghost.size(); i++)  _global_to_local_map[ghost[i]] = i;

    /* Create vector.  */
    ierr = VecCreateGhost (MPI_COMM_WORLD, petsc_n_local, petsc_n,
                           petsc_n_ghost, petsc_ghost, &_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecSetFromOptions (_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    this->_is_initialized = true;
//  this->_is_closed = true;
    if (fast == false)    this->zero ();
}

// ==============================================================
inline void PetscVectorM::init(const NumericVectorM& other, const bool fast){
    // Clear initialized vectors
    if (this->initialized())   this->clear();
    const PetscVectorM& v = libmeshM_cast_ref<const PetscVectorM&>(other);
    // Other vector should restore array.
    if (v.initialized())    {v._restore_array();}

    this->_global_to_local_map = v._global_to_local_map;
    this->_is_closed      = v._is_closed;
    this->_is_initialized = v._is_initialized;
    this->_type = v._type;

    if (v.size() != 0)   {
        int ierr = 0;
        ierr = VecDuplicate (v._vec, &this->_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    if (fast == false)   this->zero ();
}



// ======================================
inline void PetscVectorM::close (){
    this->_restore_array();
    int ierr=0;

    ierr = VecAssemblyBegin(_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecAssemblyEnd(_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    if (this->type() == GHOSTED)  {
        ierr = VecGhostUpdateBegin(_vec,INSERT_VALUES,SCATTER_FORWARD);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
        ierr = VecGhostUpdateEnd(_vec,INSERT_VALUES,SCATTER_FORWARD);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    this->_is_closed = true;
}

// ===============================================
inline void PetscVectorM::clear (){
    if (this->initialized())    this->_restore_array();

    if ((this->initialized()) && (this->_destroy_vec_on_exit))    {
        int ierr=0;
        ierr = VecDestroy(&_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    this->_is_closed = this->_is_initialized = false;
    _global_to_local_map.clear();
}

// =========================================
inline void PetscVectorM::zero (){
    this->_restore_array();
    int ierr=0;
    PetscScalar z=0.;

    if (this->type() != GHOSTED) {
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//         // 2.2.x & earlier style
//         ierr = VecSet (&z, _vec);
//         CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
        // 2.3.x & newer
        ierr = VecSet (_vec, z);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
    }
    else    {
        /* Vectors that include ghost values require a special
        handling.  */
        Vec loc_vec; ierr = VecGhostGetLocalForm (_vec,&loc_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//         // 2.2.x & earlier style
//         ierr = VecSet (&z, loc_vec);
//         CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
        // 2.3.x & newer
        ierr = VecSet (loc_vec, z);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
        ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
        CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
}

// ========================================================
inline std::auto_ptr<NumericVectorM > PetscVectorM::clone () const{
    std::auto_ptr<NumericVectorM> cloned_vector (new PetscVectorM);
    cloned_vector->init(*this, true);
    *cloned_vector = *this;
    return cloned_vector;
}

// ============================================
inline unsigned int PetscVectorM::size () const{
    assert (this->initialized());
    int ierr=0, petsc_size=0;
    if (!this->initialized())   return 0;
    ierr = VecGetSize(_vec, &petsc_size);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    return static_cast<unsigned int>(petsc_size);
}

// ================================================
inline
unsigned int PetscVectorM::local_size () const{
    assert (this->initialized());
    int ierr=0, petsc_size=0;
    ierr = VecGetLocalSize(_vec, &petsc_size);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    return static_cast<unsigned int>(petsc_size);
}

// =======================================================
inline unsigned int PetscVectorM::first_local_index () const{
    assert (this->initialized());
    int ierr=0, petsc_first=0, petsc_last=0;
    ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return static_cast<unsigned int>(petsc_first);
}

// =======================================================
inline unsigned int PetscVectorM::last_local_index () const{
    assert (this->initialized());
    int ierr=0, petsc_first=0, petsc_last=0;
    ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return static_cast<unsigned int>(petsc_last);
}

// ============================================================
inline unsigned int PetscVectorM::map_global_to_local_index (const unsigned int i) const{
    assert (this->initialized());

    int ierr=0, petsc_first=0, petsc_last=0;
    ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    const unsigned int first = static_cast<unsigned int>(petsc_first);
    const unsigned int last = static_cast<unsigned int>(petsc_last);

    if ((i>=first) && (i<last))    {        return i-first;    }

    GlobalToLocalMap::const_iterator it = _global_to_local_map.find(i);
    assert (it!=_global_to_local_map.end());
    return it->second+last-first;
}

// ====================================================
inline Real PetscVectorM::operator() (const unsigned int i) const{
    this->_get_array();
    const unsigned int local_index = this->map_global_to_local_index(i);
#ifndef NDEBUG
    if (this->type() == GHOSTED) assert(local_index<_local_size);
#endif
    return static_cast<Real>(_values[local_index]);
}

// =====================================================================
inline void PetscVectorM::get(const std::vector<unsigned int>& index, std::vector<Real>& values) const{
    this->_get_array();

    const unsigned int num = index.size();
    values.resize(num);

    for (unsigned int i=0; i<num; i++) {
        const unsigned int local_index = this->map_global_to_local_index(index[i]);
#ifndef NDEBUG
        if (this->type() == GHOSTED) assert(local_index<_local_size);
#endif
        values[i] = static_cast<Real>(_values[local_index]);
    }
}

// ============================================
inline Real PetscVectorM::min () const{
    this->_restore_array();
    int index=0, ierr=0;    PetscReal min=0.;
    ierr = VecMin (_vec, &index, &min);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // this return value is correct: VecMin returns a PetscReal
    return static_cast<Real>(min);
}

// =======================================
inline Real PetscVectorM::max() const{
    this->_restore_array();
    int index=0, ierr=0;  PetscReal max=0.;
    ierr = VecMax (_vec, &index, &max);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // this return value is correct: VecMax returns a PetscReal
    return static_cast<Real>(max);
}

// =====================================================
inline void PetscVectorM::swap (NumericVectorM &other){
    NumericVectorM::swap(other);
    PetscVectorM& v = libmeshM_cast_ref<PetscVectorM&>(other);
    std::swap(_vec, v._vec);
    std::swap(_destroy_vec_on_exit, v._destroy_vec_on_exit);
    std::swap(_global_to_local_map, v._global_to_local_map);
    std::swap(_array_is_present, v._array_is_present);
    std::swap(_local_form, v._local_form);
    std::swap(_values, v._values);
}

// =================================================
inline void PetscVectorM::_get_array(void) const{
    assert (this->initialized());
    if (!_array_is_present)    {
        int ierr=0;
        if (this->type() != GHOSTED)	{
            ierr = VecGetArray(_vec, &_values);
            CHKERRABORT(MPI_COMM_WORLD,ierr);
        }
        else	{
            ierr = VecGhostGetLocalForm (_vec,&_local_form);
            CHKERRABORT(MPI_COMM_WORLD,ierr);
            ierr = VecGetArray(_local_form, &_values);
            CHKERRABORT(MPI_COMM_WORLD,ierr);
#ifndef NDEBUG
            int local_size = 0;
            ierr = VecGetLocalSize(_local_form, &local_size);
            CHKERRABORT(MPI_COMM_WORLD,ierr);
            _local_size = static_cast<unsigned int>(local_size);
#endif
        }
        _array_is_present = true;
    }
}

// ====================================================
inline void PetscVectorM::_restore_array(void) const{
    assert (this->initialized());
    if (_array_is_present) {
        int ierr=0;
        if (this->type() != GHOSTED) {
            ierr = VecRestoreArray (_vec, &_values); CHKERRABORT(MPI_COMM_WORLD,ierr);
            _values = NULL;
        }
        else	{
            ierr = VecRestoreArray (_local_form, &_values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
            _values = NULL;
            ierr = VecGhostRestoreLocalForm (_vec,&_local_form);  CHKERRABORT(MPI_COMM_WORLD,ierr);
            _local_form = NULL;
#ifndef NDEBUG
            _local_size = 0;
#endif
        }
        _array_is_present = false;
    }
}



#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_vector_h__








