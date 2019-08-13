#include "GenCase.hpp"

// C++
#include <sstream>
#include <cassert>
#include <cmath>
#include <algorithm> 
// FEMuS
#include "FemusConfig.hpp"
#include "FemusDefault.hpp"

#include "Domain.hpp"
#include "Box.hpp"

#include "VBTypeEnum.hpp"
#include "Files.hpp"
#include "XDMFWriter.hpp"
#include "GeomElemBase.hpp"

// LibMesh
#ifdef HAVE_LIBMESH
#include "libmesh/enum_elem_type.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/boundary_info.h"

using namespace libMesh;
#endif


namespace femus {

// ========================================================
GenCase::GenCase(const unsigned nolevels, const unsigned dim, const GeomElType geomel_type, const std::string mesh_file_in)
     : MultiLevelMeshTwo(nolevels,dim,geomel_type,mesh_file_in)
{

   _feelems.resize(QL);
  for (int fe=0; fe<QL; fe++) _feelems[fe] = GeomElemBase::build(_geomelem_id[get_dim()-1].c_str(),fe);
 
}

GenCase::~GenCase() {

    _el_fm_libm.clear();  // AAA These must be destroyed HERE because they are created in PARALLEL by both processors!
    _el_fm_libm_b.clear();  //TODO is this a real destructor?

    _nd_fm_libm.clear();

//   delete _el_libm_fm; //TODO how come i cant delete this one here? when i'm serial yes, but not when i'm parallel
    // Oh well YES, because i put this INSIDE an IFPROC==0 !!!!!!!!!

    //clean
  //remember that the CHILD destructor is not called here. This is because 
  //these things are of GeomElemBase type
  //TODO Do I have to put a delete [] or only delete? Do standard vectors have an overloading of the "delete" operator?
  for (int fe=0; fe<QL; fe++) delete _feelems[fe];
  
}



// =======================================================
void GenCase::GenerateCase(const std::string output_path)   {
#ifdef HAVE_LIBMESH //i am putting this inside because there are no libmesh dependent arguments
  
  
#ifdef DEFAULT_PRINT_TIME
    std::clock_t start_timeA=std::clock();
#endif

    _msh_coarse = new libMesh::Mesh( (libMesh::Parallel::Communicator) MPI_COMM_WORLD,get_dim());

    GenerateCoarseMesh();

    _msh_all_levs = new libMesh::Mesh(*_msh_coarse);

    RefineMesh();

    _bd_msht = new  libMesh::BoundaryMesh( (libMesh::Parallel::Communicator) MPI_COMM_WORLD, _msh_all_levs->mesh_dimension()-1);

    GenerateBoundaryMesh();

#ifdef DEFAULT_PRINT_TIME
    std::clock_t start_timeC=std::clock();
#endif

    GrabMeshinfoFromLibmesh();  //only proc==0

    delete _bd_msht;
    delete _msh_all_levs;
    delete _msh_coarse;

    CreateMeshStructuresLevSubd(output_path);    //only proc==0
    
    ComputeAndPrintMGOperators(output_path);    //only proc==0

    Delete();

#ifdef DEFAULT_PRINT_TIME
    std::clock_t end_timeC = std::clock();
    std::cout << " Print and Operators time ="<< double(end_timeC- start_timeC) / CLOCKS_PER_SEC << std::endl;
#endif
#ifdef DEFAULT_PRINT_TIME
    std::cout << " +*+*+* Total time ="<< double(end_timeC- start_timeA) / CLOCKS_PER_SEC << std::endl;
#endif

    
#endif //end have_libmesh    
    
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif       
    return;
}


//===============================================================================
//================= LIBMESH coarse Mesh OBJECT from FILE or FUNCTION ============
//===============================================================================
//here,the information about the shape must be given a priori here,
//    while in the case of external mesh it should be given consistently
void GenCase::GenerateCoarseMesh() const {
#ifdef HAVE_LIBMESH

#ifdef DEFAULT_PRINT_TIME
    std::clock_t start_timeA=std::clock();
#endif
    
        std::string config_dir  = DEFAULT_INPUTDIR;
        std::string f_mesh_read = _mesh_file;

        std::ostringstream mesh_infile;
        mesh_infile << "./" << config_dir << f_mesh_read;
        std::ifstream inf(mesh_infile.str().c_str());

    if (!inf || f_mesh_read == ""  ) {

        std::cout << " Internal mesh generator at level 0 \n";

        if ( GetDomain()->GetDomainFlag() == 0 ) {

	  
       Box* box = static_cast<Box*>(GetDomain());

//---Meshing -------
            uint* ninterv = new uint[get_dim()];
	    ninterv[0] = box->_domain_rtmap.get("nintervx");
            ninterv[1] = box->_domain_rtmap.get("nintervy");
            if ( get_dim() == 3 ) ninterv[2] = box->_domain_rtmap.get("nintervz");

            // fem element definition --------------------------------
            libMesh::ElemType libmname; //convert the _geomel name into the libmesh geom el name

            if ( get_dim() == 2 ) {
            if (      _geomelem_id[get_dim()-1] == "quad") libmname = libMesh::QUAD9;
            else if ( _geomelem_id[get_dim()-1] == "tri")  libmname = libMesh::TRI6;
            libMesh::MeshTools::Generation::build_square
            (*_msh_coarse, ninterv[0], ninterv[1], box->_lb[0], box->_le[0], box->_lb[1], box->_le[1],libmname);
	    }
	    else if ( get_dim() == 3 ) {
            if (      _geomelem_id[get_dim()-1] == "hex")  libmname = libMesh::HEX27;
            else if ( _geomelem_id[get_dim()-1] == "tet")  libmname = libMesh::TET10;
            libMesh::MeshTools::Generation::build_cube
            (*_msh_coarse,  ninterv[0], ninterv[1],  ninterv[2], box->_lb[0], box->_le[0], box->_lb[1], box->_le[1], box->_lb[2], box->_le[2],libmname);
	    }
            else {         std::cout << " Dim 1 not implemented \n"; abort(); }
            
        } //box
        
        else { std::cout << " Domain shape not implemented for libmesh generation \n"; abort();  }

    }
    else {
        std::cout << " Reading Mesh File at level 0 \n";

        _msh_coarse->read(mesh_infile.str().c_str());  //is this read in parallel or only by proc=0?

    }
    
    inf.close();

    _msh_coarse->print_info();

#ifdef DEFAULT_PRINT_TIME
    std::clock_t end_timeA=std::clock();
    std::cout << " *+* Generation/Reading coarse mesh time ="
              << double(end_timeA- start_timeA) / CLOCKS_PER_SEC << std::endl;
#endif

#endif //end have_libmesh
    return;
}

//==============================================================================
//=============== GENERATE all the LEVELS for the LIBMESH Mesh OBJECT ==========
//==============================================================================
// LibMesh refinement -----------------------------------
//  well, one day we could make a DISTINCTION between the
//  LEVELS for MESH REFINEMENT and
//  LEVELS for MULTIGRID (NoLevels now)
// I mean, we could say that we refine the mesh up to some level,
// but then with the multigrid
// The point is that now the code is written WITHOUT this distinction
// and making this distinction is not trivial, one must check everywhere what to use
//   const uint mesh_refine = _utils.get_par("mesh_refine");
//   if (mesh_refine) {
void GenCase::RefineMesh() const {
#ifdef HAVE_LIBMESH

#ifdef DEFAULT_PRINT_TIME
    std::clock_t start_timeB=std::clock();
#endif

    std::cout << "\n LibMesh Mesh Refinement ---------  \n";
    libMesh::MeshRefinement mesh_refinement(*_msh_all_levs);
    mesh_refinement.uniformly_refine(_NoLevels-1);

#ifdef DEFAULT_PRINT_TIME
    std::clock_t end_timeB=std::clock();
    std::cout << " *+* Generation refined mesh time ="
              << double(end_timeB- start_timeB) / CLOCKS_PER_SEC << std::endl;
#endif

#endif //end have_libmesh
    return;
}


//==============================================================================
//=============== GENERATE BOUNDARY MESH =======================================
//==============================================================================
//  Generating Boundary Mesh from top level -----------------------
// TODO what is wrong with this sync() function that we have to redefine?
//why do we actually need to modify the sync function? 
//isnt there just a function that fills the BoundaryMesh accordingly to the input mesh,
// in such a way that it keeps the properties correctly?
// it seems like the sync function fills the boundary mesh.

//TODO can we exploit the fact that BoundaryInfo is useful for containing boundary conditions
//to READ from GAMBIT and ASSOCIATE FLAGS From Gambit to Libmesh, AND THEN from LIBMESH to FEMUS?

void GenCase::GenerateBoundaryMesh() const {
#ifdef HAVE_LIBMESH

    std::cout << " LibMesh BOUNDARY generation --------- \n";
    _msh_all_levs->boundary_info->sync(*_bd_msht);
    _bd_msht->print_info();

#endif //end have_libmesh
    return;
}

//==============================================================================
//=============== GRAB MESH INFORMATION from LIBMESH (only proc0) ==============
//==============================================================================
/// This function print the mesh, or better does all the stuff...
/// we must isolate where we call libmesh functions from where
/// we use OUR arrays

///Basically, grabbing information from libmesh means
/// filling the _el_sto[] and _el_sto_b[] arrays,
/// together with the _nod_coords[] array
///once you have this interface with libmesh, you do the rest only in FEMuS.

void  GenCase::GrabMeshinfoFromLibmesh() {
#ifdef HAVE_LIBMESH

    if (_iproc == 0)  {  //serial function

//   msht contains ALL LEVELS
//  msh0 contains THE COARSE
//how about bd_msht? i think it has ONLY THE FINE

// ===================================
// Setup
// ===================================

//here you cannot know earlier what element type will be chosen.
// All you can do is prepare a list of possible element types
//and prepare a switch that can pick one or another
//this switch must yield a child "under the cover of the father"
//and all the child data can be retrieved from the father with get functions
//no, in our case it is impossible to make a choice at RUN-TIME.
//so, one has to do the new necessarily
//we can make a class that reads things the nodes from the constructor,
// then you instantiate and initialize each object of that class

        _n_nodes = _msh_all_levs->n_nodes();                   //these are the FINE nodes
        _n_elements_sum_levs[VV] = _msh_all_levs->n_elem();                 //from mesh

        libMesh::Mesh::const_element_iterator         it_t00 = _msh_coarse->elements_begin();
        const libMesh::Mesh::const_element_iterator  end_t00 = _msh_coarse->elements_end();
        libMesh::Mesh::const_element_iterator          it_tr = _msh_all_levs->elements_begin();
        const libMesh::Mesh::const_element_iterator   end_tr = _msh_all_levs->elements_end();
        libMesh::Mesh::const_element_iterator it_el = _msh_all_levs->elements_begin();
// counting
        _n_elements_vb_lev     = new uint*[VB];
        _n_elements_vb_lev[VV] = new uint[_NoLevels];
        _n_elements_vb_lev[BB] = new uint[_NoLevels];
        _NoNodesXLev           = new uint[_NoLevels + 1];

        for (int ilev=0; ilev<_NoLevels; ilev++) {
            _n_elements_vb_lev[VV][ilev] = 0;
            _n_elements_vb_lev[BB][ilev] = 0;
        }
        for (int ilev=0; ilev<_NoLevels+1; ilev++) _NoNodesXLev[ilev] = 0;

// structures with original mesh
        _nd_coords_libm = new double[_n_nodes*3];  // why 3? because always 3 coords
        for (int i=0;i<_n_nodes*3;i++) _nd_coords_libm[i] = 0.;

//===== store VV elem info
        // this one is initialized to all -1
        _el_sto = new ElemStoVol*[_n_elements_sum_levs[VV]];
        for (int i=0;i<_n_elements_sum_levs[VV];i++)     _el_sto[i] = new ElemStoVol(_elnodes[VV][QQ],get_dim());

//===== compute the number of boundary elements
        int n_el_bdry_all_levs = 0;
        for (; it_el != end_tr; ++it_el) {
            Elem * elem = *it_el;
            for (uint s=0; s < elem->n_sides(); s++) {
                if (elem->neighbor(s) == NULL)    {
                    n_el_bdry_all_levs++;
                }
            }
        }
        _n_elements_sum_levs[BB] = n_el_bdry_all_levs;

//===== store BB elem info
        // this one is initialized to all 0
        _el_sto_b = new ElemStoBdry*[_n_elements_sum_levs[BB]];
        for (int i=0; i<_n_elements_sum_levs[BB];i++)    _el_sto_b[i]= new ElemStoBdry(_elnodes[BB][QQ],get_dim());


        //============= all the "new" done so far are deleted in the other routine, TODO check that

        //elem sto is in the order given by LIBMESH
        //we must construct a new elem ordering,
        //and always be able to go back to libmesh ordering
// ================================================
//    ELEMENTS (in elem_sto and bd_elem_sto)
// ================================================
        // Storing the mesh information
        //we fill _el_sto by count_el, which means "how the mesh iterator runs"
        //one expects that the mesh iterator goes by elem id, but it may not be true, i guess
        //instead count_eb goes with "when we find a boundary element, we add it"
	//TODO do we have, both in serial and in parallel, that count_e == elem->id(), so that 
	//the iterator goes from small to large values of elem->id() ?
	//I mean, the ELEMENT ID may have nothing to do with the ELEMENT POSITION in the list explored by the iterator!
	//I put a check and it seems like this thing is respected.
	// TODO there is another thing that must be considered!!! The VOLUME PROC PARTITIONING is done at the COARSE LEVEL,
	// while the BOUNDARY MESH PROC PARTITIONING IS DONE AFTER the REFINEMENTS!
	//This might cause MISMATCH in PROC NUMBERS!
        uint count_e=0;
        uint count_eb=0;

        for (; it_tr != end_tr; ++it_tr) {
            Elem * elem = *it_tr;
            _el_sto[count_e]->_id = elem->id();
	    /*CHECK*/ if (count_e != elem->id()) { std::cout << " The elements are not ordered by id ====== " << std::endl; abort();}
            _el_sto[count_e]->_lev  = elem->level();
//             _el_sto[count_e]->_subd = elem->processor_id();
            _n_elements_vb_lev[VV][elem->level()]++;

            for (uint inode=0; inode < _elnodes[VV][QQ];inode++) {
                int knode=elem->node(inode);      //libmesh node numbering
                _el_sto[count_e]->_elnds[inode]=knode;
            // coordinates storage
                for (int idim=0; idim<get_dim(); idim++) {
                    double xyz=  _msh_all_levs->point(knode)(idim);
                    _nd_coords_libm[knode+idim*_n_nodes]=xyz;
                }
            }

            // multigrid refinement
            Elem *parent=elem->parent();
            if (parent != NULL)   _el_sto[count_e]->_par=parent->id();
            
            int n_childs=elem->n_children();
            if (elem-> has_children()) {
                _el_sto[count_e]->_nch= n_childs;
                // children in the element
                for (int i_ch=0; i_ch< n_childs; i_ch++) {
                    Elem* child=(*it_tr)->child(i_ch);
                    _el_sto[count_e]->_elchs[i_ch]=child->id();
                }
            }
	    
            //  boundary
            for (uint s=0; s < elem->n_sides(); s++) {
                if (elem->neighbor(s) == NULL)    { //every element enters here only once
                    _el_sto_b[count_eb]->_id=count_eb;        //while for volume elems we use libmesh ordering, for boundary ones we count them "as we go"
                    _el_sto_b[count_eb]->_vol_id= elem->id();
                    _el_sto_b[count_eb]->_lev  = elem->level();
//                     _el_sto_b[count_eb]->_subd = elem->processor_id();
                    _n_elements_vb_lev[BB][elem->level()]++;
                    _el_sto_b[count_eb]->_nside= s;
                    AutoPtr<Elem> side(elem->build_side(s));
                    for (uint ns=0; ns< side->n_nodes(); ns++)  _el_sto_b[count_eb]->_elnds[ns]=side->node(ns);
                    count_eb++;
                }
            }  //side elements of each Volume Element
            
        count_e++;
       } //iterator for elems

        /*CHECK*/  assert((uint) _n_elements_sum_levs[VV] == count_e);
        /*CHECK*/  assert((uint) _n_elements_sum_levs[BB] == count_eb);
         
          
        // ================================================
        // ELEMENT GROUPING  (lev(0)>lev(n))
        // ================================================
        //TODO Notice here what happens: you are changing the processor numbers to the element list,
        // in such a way that FINER elements get the SAME PROCESSOR as the coarse ones, which is very correct!
        for (; it_t00 != end_t00; ++it_t00) {
            Elem * elem=*it_t00;
            int el_id=elem->id();
            _el_sto[el_id]->_lev =0;
            _el_sto[el_id]->_subd=elem->processor_id();
        }
        // subdomain ordering from level 0
        //for every level, pick the elements of that level
        //pick its children and set the same processor id as the father!
        //I also want to pick the BOUNDARY ELEMS and set the same proc as the father!
        //look, this process is not the best, it may be improved. but it works.
        for (int ilev=0; ilev<_NoLevels-1; ilev++) {
            for (int ielem=0;ielem<_n_elements_sum_levs[VV];ielem++) {
                if ( _el_sto[ielem]->_lev == ilev) {
                    for (int ich=0;ich < _el_sto[ielem]->_nch;ich++) {
                        _el_sto[_el_sto[ielem]->_elchs[ich]]->_subd = _el_sto[ielem]->_subd;
                    }
                }
            }
        }
        
        
        // ================================================
        //UPDATE BOUNDARY elems ==============
        // ================================================
        //The point is that I do not have the map that gives me the LIBMESH POSITION out of the count_eb position, so I'll do it again.
        //Now the processors of all levels and all 
        //these iterators are strange; for instance end_tr must be a "const const_element_iterator", otherwise the operator overloading of != does not work
        libMesh::Mesh::const_element_iterator    my_iter = _msh_all_levs->elements_begin();
        uint count_eb2=0;
            for (; my_iter != end_tr; ++my_iter) {
	               Elem* elem = *my_iter; 
	         for (uint s=0; s < elem->n_sides(); s++) {
                     if (elem->neighbor(s) == NULL)    {
                         _el_sto_b[count_eb2]->_subd = _el_sto[ _el_sto_b[count_eb2]->_vol_id]->_subd;  //PAY ATTENTION this is DIFFERENT from "elem->processor_id()"!
                         count_eb2++;
                       }
                 }
	    }
    

/*CHECK*/  for (int i=0; i < _n_elements_sum_levs[BB]; i++) {
	     if (_el_sto_b[i]->_subd !=  _el_sto[_el_sto_b[i]->_vol_id]->_subd )  { std::cout << " Subdomain mismatch between Boundary and Volume Elements " << std::endl; abort(); }
          }

 //here is the point: we are reattributing the numbers to the volume elements, so we must do the same for the boundary as well!!!         
 //Now we are first looping over volume elements, then picking the corresponding boundary elems,
 //and associating the boundary the proc number of the "father"
 //but then we CHANGE the PROC NUMBER of the FATHER and we DO NOT UPDATE the VOLUME NUMBER of the BOUNDARY ELEMENTS.
 //So, we should maybe do a FUNCTION that gives the REFERENCE to the VOLUME PROCESSOR;
 //but, since we are in a STRUCT, we must REUPDATE the NUMBERS for the BOUNDARY ELEMENTS!
 //as a matter of fact, it makes no sense to assign the processor numbers BEFORE and REASSIGN THEM AFTER.
 //let us do that ONLY AFTER!
 //Now the point is: given some boundary node we know what is the corresponding VOLUME id;
 //on the other hand, given a volume, do we know what are the corresponding BOUNDARY ID's?
 //they may be one, two, even more...
 //every volume element should have a LIST of WHAT SIDES are BOUNDARY SIDES and have their ID
 //of course the length of this list is different for each boundary side
 //so far we do not have this association from volume to boundary, let us still use from boundary to volume

 //can I get the libmesh element back from the element id? so get the element whose id is the given number,
 //without using iterator
 //I think it is maybe faster to increase the iterator rather than picking the element from the id.
 //now the point is:
 //first you loop to get the level
 //then you assign the processor to the coarse elems
 //then you swap all levels and for every level you assign the children with the same processor as the father
 
 //things are not done as nicely as I want.
 //First you loop over the elements to get the LEVEL.
 //Then, you consider the COARSE SUBDIVISION and based on that, you give the same subdivision to ALL the CHILDREN PROCESSORS, and to the BOUNDARY PROCESSORS also
 
 

    } //end proc==0
    
#endif //end have_libmesh
    return;
}


// =======================================================
// ========FOR EVERY CHILD AT ALL LEVELS, GIVE THE CORRESPONDING (UNIQUE) FATHER =======
// =====================================================
//so far I basically know HOW MANY ELEMENTS are there per SUBDOMAIN and PER LEVEL
//now it is time to associate to each boundary element the corresponding volume element,
//according to the FEMuS mesh ordering
//it restarts from zero at every Level
//so you'll make a map for every level
//actually our program is parallel, so we must also distinguish between subdomains
//The point is: first of all I loop over the volume elements
//and reorder them by LEV and SUBD
//Then, I loop over the Boundary elements
//and I reorder them by Vol and Subdomain.
//but, at that point, I lose the link between the BOUNDARY and VOLUME elements,
//they live like in two separate worlds.
//That is because eventually we associate each element,
//be it Boundary or Volume, to the NODE DOFS.
//but now, if we want to make constant FE, we must associate of course
//the boundary elements to VOLUME ELEMENT DOFS;
//so, we need this association.
//Now, the particular element can be retrieved by the list of its CONNECTIVITY
//which is printed exactly in the order given by the reordering.

//so, we want to make a map that from every boundary element in the femus ordering
// gives us back the corresponding FATHER element again in the femus ordering
//so, we must use libmesh.
//we loop on the boundary elements in femus ordering
//with el_sto_b we get the FATHER ID in libmesh
//from the id in libmesh we must get the id in femus

//qui non dobbiamo semplicemente memorizzare gli offset,
//ma dobbiamo stampare una PROPRIETA' RIFERITA ai SINGOLI ELEMENTI

//all'interno di ogni livello dobbiamo ordinare per sottodomini

//here you use the _off_el_vb in an ABSOLUTE MANNER, somehow.
// in this way they are not all together but you JUMP from processor to processor, for a given level
//so, off_el[VB] is used to explore the vector el_fm_libm[VB]
//now the point is that, while we are jumping along the el_fm_libm[VB] vector,
//we don't wanna jump along the _el_child_to_fath[lev] vector,
//so we have to put an ielem=0 outside the subdomain loop and increase it at the end of each

// when i explore the el_fm_libmb map, and I know that it is gonna be ordered by LEVELS and SUBDOMAINS
//now the point is: i wanna explore all the levels separately,
//so I know i'll have to JUMP along the el_fm_libmb map
//then, within the same level and subdomain, the order follows the el_fm_libmb order
//for the nodes it is the same: within the given level and subdomain, the order follows the v vector

//now, i can print this child_to_father map
//then, how do I read it in terms of subdomains?
//it is separated by level
//_off_el_vb is used for unique arrays that are ordered both by subdomain (external) and by level (internal)
//if an array is already separated into a given level, then
//I cannot use _off_el_vb in an ABSOLUTE MANNER
//   Again, let us stick to the fact that before we had to explore something BY JUMPS.
//    So, we know how to explore things BY JUMPS and turn them into things ATTACHED
//     when reading we have to do the inverse somehow.
// We have something which is ATTACHED and we know how to explore it when it is not attached.
//   So, from the NOT ATTACHED part we can deduce the offsets for the attached part,
//      in fact we can do:

// 	   j=0;
//       for (int isub=0;isub<_n_subdomains;isub++) {
//for (int i=0; i < _off_el_vb[BB][isub*_n_levels+lev+1] - _off_el_vb[BB][isub*_n_levels+lev]; i++) {//"i" is just a counter
//         _el_child_to_fath[lev][j] =
//          j++;
//            }
//          }

//ECCO DOV'E' il PROBLEMA! IL VOLUME MESH e il BOUNDARY MESH non sono suddivisi allo stesso modo!!!
// il boundary non appartiene allo stesso sottodominio del volume a cui appartiene!
//OVVIAMENTE IO VOGLIO CHE IL BOUNDARY APPARTENGA ALLO STESSO PROCESSORE DEL VOLUME!

//==============================================================================
//=============== CREATE FEMUS MESH, MAT, PROL, REST (only proc0) ==============
//==============================================================================
void GenCase::CreateMeshStructuresLevSubd(const std::string output_path) {

    if (_iproc == 0)   {  //serial function
//================================================
// AT THIS POINT ALL THE LIBMESH CALLS are over, we are only FEMuS
//=====================================

        // ================================================
        //      ELEMENTS
        _n_elements_sl_vb  = new int*[VB];
        for (int vb=0;vb < VB; vb++) {
            int offdim = _NoSubdom*_NoLevels+1;
            _n_elements_sl_vb[vb]= new int[offdim];
            for (int i=0;i<offdim; i++)     _n_elements_sl_vb[vb][i]=0;
        }
        ReorderElementBySubdLev_VV();
        ReorderElementBySubdLev_BB();
        ComputeElemOffsetsBySubdLevel();

        // ================================================
        //       NODES
        FillNodeSto();
        ReorderNodesBySubdLev();
        ComputeNodeOffsetsBySubdLevel();

//===== node sto no longer needed, it has filled nd_libm_fm
        for (int i=0;i<_n_nodes;i++)    delete  _nd_sto[i];
        delete [] _nd_sto;

        ElemChildToFather();

        ComputeMaxElXNode();

        ComputeNodeMapExtLevels();

        XDMFWriter::PrintMeshBiquadraticHDF5(output_path,*this);

// delete the boundary part, no more needed
        for (int i=0;i<_n_elements_sum_levs[BB];i++)      delete  _el_sto_b[i];
        delete [] _el_sto_b;

        delete [] _nd_coords_libm;


    } //end proc==0

    return;
}



void GenCase::ComputeAndPrintMGOperators(const std::string output_path) {

    if (_iproc == 0)   {  //serial function
      
            //this involves only VOLUME STUFF, no boundary stuff
            // instead, not only NODES but also ELEMENTS are used
        ComputeAndPrintMatrix(output_path); 
        ComputeAndPrintProl(output_path); 
        ComputeAndPrintRest(output_path);

    } //end proc==0

    return;
}


void GenCase::Delete() {

    if (_iproc == 0)   {  //serial function
//=====================================
//delete
//=====================================

        //==== ELEMENTS =========
        for (int i=0;i<VB;i++)     delete [] _n_elements_vb_lev[i];
        delete [] _n_elements_vb_lev;

        for (int vb=0;vb < VB; vb++)    delete [] _n_elements_sl_vb[vb];
        delete [] _n_elements_sl_vb;

        for (int i=0;i<VB;i++) delete []_off_el[i];
        delete [] _off_el;

        for (int i=0;i<_n_elements_sum_levs[VV];i++)   delete  _el_sto[i];    // delete [] el_sto[i]; with [] it doesnt work, check the static_cast thing TODO
        delete [] _el_sto;

        for (int lev=0; lev< _NoLevels; lev++) delete [] _el_child_to_fath[lev];
        delete [] _el_child_to_fath;

        //==== NODES =========
        delete [] _NoNodesXLev;

        for (int fe=0;fe < QL_NODES; fe++)   delete []_n_nodes_sl_ql[fe];
        delete [] _n_nodes_sl_ql;

        for (int i=0;i<QL_NODES;i++)    delete [] _off_nd[i];
        delete [] _off_nd;

        for (int ind=0;ind< _NoLevels+1; ind++) {
	  delete [] _Qnode_fine_Qnode_lev[ind];
          delete [] _Qnode_lev_Qnode_fine[ind];	  
	}
        delete [] _Qnode_fine_Qnode_lev;
        delete [] _Qnode_lev_Qnode_fine;


    } //end proc==0

    return;
}











// ==================================================
//Prol_ contains the POSITIONS: this is the vector of positions for each row (TODO should do it as a matrix)
//values_ contains the VALUES
//len_ off contains the OFFSETS
//here FINE and COARSE do not have anything in common with Quadratic and Linear

//FIRST, he computes the VALUES in all columns of each row.
//The values are initialized to zero. Then he fills.
//Then, if the values are nonzero, he COMPRESSES in order to obtain the SPARSITY PATTERN
//he does this for all rows
//il livello piu' basso

//In pratica in queste routine Mat Res Prol si usano per le prime volte gli ordinamenti dof 
//stabiliti dalla suddivisione in proc e livelli


void GenCase::ComputeAndPrintProl(const std::string output_path)  {

  int NegativeOneFlag = -1;
  double   PseudoZero = 1.e-8;
  
    std::string f_prol    = DEFAULT_F_PROL;
    std::string ext_h5    = DEFAULT_EXT_H5;

    std::ostringstream name;
    name << output_path << "/" << f_prol << ext_h5;
    hid_t file = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);

  for (int Level1 = 1; Level1 < _NoLevels; Level1++) {  //Level1 is the OUTPUT level (fine level) (the level of the ROWS)

    int Lev_c = Level1 - 1;
    int Lev_f = Level1;

    std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Lev_c << "_" << Lev_f;
    hid_t group = H5Gcreate(file, groupname_lev.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

//These FElevels return EXTENDED levels, to be used in the _Qnode_fine_Qnode_lev map
        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Lev_c; 
        FEXLevel_c[LL] = (Level1-1+_NoLevels)%(_NoLevels+1);
        FEXLevel_c[KK] = Lev_c;       
        int FEXLevel_f[QL];
        FEXLevel_f[QQ] = Level1;     
        FEXLevel_f[LL] = Level1-1;                                // AAA look at the symmetry, this is exactly (_n_levels + Level1 + 1)%(_n_levels + 1); ! //FINE Level for LINEAR:   Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_f[KK] = Lev_f;   

        uint** n_dofs_fe_lev;
        n_dofs_fe_lev = new uint*[QL];
        n_dofs_fe_lev[QQ] = _NoNodesXLev;            //equal the pointers
        n_dofs_fe_lev[LL] = _NoNodesXLev;            //equal the pointers
        n_dofs_fe_lev[KK] = _n_elements_vb_lev[VV];  //equal the pointers


      int **     Prol_pos = new int*[QL];
      double **  Prol_val = new double*[QL];
      int **          len = new int*[QL];
      int **       lenoff = new int*[QL];
      int    count[QL];
      int countoff[QL];
     
       
          for (int fe=0;fe<QL;fe++) {
               len[fe] = new int [n_dofs_fe_lev[fe][FEXLevel_f[fe]]+1];
            lenoff[fe] = new int [n_dofs_fe_lev[fe][FEXLevel_f[fe]]+1];
             count[fe] = 0;
          countoff[fe] = 0;
               len[fe][0] = 0;
            lenoff[fe][0] = 0;      
          }     
        
          for (int fe=0;fe<QL;fe++) {
	      int dim_large = n_dofs_fe_lev[fe][FEXLevel_f[fe]]*_elnodes[VV][fe];  //TODO the length of this array depends on HOW you EXPLORE things:
	                                                                                //here for every node we loop over its father nodes, so here's why
	                                                                                // in the RESTRICTOR the way of exploring is different
	                                                                                //actually I think that for the elements you don't need the -1, you can make it compact from the beginning...
	                                                                                //for the elements I could just allocate the exact number, because every fine element is influenced by only one coarse element
	                                                                                //so far I will keep the same scheme
            Prol_pos[fe] = new    int [ dim_large ];
            Prol_val[fe] = new double [ dim_large ];
            for (uint i=0; i < dim_large; i++) {
                Prol_pos[fe][i] = NegativeOneFlag;
                Prol_val[fe][i] =  0.;
            }    
        }

        // Prolongation Operator
        // For any FINE ELEMENT, find its own FATHER, then find WHAT CHILD our element is
        // Then loop over the child nodes
        //For any child node, loop over the SURROUNDING FATHER NODES
        // basically for every FINE DOF you have to find all the RELATED COARSE DOFS
        //for the element the related coarse dof is just the FATHER!
        for (int pr=0; pr<_NoSubdom; pr++) {
            for (int iel = _off_el[VV][pr*_NoLevels + Lev_f]; iel< _off_el[VV][pr*_NoLevels + Lev_f+1]; iel++) {
                int el_lib = _el_fm_libm[iel].second;
                int par_el = _el_sto[el_lib]->_par;
                int ch=0;
                int ch_id=0;
                while (el_lib != _el_sto[par_el]->_elchs[ch]) ch_id=++ch;

        for (int fe=0;fe<QL;fe++) {
		if (fe < KK) {
		for (uint node_of_child_elem=0; node_of_child_elem < _elnodes[VV][fe]; node_of_child_elem++) {
                    int el_n = _nd_libm_fm[ _el_sto[el_lib]->_elnds[node_of_child_elem] ];
                    for (uint node_of_par_elem=0;node_of_par_elem < _elnodes[VV][fe]; node_of_par_elem++) {
                        int par_n = _nd_libm_fm[ _el_sto[par_el]->_elnds[node_of_par_elem] ];
			  int mypp_node = _Qnode_fine_Qnode_lev[ FEXLevel_f[fe] ][el_n]*_elnodes[VV][fe]+node_of_par_elem;
                        Prol_pos[fe][ mypp_node ] = _Qnode_fine_Qnode_lev[ FEXLevel_c[fe] ][par_n];
                        Prol_val[fe][ mypp_node ] = _feelems[fe]->get_embedding_matrix(ch_id, node_of_child_elem, node_of_par_elem);
                    }
                }
	     }
		else if (fe == KK)  {
		  int sum_previous_sd=0;      //sum of the dofs AT THAT LEVEL coming from the PREVIOUSLY CONSIDERED SUBDOMAINS
                   for (int sd=0;sd<pr;sd++) sum_previous_sd += _off_el[VV][ sd*_NoLevels + Lev_f+1 ] - _off_el[VV][ sd*_NoLevels+Lev_f ];
                        int  dof_r = iel - _off_el[VV][pr*_NoLevels+Lev_f] + sum_previous_sd;

		   //what is the number in femus indexing of the parent element? par_el is in libmesh indexing!
		  int iel_femus_c = _el_libm_fm[par_el];
		  int sum_previous_sd_coarse = 0;
                   for (int sd=0;sd<pr;sd++) sum_previous_sd_coarse += _off_el[VV][ sd*_NoLevels + Lev_c+1 ] - _off_el[VV][ sd*_NoLevels+Lev_c ];
                        int  dof_c = iel_femus_c - _off_el[VV][pr*_NoLevels+Lev_c] + sum_previous_sd_coarse;
		  
		  Prol_pos[fe][ dof_r ] = dof_c;
		  Prol_val[fe][ dof_r ] = 1;
		}
		
	    }  //end fe
          } //iel
        }  //subd
        
        //========== COMPRESSION ================
        //this is needed in order to know what is in-diagonal and what is off-diagonal in the columns.
        //therefore, columns = coarse level
        int **dof_extreme_at_coarse_lev_sd = new int*[QL];
        for (int fe=0; fe < QL; fe++) {
            dof_extreme_at_coarse_lev_sd[fe] = new int[_NoSubdom+1];
            dof_extreme_at_coarse_lev_sd[fe][0] = 0;
            for (int idom=0;idom<_NoSubdom;idom++){
	           if (fe  < KK)  dof_extreme_at_coarse_lev_sd[fe][idom+1] =  dof_extreme_at_coarse_lev_sd[fe][idom] + ( _off_nd[fe][idom*_NoLevels + Lev_c+1] - _off_nd[fe][idom*_NoLevels]); 
              else if (fe == KK)  dof_extreme_at_coarse_lev_sd[fe][idom+1] =  dof_extreme_at_coarse_lev_sd[fe][idom] + ( _off_el[VV][idom*_NoLevels + Lev_c+1] - _off_el[VV][idom*_NoLevels + Lev_c]);
             }
        }        


        for (int fe=0; fe < QL; fe++) {
	  
        int ki=0;
        for (int kdom=0;kdom <_NoSubdom; kdom++) {
                    int delta_dof_at_fine_lev_sd;
                    if       (fe < KK)  delta_dof_at_fine_lev_sd = _off_nd[fe][kdom*_NoLevels + Lev_f+1] - _off_nd[fe][kdom*_NoLevels];
                    else if (fe == KK)  delta_dof_at_fine_lev_sd = _off_el[VV][kdom*_NoLevels + Lev_f+1] - _off_el[VV][kdom*_NoLevels+Lev_f];
  
            for (int iki = 0; iki < delta_dof_at_fine_lev_sd; iki++) {

                for (uint kj=0;kj < _elnodes[VV][fe]; kj++) {
                    int kk= _elnodes[VV][fe]*ki+kj;
                    if (fabs(Prol_val[fe][kk])> PseudoZero) {
                        Prol_pos[fe][ count[fe] ] = Prol_pos[fe][kk];
                        Prol_val[fe][ count[fe] ] = Prol_val[fe][kk];
                        count[fe]++;
                        if (Prol_pos[fe][kk] >=(int) dof_extreme_at_coarse_lev_sd[fe][kdom+1] ||
			    Prol_pos[fe][kk] <(int)  dof_extreme_at_coarse_lev_sd[fe][kdom]) countoff[fe]++;
                    }
                }
                   len[fe][ki+1] = count[fe];
                lenoff[fe][ki+1] = countoff[fe];
                ki++;
            }
        }
            XDMFWriter::PrintOneVarMGOperatorHDF5(name.str(),groupname_lev.str(),n_dofs_fe_lev[fe],count[fe],Prol_pos[fe],Prol_val[fe],len[fe],lenoff[fe],FEXLevel_f[fe],FEXLevel_c[fe],fe);

          delete []  dof_extreme_at_coarse_lev_sd[fe];

    } //end fe

//=======================================

     for (int fe=0;fe<QL;fe++) {
       delete [] Prol_pos[fe];
       delete [] Prol_val[fe];
       delete [] len[fe];
       delete [] lenoff[fe];
     }
     delete [] Prol_pos;
     delete [] Prol_val;
     delete [] len;
     delete [] lenoff;
     delete []  dof_extreme_at_coarse_lev_sd;

     H5Gclose(group);
   
    }
    //end Level1

    H5Fclose(file);

#ifdef DEFAULT_PRINT_INFO
    std::cout<< " GenCase::compute_and_print_MGOps: compute_prol end  \n";
#endif

    return;
}

// ===============================================================
/// This function computes the sparse matrix pattern
//from el_fm_libm you can get libmesh elem numbering
//what about femus elem numbering? you can do it just by running within _off_el
//from nd_libm_fm you can get FEMuS node numbering from libmesh node numbering
//if you want to get libmesh node numbering also, you can get it through an element loop,
//if you are doing an element loop, or you can get it with v.second if you are doing
// a direct node loop: here all the nodes we need are taken from element loops, so we dont need v.second
// in fact here we are not picking nodes just because they are nodes,
// but because they are nodes of some element in some subdomain and level
// once we distinguish things by level, in order to get the rest we must loop over subdomain,
//get the elements of that lev and subd, and from them you get whatever you want
// the alternative would be to have a structure with the nodes of each level and subdomain,
// that we actually have... _off_nd_ql
//so getting the nodes from there should be equivalent
//Mat is the sparsity pattern
// we want to do a unique function that computes the sparsity pattern of a matrix
//with rows of ONE FEM (quadratic, linear, const) and columns of ANOTHER FEM (quadratic, linear const)
// Mat[FEM1][FEM2][row][col] ... or at least Mat[FEM1][FEM2][row X MAXNONZEROCOL]

//For every level, compute the MATRIX STRUCTURE and print it to file
//So for each level we allocate, then we print, then we destroy.
//The submatrices are all separate so it seems reasonable to obtain a unique routine

//- first allocate the structures with the highest dimension possible (quadratic, finest level)
//- then, compress them


//we have to count how many nonzero numbers we have per row of the submatrix LQ
//also, if we are in parallel, we have to count how many of those columns
//belong to MYPROC and how many to ALL THE OTHERS
//consider the number of linear nodes for that level and that subdomain
//it seems like the problem is HOW the matrix Mat_lq
//Mat_lq contains the POS
// it must contain positive values where we have a node, otherwise -1

//here in practice we have a Mat_lq which is very large,
//and we rewrite it from the beginning, and then we pass the pointer
// with only the required length and we only print THAT PART of the big vector

//It looks to me like there is some LACK OF SYMMETRY in aux: what is it for?
//Ok, now we go home and we think of what happens in this spot
// The problem is only here, because it's here that  Mat_lq is used
// and len_lq is filled

// when you go from serial to parallel, the number of rows and columns is the same of course;
// the length of each row should be the same, but MAYBE the NODE NUMBERING in PARALLEL is DIFFERENT.
// IN FACT THAT'S EXACTLY WHAT HAPPENS!
//the node numbering changes because it must reflect the structure in Levels and Subdomains
//remember that this function is performed only by PROCESSOR 0!
//I have to put a breakpoint for this processor and debug only this process
//and see the difference in debugging the SERIAL code

//now, the point is: how come that we have to allocate such big vectors for the sparsity pattern?
// 	n_nodes*_elnodes[VV][QQ]*_maxelxnode
// in principle, we loop over nodes, and for each node we loop over the elements it belongs to,
//and for every element it belongs to we loop over its element nodes

//Remember that here we are building MODULAR sparsity patterns.
//Now the mesh corresponds to the FINE LEVEL.
// In the SERIAL CASE, the mesh is ordered with FIRST the LINEAR NODES and then all the others
//So the mesh allows us to distinguish FINE LINEAR and FINE QUADRATIC nodes, and so FINE LINEAR and FINE QUADRATIC dofs
//But we must do this AT ALL LEVELS.
// We can use one trick: nodes that are FINE LINEAR at the finest level correspond to the COARSE QUADRATIC at the coarser level!
//In this way we can ORDER the MATRIX SPARSITY PATTERNS at ALL LEVELS
//Of course, the order of QUADRATIC ROWS is the ORDER of the NODES in the MESH
//The FIRST NODES in the MESH are the LINEAR NODES

//So, how come that when we are in parallel we DO NOT INCREMENT THE COUNT?
// we know that the function that compute the sparsity pattern is performed only
//by the FIRST PROCESSOR

// We have to understand those few lines where the count is incremented
// First of all, what does the count give?
//For every [ROW][COL], count gives the OVERALL number of nonzero elements in the MATRIX[ROW][COL]
//So that is the length of the POS vector
//POS is the REAL SPARSITY PATTERN vector
//so somehow when we go parallel, we do not count the NONZERO ELEMENTS in the OTHER SUBDOMAINS

// Wait, i think that inside the gencase we can also work with ONE element
// Of course in that case you cannot do any PARALLEL thing, because you don't do
// any DOMAIN DECOMPOSITION of ONE ELEMENT!!!

//Now, POS gives all the positions for each row of the sparsity pattern
//     LEN gives the LENGTH of each row
//  OFFLEN gives the number of columns of the NODES (dofs) that BELONG to OTHER PROCESSORS

//Now, the first time the POS vector is filled it is actually MUCH BIGGER THAN IT SHOULD BE,
//it is full of -1, but it must have the BLOCKS of the LISTS OF COLUMNS, and these blocks must be NON-DISJOINT
//Actually these blocks MAY OR MAY NOT BE SEPARATED by -1 values, in some cases they may be CONTIGUOUS;
//so we must individuate a separation in two ways:
// either we find a -1, or we find reached the MAXIMUM possible value in the LinearToQuadratic correspondance.

//Remember that the sparsity patterns we give here are MODULAR,
// in the sense that they are referred to a COUPLE of VARIABLES (e.g. ROW linear, COL quadratic)
//where the indices start from zero in both cases;
//therefore, we refer directly to the MESH ENTITIES (LINEAR nodes, QUADRATIC nodes, ELEMENTS)
// as DofObjects.

//====== PARALLEL vs SERIAL ============
//When you go parallel, every MODULE of SPARSITY PATTERN still gives you a GLOBAL MATRIX,
// so for instance you have ALL the LINEAR ROWS against ALL the QUADRATIC COLUMNS,
//but what changes is basically THE ORDER IN WHICH YOU PUT the ROWS (and correspondingly the COLUMNS).
// The DOF OBJECT ORDERING reflects the LEVEL and SUBDOMAIN division

//============ LEVEL LOOP =================
//In this routine we have a loop over level, so we build the sparsity pattern for all COUPLE of DOF OBJECTS, we print and we start again.


//======= CAN WE MAKE THINGS MORE APPROPRIATE? ==============
//Now, the point is: we want to build the sparsity pattern of the matrix
// for each level
// for each FE COUPLE
//Now, of course we have to base ourselves on the fine quadratic level of the mesh,
//so eventually we'll have to deal with it,
//but at least can we construct something that stays on the level?

//now, the point is: we pick the numbers from the _Qnode_fine_Qnode_lev thing.
//So of course we have to have the CORRESPONDANCE between the MESH FINE QUADRATIC NODES
// and the DOFS, AT EVERY LEVEL

//_Qnode_fine_Qnode_lev[][] goes from the FEMUS NODE at the FINE LEVEL
// to the DOF value of the CURRENT LEVEL for one scalar value of that FE Family
// Here the construction of the dof map starts from a FINE and QUADRATIC MESH.
//If for instance we had a FINE and LINEAR mesh, we should turn it into a
// FINE and QUADRATIC (Libmesh has functions to do that)
// and then build everything on top of that.

//here we are preparing the sparsity pattern in a one-dimensional array first,
// in a sort of "Uncompressed row format"
// with a length still related to the fine world

//What is the maximum length of a line, in the worst case?
//It would be a line where all the nodes of the adjacent elements are involved.
//Every adjacent element may have a NUMBER of DOFS that DEPENDS on the FAMILY:
//1 dof for constants
//4 dof for linear (2D)
//9 dof for quadratic (2D)
//so the maximum number of involved dofs is _elnodes[VV][FE_COL]*_maxelxnode
//of course it is not as big as this one because there are REPEATED DOFS, like the DIAGONAL dof for sure.
//and of course this is multiplied by the total number of ROW DOFS.

//Praticamente con _gindexL sono i LIVELLI che distinguono tra i FE quadratici e i FE lineari

//TODO: We should use HDF5 matrices

//TODO: guardare tutti i loop (r c) che vengono aperti e poi chiusi
//============
//remember that we have to REORDER the NODE NUMBERING and ELEMENT NUMBERING
// so that we can distinguish
// the FINEST LEVEL they belong to
// and the SUBDOMAIN they belong to.

//If we didnt do the REORDERING, we could not distinguish SUBDOMAINS and LEVELS
//The reordering is based on the MESH LEVEL, which is the FINEST LEVEL, the FINE QUADRATIC LEVEL.
//based on this finest quadratic level,
//we renumber the nodes by distinguishing them in the "EXTENDED LEVEL" format, i.e.
//in every SUBDOMAIN:
//                0 LINEAR
// 0 QUADRATIC == 1 LINEAR
// 1 QUADRATIC == 2 LINEAR
// (n-1) QUADRATIC == (n) LINEAR
// (n) QUADRATIC
//============

//=========
// Qui abbiamo dei loop per
//  LIVELLO
//  SUBDOMAIN
// ELEMENTO
// righe QL
// colonne QL
//==========

//here we compute the matrix sparsity pattern for EVERY LEVEL,
//where LEVEL is associated to QUADRATIC LEVEL.

//==============================================================
//============ LEVEL LOOP ======================================
//==============================================================
//For every Level1 ("quadratic" level, or "reference" level Level1):
// - setup the file;
// - for every subdomain in that Level1,
//       and for every couple (FE_r, FE_c),
//          gather the list of column dofs that affect a row dof (by looping over elements; of course you'll have double occurrences)
// - remove all double occurrences;
// - print to file
// - delete the arrays

//===========================
// Pay attention because  _n_elements_vb_lev[VV]
//     and _n_nodes_lev
//     are defined on two different levels
//   the nodes are on n_levels+1
//   the dofs are on n_levels

//   so n_dofs_lev must be of dimension n_levels+1 for linear, but n_levels for Quadratic and n_levels for Konstants

// memG[][] gives me, for every row dof number, basically the number of elements to which that node belongs ( multiplied by a factor _elnodes[VV][QQ/*c*/] )

//==== Level1 is the REFERENCE GEOMETRICAL LEVEL:
//==== it goes from (0) to (NoLevels-1),
//==== it is used to get infos from the geometrical entities off_nd_ql, off_el

//===================
//_elnodes should be the "number of dofs in that element",
//so that inside an element you pick all the possible dofs that couple with your row dof.
//remember that you do not have to loop over the "geometrical entities
//to which there CAN BE ASSOCIATED DOFS",
//but we want to directly loop over the dofs that the element can have.
//now the point is all about:
//HOW DO WE EXTRACT DOF NUMBERS FROM GEOMETRICAL ENTITIES?
//at this moment, we extract DOF NUMBERS for the NODES from NODE GEOMETRICAL ENTITIES INFORMATIONS
//and we extract DOF NUMBERS for the ELEMENTS from ELEMENT GEOMETRICAL ENTITIES
//so, every dof should have the information on WHAT GEOMETRICAL ENTITY is associated to it.
//in fact, first of all the GEOMETRICAL ENTITIES in our code have been
//NUMERATED and DIVIDED by LEVELS and SUBDOMAINS.
//         from the DIVISION by LEVELS and SUBDOMAINS of the GEOMETRICAL ENTITIES
//we can obtain the DIVISION by LEVELS and SUBDOMAINS of the DOFS

//also, not always we have DOUBLE OCCURRENCES of DOFS in a given ROW.
//for instance, if the COLUMNS are KK, then that ELEMENT DOF influences
//only the CURRENT ELEMENT, so it is NOT picked MORE THAN ONCE.
//if we have an ELEMENT DOF ROW, again it is picked only once,
//so in its row the influences of the other dofs in that element cannot be repeated.
//so when the (KK,*) or  (*,KK) couples are involved,
//we do not need to REMOVE MULTIPLE COLUMN OCCURRENCES.

//So for the KK fe,instead of
//      int nd_r  = nd_libm_fm[ _el_sto[el_lib]->_elnds[k_r] ]; //node corresponding to the row k_r IN FEMUS ORDERING
//      int dof_r = _Qnode_fine_Qnode_lev[FELevel[r]][nd_r];              //dof corresponding to the node for that level
//where "nd_r" has the meaning of "GEOMETRICAL ENTITY TO WHICH ONE OR SOME DOF IS ASSOCIATED"
//we'll have something like
//         int nd_r  = iel;
//         int dof_r = nd_r;

// NO, PAY ATTENTION!!!
//always remember that the arrays _off_nd_ql and _off_el_vb
//tell us the ENDS of the INTERVALS in the array PER SUBDOMAIN and PER LEVEL.
//NOW, the point is that IN THE SAME ARRAY we have GEOMETRICAL ENTITIES (Elements and Nodes)
//belonging to DIFFERENT LEVELS and DIFFERENT SUBDOMAINS,
//but in the arrays they LIVE TOGETHER IN A UNIQUE LIST!
//therefore, we have to ALWAYS START FROM THE CORRECT POINT!
// So we should perhaps not use nd_r and directly do "FROM EL to DOF" like we do "FROM NODE TO DOF":
// so we should have a general map called "FROM GEOMETRICAL ENTITY TO DOF",
// or say: OK WE HAVE A DOF, what is the geometrical entity associated to it? Pick it and use the corresponding map
// yeah nd_r is not important in this case
//so we need the dofs from 0 to the number in our level (the sum of those between angle brackets)
//         <-->            <-->            <-->
//       || -- | -- | -- || -- | -- | -- || -- | -- | -- ||

// The number of dofs from previous subdomains is for (sd < current_pr) sum_previous_sd += _off_el_vb[VV][sd*_n_levels+Level1+1] - _off_el_vb[VV][sd*_n_levels+Level1];
//So in order to compute dof_r, you have to do sthg like
//        ...
//         int sum_previous_sd=0;      //sum of the dofs AT THAT LEVEL coming from the PREVIOUSLY CONSIDERED SUBDOMAINS
//         for (int sd=0;sd<pr;sd++) sum_previous_sd += _off_el_vb[VV][sd*_n_levels+Level1+1] - _off_el_vb[VV][sd*_n_levels+Level1];
//         int dof_r = iel - _off_el_vb[VV][pr*_n_levels+Level1] + sum_previous_sd;



//Then _off_nd_ql[fe][nsub*nlevs+1] = _off_el_vb[VV][nsub*nlevs+1]
//then the only other thing that changes is the way to go from GEOMETRICAL ENTITY to DOF,
// so the way you deal with (nd_r, dof_r) or (nd_c, dof_c)

//Since for constant element you do not have MULTIPLE OCCURRENCES,
//then you do not even need the compression stage, but i think that you can leave it there and it should work in principle.
//If you do not want to use it, then you may set _maxelxnode=1 AND multiply memG by ZERO (you dont have to use memG).
//As a matter of fact, _maxelxnode tells you "how many times one Node can influence other Nodes",
// so "how many times one Dof can INFLUENCE another Dof",
// and at the same time "How many times one Dof can BE INFLUENCED by another Dof".
//The number of times is with respect to an ELEMENT LOOP OVER THE WHOLE DOMAIN.
//============

// ========== DOF_BEGIN_SD ====================================
//we need to know if the node we are considering BELONGS or NOT to the CURRENT SUBDOMAIN,
//whether it's a quadratic or a linear node
//so we can explore the map per subdomain and per level,
//we simply have to construct the RANGE per SUBDOMAIN,
// and so we simply add the delta offsets which give us the AMOUNT of nodes in that range.
//Not only we know the AMOUNT of NODES in that range, but also we REORDERED them
//so that their numbers are CONTIGUOUS.
//dof_extreme_at_lev_sd holds the value of the Q or L NODES that SEPARATE one SUBDOMAIN to THE FOLLOWING.
//Since all the nodes are ordered SUBDOMAIN after SUBDOMAIN,
//and, inside each subdomain,EXTENDED LEVEL by EXTENDED LEVEL

//dof_extreme_at_lev_sd is needed because we need to know if a certain column dof is within or outside the range of the current subdomain.
// for a given Level1, the dofs are counted consecutively from zero by increasing the subdomain number.
//you need to do by increments because you would not have consecutive numbers for the matrix rows

//now pay attention to this guy: the offset for nodes is computed at given level Level1 is computed starting from level 0 to Level1.
// of course, in fact coarse nodes are also fine nodes
//on the other hand, elements only belong to ONE LEVEL
//so it doesnt seem to me to be appropriate to do a unique off_dof thing, or at least i must redefine it accordingly

// ==============================================
// ====== Mesh loop =============================
// ==============================================
//For every element, you loop over its dofs, either quadratic, or linear, or constant, ...
//For every dof, you get a matrix row, nd_r,
//Now, you stay inside that element, and you loop over the dofs inside that element, of any FE Family
//So, what are the "column_dofs" inside that element that may be involved with that "row_dof"?
//Every row dof is visited more than once, and so are the column_dofs for that node
//to be precise, we know that that dof_row is visited "_elxnode" times,
//because we know for every dof to how many elements it belongs
//so for every dof we actually have "_elxnode" connectivities,
//so we could dimension exactly the MatG[][] vector and we could have NO -1 AT ALL.
//the only thing is to see what happens next when we have to GATHER the DOUBLE/TRIPLE/QUADRUPLE/... COLUMNS CONNECTIVITIES corresponding to a certain ROW DOF
//So then, in each row we'll have to COMPRESS the MULTIPLE OCCURRENCES of the COLUMNS
//Even if we visit a row multiple times, we do not have multiple occurrences of rows, because every MatG is n_nodes x [_elnodes[VV][FE_COL]*_maxelxnode]
//So, if the column FE is LINEAR, then the maximum number of nodes in multiple connectivities is "_elxnode" times _elnodes[LL]

//So, one variation we could do first is to substitute the _elnodes[VV][QQ] with  _elnodes[VV][LL] in the increments, and that should bring no harm

//It is a good thing to base everything on the number of subdomains first so that it should be better prepared for a parallelization
//The idea would be that once you visit an element you fill the "expanded" sparsity pattern for all the possible FE couples
//So the filling of the sparsity pattern is "ELEMENT BASED", just like the ASSEMBLY of the MATRIX,
//so it means that all the blocks are filled at the same time, so only at the very end you will have all of them completed.
//Now, the point is this: with the ASSEMBLY you do not have to worry about passing MORE THAN ONCE on the SAME MATRIX POSITION,
// because you have to SUM ALL THE CONTRIBUTIONS.
//Instead with the sparsity pattern you DO NOT SUM... well, you could also SUM, and count the number of Elements per node, and divide by that number!
//So here for every FErow we fill the sparsity pattern of EVERY COLUMN,
// and then we go ahead with the next row (all the columns) and that's it.

// #if FM_DEVEL_GENCASE==1
//for the current Level1, this dof_r must be in the range for this level and this subdomain
//no wait a second... when you explore one element you may actually encounter a node that belongs to the other processor
//                         if (  dof_r < dof_extreme_at_lev_sd[r][pr] || dof_r >= dof_extreme_at_lev_sd[r][pr+1] )  {std::cout << " ====== Row mistake in matrix sparsity pattern construction ======== " << dof_r << " " << dof_extreme_at_lev_sd[r][pr] << " " << dof_extreme_at_lev_sd[r][pr+1] << " " << std::endl; abort();}
//this test is not correct
//as a matter of fact, in this moment I may be filling a row belonging to another processor
//by the way, if I were in parallel, I would need to communicate the "partial sparsity pattern" coming from this processor
//so is there some other way that I can check if this dof is well chosen?
//it must belong to all the possible ranges of the current level,
//because we are living in only one level of course,
//but it may or may not belong to the current processor
//we must find a way to verify it on-the-go.
//so the test for the "LEVEL RANGE" is
//OK			      if ( (dof_r < dof_extreme_at_lev_sd[r][0] || dof_r >= dof_extreme_at_lev_sd[r][_n_subdomains]))  { std::cout << " ====== Row mistake in matrix sparsity pattern construction ======== " << dof_r << " " << dof_extreme_at_lev_sd[r][pr] << " " << dof_extreme_at_lev_sd[r][pr+1] << " " << std::endl; abort();}
//    la condizione giusta e' che ALMENO UN INTERVALLO (se e' uno, e' uno solo) contenga il nostro dof_r.
//    la condizione da bloccare e' che NESSUN INTERVALLO contenga dof_r
//    quindi il controllo della condizione da bloccare bisogna farlo alla fine dell'ispezione;
//    invece il controllo della giustezza bisognerebbe farlo dentro il ciclo for, dicendo di interrompere il ciclo for.
//    Per fare questo si usa break.
//    Inoltre, questo e' proprio un classico esempio di ciclo "while", perche' tu non sai a priori il numero di iterazioni
//    che devi fare
// // //                           int truth=0;
// // //                           for (int sd=0; sd<_n_subdomains; sd++) {
// // // 			      if ( dof_r >= dof_extreme_at_lev_sd[r][sd] || dof_r < dof_extreme_at_lev_sd[r][sd+1] )  { truth += 1; break; }
// // // 			  }
// // // 			    if (truth==0) { std::cout << " ====== Row mistake in matrix sparsity pattern construction ======== " << dof_r << " " << dof_extreme_at_lev_sd[r][pr] << " " << dof_extreme_at_lev_sd[r][pr+1] << " " << std::endl; abort();}
// // // //                         if (r<KK) {
// // // //                             int sd=0;
// // // //                             while ( !(dof_r >= dof_extreme_at_lev_sd[r][sd] || dof_r < dof_extreme_at_lev_sd[r][sd+1]) && sd<_n_subdomains)  {
// // // //                                 sd++;
// // // //                             }
// // // //                             if (sd==_n_subdomains) {
// // // //                                 std::cout << " ====== Row mistake in matrix sparsity pattern construction ======== " << dof_r << " " << dof_extreme_at_lev_sd[r][pr] << " " << dof_extreme_at_lev_sd[r][pr+1] << " " << std::endl;
// // // //                                 abort();
// // // //                             }
// // // //                         }
// // // //                         else if (r==KK) {
// // // // // the control for the KK case is: dof_r is smaller than the sum of all subdomains
// // // //                             int sum_all_sd=0;
// // // //                             for (int sd=0;sd<_n_subdomains;sd++)  sum_all_sd += _off_el_vb[VV][sd*_n_levels+Level1+1] - _off_el_vb[VV][sd*_n_levels+Level1];
// // // //                             if (dof_r>sum_all_sd ) {
// // // //                                 std::cout << " ====== Row mistake in matrix sparsity pattern construction ======== " << std::endl;
// // // //                                 abort();
// // // //                             }
// // // //                         }
// #endif
//==============================================
//====== Begin the Big (r,c) FECouple loop =========
//==============================================
//Nota che il loop precedente aveva n_subdomains all'esterno
//e il loop sui vari tipi di elementi finiti in riga e colonna era all'interno
// Qui invece il loop FE r,c e' esterno mentre il loop sui subdomain e' interno
// in entrambi i casi la parallelizzazione avviene al livello del loop dei subdomains.
//Credo che tale loop debba essere eliminato(ogni processore lavora sul suo subdomain)
// e poi ad un certo punto bisogna fare in modo di comunicare tutto al processore di riferimento
//che si occupa di eliminare le doppie occorrenze
//e di fare la stampa

//=============================================
//============== COMPRESS MatG ================
//=============================================
//for every row (dof) we first build the "uncompressed" sparsity pattern.
//then we compress it here. How do we do?
//We have to:
// - eliminate all the double occurrences,
// - eliminate all the -1,
// - then establish the LENGTH and OFFLENGTH
// Taking into account a parallel structure, it seems to be better
// to loop over subdomains first, rather than on dofs directly

//remember that for the KK elements I dont need eliminating multiple occurrences


void GenCase::ComputeAndPrintMatrix(const std::string output_path) {

#ifdef DEFAULT_PRINT_INFO
    std::cout << " GenCase::compute_matrix:  start \n";
#endif

      const int NegativeOneFlag = -1;

//==============================================================
//============ ALLOCATIONS =====================================
//==============================================================

    int ** off_dof = new int*[QL];  //TODO do the delete for this
    off_dof[QQ] = _off_nd[QQ];  //You set the pointers to be equal
    off_dof[LL] = _off_nd[LL];  //You set the pointers to be equal
    off_dof[KK] = _off_el[VV];  //You set the pointers to be equal

    int *** lenG;
    int *** lenoffG;
    int *** MatG;
    int *** memG;

//========= CREATE THE FILE ============================
    std::string f_matrix  = DEFAULT_F_MATRIX;
    std::string ext_h5    = DEFAULT_EXT_H5;

        std::ostringstream name;
        name << output_path << "/" << f_matrix << ext_h5;
        hid_t file = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);

//==============================================================
//============ LEVEL LOOP ======================================
//==============================================================

    for (int Level1 = 0; Level1 < _NoLevels; Level1++) {

#ifdef DEFAULT_PRINT_TIME
        std::clock_t   start_time=std::clock();
#endif

    std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Level1;
    hid_t group = H5Gcreate(file, groupname_lev.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

//============== SET LEVELS and DOFS =======================
        int FELevel[QL];
        FELevel[QQ] = Level1;
        FELevel[LL] = (Level1+_NoLevels)%(_NoLevels+1);
        FELevel[KK] = Level1;
        uint** n_dofs_lev_fe;
        n_dofs_lev_fe = new uint*[QL];
        n_dofs_lev_fe[QQ] = _NoNodesXLev;            //equal the pointers
        n_dofs_lev_fe[LL] = _NoNodesXLev;            //equal the pointers
        n_dofs_lev_fe[KK] = _n_elements_vb_lev[VV];  //equal the pointers

        lenG       = new int**[QL];
        lenoffG    = new int**[QL];
        MatG       = new int**[QL];
        memG       = new int**[QL];

        //initial setups
        for (int r=0; r<QL;r++) {
            MatG[r]       = new int*[QL];
            memG[r]       = new int*[QL];

            for (int c=0; c<QL;c++) {

                memG[r][c]       = new int[ n_dofs_lev_fe[r][FELevel[r]] ];   //for every row, it memorizes the offset for multiple occurrences
                for (int im = 0; im < n_dofs_lev_fe[r][ FELevel[r] ]; im++)     memG[r][c][im] = 0;

                MatG[r][c]       = new int[n_dofs_lev_fe[r][FELevel[r]]*_elnodes[VV][c]*_maxelxnode];       // or even better (n_dofs_lev_fe[level][FE_ROW]) * ( _elnodes[VV][FE_COL]*_elxnode[ROW_DOF])
                for (uint im = 0; im <     n_dofs_lev_fe[r][FELevel[r]]*_elnodes[VV][c]*_maxelxnode; im++)  MatG[r][c][im] = NegativeOneFlag;
            }
        }

// ========== DOF EXTREMES AT LEVEL ====================================

        int ** dof_extreme_at_lev_sd;
        dof_extreme_at_lev_sd= new int*[QL];
        for (int fe=0; fe<QL; fe++)   dof_extreme_at_lev_sd[fe] = new int [_NoSubdom+1];

        for (int fe=0; fe<QL; fe++)  {
            dof_extreme_at_lev_sd[fe][0]=0;
            for (int idom=0;idom<_NoSubdom;idom++) {
                if (fe < KK)       dof_extreme_at_lev_sd[fe][idom+1] =  dof_extreme_at_lev_sd[fe][idom] + ( _off_nd[fe][idom*_NoLevels+Level1+1] - _off_nd[fe][idom*_NoLevels]);
                else if (fe == KK) dof_extreme_at_lev_sd[fe][idom+1] =  dof_extreme_at_lev_sd[fe][idom] + ( _off_el[VV][idom*_NoLevels+Level1+1] - _off_el[VV][idom*_NoLevels+Level1]);
                //from this you see that on a given level the dofs are of course numbered CONTIGUOUSLY; in practice this is the construction of the dof map
            }
        }

// ==============================================
// ====== Mesh loop =============================
// ==============================================

        for (int pr=0;pr<_NoSubdom;pr++) {
            for (int iel = _off_el[VV][pr*_NoLevels+Level1]; iel < _off_el[VV][pr*_NoLevels+Level1+1]; iel++) {
                int el_lib = _el_fm_libm[iel].second;

                for (int r=0; r<QL;r++) {
                    for (uint k_r=0; k_r < _elnodes[VV][r]; k_r++) { //this loop can be optimized because some operations are not done for KK
                        int nd_libm_r = _el_sto[el_lib]->_elnds[k_r];
                        int nd_r  = _nd_libm_fm[ nd_libm_r ]; //node corresponding to the row k_r
                        int dof_r = _Qnode_fine_Qnode_lev[FELevel[r]][nd_r];              //dof corresponding to the node for that level

                        if (r==KK) {
                            int sum_previous_sd=0;      //sum of the dofs AT THAT LEVEL coming from the PREVIOUSLY CONSIDERED SUBDOMAINS
                            for (int sd=0;sd<pr;sd++) sum_previous_sd += _off_el[VV][sd*_NoLevels+Level1+1] - _off_el[VV][sd*_NoLevels+Level1];
                            dof_r = iel - _off_el[VV][pr*_NoLevels+Level1] + sum_previous_sd;
                        }
// put debug control here on dof_r
                        for (int c=0; c<QL;c++) {
                            for (uint k_c=0;k_c < _elnodes[VV][c]; k_c++) {
                                int nd_c  = _nd_libm_fm[ _el_sto[el_lib]->_elnds[k_c]];  //node corresponding to the column k_c
                                int dof_c = _Qnode_fine_Qnode_lev[FELevel[c]][nd_c];

                                if (c==KK) {
                                    int sum_previous_sd=0;      //sum of the dofs AT THAT LEVEL coming from the PREVIOUSLY CONSIDERED SUBDOMAINS
                                    for (int sd=0;sd<pr;sd++) sum_previous_sd += _off_el[VV][sd*_NoLevels+Level1+1] - _off_el[VV][sd*_NoLevels+Level1];
                                    dof_c = iel - _off_el[VV][pr*_NoLevels+Level1] + sum_previous_sd;
                                }

                                MatG[r][c][ dof_r*_elnodes[VV][c]*_maxelxnode + memG[r][c][dof_r]*_elnodes[VV][c] + k_c ] = dof_c;

                            }   //>>>>>>> end k_c

                            //for the current row "r", we hit this column "c" with the current element, so we must remember it when we hit it with another element
                            //for every row dof nd_r... memG is basically the "number of involved elements"   //jump by another element in that row
                            memG[r][c][dof_r] +=  1;

                        } //>>>>>>>end c
                    }   //>>>>>>> end k_r
                }  //>>>>>>> end r
            } //>>>>>>> end element
        }  //>>>>>>> end subdomain


//here you don't use memG anymore

//  debug      for (uint im = 0; im < n_nodes_f*_elnodes[VV][QQ]*_maxelxnode; im++) {   std::cout << MatG[LL][QQ][im] << std::endl; }

//================================================================
//====== COMPRESS, COUNT LENGTHS, then PRINT, then DESTROY ===============
//===============================================================

        int** countG    = new int*[QL];      //[FE_R][FE_C] this is on the fly, for every row involving FE_R and FE_C; these are stored in lenG and lenoffG
        int** countoffG = new int*[QL];

//==============================================
//====== Begin the Big (r,c) FECouple loop =========
//==============================================
        for (int r=0; r<QL;r++) {

            countG[r]     = new int[QL];
            countoffG[r]  = new int[QL];
            lenG[r]       = new int*[QL];
            lenoffG[r]    = new int*[QL];

            for (int c=0; c<QL;c++) {

//=============================================
//============== COMPRESS MatG ================
//=============================================
                int *  aux;
                lenG[r][c]       = new int[ n_dofs_lev_fe[r][FELevel[r]] + 1 ];   //notice how actually these do not depend on "c"
                lenoffG[r][c]    = new int[ n_dofs_lev_fe[r][FELevel[r]] + 1 ];   //notice how actually these do not depend on "c"
                aux = new int[ n_dofs_lev_fe[r][FELevel[r]]*_elnodes[VV][c]*_maxelxnode];  //This has the same length as MatG

                for (uint i=0; i < n_dofs_lev_fe[r][FELevel[r]]*_elnodes[VV][c]*_maxelxnode; i++) aux[i] = 0;

                countG[r][c] = 0;
                countoffG[r][c] = 0;
                lenG[r][c][0] = 0;
                lenoffG[r][c][0] = 0;

                int ki=0; //the maximum value for ki is what? The sum of all the _off_nd_ql for all _nsubdomains, so it is exactly the number of nodes of the current level Level1
                for (int kdom=0;kdom <_NoSubdom; kdom++) {
                    int delta_dof_at_lev_sd;
                    if       (r < KK)  delta_dof_at_lev_sd =  _off_nd[r][kdom*_NoLevels+Level1+1] -  _off_nd[r][kdom*_NoLevels];
                    else if (r == KK)  delta_dof_at_lev_sd = _off_el[VV][kdom*_NoLevels+Level1+1] - _off_el[VV][kdom*_NoLevels+Level1];
                    for (int iki=0; iki < delta_dof_at_lev_sd; iki++) {
                        //is this Level1 or FELevel[r]? No, this must be Level1 because you do not have "nsubdomains*(nlevels+1)" groups.
                        //this is the "OFFSET of the DOFS (q,l,constant,...) WITH RESPECT TO the SUBDOMAINS and the "NUMBER of LEVELS of the QUADRATIC MESH" (Level1)
                        //you see that nothing depends on iki, so it is just a counter
                        //also notice that iki is DIFFERENT from ki !!!
                        //iki goes from zero to the sum of dofs PER SUBDOMAIN and PER LEVEL;
                        // ki goes from zero to the sum of dofs PER LEVEL, so on ALL subdomains
                        for (uint kj=0;kj < _elnodes[VV][c]*_maxelxnode; kj++) {
                            int kdofG = MatG[r][c][_elnodes[VV][c]*_maxelxnode*ki+kj];
                            if (kdofG >= 0)    aux[kdofG]=1;  //aux=1 means that there is a node of the sparsity pattern
                        }

                        for (uint kj=0;kj < _elnodes[VV][c]*_maxelxnode; kj++) {
                            int kdofG = MatG[r][c][_elnodes[VV][c]*_maxelxnode*ki+kj];
                            if (kdofG>=0 &&   aux[kdofG]==1) {
                                aux[kdofG]=0;                      //if you found a node of the sparsity pattern, you do not have to count it anymore the next times
                                MatG[r][c][countG[r][c]]=kdofG;    //WE ARE REWRITING ON TOP OF THE SAME ARRAY!
                                //yes, because countG[r][c] is never gonna be larger than "_elnodes[VV][QQ/*c*/]*_maxelxnode"
                                //then on the next row you will go ahead ATTACHED to the PREVIOUS, and of course you're gonna have
                                //some DOUBLE OCCURRENCES so for sure you will not superimpose.
                                //TODO: this is compatible also with having _elxnode instead of _maxelxnode
                                countG[r][c]++;                    //increase the count of the length of that row
                                if ( kdofG  < (int) dof_extreme_at_lev_sd[c][kdom] || kdofG >= (int) dof_extreme_at_lev_sd[c][kdom+1] )  countoffG[r][c]++;  //increase the count of the off-processor columns
                                //if the dof number is outside the range, then increase
                            }
                        }
                        //here you finish the contributions of ONE ROW: so you summarize the resulting counts
                           lenG[r][c][ki+1] =    countG[r][c];
                        lenoffG[r][c][ki+1] = countoffG[r][c];

                        ki++; //increase the row index for that level (you are inside "r" and "c", so you re-use it everytime from scratch)
                    }
                }
                //=============
                //at this point countG and countoffG give you the length of Mat
                //the arrays "len" have length equal to the number of row dofs
//============================================================
//============== PRINT MatG,lenG,offG TO FILE ================

                XDMFWriter::PrintOneVarMatrixHDF5(name.str(),groupname_lev.str(),n_dofs_lev_fe,countG[r][c],MatG[r][c],lenG[r][c],lenoffG[r][c],r,c,FELevel);

//============================================================
//============ DELETE THE CURRENT c ==========================

                delete [] aux;             //valgrind wants the []
                delete [] lenoffG[r][c];   //valgrind wants the []
                delete [] lenG[r][c];      //valgrind wants the []


            }  //*********************************
            //*************** end c *************
            //*********************************

//============================================================
//============ DELETE THE CURRENT r ==========================
            delete [] lenG[r];
            delete [] lenoffG[r];
            delete [] countG[r];
            delete [] countoffG[r];

        }   //**********************************
        //************** end r ***************
        //**********************************

        for (int fe=0; fe<QL; fe++)   delete [] dof_extreme_at_lev_sd[fe];

        delete [] dof_extreme_at_lev_sd;


        delete [] countG;
        delete [] countoffG;


        for (int r=0; r<QL;r++) {

            for (int c=0; c<QL;c++) {
                delete [] memG[r][c];
                delete [] MatG[r][c];
            }

            delete [] MatG[r]  ;
            delete [] memG[r]  ;

        }

        delete [] MatG       ;
        delete [] memG       ;

        delete [] lenG       ;
        delete [] lenoffG    ;

	H5Gclose(group);

#ifdef DEFAULT_PRINT_TIME
        std::clock_t end_time=std::clock();
        std::cout << " Reference level = " <<  Level1 << " Matrix compute time ="<< double(end_time- start_time) / CLOCKS_PER_SEC << std::endl;
#endif

    }

//==============================================================
//================= END LEVEL LOOP =============================
//==============================================================

// CLOSE THE FILE ============================
      H5Fclose(file);

#ifdef DEFAULT_PRINT_INFO
    std::cout << " GenCase::compute_and_print_MGOps: compute_matrix end \n";
#endif

    return;
}


// ========================================
// In the restrictor and prolongator there is no FE couple,
//there is just the restriction of the same FE type

//FINE LEVEL = COLUMN LEVEL
//COARSE LEVEL = ROW LEVEL

//use _elxnode instead of _maxelxnode would be better. But pay attention,
// because in that case you always have to perform ELEMENT-BASED LOOPS
//and not directly NODE-BASED LOOPS, because if you refer directly to NODES
//you would not be able to DETERMINE to WHAT PROCESSOR THAT NODE BELONGS

//loop on what has been visited so far
//memG massimo e' _elnodes * _maxelxnode
//loop over all the previous times you picked that  //for all the times that we picked that COARSE ROW, we do not have to set again Rest, values, memG
//if we already picked the fine node jcln, then do not set it again
//fine_node_is_set is inside the child nodes loop k_ch
//so, for every child element and child node in that element, you check if the column for that restrictor row
//was already set. If so, do not go ahead setting again
//actually, the range of fine_node_is_set is only INSIDE the if>1.e-6,
//so, only for the fine nodes that contribute to the coarse row,
//you check if they were already set.
//If they were not set, set them
//the point is: being like this, it seems like memG cannot increase to more than the value 1 for each "irow",
//and so also indx_mem cannot increase to more than one
//also jcln, indx0, indx_mem live inside if>1.e-6
//now, the point is: what is the maximum number of nonzero values for the imbedding matrix?
//that value is smething geometrical
//then, since memG[QQ][irow]++ is only inside the  if>1.e-6,
//the maximum increase of memG[QQ][irow] is not dependent on the loops containing irow,
// which would be _elnodes[QQ][VV] * _elnodes[QQ][VV],
//but it is rather (_elnodes[QQ][VV]_of_the_father) * (max number of child nodes that may influence a given coarse node)
//This last number is  smaller than _elnodes[QQ][VV]* _maxelxnode,
//(or maybe we can make it equal to the number of dofs in a QUADRATIC ELEMENT (if we are linear...))...
//no, let us stick around _elnodes[QQ][VV]* _maxelxnode, in an element-based logic.

//------
// if (fabs(val)>1.e-6) {  // "if your fine node (QQ,LL,KK) contributes to the coarse node (QQ,LL,KK)"
// if the coarse shape function is not zero at the fine point;
//so, if the fine node has a nonzero value of the coarse shape function
// then the fine node INFLUENCES the COARSE NODE
// otherwise, you would not just need to update the sparsity pattern


void GenCase::ComputeAndPrintRest(const std::string output_path) {

  int NegativeOneFlag = -1;
  double   PseudoZero = 1.e-8;
  
        std::string f_rest    = DEFAULT_F_REST;
        std::string ext_h5    = DEFAULT_EXT_H5;

        std::ostringstream filename;
        filename << output_path << "/" << f_rest << ext_h5;
        hid_t file = H5Fcreate(filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);

  
    for (int Level1 = 0; Level1 < _NoLevels-1; Level1++) {  //Level1 is the COARSE level (OUTPUT level)

        int Lev_c = Level1;
        int Lev_f = Level1+1;
	
    std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Lev_f << "_" << Lev_c;
    hid_t group = H5Gcreate(file, groupname_lev.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Level1;                                 //COARSE Level for QUADRATIC
        FEXLevel_c[LL] = (_NoLevels + Level1)%(_NoLevels + 1);   //COARSE Level for LINEAR: Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_c[KK] = Level1;                                 //COARSE Level for CONSTANT
        int FEXLevel_f[QL];
        FEXLevel_f[QQ] = Level1+1;                                  //FINE Level for QUADRATIC
        FEXLevel_f[LL] = Level1;                                // AAA look at the symmetry, this is exactly (_n_levels + Level1 + 1)%(_n_levels + 1); ! //FINE Level for LINEAR:   Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_f[KK] = Level1+1;                                  //FINE Level for CONSTANT

        uint** n_dofs_fe_lev;
        n_dofs_fe_lev = new uint*[QL];
        n_dofs_fe_lev[QQ] = _NoNodesXLev;            //equal the pointers
        n_dofs_fe_lev[LL] = _NoNodesXLev;            //equal the pointers
        n_dofs_fe_lev[KK] = _n_elements_vb_lev[VV];  //equal the pointers

        int**     Rest_pos = new    int*[QL];
        double**  Rest_val = new double*[QL];
        int**    mult_cols = new    int*[QL];   //here mult_cols counts, for every dof row (coarse row), the number of FINE NODES that INFLUENCE SOME GIVEN COARSE NODE  //it counts only the nonzero's, doesn't it? yes

        for (int fe=0;fe<QL;fe++) {
	    if (fe < KK)  {
            Rest_pos[fe]   = new    int [n_dofs_fe_lev[fe][FEXLevel_c[fe]]*_elnodes[VV][fe]*_maxelxnode];
            Rest_val[fe]   = new double [n_dofs_fe_lev[fe][FEXLevel_c[fe]]*_elnodes[VV][fe]*_maxelxnode];
            mult_cols[fe]  = new    int [n_dofs_fe_lev[fe][FEXLevel_c[fe]]];
            for (uint i=0; i < n_dofs_fe_lev[fe][FEXLevel_c[fe]]*_elnodes[VV][fe]*_maxelxnode; i++) {
                Rest_pos[fe][i] = NegativeOneFlag;
                Rest_val[fe][i] = 0.;
               }
            for (int im = 0; im < n_dofs_fe_lev[fe][FEXLevel_c[fe]]; im++)     mult_cols[fe][im]=0;

	    }
	    else if (fe == KK)  {  //you dont need excess of space for L2 elements
            Rest_pos[fe]   = new    int [ n_dofs_fe_lev[fe][FEXLevel_c[fe]]*4*(get_dim()-1) /**_elnodes[VV][fe]*/]; //here I have to put the number of childs
            Rest_val[fe]   = new double [ n_dofs_fe_lev[fe][FEXLevel_c[fe]]*4*(get_dim()-1) /**_elnodes[VV][fe]*/];
            mult_cols[fe]  = new    int [ n_dofs_fe_lev[fe][FEXLevel_c[fe]] ];  //TODO unused
            for (uint i=0; i < n_dofs_fe_lev[fe][FEXLevel_c[fe]]*4*(get_dim()-1) /**_elnodes[VV][fe]*/; i++) {
                Rest_pos[fe][i] = NegativeOneFlag;
                Rest_val[fe][i] = 0.;
               }
            for (int im = 0; im < n_dofs_fe_lev[fe][FEXLevel_c[fe]]; im++)     mult_cols[fe][im]=0;
	    }
	    
          }


        for (int pr=0; pr<_NoSubdom; pr++) {
            for (int iel = _off_el[VV][pr*_NoLevels + Lev_c]; iel<_off_el[VV][pr*_NoLevels + Lev_c + 1]; iel++) {  //we are looping over the COARSE ELEMENTS,
                int   el_lib = _el_fm_libm[iel].second;
                int n_childs = _el_sto[el_lib]->_nch;
                for (int i_ch=0; i_ch< n_childs; i_ch++) {
                    int ch_el = _el_sto[el_lib]->_elchs[i_ch];

                    for (int fe=0;fe<QL;fe++) {
                        if (fe < KK)  {
                        for (uint k_fath=0; k_fath < _elnodes[VV][fe]; k_fath++) {

                            int op   = _nd_libm_fm[_el_sto[el_lib]->_elnds[k_fath]]; //coarse node (for rows)
                            int irow = _Qnode_fine_Qnode_lev[FEXLevel_c[fe]][op]; //coarse dof (for rows)    //the COARSE DOF determines the ROW

                            for (uint k_ch = 0; k_ch < _elnodes[VV][fe]; k_ch++) {
                                int ch_op = _nd_libm_fm[_el_sto[ch_el]->_elnds[k_ch]];   //child node (fine node) (for columns)

#if DEFAULT_REST_SIMPLE == 1
//the only contributions from the fine grid to the coarse grid come from the fine nodes which also belong to the coarse grid, very easy.
                                if (ch_op == op) {   //if the fine dof is also a coarse dof
                                    std::cout << "=============== SIMPLE RESTRICTOR, CHECK HERE AND IN THE READING =============" << std::endl;
                                    Rest_pos[fe][irow*_elnodes[VV][fe]*_maxelxnode] = _Qnode_fine_Qnode_lev[FEXLevel_f[fe]][ch_op];   //put the fine level
                                    Rest_val[fe][irow*_elnodes[VV][fe]*_maxelxnode] = 1.;
                                    mult_cols[fe][irow] += 1;
//How it was before
//quadratic simple
//                              Rest_val[fe][irow*_elnodes[VV][fe]*max_elnd] = _elnodes[VV][LL]; //TODO: what is this?
//                                 mult_cols[fe][irow]==1;                                              //TODO: what is this? double equal?
//linear simple
//                                 Rest_val[LL][ g_indexL[ FELevel_c[LL] ][op]*_elnodes[VV][/*LL*/QQ]*max_elnd ] = 1.;    //TODO: what is this?
//                                mult_cols[LL][ g_indexL[ FELevel_c[LL] ][op]]=_elnodes[VV][LL];                   //TODO: what is this?
                                }
#else
                                double val = _feelems[fe]->get_embedding_matrix(i_ch,k_ch,k_fath);   //the value AT THE FINE (CHILD) NODE of the SHAPE FUNCTION of the COARSE (FATHER) NODE
                                if ( fabs(val) > PseudoZero ) {

                                    int jcln     = _Qnode_fine_Qnode_lev[FEXLevel_f[fe]][ch_op];  //put the fine level, because it is the column level
                                    int indx0    = irow*_elnodes[VV][fe]*_maxelxnode;
                                    int indx_mem = mult_cols[fe][irow];   //this is to keep track of how many times we picked that COARSE ROW SO FAR
                                    int fine_node_is_set = 0;
                                    for (int ks=0;ks<indx_mem;ks++) if (Rest_pos[fe][indx0+ks] == jcln) fine_node_is_set=1;

                                    if (fine_node_is_set == 0) {  //if that fine node has not been set in the current row "irow" yet
                                        Rest_pos[fe][indx0+indx_mem] = jcln;
                                        Rest_val[fe][indx0+indx_mem] = val;
                                        mult_cols[fe][irow]++;
                                    }
                                }

#endif

                            } //end k_ch

                        }  //k_fath  //irow changes inside here, so mult_cols[irow] changes inside here, so the maximum value of mult_cols is (range of k_fath) * (range of k_ch) = _elnodes[QQ][VV] * _elnodes[QQ][VV]
                        //no, it is not true!! The maximum value of mult_cols is _elnodes[fe][VV]* _maxelxnode, because that "irow" is picked _maxelxnode times (those are coarse elements), and each of these times
                        //the maximum number of nodes that can influence is _elnodes[fe][VV]
                        //(_elnodes[fe][VV] is here the number of NONZERO EVALUATIONS of the COARSE TEST_ FUNCTION of irow at the FINE NODES belonging to the COARSE ELEMENT that is in the support of irow)
                    
		      
		      
			} //end temp fe
		       else if (fe == KK)  {
			 
			 int sum_previous_sd_c = 0;
                            for (int sd=0;sd<pr;sd++) sum_previous_sd_c += _off_el[VV][sd*_NoLevels+Lev_c+1] - _off_el[VV][sd*_NoLevels+Lev_c];
                            int dof_pos_c = iel - _off_el[VV][pr*_NoLevels+Lev_c] + sum_previous_sd_c;
			 
			 int el_child_fm = _el_libm_fm[ch_el];
			 int sum_previous_sd = 0;
                            for (int sd=0;sd<pr;sd++) sum_previous_sd += _off_el[VV][sd*_NoLevels+Lev_f+1] - _off_el[VV][sd*_NoLevels+Lev_f];
                            int dof_pos_f = el_child_fm - _off_el[VV][pr*_NoLevels+Lev_f] + sum_previous_sd;

                             Rest_pos[fe][ dof_pos_c*4*(get_dim()-1) + i_ch ] = dof_pos_f;
                             Rest_val[fe][ dof_pos_c*4*(get_dim()-1) + i_ch ] = 1;
                            mult_cols[fe][ dof_pos_c ]++;
			 
		       }
		      
		    }  //end fe

                }   //child loop

            }   //iel loop
        }  //end loop pr subdomains

//Dof extremes
//I am only using the QUADRATIC Extremes so far
//if you change the way you PRINT, you also have to change the way you READ
//Now the point is: with dof_extremes we

        int **dof_extreme_at_fine_lev_sd = new int*[QL];
        for (int fe=0;fe<QL;fe++) {
            dof_extreme_at_fine_lev_sd[fe] = new int[_NoSubdom+1];
            dof_extreme_at_fine_lev_sd[fe][0] = 0;
            for (int idom=0;idom<_NoSubdom;idom++)  { 
	           if (fe  < KK)  dof_extreme_at_fine_lev_sd[fe][idom+1] = dof_extreme_at_fine_lev_sd[fe][idom] + ( _off_nd[fe][idom*_NoLevels + Lev_f +1] - _off_nd[fe][idom*_NoLevels]); 
              else if (fe == KK)  dof_extreme_at_fine_lev_sd[fe][idom+1] = dof_extreme_at_fine_lev_sd[fe][idom] + ( _off_el[VV][idom*_NoLevels + Lev_f +1] - _off_el[VV][idom*_NoLevels + Lev_f]);
           }
	}

//==============================
//============= COMPRESSION =================
//==============================

        int ** lenG    = new int*[QL];
        int ** lenoffG = new int*[QL];
        int count[QL];
        int countoff[QL];

        for (int fe=0;fe<QL;fe++) {

            lenG[fe]    = new int [n_dofs_fe_lev[fe][FEXLevel_c[fe]]+1];
            lenoffG[fe] = new int [n_dofs_fe_lev[fe][FEXLevel_c[fe]]+1];
            count[fe]=0;
            lenG[fe][0]=0;
            countoff[fe]=0;
            lenoffG[fe][0]=0;

            int dof_row=0;
            for (int kdom=0; kdom<_NoSubdom; kdom++) {
	      int delta_dof_at_coarse_lev_sd;
	           if (fe  < KK) delta_dof_at_coarse_lev_sd = (int) (_off_nd[fe][kdom*_NoLevels + Lev_c+1] - _off_nd[fe][kdom*_NoLevels]);
              else if (fe == KK) delta_dof_at_coarse_lev_sd =        _off_el[VV][kdom*_NoLevels + Lev_c+1] - _off_el[VV][kdom*_NoLevels + Lev_c];
	    
		for (int iki=0;iki < delta_dof_at_coarse_lev_sd; iki++) {
	           if (fe  < KK) {
                    for (uint k_ch=0;k_ch < mult_cols[fe][dof_row]; k_ch++) {
                        int k_fath=_elnodes[VV][fe]*_maxelxnode*dof_row+k_ch;
                        if (fabs(Rest_val[fe][k_fath])> PseudoZero) {
                            Rest_pos[fe][count[fe]] = Rest_pos[fe][k_fath];
                            Rest_val[fe][count[fe]] = Rest_val[fe][k_fath];
                            count[fe]++;
                            if (Rest_pos[fe][k_fath] >= (int)  dof_extreme_at_fine_lev_sd[fe][kdom+1] ||  Rest_pos[fe][k_fath] < (int) dof_extreme_at_fine_lev_sd[fe][kdom])   countoff[fe]++;
                        }
                    }
		} //end temporary fe
		else if (fe == KK) {
		   for (uint k_ch = 0; k_ch < 4*(get_dim()-1); k_ch++) {
//                   Rest_pos  //dont need compression i think
//                  Rest_val		  
                  // you also dont need mult_cols
                   count[fe]++;
		  //you dont need the countoff
		   }
		}
                       lenG[fe][dof_row+1] =    count[fe];
                    lenoffG[fe][dof_row+1] = countoff[fe];
                    dof_row++;
                }
            }

            XDMFWriter::PrintOneVarMGOperatorHDF5(filename.str(),groupname_lev.str(),n_dofs_fe_lev[fe],count[fe],Rest_pos[fe],Rest_val[fe],lenG[fe],lenoffG[fe],FEXLevel_c[fe],FEXLevel_f[fe],fe);

            delete [] lenG[fe];
            delete [] lenoffG[fe];
            delete [] Rest_pos[fe];
            delete [] Rest_val[fe];
            delete [] mult_cols[fe];
	    delete [] dof_extreme_at_fine_lev_sd[fe];
// 	delete [] n_dofs_lev_fe[fe]; //DO NOT DO THIS DELETE BECAUSE IT IS A COPIED POINTER

        } //end fe


// =======================================
        delete []dof_extreme_at_fine_lev_sd;
        delete [] lenG;
        delete [] lenoffG;
        delete [] Rest_pos;
        delete [] Rest_val;
        delete [] n_dofs_fe_lev;
        delete [] mult_cols;

        H5Gclose(group);

    }
//end Level1

     H5Fclose(file);
	
#ifdef DEFAULT_PRINT_INFO
    std::cout<< " GenCase::compute_and_print_MGOps: compute_rest  \n";
#endif

    return;
}





// ==================================================================
//ok, this doesnt print the femus ordering like i wish, when i am in multiprocessor
//now, let us see.
//The point is that there is some mismatch in the construction of this function.
// In fact, this function receives a BOUNDARY ELEM in the LEVEL NUMERATION
// and it yields a VOLUME ELEMENT in the GLOBAL NUMERATION.
// Since I would prefer doing things by SEPARATE LEVEL, let us see.

// Now my point is: i want to convert from the ABSOLUTE VOLUME NUMBER to the "RELATIVE per LEVEL" VOLUME NUMBER.
//So for every level, I still want to pick the elements of only THAT level but STARTING from ZERO
// 
// 
     //ok now I modified this function consistently, in such a way that for every LEVEL i have an el_bdry_to_vol map
     //that deals with BOUNDARY ELEM numbers and VOLUME ELEM NUMBERS *at THAT LEVEL*
     // so now I already have the  element number IN THAT LEVEL, so that is what I want.

void GenCase::ElemChildToFather() {

    _el_child_to_fath = new int*[_NoLevels];
    for (int lev=0; lev<_NoLevels; lev++)  {
        std::cout << "Level=== " << lev << std::endl;
        _el_child_to_fath[lev] =  new int[_n_elements_vb_lev[BB][lev]];
        int ielem = 0;
        for (int isub=0;isub<_NoSubdom;isub++) {
            std::cout << "Subdomain= " << isub << std::endl;
            for (int iel = _off_el[BB][isub*_NoLevels+lev];
                     iel < _off_el[BB][isub*_NoLevels+lev+1]; iel++) {
                int   el_libm_b = _el_fm_libm_b[iel].second;
//                 std::cout << "boundary elem libmesh" << el_libm_b << std::endl;
                int father_id_libm = _el_sto_b[el_libm_b]->_vol_id;
//                 std::cout << "Volume elem libmesh " << father_id_libm << std::endl;
                int father_id_fm = _el_libm_fm[father_id_libm];
//                 std::cout << "Volume elem fm " << father_id_fm << std::endl;
                int boundary_proc = _el_sto_b[el_libm_b]->_subd;
//                 std::cout << "BOUNDARY PROC " << boundary_proc << std::endl;
                int volume_proc   = _el_sto[father_id_libm]->_subd;
		// father_id_fm is the ABSOLUTE number
		// I want to have the number SO THAT EACH LEVEL RESTARTS FROM ZERO
		// therefore, i want to remove 
		int el_sum_prev_sd_at_lev = 0;
		    for (int sd=0; sd < isub; sd++)  el_sum_prev_sd_at_lev += _off_el[VV][sd*_NoLevels + lev+1] - _off_el[VV][sd*_NoLevels + lev];
		int vol_elem_per_lev = father_id_fm -_off_el[VV][isub*_NoLevels + lev] + el_sum_prev_sd_at_lev;
//                 std::cout << "Volume PROC " << volume_proc  << std::endl;
                _el_child_to_fath[lev][ielem] = vol_elem_per_lev;  //father_id_fm;
//                 std::cout << "ielem pos " << ielem << std::endl;
                //TODO of course _el_sto and _el_sto_b have the Node and Elem numberings still in LIBMESH order
                //given the libmesh order, we must get the FEMuS order
                //TODO is it possible that when you have ONE PROC and ONE SUBDOMAIN
                // the LIBMESH and FEMUS orderings COINCIDE?
                //I guess I will build the INVERSE of the ELEMENT ORDERING
                /*CHECK*/
                if (ielem >= _n_elements_vb_lev[BB][lev]  ) {
                    std::cout << ielem << " Something wrong with number of elems" << std::endl;
                    abort();
                }

                ielem++;

            } //end iel
        } //end subdomain
    }//end lev


    return;
}

// ================================================
// REORDERING ELEMENT (lev>subdomain)
// ================================================
// the goal is to put the elements first by subdomain, and by level inside each subdomain
// in this way you can distinguish them easily if you want to do multigrid AND/OR domain decomposition
// WHAT IS THE WAY TO DISTINGUISH THEM? You can use the sum (L+P*N_L)
//rank_sl 0, rank_sl1: you regroup them: if they have the SAME SUBDOMAIN and SAME LEVEL, you have the SAME SUM (L+P*N_L)!
// Then you SORT them and so you have the FEMUS ordering! (This would be the ordering for CELL dofs!)
// in order to get back to libmesh ordering, just do el_fm_libm[kel].second
//the element info is always taken from the elem_sto file which is in libmesh ordering
//reordering the ELEMENTS is easier than reordering the NODES
// the point is that we have a list of elements of all levels, and a list of the FINE nodes.
// understanding what subd and lev the element have is immediate (->subd ->_lev), so you can group them.
// for the nodes you have to start from the ELEMENT info.
//for every node you have to figure out: what subdomain it belongs to, what level it belongs to, and also ...
// in practice for nodes you have to reason with an auxiliary level.
// if you have N mesh levels, i.e. N quadratic levels,
// then for the dofs you have to consider N levels for the LINEAR part and N levels for the quadratic part.
// Now, since QQ and LL can be interlaced, you can consider N+1 overall levels
// the coarsest one contains only the COARSE LINEAR nodes
//then you have the COARSE QUADRATIC NODES = FINER LINEAR NODES ... and so on
//so, for the ordering we could add a level with a +1...
// PARAVIEW ORDERING: it is the ordering we fix here. For 1 processor you first see the coarse linear nodes,
// but for more processors
//el_fm_libm is a vector of pairs
//we want to make it VB, the only point is El_Sto vs El_Sto_b
// ==================================================================
void GenCase::ReorderElementBySubdLev_VV() {

    _el_fm_libm.resize(_n_elements_sum_levs[VV]);  //every component of this vector is a pair of integers  ==> 3 things
    _el_libm_fm.resize(_n_elements_sum_levs[VV]);     //so far I need the INVERSE for the VOLUME ELEMENT LIST, not for the boundary so far
    //==== Volume ==========
    for (int kel=0;kel<_n_elements_sum_levs[VV];kel++) {
        // counting (level and subdom)
        int isubd = _el_sto[kel]->_subd;
        int  ilev = _el_sto[kel]->_lev;
        int rank_sl = ilev + isubd*_NoLevels;
        _n_elements_sl_vb[VV][rank_sl]++;
        _el_fm_libm[kel].first  = rank_sl ;
        _el_fm_libm[kel].second = kel;                                 //this is the libmesh ordering
    }

    // std reordering on pairs //pairs are reordered wrt their .FIRST element
    std::sort(_el_fm_libm.begin(),_el_fm_libm.end());

    //====== INVERSE Volume Elements ===========
    for (uint i = 0; i < _el_fm_libm.size(); i++)  _el_libm_fm[_el_fm_libm[i].second] = i;


    return;

}



// ==================================================================
void GenCase::ReorderElementBySubdLev_BB() {

    _el_fm_libm_b.resize(_n_elements_sum_levs[BB]);

    for (int kel=0;kel<_n_elements_sum_levs[BB];kel++) {
        int isubdom = _el_sto_b[kel]->_subd;  // _el_sto[_el_sto_b[kel]->_vol_id]->_subd;  //this is the only difference with VV: you get the SUBDOMAIN from the VOLUME!
        int ilev    = _el_sto_b[kel]->_lev;
        int rank_sl = ilev+isubdom*_NoLevels;
        _n_elements_sl_vb[BB][rank_sl]++;
        _el_fm_libm_b[kel].first  = rank_sl;
        _el_fm_libm_b[kel].second = kel;
    }
    // std reordering on pairs
    std::sort(_el_fm_libm_b.begin(),_el_fm_libm_b.end());

    return;

}



// ==================================================================
void GenCase::ReorderNodesBySubdLev() {

    _n_nodes_sl_ql     = new int*[QL_NODES];
    for (int fe=0;fe < QL_NODES; fe++) {
        int offdim2 = _NoSubdom*_NoLevels+1;
        _n_nodes_sl_ql[fe] = new int[offdim2];
        for (int i=0; i<offdim2; i++)      _n_nodes_sl_ql[fe][i]=0;
    }

    _nd_fm_libm.resize(_n_nodes);                        //every component of this vector is a pair of integers (some_rank_to_give_an_order, libmesh_index)
    _nd_libm_fm.resize(_n_nodes);                       //the new numbering in the desired P-L ordering

    // -----------------------------------------------
    // REORDERING NODES (subdom>lev>levP>var)
    // -----------------------------------------------
    //since the mesh is QUADRATIC and we also want to track COARSE LINEAR nodes, it's like thinking of an ADDITIONAL LEVEL
    // 1 + LP + (NL+1)(L + NL*SUBD)
    int n_variables=1;
    int n_nodes_0p=0;
    for (int knode=0;knode<_n_nodes;knode++) {
        _n_nodes_sl_ql[QQ][_NoLevels*_nd_sto[knode]->_subd + _nd_sto[knode]->_lev]++;   //order by quadratic rank
        if (_nd_sto[knode]->_levp <_NoLevels)    _n_nodes_sl_ql[LL][_NoLevels*_nd_sto[knode]->_subd + _nd_sto[knode]->_levp]++;
        if (_nd_sto[knode]->_levp == 0)    n_nodes_0p++;
        for (int klev=_nd_sto[knode]->_lev;klev<_NoLevels;klev++) _NoNodesXLev[klev]++;  // a node contributes to its level and all the finer
        _nd_fm_libm[knode].first= _nd_sto[knode]->_var
                                  + n_variables*(_nd_sto[knode]->_levp + (_NoLevels+1)*(_nd_sto[knode]->_lev + _NoLevels*_nd_sto[knode]->_subd)); //rank
        _nd_fm_libm[knode].second=knode;         //libmesh ordering value
    }
    _NoNodesXLev[_NoLevels]= n_nodes_0p;

    // std reordering on pairs
    sort(_nd_fm_libm.begin(),_nd_fm_libm.end());   //FROM NOW ON, nodes are ordered in this way
    //if i want to retrieve info about a node,  v[knode].second gives me the libmesh node ordering... BUT,
    //the point is that _nd_sto is not used later on
    //so, nd_libm_fm RECEIVES the libmesh ordering and YIELDS the FEMUS ordering
    //while for elements it is not important to store the femus index, with nodes we have to
    // nodes must be ordered some way, they are the dofs!
    //element connectivities can be order whatever you want

    for (uint i = 0; i < _nd_fm_libm.size(); i++)  _nd_libm_fm[_nd_fm_libm[i].second] = i;  //you start numbering from zero   //TODO we could overwrite v[i].first = i, without using nd_libm_fm, we dont need the rank later on


//given this node numbering, you can later fill the coords also

// el_fm_libm[kel].second receives the femus index and yields the libmesh index
// if you want the opposite just do a v_inv_el exactly like nd_libm_fm: this might be good for cell dofs!

//nd_libm_fm[i]       receives the libmesh index and yields the femus index
// v[i].second      receives the femus index and yields the libmesh index


    return;

}

// ============================================
// COMPUTE max number of elements on each quadratic node
// ============================================
// ==================================================================
void GenCase::ComputeMaxElXNode() {

    _elxnode = new int[_n_nodes];

    for (int inode=0;inode<_n_nodes;inode++) _elxnode[inode]=0;
    for (int isub=0;isub<_NoSubdom;isub++)
        for (int i = _off_el[VV][isub*_NoLevels];
                 i < _off_el[VV][isub*_NoLevels + 1];i++)
            for (uint k=0;k<_elnodes[VV][QQ];k++)   {
                _elxnode[_el_sto[_el_fm_libm[i].second]->_elnds[k]]++;
            }

    _maxelxnode=0;
    for (int inode=0;inode<_n_nodes;inode++) if (_elxnode[inode]>_maxelxnode) _maxelxnode=_elxnode[inode];
    std::cout<< " Max number of element on each node = " << _maxelxnode << " \n";

    delete []_elxnode;

    return;
}


// =========================================
void GenCase::ReadMGOps(const std::string output_path, SystemTwo * mysys) {

    std::string     f_matrix = DEFAULT_F_MATRIX;
    std::string       f_rest = DEFAULT_F_REST;
    std::string       f_prol = DEFAULT_F_PROL;
    std::string       ext_h5 = DEFAULT_EXT_H5;

    std::ostringstream filename;
    std::string filename_base;
    filename_base = output_path + "/";
    
        filename.str("");     filename << filename_base << f_rest << ext_h5;
        ReadRest(filename.str(),mysys);
    
        filename.str("");   filename << filename_base << f_matrix << ext_h5;
        ReadMatrix(filename.str(),mysys);
  
        filename.str("");    filename << filename_base << f_prol << ext_h5;
        ReadProl(filename.str(),mysys);

    return;
}



//=============================
//This function depends on _iproc
//remember that a new does not vanish when inside a for loop,
//because the variable was declared outside
// it is NOT like a DECLARATION INSIDE

//The point is that now we must be ready to READ the matrix,
//also the part with CONSTANT ELEMENTS.
//so we have to figure out how to put all the things together.
//basically we need to remember that we have a bunch of
// quantities that may be QQ, LL, KK,
// so eventually we'll have a certain number of QQ, LL and KK
//now, we simply have to decide HOW to ADD THEM in the MATRIX.
// so far, I think we can decide to put
//FIRST all the QQ, then LL, then KK.

// The DofMap was already initialized much earlier than now.
// As a matter of fact, it was needed for setting the BC flags
// It would be nice not to have a lot of places where _Dim is recomputed...

// initDofDim
// initDof
// initVectors

//The GRAPH, or SPARSITY PATTERN, with PETSC functions can be filled JUST WITH ZEROS;
//but with LASPACK  you need do put the DOF VALUES INSIDE.

// so, I got an error that seems to tell me that I basically didn't fill the graph
//correctly. It seems like one processor fills it and the other ones do not,
//half of the lines seem to be unfilled.

///@todo the std::cerr is not redirected to file, think of redirecting it also to file.
// The idea is that every output to terminal must be redirected to file.
// so printf must not be used 

// Every processor will call this function

//Now, the dimension of "graph" is global,
// so it is the dimension of the whole

//TODO: idea about SYNCHRONIZING OUTPUTS.
//I would really like all the outputs OF ALL THE LIBRARIES (HFD5, PETSC, MPI, LIBMESH) involved in FEMuS
//to be REDIRECTED to run_log, altogether.
//What is the best way to make it possible?
//I guess from shell, in principle. The point is how to generate the string for the time for the folder!
// I should use bash functions, make an environment variable, and let all the processes read that environment variable.
// And also, i would like to be able to choose whether to DUPLICATE them, 
// so both to file and to std output, or only to file, or the other cases.

//I need to have the TIME STRING. Bash has the date() function, I want to have it my way
//OK, THAT'S IT:  date +%Y_%m_%d-%H_%M_%S


  //TODO what is this mlinit used for?!?
                      // it is used to determine with WHAT DOF NUMBER the CURRENT PROCESSOR _iproc BEGINS.
                      // NOTICE THAT THE _NODE_DOF MAP IS FILLED COMPLETELY BY ALL PROCESSORS, 
                      // BUT HERE IT IS READ IN PARALLEL BY EACH SEPARATE PROCESSOR !!!
                      //THEREFORE EACH PROCESSOR WILL START WITH ITS OWN INITIAL_DOF AND FINISH WITH ITS OWN FINAL_DOF.
                      //The difference between NODES and ELEMENTS is that for nodes
                      // we keep the offset as the FINEST NODES, because somehow 
                      //we have to keep the correspondence with the NODES OF THE MESH, 
                      //especially when we have to PRINT ON THE MESH
                      //as a matter of fact, we have to bring everything to the FINE MESH
                      //so we can view what happens on the FINE MESH
                      //every processor will fill QUADRATIC, LINEAR and CONSTANT DOFS
                      //being that these dofs belong to a given processor, they will be CONTIGUOUS in the matrix
                      // So you'll have
                      
                      // QQ LL KK || QQ LL KK || QQ LL KK 
                      //< proc0  >||< proc1  >||< proc2 > 
                      
                      // Now, consider the level to be FROZEN, we concentrate on the SUBDOMAINS
                      // Every processor will begin filling the rows at some point.
                      // Every processor has a range of rows, which is mrow_lev_proc_t
                      // Within this range you'll have to put QQ, LL and KK dofs
                      //so within this range you'll have a number of dofs for every fe given by mrow_lev_proc[r],
                      // of course each of them being multiplied by the number of variables of that type
                      // Now, you have to compare how you EXPLORE the NODE_DOF MAP and the MATRIX ROWS
                      //the node_dof map is basically separated into OFFSETS of GEOMETRICAL ENTITIES,
                      //so these ranges are BASED ON THE MESH GEOMETRICAL ENTITIES.
                      // every GEOMETRICAL ENTITY (NODE or ELEMENT) can have ONE OR MORE ASSOCIATED DOFS.
                      //NODES CAN HAVE QUADRATIC DOFS or LINEAR DOFS associated
                      //ELEMENTS CAN HAVE CONSTANT DOFS associated.
                      //so, the node dof is constructed by picking
                      // FIRST NODE QUADRATIC DOFS, for all the quadratic unknowns of the system
                      //  THEN NODE LINEAR DOFS, for all the linear unknowns of the system
                      //  THEN ELEMENT CONSTANT DOFS, for all the constant unknowns of the system
                      // so "QUADRATIC NODES" can be defined as "NODES to which QUADRATIC VARIABLES are associated",
                      //  "LINEAR NODES" means                  "NODES to which LINEAR variables are associated"
                      // ELEMENTS are "GEOMETRICAL ENTITIES to which CONSTANT variables are associated"
                      // So the node dof is divided by 
                      // "GEOMETRICAL ENTITIES for QUADRATIC DOFS" = "QUADRATIC NODES", the length of which must be multiplied by the "NUMBER OF QUADRATIC DOF VARIABLES"
                      // "GEOMETRICAL ENTITIES for LINEAR DOFS"    = "LINEAR NODES",    the length of which must be multiplied by the "NUMBER OF LINEAR DOF VARIABLES"
                      // "GEOMETRICAL ENTITIES for CONSTANT DOFS"  = "ELEMENTS"         the length of which must be multiplied by the "NUMBER OF CONSTANT DOF VARIABLES"
                      //This of course assumes that every GEOMETRICAL ENTITY has the SAME NUMBER of ASSOCIATED DOFS for every FE type.  
                      //So you may say that the node_dof is ordered by GEOMETRICAL ENTITIES first, and WITHIN EACH GEOMETRICAL ENTITY you find ALL THE VARIABLES of THAT FE TYPE.
                      // This is because actually we have a ONE-TO-ONE CORRESPONDENCE between GEOMETRICAL ENTITY and ASSOCIATED DOF.
                      // In a sense, we start with a given DOF family, and, going backwards, we ask ourselves "WHAT IS THE SET OF GEOMETRICAL ENTITIES NATURALLY ASSOCIATED to THIS DOF FE Family?"
                      //For instance, we might have QUAD8 and QUAD9, and associate BOTH OF THEM to a "QUAD9 NODES".In that case we would have EMPTY points,
                      // and it is actually the same that happens when we have the LINEAR DOF FAMILY and we decide to associate the "QUADRATIC MESH NODES" to it for various reasons: multigrid, parallel computing, etc...)
                      //So our current situation is:
                      // QUADRATIC FE DOFS ---> QUADRATIC MESH NODES (at the FINE level)
                      // LINEAR FE DOFS    ---> QUADRATIC MESH NODES (at the FINE level)
                      // CONSTANT FE DOFS  ---> MESH ELEMENTS  (at the CURRENT level)
                      //So, we see that for various reasons we may actually have EMPTY SPACES for two reasons:
                      // 1) BECAUSE we are at some level which is NOT THE FINE
                      // 2) even at the FINE Level, we have LINEAR FINE DOFS built on top of  QUADRATIC FINE NODES
                      
                      //so, _nvars[QQ] is the number of QUADRATIC FE variables
                      
  //NOW, THE POINT IS THIS: FOR EVERY SUBDOMAIN, Do we want to have ALL DOFS CONTIGUOUS or NOT?
  //Yes, of course we want, so in the same subdomain we would have QUADRATIC DOFS, THEN LINEAR DOFS, THEN CONSTANT DOFS.
  //Otherwise, we would have that dofs of the SAME SUBDOMAIN are DISJOINT.
  //Keeping all the dofs of the same subdomain CONTIGUOUS, the equations do not change, they will be the same, 
  // but in this way we will have AS FEW OFF-PROCESSOR COLUMNS as POSSIBLE!
  // So, we do need to BUILD the NODE-DOF MAP APPROPRIATELY for the PARALLEL CASE!
  //Of course, in a SERIAL ENVIRONMENT everything works fine
  
//=======================
//The dimension of the Graph is the total number of rows
// This Graph is the TOTAL GRAPH for ALL THE MATRIX
//The Graph is basically a std::VECTOR of std::VECTORS,
//   therefore you can use the resize() function both on graph and on graph[].
//Also, by the way, with std::vectors you have the overloading of the operator[]
//=====================

// Now, the point is: What is mlinit[]?
// Now, based on the mesh Node and Element numbering, we built the node_dof map.
//Now we have to Read the ONE-VARIABLE-MATRICES for EVERY LEVEL, 
//and USE THEM to BUILD THE WHOLE MATRIX WITH ALL THE VARIABLES and SO ON.
//Now, we already know how the dof map is. So we gave a number to every dof.
// So after doing the NODE ORDERING and ELEMENT ORDERING (VV and BB) in the GENCASE,
// now it is time to use those orderings consequently and have the 
//DOF ORDERING.
//The node map was built in such a way that ALL the DOFS associated to DofObjects (Elems or Nodes)
//of a given subdomain be CONTIGUOUS in the DOF ROW INDEXING.
//Now, given this DOF ORDERING, every PROCESSOR WILL TAKE CARE of ITS OWN DOFS,
// which are contiguous by construction of the node_dof map;
// therefore every processor will pick that block of contiguous row numbers
// and to each of them it will associate the corresponding LENGTH and OFFLENGTH

// Notice how GRAPH has GLOBAL dimension, but SINCE THIS READ MATRIX is PARALLEL,
//Then each processor will fill only PART OF IT. So, actually it makes no sense to make graph be GLOBAL
// and then have each processor only fill one part of it.

//ReadMatrix NEEDS the NODE_DOF to be filled BEFORE
  //Now, we have to understand these steps:
  //1) How we GENERATE the SINGLE DOF VARIABLE sparsity pattern in Gencase;
  //2) How we BUILD the NODE_DOF MAP in PARALLEL based on these SINGLE VARIABLE SPARSITY PATTERNS;
  //3) How we READ the MATRIX SPARSITY PATTERNS and COMBINE THEM WITH the NODE_DOF.
  
  //Basically, you have to COMBINE the "DOFS" (Matrix.h5) with the "MESH" (mesh.h5).
                      
// The node_dof of course takes into account the NUMBER of VARIABLES for each FE type
// So, if the node_dof was not filled correctly, then when you 


void GenCase::ReadMatrix(const  std::string& namefile, SystemTwo * mysys) {

    for (uint Level = 0; Level< mysys->GetGridn(); Level++) {

      std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Level;
  
    uint*** rowcln;    //number of rows and columns for one variable couple
    uint*** length_row; //for every FE row and FE column
    uint*** length_offrow;
    uint*** pos_row;

    rowcln        = new uint**[QL];
    length_row    = new uint**[QL];
    length_offrow = new uint**[QL];
    pos_row       = new uint**[QL];

    hid_t  file = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

    for (int r=0;r<QL;r++) {  //row FE type

        rowcln[r]        = new uint*[QL];
        length_row[r]    = new uint*[QL];
        length_offrow[r] = new uint*[QL];
        pos_row[r]       = new uint*[QL];

        for (int c=0;c<QL;c++) {   //col FE type

            std::ostringstream fe_couple;
            fe_couple << "_F" << r << "_F" << c;

            // dimension
            std::ostringstream dim_name;
            dim_name << groupname_lev.str() << "/" << "DIM" << fe_couple.str();
            rowcln[r][c] = new uint[2];
            XDMFWriter::read_UIhdf5(file,dim_name.str().c_str(),rowcln[r][c]);

            //row length
            std::ostringstream len_name;
            len_name << groupname_lev.str() << "/" << "LEN" << fe_couple.str();
            length_row[r][c]=new uint[ rowcln[r][c][0]+1 ];
            XDMFWriter::read_UIhdf5(file,len_name.str().c_str(),length_row[r][c]);

            // matrix off diagonal
            std::ostringstream offlen_name;
            offlen_name << groupname_lev.str() << "/" << "OFFLEN" << fe_couple.str();
            length_offrow[r][c]=new uint[ rowcln[r][c][0]+1 ];
            XDMFWriter::read_UIhdf5(file,offlen_name.str().c_str(),length_offrow[r][c]);

            // matrix pos //must stay AFTER reading length_offrow
            std::ostringstream pos_name;
            pos_name << groupname_lev.str() << "/" << "POS" << fe_couple.str();
            pos_row[r][c]=new uint[ length_row[r][c][rowcln[r][c][0]] ];
            XDMFWriter::read_UIhdf5(file,pos_name.str().c_str(),pos_row[r][c]);

        } //end col
    } //end row
    
    H5Fclose(file);

//============================================================================
//============ compute things for the sparsity pattern =======================
//============================================================================
    
    int NoLevels  = mysys->GetMLProb().GetMeshTwo()._NoLevels;
    uint off_proc = NoLevels*mysys->GetMLProb().GetMeshTwo()._iproc;

    uint mrow_glob_t = 0;
    for (int fe=0; fe<QL; fe++) mrow_glob_t += mysys->_dofmap._nvars[fe]*rowcln[fe][fe][0];
    uint ncol_glob_t = mrow_glob_t;

    uint mrow_lev_proc_t  = 0;
    for (int fe=0; fe<QL; fe++)  mrow_lev_proc_t +=  mysys->_dofmap._DofLocLevProcFE[Level][mysys->GetMLProb().GetMeshTwo()._iproc][fe]*mysys->_dofmap._nvars[fe];
    uint ncol_lev_proc_t  = mrow_lev_proc_t;

    uint DofObjInit_lev_PrevProcs[QL];  // what is this? it is the ROW INDEX at which to begin for every processor
                      
     for (int r=0; r<QL; r++)     DofObjInit_lev_PrevProcs[r] = 0;
         
    for (uint isubd=0; isubd<mysys->GetMLProb().GetMeshTwo()._iproc; isubd++) {
        DofObjInit_lev_PrevProcs[QQ] += mysys->_dofmap._DofLocLevProcFE[Level][isubd][QQ];
        DofObjInit_lev_PrevProcs[LL] += mysys->_dofmap._DofLocLevProcFE[Level][isubd][LL];
        DofObjInit_lev_PrevProcs[KK] += mysys->_dofmap._DofLocLevProcFE[Level][isubd][KK];
    }    
    
    
    mysys->_LinSolver[Level]->_KK = SparseMatrix::build().release();
// //     mysys->_LinSolver[Level]->_KK->init(_Dim[Level],_Dim[Level], mrow_lev_proc_t, mrow_lev_proc_t); ///@todo BACK TO a REASONABLE INIT

    Graph graph;
    graph.resize(mrow_glob_t);
    graph._m  = mrow_glob_t;
    graph._n  = ncol_glob_t;
    graph._ml = mrow_lev_proc_t;
    graph._nl = ncol_lev_proc_t;
    graph._ml_start = mysys->_dofmap.GetStartDof(Level,off_proc);
    // TODO is this used? I guess it is used by update_sparsity_pattern !
    // Every subdomain has a local set of dofs, and these dofs start at a specific point.
    // Now, remember that mysys->GetMLProb().GetMeshTwo()._off_nd[QQ] should only be used for computing offsets, so differences.
    // Here, it is used ALONE, because it gives you the NODE (in femus ordering) AT WHICH THE CURRENT SUBDOMAIN BEGINS,
    //and then from _node_dof[Level] (which was already constructed) you get THE LOCAL DOF AT THE CURRENT LEVEL TO START FROM.
    // Clearly, pay attention when you add elements, because in that case you would need to REDO the _node_dof map !!!
    //TODO: also, what happens if you have a system with ONLY ELEMENT BASED DOFS?!?
    
    
    int FELevel[QL];
    FELevel[QQ] = Level;
    FELevel[LL] = (Level+mysys->GetGridn())%(mysys->GetGridn()+1); //This is the map for the level of the LINEAR DOFS
    FELevel[KK] = Level;

    int off_onevar[QL];
    off_onevar[QQ] = mysys->GetMLProb().GetMeshTwo()._NoNodesXLev[mysys->GetGridn()-1];
    off_onevar[LL] = mysys->GetMLProb().GetMeshTwo()._NoNodesXLev[mysys->GetGridn()-1];
    off_onevar[KK] = mysys->GetMLProb().GetMeshTwo()._n_elements_vb_lev[VV][Level];
    
    uint  off_EachFEFromStart[QL];
    off_EachFEFromStart[QQ] = 0;
    off_EachFEFromStart[LL] = mysys->_dofmap._nvars[QQ]*off_onevar[QQ];
    off_EachFEFromStart[KK] = mysys->_dofmap._nvars[QQ]*off_onevar[QQ] + mysys->_dofmap._nvars[LL]*off_onevar[LL];

 //==============   
     for (int r=0;r<QL;r++) {
      for (uint ivar=0; ivar < mysys->_dofmap._nvars[r]; ivar++) {

        for (uint DofObj_lev = DofObjInit_lev_PrevProcs[r]; DofObj_lev < DofObjInit_lev_PrevProcs[r] + mysys->_dofmap._DofLocLevProcFE[Level][mysys->GetMLProb().GetMeshTwo()._iproc][r]; DofObj_lev++) {

            int dof_pos, irow;
	         if  (r<KK) {  dof_pos = mysys->GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[FELevel[r]][ DofObj_lev ];  }
            else if (r==KK) {  dof_pos = DofObj_lev; }
                    irow = mysys->_dofmap.GetDof(Level,r,ivar,dof_pos); 

	    int len[QL];   for (int c=0;c<QL;c++) len[c] = 0;
	    for (int c=0;c<QL;c++) len[c] = (length_row[r][c][ DofObj_lev+1 ]-length_row[r][c][ DofObj_lev ]);
            int rowsize = 0; 
	    for (int c=0;c<QL;c++) rowsize +=mysys->_dofmap._nvars[c]*len[c];
	      graph[irow].resize(rowsize + 1);  //There is a +1 because in the last position you memorize the number of offset dofs in that row

            int lenoff[QL];  for (int c=0;c<QL;c++) lenoff[c] = 0;
	    for (int c=0;c<QL;c++)  lenoff[c] = length_offrow[r][c][DofObj_lev+1] - length_offrow[r][c][ DofObj_lev ];
	    int lenoff_size = 0;
	    for (int c=0;c<QL;c++) lenoff_size += mysys->_dofmap._nvars[c]*lenoff[c];
            graph[irow][rowsize] = lenoff_size;      // last stored value is the number of in-matrix nonzero off-diagonal values

        }
    }
  } //end r
     
//===========================
        std::cout << " Matrix \n";
	graph.print();
//===========================

    mysys->_LinSolver[Level]->_KK->update_sparsity_pattern_old(graph);

    //  clean ===============
    graph.clear();

    for (int r=0;r<QL;r++) {  //row FE type

        for (int c=0;c<QL;c++) {   //col

            delete []  rowcln[r][c];
            delete []  length_row[r][c];
            delete []  length_offrow[r][c];
            delete []  pos_row[r][c];

        }

        delete []  rowcln[r];
        delete []  length_row[r];
        delete []  length_offrow[r];
        delete []  pos_row[r];

    }

    delete  []  rowcln;
    delete  []  length_row;
    delete  []  length_offrow;
    delete  []  pos_row;

    
    } //end levels   

#ifdef DEFAULT_PRINT_INFO
    std::cout << " ReadMatrix: matrix reading "  << std::endl;
#endif

    return;
}



//=============================
//This function depends on _iproc
void GenCase::ReadProl(const std::string& name, SystemTwo * mysys) {

  
    vector < SparseMatrix* > &_PP = mysys->GetProjectionMatrix(); //added by Eugenio TO BE TESTED
    
    for (uint Level = 1; Level< mysys->GetGridn(); Level++) {
  
    uint Lev_c = Level-1;
    uint Lev_f = Level;    

        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Level-1;                                 //COARSE Level for QUADRATIC
        FEXLevel_c[LL] = (Level-1+mysys->GetGridn())%(mysys->GetGridn()+1);
        FEXLevel_c[KK] = Level-1;                                 //COARSE Level for CONSTANT //TODO is this used?
        int FEXLevel_f[QL];
        FEXLevel_f[QQ] = Level;                                  //FINE Level for QUADRATIC
        FEXLevel_f[LL] = Level-1;                                // AAA look at the symmetry, this is exactly (_n_levels + Level1 + 1)%(_n_levels + 1); ! //FINE Level for LINEAR:   Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_f[KK] = Level;                                  //FINE Level for CONSTANT //TODO is this used?

    std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Lev_c << "_" << Lev_f;
  
    int** rowcln;
    int** len;
    int** lenoff;
    int** Prol_pos;
    double** Prol_val;

    rowcln     = new int*[QL];
    len        = new int*[QL];
    lenoff     = new int*[QL];
    Prol_pos   = new int*[QL];
    Prol_val   = new double*[QL];

    hid_t  file = H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT); //TODO do I need to open it here for every level?!?

    for (int fe=0; fe<QL; fe++) {

        std::ostringstream fe_family;
        fe_family <<  "_F" << fe;

        //==== DIM ========
        std::ostringstream name0;
        name0 << groupname_lev.str() << "/" << "DIM" << fe_family.str();
        rowcln[fe]= new int[2];
        hid_t dataset = H5Dopen(file, name0.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,rowcln[fe]);
        int nrow=rowcln[fe][0];
        //===== LEN =======
        std::ostringstream name2;
        name2 << groupname_lev.str() << "/" << "LEN" << fe_family.str();
        len[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name2.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len[fe]);
        //==== OFFLEN ========
        std::ostringstream name3;
        name3 << groupname_lev.str() << "/" << "OFFLEN" << fe_family.str();
        lenoff[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name3.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,lenoff[fe]);
        //===== POS =======
        std::ostringstream name1;
        name1 << groupname_lev.str() << "/" << "POS" << fe_family.str();
        uint count3 = len[fe][nrow];
        Prol_pos[fe] = new int[count3];
        dataset = H5Dopen(file,name1.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,Prol_pos[fe]);
        //===== VAL =======
        std::ostringstream name1b;
        name1b << groupname_lev.str() << "/" << "VAL" << fe_family.str();
        Prol_val[fe] = new double[count3];
        dataset = H5Dopen(file,name1b.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,Prol_val[fe]);

    }
   
    H5Fclose(file);

//=================== end reading =====
    
//======= From here on the EQUATION comes into play, because we have to take into account
// the number of variables of every FE type
//Level goes from 1 to NoLevels-1
   
    uint off_proc = mysys->GetMLProb().GetMeshTwo()._iproc*mysys->GetGridn();

    //mysys->
    _PP[Lev_f] = SparseMatrix::build().release();
// // //     _Prl[ Lev_f ]->init(0,0,0,0); //TODO BACK TO A REASONABLE INIT

    // local matrix dimension
    uint ml[QL]; uint nl[QL];
     for (int fe=0; fe<QL; fe++) { 
       if (fe < KK) {
       ml[fe] = mysys->GetMLProb().GetMeshTwo()._off_nd[fe][off_proc + Lev_f +1] - mysys->GetMLProb().GetMeshTwo()._off_nd[fe][off_proc];    //  local quadratic    //COARSE (rows)
       nl[fe] = mysys->GetMLProb().GetMeshTwo()._off_nd[fe][off_proc + Lev_c +1] - mysys->GetMLProb().GetMeshTwo()._off_nd[fe][off_proc];    // global quadratic  //FINE QUADRATIC (cols)
       }
       else if (fe == KK) { 
       ml[fe] = mysys->GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_f +1] - mysys->GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_f];
       nl[fe] = mysys->GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_c +1] - mysys->GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_c];
      }
     }
    // pattern dimension
        int nrowt=0;int nclnt=0;
        for (int fe=0;fe<QL;fe++) {
	  nrowt += mysys->_dofmap._nvars[fe]*rowcln[fe][0];
          nclnt += mysys->_dofmap._nvars[fe]*rowcln[fe][1];
	}
    
    Graph pattern;
    pattern.resize(nrowt);
    pattern._m=nrowt;
    pattern._n=nclnt;
    pattern._ml = 0;
    pattern._nl = 0;
     for (int fe=0;fe<QL;fe++) { 
        pattern._ml += mysys->_dofmap._nvars[fe]*ml[fe];  //  local _m
        pattern._nl += mysys->_dofmap._nvars[fe]*nl[fe];  //  local _n
     }
    uint ml_start = mysys->_dofmap.GetStartDof(Level,off_proc);
    pattern._ml_start = ml_start;

    uint ml_init[QL]; //up to the current processor
      for (int fe=0;fe<QL;fe++) { 
           ml_init[fe]=0;
        for (uint isubd=0;isubd<mysys->GetMLProb().GetMeshTwo()._iproc; isubd++) {
       if (fe < KK)       ml_init[fe] += mysys->GetMLProb().GetMeshTwo()._off_nd[fe][isubd*mysys->GetGridn() + Lev_f +1] - mysys->GetMLProb().GetMeshTwo()._off_nd[fe][isubd*mysys->GetGridn()];
       else if (fe == KK) ml_init[fe] += mysys->GetMLProb().GetMeshTwo()._off_el[VV][isubd*mysys->GetGridn() + Lev_f +1] - mysys->GetMLProb().GetMeshTwo()._off_el[VV][isubd*mysys->GetGridn() + Lev_f];
	}
      }
  
    //============= POSITION =========
   for (int fe=0;fe<QL;fe++) {

     for (uint ivar=0;ivar<mysys->_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = ml_init[fe]; i < ml_init[fe]+ml[fe]; i++) {
          int dof_pos_f;
          if (fe < KK)         dof_pos_f = mysys->GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
          else if (fe == KK)   dof_pos_f = i;
          
            int irow  = mysys->_dofmap.GetDof(Lev_f,fe,ivar,dof_pos_f);

	    uint ncol =    len[fe][i+1] -    len[fe][i];
            uint noff = lenoff[fe][i+1] - lenoff[fe][i];
            pattern[irow].resize(ncol+1);
            pattern[irow][ncol] = noff;
// #ifdef FEMUS_HAVE_LASPACK
            for (uint j=0; j<ncol; j++) {
	      int dof_pos_lev_c = Prol_pos[fe][j+len[fe][i]];
	      int dof_pos_c;
	      if      (fe  < KK) dof_pos_c = mysys->GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ dof_pos_lev_c ];
              else if (fe == KK) dof_pos_c = dof_pos_lev_c; 
	      
	      pattern[irow][j] = mysys->_dofmap.GetDof(Lev_c,fe,ivar,dof_pos_c);
            }
// #endif

        }
	    }//end fe < KK
   } //end fe
   

    std::cout << "Printing Prolongator ===========" << std::endl;
    pattern.print();
    //mysys->
    _PP[Lev_f]->update_sparsity_pattern_old(pattern);

//=========== VALUES ===================
    DenseMatrix *valmat;
    std::vector<uint> tmp(1);
   for (int fe=0; fe<QL; fe++) { 
       for (uint ivar=0;ivar < mysys->_dofmap._nvars[fe];ivar++) {
          for (unsigned int i=ml_init[fe]; i<ml_init[fe] + ml[fe]; i++) {
          int dof_pos_f;
          if (fe < KK)        dof_pos_f = mysys->GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
          else if (fe == KK)  dof_pos_f = i;

          int irow  = mysys->_dofmap.GetDof(Lev_f,fe,ivar,dof_pos_f);

            uint ncol = len[fe][i+1]-len[fe][i];
            tmp[0] = irow;
            std::vector< uint> ind(pattern[irow].size()-1);
            for (uint j=0; j<ind.size(); j++) ind[j] = pattern[irow][j];
            valmat = new DenseMatrix(1,ncol);
            for (uint j=0; j<ncol; j++)(*valmat)(0,j) = Prol_val[fe][j+len[fe][i]];
            //mysys->
	    _PP[Lev_f]->add_matrix(*valmat,tmp,ind);
            delete  valmat;
        }
     }
   }  //end fe


    for (int fe=0;fe<QL;fe++) {
        delete [] Prol_val[fe];
        delete [] Prol_pos[fe];
        delete [] len[fe];
        delete [] lenoff[fe];
    }

    delete [] Prol_val;
    delete [] Prol_pos;
    delete [] len;
    delete [] lenoff;

    pattern.clear();

    //mysys->
    _PP[Lev_f]->close();
//     if (mysys->GetMLProb().GetMeshTwo()._iproc==0) _Prl[  Lev_f ]->print_personal();
//     _Prl[  Lev_f ]->print_graphic(false); //TODO should pass this true or false as a parameter
   } //end levels
    
#ifdef DEFAULT_PRINT_INFO
    std::cout << " ReadProl(B): read Op " << name.c_str() << std::endl;
#endif

    return;
}



// =================================================================
//This function depends on _iproc
//TODO AAA This operator here DEPENDS on the BOUNDARY CONDITIONS
// while the prolongator does not depend on the bc's... CHECK THAT
//Also notice that the Matrix, Restriction and Prolongation sparsity patterns
//are assembled by looping over NODES.
// Then, the values are set with NODE LOOPS for Rest and Prol,
// and with ELEMENT LOOP for the MATRIX (the Assemble Function)
//Level goes from 0 to < GetGridn() - 1 ==> Level is COARSE here
  
//    uint Lev_c = Level;
//   uint Lev_f = Level+1;
	//with these you explore arrays that go from 0  to GetGridn() - 1
                           //so where the distinction between QQ and LL is already made
                           // with the EXTENDED levels you explore things that have an additional level,
                           // and so can work both with QQ and with LL
    //the point is: there are parts where you cannot use extended levels, and parts where you can
    //for instance, in this routine the FINE LEVELS and the FINE EXTENDED LEVELS will both be ok,
    //so we can use them in both cases, but we cannot say the same for the COARSE levels and ext levels
    
    //AAA fai molta attenzione: per esplorare la node_dof devi usare Lev_c e Lev_f,
    //perche' sono legati ai DOF (devi pensare che la questione del mesh e' gia' risolta)
void GenCase::ReadRest(const std::string& name, SystemTwo * mysys) {
 
  vector < SparseMatrix* > &_RR = mysys->GetRestrictionMatrix(); //added by Eugenio TO BE TESTED
  
  for (uint Level = 0; Level< mysys->GetGridn() - 1; Level++) {
    
    uint Lev_c = Level;
    uint Lev_f = Level+1;
    
    std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Lev_f << "_" << Lev_c;
   
    int** rowcln;
    int** len;
    int** lenoff;
    int** Rest_pos;
    double** Rest_val;

    rowcln   = new int*[QL];
    len      = new int*[QL];
    lenoff   = new int*[QL];
    Rest_pos = new int*[QL];
    Rest_val = new double*[QL];

        hid_t  file = H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

	for (int fe=0;fe<QL;fe++) {

        std::ostringstream fe_family;
        fe_family <<  "_F" << fe;

        //==== DIM ========
        std::ostringstream name0;
        name0 << groupname_lev.str() << "/" << "DIM" << fe_family.str();
        rowcln[fe]= new int[2];
        hid_t dataset = H5Dopen(file, name0.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,rowcln[fe]);
        int nrow=rowcln[fe][0];
        //===== LEN =======
        std::ostringstream name2;
        name2 << groupname_lev.str() << "/" << "LEN" << fe_family.str();
        len[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name2.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len[fe]);
        //==== OFFLEN ========
        std::ostringstream name3;
        name3 << groupname_lev.str() << "/" << "OFFLEN" << fe_family.str();
        lenoff[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name3.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,lenoff[fe]);
        //===== POS =======
        std::ostringstream name1;
        name1 << groupname_lev.str() << "/" << "POS" << fe_family.str();
        uint count3 = len[fe][nrow];
        Rest_pos[fe] = new int[count3];
        dataset = H5Dopen(file,name1.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,Rest_pos[fe]);
        //===== VAL =======
        std::ostringstream name1b;
        name1b << groupname_lev.str() << "/" << "VAL" << fe_family.str();
        Rest_val[fe] = new double[count3];
        dataset = H5Dopen(file,name1b.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,Rest_val[fe]);

    }
    
        H5Fclose(file);

//======= From here on the EQUATION comes into play, because we have to take into account
// the number of variables of every FE type
//TODO if you watch carefully, you see that the COARSE levels are EXTENDED,
// while the FINE LEVELS never risk to be "extended", the maximum for them is (mysys->GetGridn()-1) !
    
        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Level;                                 //COARSE Level for QUADRATIC
        FEXLevel_c[LL] = (mysys->GetGridn() + Level)%(mysys->GetGridn() + 1);   //COARSE Level for LINEAR: Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_c[KK] = Level;                                 //COARSE Level for CONSTANT
        int FEXLevel_f[QL];
        FEXLevel_f[QQ] = Level+1;                                  //FINE Level for QUADRATIC
        FEXLevel_f[LL] = Level;                                // AAA look at the symmetry, this is exactly (_n_levels + Level1 + 1)%(_n_levels + 1); ! //FINE Level for LINEAR:   Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_f[KK] = Level+1;                                  //FINE Level for CONSTANT

    uint off_proc=mysys->GetGridn()*mysys->GetMLProb().GetMeshTwo()._iproc;

    //mysys->
    _RR[Lev_c] = SparseMatrix::build().release();
// // //     _Rst[Lev_c]->init(0,0,0,0);   //TODO BACK TO A REASONABLE INIT  //we have to do this before appropriately!!!

    int nrowt=0;int nclnt=0;
        for (int fe=0;fe<QL;fe++) {
	  nrowt += mysys->_dofmap._nvars[fe]*rowcln[fe][0];
          nclnt += mysys->_dofmap._nvars[fe]*rowcln[fe][1];
	}


    Graph pattern;
    pattern.resize(nrowt);
    pattern._m  = nrowt;
    pattern._n  = nclnt;        // global dim _m x _n
    pattern._ml = 0;            //  local _m
    pattern._nl = 0;            //  local _n
     for (int fe=0;fe<QL;fe++) { 
        pattern._ml += mysys->_dofmap._nvars[fe]*mysys->_dofmap._DofLocLevProcFE[Lev_c][mysys->GetMLProb().GetMeshTwo()._iproc][fe]; 
        pattern._nl += mysys->_dofmap._nvars[fe]*mysys->_dofmap._DofLocLevProcFE[Lev_f][mysys->GetMLProb().GetMeshTwo()._iproc][fe];
     } 
     
    // starting indices for local matrix
   uint ml_start = mysys->_dofmap.GetStartDof(Lev_c,off_proc);    //  offset proc nodes
   pattern._ml_start = ml_start;
   uint DofObjInit_lev_PrevProcs_c[QL];
        for (int fe=0;fe<QL;fe++) { DofObjInit_lev_PrevProcs_c[fe] = 0;  }
        for (int fe=0;fe<QL;fe++) { 
    for (uint isubd=0;isubd<mysys->GetMLProb().GetMeshTwo()._iproc; isubd++) { //up to the current processor
       if (fe < KK)       /*mlinit*/ DofObjInit_lev_PrevProcs_c[fe] += mysys->GetMLProb().GetMeshTwo()._off_nd[fe][ isubd*mysys->GetGridn() + Lev_c +1 ]  - mysys->GetMLProb().GetMeshTwo()._off_nd[fe][isubd*mysys->GetGridn()];  
       else if (fe == KK) /*mlinit*/ DofObjInit_lev_PrevProcs_c[fe] += mysys->GetMLProb().GetMeshTwo()._off_el[VV][ isubd*mysys->GetGridn() + Lev_c +1 ]  - mysys->GetMLProb().GetMeshTwo()._off_el[VV][isubd*mysys->GetGridn() + Lev_c];   
         }
     }

    //============= POSITION =========
        for (int fe=0; fe<QL; fe++) { 
    for (uint ivar=0;ivar<mysys->_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + mysys->_dofmap._DofLocLevProcFE[Lev_c][mysys->GetMLProb().GetMeshTwo()._iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = mysys->GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
          else if (fe == KK)   dof_pos_c = i;
	  
            int irow  = mysys->_dofmap.GetDof(Lev_c,fe,ivar,dof_pos_c);

	    uint ncol = len[fe][i+1]-len[fe][i];
            uint noff = lenoff[fe][i+1] - lenoff[fe][i];
            pattern[irow].resize(ncol+1);  //when you do resize, it puts a zero in all positions
	    pattern[irow][ncol] = noff;
// pattern structure (was for laspack only)
            for (uint j=0; j<ncol; j++) {
	      int dof_pos_lev_f = Rest_pos[fe][ j+len[fe][i] ];
	      int dof_pos_f;
    
	      if      (fe  < KK)  dof_pos_f = mysys->GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ dof_pos_lev_f ];
              else if (fe == KK)  dof_pos_f = dof_pos_lev_f;
              pattern[irow][j] = mysys->_dofmap.GetDof(Lev_f,fe,ivar,dof_pos_f);
                }

              }
           }
	} //end fe loop


    std::cout << "Printing Restrictor ===========" << std::endl;
    pattern.print();
    //mysys->
    _RR[Lev_c]->update_sparsity_pattern_old(pattern);  //TODO see 
//         _Rst[Lev_c]->close();
//     if (mysys->GetMLProb().GetMeshTwo()._iproc==0) _Rst[Lev_c]->print_personal(); //there is no print function for rectangular matrices, and print_personal doesnt seem to be working...
// la print stampa il contenuto, ma io voglio solo stampare lo sparsity pattern!
     //Allora cosa faccio: riempio di zeri e poi la stampo! No, di zeri no!!! devi riempirla con qualcos'altro!
//TODO how can I print the sparsity pattern in Petsc BEFORE FILLING the MATRIX ?!?!
    //==============================
//========= SET VALUES =========
//==============================

    DenseMatrix *valmat;
    std::vector<uint> tmp(1);
        for (int fe=0;fe<QL;fe++) {

    for (uint ivar=0;ivar<mysys->_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + mysys->_dofmap._DofLocLevProcFE[Lev_c][mysys->GetMLProb().GetMeshTwo()._iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = mysys->GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
          else if (fe == KK)   dof_pos_c = i;
            int irow     = mysys->_dofmap.GetDof(Lev_c,fe,ivar,dof_pos_c);
            int irow_top = mysys->_dofmap.GetDof(mysys->GetGridn()-1,fe,ivar,dof_pos_c);
            uint ncol = len[fe][i+1]-len[fe][i];
            tmp[0]=irow;
            std::vector< uint> ind(pattern[irow].size()-1);
// 	    std::cout << "\n ==== " << irow << ": ";
            for (uint i1=0;i1<ind.size();i1++) { ind[i1] = pattern[irow][i1]; /*std::cout << " " << ind[i1] << " ";*/}
            valmat = new DenseMatrix(1,ncol);  //TODO add a matrix row by row...
            for (uint j=0; j<ncol; j++) (*valmat)(0,j) = mysys->_bcond._bc[irow_top]*Rest_val[fe][ j+len[fe][i] ];
            //mysys->
	    _RR[Lev_c]->add_matrix(*valmat,tmp,ind);
            delete  valmat;
            }// end dof loop
         } // end var loop
       } //end fe

    for (int fe=0;fe<QL;fe++) {
        delete [] Rest_val[fe];
        delete [] Rest_pos[fe];
        delete [] len[fe];
        delete [] lenoff[fe];
    }

    delete [] Rest_val;
    delete [] Rest_pos;
    delete [] len;
    delete [] lenoff;

    pattern.clear();

    //mysys->
    _RR[Lev_c]->close();
//     if (mysys->GetMLProb().GetMeshTwo()._iproc==0)  _Rst[Lev_c]->print_personal(std::cout);
//     _Rst[Lev_c]->print_graphic(false); // TODO should pass this true or false as a parameter

  } //end levels
  
#ifdef DEFAULT_PRINT_INFO
    std::cout << " ReadRest(B): read Op " << name.c_str() << std::endl;
#endif
    return;
}






} //end namespace femus



