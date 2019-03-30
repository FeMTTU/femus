#include "DofMap.hpp"

#include "MultiLevelMeshTwo.hpp"
#include "SystemTwo.hpp"
#include "Quantity.hpp"

namespace femus {
  
  
 DofMap::DofMap(const SystemTwo* eqn_in, const MultiLevelMeshTwo& mesh_in) : _eqn(eqn_in),_mesh(mesh_in) { }

   int DofMap::GetDofQuantityComponent(const uint Level, const Quantity* quantity_in, const uint quantity_ivar,const uint dofobj) const {
    
   int off_previous = 0;
     for (uint i = 0; i < quantity_in->_pos; i++) off_previous += _eqn->GetUnknownQuantitiesVector()[i]->_dim *_DofNumLevFE[ Level ][ _eqn->GetUnknownQuantitiesVector()[i]->_FEord ];

         return GetDofPosIn(Level, dofobj + quantity_ivar*_DofNumLevFE[ Level ][quantity_in->_FEord] + off_previous);
    
  }
  

//=========================
//based on the Qty Internal Vector, i can define HOW MANY variables this equation has
//then compute the total number of variables
void DofMap::initNVars()  {

    for (uint fe = 0; fe < QL; fe++)  _nvars[fe]= 0;

    for (uint i=0;i< _eqn->GetUnknownQuantitiesVector().size(); i++) {
        for (uint fe = 0; fe < QL; fe++)  if ( _eqn->GetUnknownQuantitiesVector()[i]->_FEord == fe) {
                _nvars[fe] +=  _eqn->GetUnknownQuantitiesVector()[i]->_dim;
            }
        }

    
    _n_vars = 0;
    for (uint fe = 0; fe < QL; fe++) _n_vars += _nvars[fe];

    
  for (uint fe = 0; fe < QL; fe++) {
    _VarOff[fe] = 0; 
  for (uint fe2 = 0; fe2 < fe; fe2++) _VarOff[fe] += _nvars[fe2];
  }
  
    
    if (MESH_ORDER == 1 && _nvars[QQ]>0) {
        std::cout << "Can't handle linear mesh with quadratic variables, do a quadratic mesh" << std::endl;
        abort();
    }

    return;
}




 DofMap::~DofMap()  {
   
 //========= DOF MAP ==============================
    for (uint Level = 0; Level < _mesh._NoLevels; Level++)  delete [] _node_dof[ Level];
    delete [] _node_dof;
    for (uint Level = 0; Level < _mesh._NoLevels; Level++)  { delete [] _DofNumLevFE[Level]; delete [] _DofOffLevFE[Level]; } 
    delete [] _DofNumLevFE;
    delete [] _DofOffLevFE;
    delete [] _Dim;           // dimension system Ax=b

 }
 
 


// ============================================================
/// This function initializes the system degrees of freedom (dof)
/// for every Level you have a different dof map,
/// but this is allocated as the fine level, therefore
/// except the fine level the _nodedof is "full of holes"
/// Why dont we allocate the dof map "without holes" ?
/// because we want to refer to the FINE MESH node numbering
/// which is unique.
/// In fact, off_nd_q and off_nd_l give us the RANGES in the FINE NODES
/// The dof map IS NOT THE IDENTITY !
/// node numbering is what is given by the GENCASE;
/// dof numbering instead starts from ZERO for every level!
/// i think it is the identity only on the COARSE LEVEL
/// For every level, you construct a _node_dof[Level]
/// here you start giving an index to the dofs
/// starting SUBDOMAIN BY SUBDOMAIN
/// For each subdomain, you pick FIRST the QUADRATIC then the LINEAR dofs
// void DofMap::InitMeshToDof(const uint Level) {
//let us do the level loop inside here,
//so we can allocate and fill the nodedof
//clearly, if we put a new inside here, we cannot call this function twice,
//because we have only one delete.
//so, we should either call this function in the base constructor or something
//the point is that here the loop is over the SCALAR components already,
//but first it should be over the ordered Quantities
//Pick a quantity, pick its FE order, its number of scalar components
// _off_nd[QQ] gives the node indices corresponding to the quadratic dofs for that lev and subd
//_off_nd[LL] is different: it gives how many linear nodes you have after one quadratic node
//so, to get the node index, you have to start from the quadratic node and add the offset of the linear
// so, _off_nd[QQ] gives ABSOLUTE NODE values
// instead, _off_nd[LL] gives "RELATIVE" node values wrt the quadratic nodes:
// in order to obtain ABSOLUTE nodes you have to START FROM THE ABSOLUTE QUADRATIC and move
// according to the relative offset
// "quadratic nodes" = "nodes which correspond to quadratic dofs"
// "linear nodes" = "nodes which correspond to linear dofs"
// so those that are ABSOLUTE may be used also as RELATIVE when you do DIFFERENCES
// but those that are RELATIVE may only be used as RELATIVE
// notice that for each level, you have to SUBDIVIDE into SUBDOMAINS,
//and within each subdomain you'll have to put ux,uy,uz and p
// la versione "ASSOLUTA" e' in pratica la _node_map
// l'unica differenza e' che la node_map e' su UN LIVELLO IN PIU'
//noi praticamente dobbiamo prendere DUE LIVELLI CONTIGUI e metterli UNO IN FILA ALL'ALTRO
//infatti io guarderei anche alla node_map per riempire questa node_dof,
//perche' secondo me sono la stessa cosa!

//here we have something to consider carefully
//we must understand how to build the nodedof map
//here, we first put all the QUADRATIC VARIABLES
// and then all the LINEAR variables
// - so we see that things here are based on nvars
// all divided. Suppose that you have quadratic
//linear and quadratic, it first puts the
// two quadratic and then the linear one
//this seems to be alright

//Now we have to add the elements as dofs
//InitMeshToDof basically fills the _node_dof map
//now we are gonna do an _elem_dof thing

//From here we see that QUADRATIC DOFS are put first,
// then LINEAR DOFS,
//and then we'll add CONSTANT DOFS.

//Devo pensare se fare elem_dof oppure una unica node_dof,
//che diventera' piuttosto "_geometrical_entity_to_dof"

//count starts from zero for every level and then it grows.

//node dof e' una mappa che tiene conto di TUTTI i PROCESSORI

//Now the point is: Once you introduce ELEMENTS,
//then you're gonna have NODES to which correspond
// BOTH QUADRATIC DOFS
// AND LINEAR DOFS
// and then you have ELEMENTS to which correspond only CONSTANT DOFS.
//So, do we need to INCREASE the OFFSET for the "GEOMETRICAL ENTITIES"?!
//We'll see that.
//Now, actually the biggest number you could have of GEOMETRICAL ENTITIES
//is the NUMBER OF QUADRATIC NODES.
//But let us not think of QUADRATIC NODES as "quadratic".
//Let us think of NODES as MESH GEOMETRICAL ENTITIES.
//We have NODES, we have ELEMENTS,
//those are GEOMETRICAL ENTITIES.
//On the other hand, we have FINITE ELEMENT FAMILIES,
//each of which is characterized by a certain SET OF DEGREES OF FREEDOM,
//so, HOW DO WE MATCH THE TWO?

//===========
// ogni processore inizia ad occupare la node_dof per un pezzo ben preciso,
//un blocco contiguo che pero' e' diverso per ogni processore
//     OFF_ND
//     || ---|---|---  || ---|---|---  ||
//       <--SUBD0----->  <---SUBD1----->
//       <--->           <--->              Level0 - NODES
//       <------->       <------>           Level1 - NODES
//       <------------>  <------------>     Level2 - NODES

//     This was for the NODES, for the ELEMENTS we have sthg like
//     OFF_EL
//     || ---|---|---  || ---|---|---  ||
//       <--SUBD0----->  <---SUBD1---->
//       <--->           <--->              Level0 - ELEMENTS
//           <-->             <-->          Level1 - ELEMENTS
//                <--->           <--->     Level2 - ELEMENTS

//============
// linear -----------------------------------
//why do we put the linears in those positions?
//maybe because we want some correspondence with the mesh
//anyway, the point is that this node_dof is actually full of -1,
//at all levels, both QQ and LL
//the only level that is completely WITHOUT -1 is the QUADRATIC FINEST level.
//now, i dont understand why we don't do something like
//       for (int k1= 0;
//                k1< _mesh._off_nd[LL][off_proc + Level+1]  - _mesh._off_nd[LL][off_proc]; k1++)
//So, why don't we start from the ZERO position? Because we are supposed to put the dof number
//at the CORRECT NODE POSITION; so, if it is a quadratic node, we put it at the quadratic position,
//if it is a linear node, we put it at the linear position,
//if it is an element dof, we put it at the element position, and so on and so on.
//Now, for the elements things may be simpler, because one element belongs to ONLY ONE Level,
//as a matter of fact, you dont really have an ABSOLUTE NUMBERING of ELEMENTS.
//Well, the connectivity gives the numbering of the Elements at each level
//Since we do not have to rely on some kind of Absolute Numbering for the Elements,
//while we have to consider that for the NODES (because we have a well established list of nodes
//on which we have to print all our fields),
//then for the elements we will simply start from zero for every level
// (of course we will go ahead after the previous quadratic and linear dofs)

//So, remember that you have to put the DOF in the position of the corresponding NODE,
//so that you can use this map!
//So since the positions are always with respect to the QUADRATIC MESH,
// you have to set the quadratic nodes, and as a matter of fact, linear nodes 
// are a subset of quadratic nodes

//
// From the fact that we add this guy for the linear dofs _mesh._off_nd[QQ][off_proc]
//I think that it means that the NODES are numbered FIRST QUADRATIC, THEN LINEAR... is that true?
//Wait, some quadratic nodes are also linear already...

//The count variable will change from Level to Level.
//We have that "count" must be equal to _Dim[Level]. 
//If that is not the case, we are missing something


//Now, the point is: Every subdomain fills different pieces of the _node_dof map

//TODO I want to write the three loops for QUADRATIC, LINEAR and CONSTANT ALTOGETHER

//Now, what happens here... 
//NODES are ordered in a unique set, separated by subdomain and by Level
//Somewhat it makes sense to have nodes in a unique set because the same node belongs to DIFFERENT LEVELS.
//It happens that elements are ordered as a unique set too,
// but in this case each Element belongs to ONE and ONLY ONE LEVEL.
// So there is a unique list of element, which we will call ABSOLUTE list,
// and the elements in there have an ABSOLUTE number.
// The ABSOLUTE element number is used to retrieve the connectivity
// Now, for every Level, we want to have the list of elements 
// starting from ZERO. In this way at each level we can fill the node_dof appropriately.
// Of course, you want each subdomain to start at a DIFFERENT POINT,
//So you dont want k1 to start from zero for all processors
// Maybe the best way to do that is:
// 1) to SUM the elements of the previous domains,
// and put that as an OFFSET for starting;
// 2) also, it could be good to read the elem connectivities at SEPARATE levels.
// so, for every level we have to subtract the sum of the elements at the previous levels!

// So, so far, we decide we'll explore the ELEM CONNECTIVITY with the ABSOLUTE Element number,
// and the node_dof map with the "LEVEL-BASED" Element number.
//so we must subtract the sum of the previous levels, 
//which is a function of the current subdomain.

//TODO we must make a CLEAR DISTINCTION between DOF_OBJECTS and DOF_INDICES!!!
// The DOF_OBJECTS are the INPUTS in the node_dof map
// The DOF_INDICES are the OUTPUTS
// Every GEOMETRICAL ENTITY (NODE or ELEMENT) is characterized by MORE DOFS;
// so, in order to make a FUNCTION, we must do the association
//    GEOMETRICAL ENTITY --> DOF_OBJECTS --> DOF_INDICES
// Now the idea is this: for the NODES, everything was built in such a way that
// the DOF_OBJECTS are the SAME AT ALL LEVELS.
//So, if you have a dof_object of some node at some level L, that is going to be the SAME dof_object 
// also at the FINE LEVEL.
// So, dof_object means "INPUT POSITION in the node_dof map"
// Now, for ELEMENTS the situation is different.
// We can choose two possibilities:
// either build a _bc array FOR EVERY LEVEL for the element;
// or building a _bc array ONLY AT THE FINE LEVEL,
// and then know how to convert from the fine level to the CURRENT LEVEL.
// In this second case you would say: for every level, what is the boundary condition
// at my level given the boundary conditions on the fine level?

// ok, i need to print some FIELD at ALL LEVELS
// quindi mi ci vuole una routine che mi LINEARIZZA le CONNETTIVITA' a TUTTI I LIVELLI,
// in modo che io posso sempre vedere, siccome paraview non mi stampa quadratico in 3D.
// dopodiche' la stampa verra' fatta a ciascun livello su ciascun mesh LINEARIZZATO.
// Allora, il punto e' questo: se usi un multigriglia e fai su e giu', allora i vettori intermedi
//non hanno significato perche' contengono valori temporanei, residui, ecc...
//ma se fai la stessa risoluzione a diversi livelli, allora ti puo' interessare stampare indipendentemente 
// a diversi livelli! sia per le bc sia per x_old. Facciamolo.

//=====================================
// The DofObjects here, i.e. the "dof carriers", are NODES and ELEMENTS.
// Assume a QUADRATIC MESH.
// - For QQ FE, all Nodes are dof carriers
// - For LL FE, SOME Nodes are dof carriers (and the others are -1)
// - For KK FE, all Elements are dof carriers

// - For the NODES, we always refer to the FINE numbering (Then, you may always want to use the maps 
// to obtain the LEVEL numberings)
// - For the ELEMENTS, it is different, because any Element belongs to ONLY ONE level.
// So, what you do is you count the number of elements of the previous processors at THAT level,
// and then you go ahead up to the following processor
//======================================

// When do you have a -1 in the node_dof? For the LINEAR nodes, and also for the QUADRATIC nodes at the non-fine level.

//TODO this function is based on the assumption that we order
// FIRST QUADRATIC, THEN LINEAR, THEN CONSTANT DOFS
//I would like to see what happens if we remove this ordering
// (of course all the matrices and MG operators would have consistent ordering)
// also, it affects the way the ELEM MATRIX and RHS are FILLED

// The node dof is used by MATRICES, RECTANGULAR or NOT, and VECTORS,
//so it is basically the fundamental part of the equation

void DofMap::ComputeMeshToDof() {
  
  
    assert(  _n_vars > 0);
    assert(_mesh._NoLevels > 0);

    _Dim      = new uint [ _mesh._NoLevels];

    for (uint Level=0; Level <  _mesh._NoLevels; Level++)   {
      
      int ndof_onevar[QL];
      ndof_onevar[QQ] = _mesh._NoNodesXLev[Level];
      ndof_onevar[LL] = _mesh._NoNodesXLev[ _mesh._NoLevels];
      if (Level>0) ndof_onevar[LL] = _mesh._NoNodesXLev[Level-1];
      ndof_onevar[KK] =  _mesh._n_elements_vb_lev[VV][Level];

      _Dim[Level] = 0;
      for (uint fe=0;fe<QL;fe++) _Dim[Level] += _nvars[fe]*ndof_onevar[fe];

        std::cout << " Level "<< Level
                  << ", Mesh Quadratic Level Nodes " << _mesh._NoNodesXLev[Level]
                  << ", Quadratic DOFs " << _nvars[QQ]*ndof_onevar[QQ]
                  << ", Linear DOFs "    << _nvars[LL]*ndof_onevar[LL]
                  << ", Constant DOFs "  << _nvars[KK]*ndof_onevar[KK]  << std::endl;
		  
    }  
  
  
  _DofLocLevProcFE = new uint**[ _mesh._NoLevels];     
   for (uint Level = 0; Level <  _mesh._NoLevels; Level++) {
       _DofLocLevProcFE[Level] = new uint*[_mesh._NoSubdom];
        for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
            _DofLocLevProcFE[Level][isubdom] = new uint[QL];
            _DofLocLevProcFE[Level][isubdom][QQ] = _mesh._off_nd[QQ][isubdom* _mesh._NoLevels+Level+1] - _mesh._off_nd[QQ][isubdom* _mesh._NoLevels];
            _DofLocLevProcFE[Level][isubdom][LL] = _mesh._off_nd[LL][isubdom* _mesh._NoLevels+Level+1] - _mesh._off_nd[LL][isubdom* _mesh._NoLevels];
            _DofLocLevProcFE[Level][isubdom][KK] = _mesh._off_el[VV][isubdom* _mesh._NoLevels+Level+1] - _mesh._off_el[VV][isubdom* _mesh._NoLevels+Level];
       }
    }
  
  //SETUP DOF OFFSETS for NODE DOF MAP ********************
  //I need to set up the dof offsets for every level
  //Since for nodes they do not depend on fe, but for elements they do, in general they depend on fe
  _DofNumLevFE = new uint*[ _mesh._NoLevels];
  for (uint Level = 0; Level <  _mesh._NoLevels; Level++) {
  _DofNumLevFE[Level] = new uint[QL];
  _DofNumLevFE[Level][QQ] = _mesh._NoNodesXLev[ _mesh._NoLevels-1];
  _DofNumLevFE[Level][LL] = _mesh._NoNodesXLev[ _mesh._NoLevels-1];
  _DofNumLevFE[Level][KK] = _mesh._n_elements_vb_lev[VV][Level];
   }
      
  _DofOffLevFE = new uint*[ _mesh._NoLevels];
  for (uint Level = 0; Level <  _mesh._NoLevels; Level++) {
  _DofOffLevFE[Level] = new uint[QL];
  
 for (int fe=0; fe<QL; fe++) { 
  _DofOffLevFE[Level][fe] = 0;
  for (int fe2=0; fe2 < fe; fe2++) _DofOffLevFE[Level][fe] += _nvars[fe2]*_DofNumLevFE[Level][fe2];
    }
  }

  //fill node dof ******************************************
    _node_dof = new int*[ _mesh._NoLevels];

    for (uint Level = 0; Level <  _mesh._NoLevels; Level++) {
      
       int nodedof_size = 0;
            for (uint fe=0;fe<QL;fe++)      nodedof_size += _nvars[fe]*_DofNumLevFE[Level][fe];
   
        _node_dof[Level] = new int[nodedof_size];
            for (int k1=0;k1< nodedof_size;k1++) _node_dof[Level][k1] = -1;

        uint dof_count_lev=0;
        for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
            uint off_proc=isubdom* _mesh._NoLevels;

            // quadratic -----------------------------------
            for (uint ivar=0;ivar<_nvars[QQ];ivar++) {
                for (int fine_node = _mesh._off_nd[QQ][off_proc];
                         fine_node < _mesh._off_nd[QQ][off_proc+Level+1];fine_node++) {
 #ifndef NDEBUG
 		  if ( fine_node >= (int) _DofNumLevFE[Level][QQ] ) {  std::cout << " Wrong QQ " << std::endl; abort(); }
 #endif	    
                    _node_dof[Level][ fine_node + ivar*_DofNumLevFE[Level][QQ] + _DofOffLevFE[Level][QQ] ] = dof_count_lev;
                    dof_count_lev++;
                }
            }

            // linear -----------------------------------
            for (uint ivar=0;ivar<_nvars[LL];ivar++) {
                for (int fine_node = _mesh._off_nd[QQ][off_proc];
                         fine_node < _mesh._off_nd[QQ][off_proc]
                        + _mesh._off_nd[LL][off_proc + Level+1]
                        - _mesh._off_nd[LL][off_proc]; fine_node++) {
#ifndef NDEBUG
  		  if ( fine_node >= (int) _DofNumLevFE[Level][LL] ) {  std::cout << " Wrong LL " << std::endl; abort(); }
#endif
                   _node_dof[Level][ fine_node + ivar*_DofNumLevFE[Level][LL] + _DofOffLevFE[Level][LL] ] = dof_count_lev;
                    dof_count_lev++;
                }
            }
           
            // constant -----------------------------------
            int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < isubdom; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[VV][pr* _mesh._NoLevels + Level + 1] - _mesh._off_el[VV][pr* _mesh._NoLevels + Level]; }
           
            for (uint ivar = 0; ivar<_nvars[KK]; ivar++) {
                for (uint elem =  sum_elems_prev_sd_at_lev;
		                 elem <  sum_elems_prev_sd_at_lev
		                     + _mesh._off_el[VV][off_proc + Level+1]
                                     - _mesh._off_el[VV][off_proc + Level]; elem++) {
#ifndef NDEBUG
		  if ( elem >= _DofNumLevFE[Level][KK] ) {  std::cout << " Wrong KK " << std::endl; abort(); }
#endif
                    _node_dof[Level][elem + ivar*_DofNumLevFE[Level][KK] + _DofOffLevFE[Level][KK] ] = dof_count_lev;
                    dof_count_lev++;
                }
            }

        }  //end subdomain
        
#ifndef NDEBUG
       if ( dof_count_lev  != _Dim[Level] ) { std::cout << "There is a mismatch between count and _Dim[Level] " << dof_count_lev << " " << _Dim[Level] << std::endl; abort(); }
#endif

        
#ifdef DEFAULT_PRINT_INFO
        std::cout << "DofMap::InitMeshToDof(D)   Level= " << Level  << std::endl;
#endif

    } //end Level
    
    PrintMeshToDof();

    return;
}



// For every Level, we loop over the DOF FE FAMILIES.
// For every DOF FE FAMILY, we associate the GEOMETRICAL ENTITIES on top of which that FE FAMILY is built.
// So, QUAD9 FE is associated to QUADRATIC QUAD9 NODES
// You could also associate QUAD8 FE to QUADRATIC QUAD9 NODES,

//On the other hand, we could think first of the GEOMETRICAL ENTITIES: "QUADRATIC NODES", "LINEAR NODES", "ELEMENTS",
//And loop over the GEOMETRICAL ENTITIES, and ask: WHAT ARE THE DOFS ASSOCIATED TO THESE ENTITIES?
//So, I would better just consider "QUADRATIC NODES" and "ELEMENTS", or "QUAD9 NODES" and "ELEMENTS";
// or, if you have a "QUAD8 MESH", you would have "QUAD8 NODES" and "ELEMENTS".
//Now, what do you do if you have a QUAD9 FE Family? No problem, in your computations you use another DOF,
//which is not directly associated to any Geometrical Entity (in the "computer science" sense that the 9th node does not exists in the file..).
//But, eventually, when you PRINT the QUAD9 ONTO the QUAD8 MESH, you somehow should take into account the presence of the 9TH DOF,
// WHICH SHOULD CONTRIBUTE TO the DOF COEFFICIENTS in the PROJECTED SPACE.

//So, putting the GEOMETRICAL ENTITIES FIRST, "NODES" and "ELEMENTS",
// you may have:
//NODES with QUADRATIC AND LINEAR FE DOFS
//NODES with QUADRATIC FE DOFS
//ELEMENTS with CONSTANT DOFS

//What if you have "GEOMETRICAL ENTITIES" WITH "NO ASSOCIATED DOFS"? They just do not contribute.
// Or "DOFS" with "NO GEOMETRICAL ENTITY", viceversa? You just PROJECT their contribution to the OTHER EXISTING GEOMETRICAL ENTITIES that HAVE DOFS of THAT SAME FE FAMILY.
//So, another example, if you want to PRINT a QUAD4 FE SPACE onto a QUAD9 NODE GEOMETRY, 
//you would be going from a SMALLER to a LARGER space,
// so you would do some PROLONGATION.
//If you went from a LARGER to a SMALLER space,
// you would do some sort of PROJECTION ONTO SUBSPACE.
// Somehow you may see the MG operators in the same way:

// Restriction  = FROM LARGER TO SMALLER SPACE (PROJECTION ONTO SUBSPACE)
// Prolongation = FROM SMALLER TO LARGER SPACE

// Now, let us go back, and let us consider being in the case of looping over the DOF FE FAMILIES, 
// and for each FE Family we pick the CORRESPONDING "DOF GENERATING" GEOMETRICAL ENTITIES.
// The point is: the existence of the DOF should actually INTRINSIC, in the sense that it does not depend
// on the UNDERLYING "DOF GENERATING" GEOMETRICAL ENTITIES.
// You might wanna have QUAD9 FE on a linear mesh, for instance, and do all the computations as quadratic
// and only at the end you project to the LINEAR MESH.

// "GEOMETRICAL ENTITY" (NODES or ELEMENTS) may be also called "DOF GENERATING GEOMETRICAL OBJECT" or "DOF OBJECT"

//The only thing that we are not printing here is the SUBDOMAIN SUBDIVISION

//You should print this map while you are building it, and show precisely
//where you SEPARATE things between one SUBDOMAIN and the other

void DofMap::PrintMeshToDof() const {

  
      for (uint Level = 0; Level<  _mesh._NoLevels; Level++) { 
	
	std::cout << "========== Level " << Level << " =========="  << std::endl;

	for (int fe=0; fe<QL; fe++) {
	  
	  	std::cout << "========== FE " << fe << " =========="  << std::endl;
	  
	for (uint ivar=0; ivar<_nvars[fe]; ivar++) {
	  
  	  	std::cout << "========== VAR " << ivar << " =========="  << std::endl;
		
	     for (uint i=0; i < _DofNumLevFE[Level][fe]; i++) {
  
	  std::cout <<  _node_dof[Level][ i + ivar*_DofNumLevFE[Level][fe] + _DofOffLevFE[Level][fe] ]   << std::endl;
	  
	  //inside here I can also find a way to print at the end of each subdomain, depending on "i"
	  
	     } //end i
	     
	  }
	}
	
      }
  
  
 return; 
}





} //end namespace femus