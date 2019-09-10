#include "BoundaryConditions.hpp"


#include "MultiLevelProblem.hpp"
#include "NormTangEnum.hpp"
#include "DofMap.hpp"
#include "SystemTwo.hpp"
#include "NumericVector.hpp"


#include "Box.hpp"


namespace femus {

// ========= ELEM BC AUX ==============
// const int BoundaryConditions::_number_tang_comps[3] = {0,1,3};

     BoundaryConditions::BoundaryConditions(const DofMap* dofmap_in) : _dofmap(dofmap_in)/*,_Dir_pen_fl(0)*/  {   }

         BoundaryConditions::~BoundaryConditions() {

 //=== node
    delete[] _bc;
 //===penalty
//     clearElBc();   /*if (_Dir_pen_fl==1)*/ //DO IT ALWAYS!

};


//=== every processor does ALL because bc is SERIAL
 //Now here we have to think how to impose the boundary conditions for KK
 //Now, the bc_read function orders the boundary condition flags of a QUADRATIC or LINEAR DOF OBJECT.
 //=== The point is that you already have to know the ORDER in which the UNKNOWNS are settled in your system,
 // FIRST QUADRATIC, THEN LINEAR, THEN CONSTANT.
 // as a matter of fact, the bc_field should be ordered like a DOF vector,
 //so have first quadratic, then linear, then constant variables.
 //The fact is that by now we only have ONE bc FIELD, ordered like dofs yes, but ONLY FOR THE *FINE* LEVEL.
 // This is because every coarser NODE is always also FINE, and so that is ok.
 //But, for constant elements I want to make a LEVEL DEPENDENT bc_map,
 //  whose length changes with Level, just like the node_dof map.
 //Then you will have to loop over the BOUNDARY ELEMENT at EACH LEVEL,
 // and assign a FLAG to the VOLUME DOF ELEMENTS.
 // climbing back from Boundary to Volume (we know how to do that)
 // In order to pick the boundary element, we pick its MIDDLE POINT.

 //=== why don't we better impose the boundary conditions
 // based on QUANTITIES and then SCALAR COMPONENTS of QUANTITIES?
 // Then depending on the type of dof you would have to
 // loop over the DOF OBJECTS generating those dofs!
 // So you would loop either over NODES or over ELEMENTS and so on.

 //Wait, we did two functions for QQ-LL and KK, but actually one function seems to be enough!!!
 //The only thing is that you need to distinguish between passing the NODE COORDINATE or the MIDDLE POINT of an ELEMENT,
 //but we can MODIFY the bc_read INTERFACE and MAKE A UNIQUE FUNCTION that ENCOMPASSES EVERYTHING.
 // Yes we can do that.

 //Now, can we make a UNIQUE Function for filling the _bc fields, both for Nodes and for Elements?
 //One of the two is defined only on the FINE level, the other on ALL LEVELS

 //BC is a property of the DOFS. Now, since to every DOF at any level has a corresponding FINE DOF,
 // and for consistency the bc_flag must be THE SAME at ALL LEVELS, then for the nodes
 // we have a UNIQUE FIELD for all levels.
 // Now, of course the DOF NUMBER CHANGES from one level to another.
 // But, the POSITION of the DofObject (= "source" of the dof number) in the _node_dof[Level]
 // is THE SAME!
 //That is why we have a lot of (-1) in the _node_dof map.

 //Ok, now i did a UNIQUE FUNCTION GenerateBdc()
 // The point here is that for every equation I have to set TWO FUNCTIONS, bc_read and BCElemKKRead.
 //One is based on the NODE COORDINATE, one on the MIDDLE POINT of the BOUNDARY ELEMENT
 // So, for every equation i should implement two separate functions bc_read and BCElemKKRead,
 // which i don't really like
 // The problem is that bc_read is called only on the fine level... no problem, I can call it also at the other levels.
 // but, in one case i will only use xp, in the other case i will use only el_xm, so i should pass fake variables in one case or the other
 //As far as i use two SEPARATE FIELDS for the boundary conditions, _bc AND _bc_fe_kk,
 // it is obvious that I have to use two separate functions
 //Then the problem is when you try to switch one Quantity from a NODE BASED FE Family to an ELEM BASED FE FAMILY
 //You have to change the way you  ENFORCE the BOUNDARY CONDITIONS.
 //=== the problem is that here we should write things better.
 // We should write things in such a way that if we set an equation WITHOUT VARIABLES the code still runs.
 // Here this does not happen.
 // if you allocate bc_flag with zero components, then you pass the pointer to the bc_read function which will set the components.
 // So you have to call the bc_read function only if the number of variables is not zero, and in the same way you have to call the function
 // BCElemKKRead only if the number of ElemBased variables is not zero!


 //SWITCHING one VARIABLE from KK to LL.
 // since the length of bc_value and bc_value_kk change, in the routine i cannot leave both selections fixed!

 //Now the point is: I am sorry but we have to modify the structure of the boundary conditions...
 //I think we have to enforce the boundary conditions SEPARATELY for EACH QUANTITY
 //If the quantity is NODE-BASED, then we will use bc_read;
 //If the quantity is ELEMENT BASED, we will use bcelemKKread;
 //Now the point is only this: we want to do this in a UNIQUE routine,
 //so that we write the if's of the geometry only once.
 //I guess first of all I should base everything on the BOUNDARY MIDDLE POINT.

//As a matter of fact I am taking the nodes but these are taken in an element loop,
// so actually i am not looping straight away on the nodes, so i have superimpositions anyway.
//so YES, i'll base anything on the ELEMENT MIDDLE POINT!
//Ok, I did that. Now the list of the bc for each scalar variable is only in one place, under bc_read

//Now that routine is done in such a way that you have to order FIRST QQ, THEN LL, THEN KK;
//In this way you will read correctly

//=== one thing that may be optimized here is that there are computations that could be avoided when you DO NOT HAVE QQ variables, or LL variables, or KK variables.

// The good thing here would be to have a unique bc array at EACH LEVEL.
// The distinction between level is especially good for the elements, because the elements
//only belong to ONE level, it's not like the nodes...

//=== now i am looping over elements first, and nodes inside each element next.
//therefore, if an element puts a zero and a following element puts a one, the one wins.
//so we must find a way which is independent of the element order
//if you loop over the nodes you involve them only once...
//well, there are two reasons for which you superimpose elements:
// 1 - because you loop over the boundary elements
// 2 - because you loop over all faces which have "intersections"

// I think we should do two routines for the flags: one that loops
// over the BOUNDARY ELEMENTS,
// and one that loops over the VOLUME ELEMENTS.

//So far, let us concentrate on the boundary loop. What is the way to check
//not to superimpose the ones on the zeros?
//just put a check before checking

//=== ricorda che i bc sono praticamente dei dof fields,
// quindi in principio possono essere fatti in PARALLELO!
// anche x_old e' un vettore di dof. Alla fine non sono cose
//molto diverse in fondo!

//Quindi bisogna pensare ad una cosa piu' generale.
//Faro' il bc come un Numeric Vector di INTERI,
//cosi' come x_old e' un NumericVector di DOUBLE.

//in questo modo il loop per riempire le bc
// o il loop per riempire x_old
// non saranno molto diversi, anzi molto simili!

//per come e' fatto il mio sistema,
//direi che la cosa giusta e' fare le BC a TUTTI I LIVELLI.
// in questo modo e' tutto molto piu' chiaro.
//la cosa che devo capire un attimo e' come avere
// il valore fisico A TUTTI I LIVELLI.
//Per i nodi e' immediato, perche' un nodo coarser e' anche
// un nodo fine.
// Per un elemento coarser, bisogna fare cosi': ottieni la soluzione
// FINE per gli elementi. Poi la ridistribuisci sugli elementi coarse,
// a tutti i livelli.

//Se volessimo fare le cose ottimizzate,
//dovremmo trattare cosi' lo STORAGE della SOLUZIONE:
//Per i dof il cui dof object e' un NODO,
// mettiamo in un vettore SOLO FINE.

//per i dof il cui dof object e' un ELEMENTO,
//mettiamo in un vettore DI TUTTI I LIVELLI.

//in questo modo minimizzeremmo lo spazio di memoria, che mi sembra
//una cosa giusta, anche se si devono splittare un attimo queste cose.

//il punto e' che se magari facciamo prima QQ, poi KK, poi LL,
// avremmo dei pezzi separati.
//Ma non c'e' problema, perche' noi proprio all'atto della stampa
//facciamo la separazione!
//ora, dobbiamo farla anche durante l'ALGORITMO DI SOLUZIONE.
//Pertanto, definiro' nell'equazione
// un VETTORE di STORAGE NODALE per le VARIABILI NODALI
// e un VETTORE di STORAGE di ELEMENTO per le VARIABILI DI ELEMENTO.

//Questo vettore servira' per RETRIEVE THE DOFS FOR THE TRUE VALUES FOR ANY EQUATION,
// and to PRINT ALL THE VARIABLES, node or cell.

//
// here the boundary conditions are STEADY-STATE right now
void BoundaryConditions::GenerateBdc() {

    const uint Lev_pick_bc_NODE_dof = _dofmap->_mesh._NoLevels - 1;  //we use the FINE Level as reference

 //************************************************
 //******** NODE BASED ****************************
    const uint offset      = _dofmap->_mesh._NoNodesXLev[_dofmap->_mesh._NoLevels-1];

    std::vector<int> bc_flag(_dofmap->_n_vars);

    _bc = new int[_dofmap->_Dim[Lev_pick_bc_NODE_dof]];
    for (uint i1=0;i1< _dofmap->_Dim[Lev_pick_bc_NODE_dof];i1++) _bc[i1] = DEFAULT_BC_FLAG;    // set 1 all the points for  bc (boundary condition)

    //both here and in bc_read you have to put the value bc=0
    //so that you get the IDENTITY OPERATOR and you only have to provide the
    //function for the RHS
 //******** NODE BASED ****************************
 //***************************************************

 //**************************************************
 //******** ELEM BASED ******************************

    int* DofOff_Lev_kk    =  new int[_dofmap->_mesh._NoLevels];

    for (uint Level=0; Level <_dofmap->_mesh._NoLevels; Level++)   { //loop over the levels

          DofOff_Lev_kk[Level] = _dofmap->_nvars[KK]*_dofmap->_DofNumLevFE[Level][KK];

    }
 //******** ELEM BASED ******************************
 //************************************************

 // Both the NODE BASED and the ELEM BASED are initialized to 1 over THE WHOLE VOLUME,
// so over ALL THE DOFS.
// Then you call the BOUNDARY FUNCTIONS, and these act only on boundary nodes or boundary elements,
// but then of course their contribution goes into the volume,
// because every Boundary NODE belongs to the VOLUME
// and every Boundary ELEMENT communicates its value ("as a TRACE") to the corresponding VOLUME ELEM DOF.


    for (uint Level=0; Level <_dofmap->_mesh._NoLevels;Level++)   { //loop over the levels

  Mesh		*mymsh		=  _dofmap->_eqn->GetMLProb()._ml_msh->GetLevel(Level);
  
   unsigned iel0 = mymsh->_elementOffset[ mymsh->processor_id() ]; //WARNING CHANGE MADE WITH NO TESTING!!!!!!
   const uint el_nnodes_b = mymsh->GetElementFaceDofNumber(iel0,ZERO_FACE,BIQUADR_FE);
   
        for (uint isubd=0; isubd<_dofmap->_mesh._NoSubdom; ++isubd) {
            uint iel_b = _dofmap->_mesh._off_el[BB][ _dofmap->_mesh._NoLevels*isubd + Level];
            uint iel_e = _dofmap->_mesh._off_el[BB][ _dofmap->_mesh._NoLevels*isubd + Level+1];

            for (uint iel=0; iel < (iel_e - iel_b); iel++) {

                CurrentElem<double>       currelem(iel,isubd,Level,BB,_dofmap->_eqn,_dofmap->_mesh,_dofmap->_eqn->GetMLProb().GetElemType(),mymsh);

	        currelem.SetDofobjConnCoords();
                currelem.set_elem_center(iel,BIQUADR_FE);

 	    for (uint ivar=0; ivar< _dofmap->_n_vars; ivar++)  bc_flag[ivar] = DEFAULT_BC_FLAG; //this is necessary here to re-clean!

      uint count = 0;
        for (uint i = 0; i < _dofmap->_eqn->GetUnknownQuantitiesVector().size(); i++) {
	  std::vector<int>  bc_temp(_dofmap->_eqn->GetUnknownQuantitiesVector()[i]->_dim,DEFAULT_BC_FLAG);
// 	  _dofmap->_eqn->GetUnknownQuantitiesVector()[i]->bc_flag_txyz(0.,currelem.get_elem_center(),bc_temp);
	  for (uint j = 0; j < _dofmap->_eqn->GetUnknownQuantitiesVector()[i]->_dim; j++) {
	    bc_flag[count] = bc_temp[j];
	    count++;
	  }
	}

  //******************* ONLY FINE LEVEL, NODE VARS *****************
   if (Level == Lev_pick_bc_NODE_dof)  {
                for (uint i=0; i<  el_nnodes_b; i++)  {
                        const uint fine_node = _dofmap->_mesh._el_map[BB][(iel+iel_b)*el_nnodes_b+i];

                    //Set the quadratic fields
                    if ( i < currelem.GetElemType(QQ)->GetNDofs() )
		      for (uint ivar=0; ivar<_dofmap->_nvars[QQ]; ivar++) {
                            int kdof = _dofmap->GetDof(Lev_pick_bc_NODE_dof,QQ,ivar,fine_node);
                           if (_bc[kdof] != 0) _bc[kdof] = bc_flag[ ivar + _dofmap->_VarOff[QQ]];
                        }
                    // Set the linear fields
                    if ( i < currelem.GetElemType(LL)->GetNDofs() ) {
                        for (uint ivar = 0; ivar < _dofmap->_nvars[LL]; ivar++) {
                            int kdof = _dofmap->GetDof(Lev_pick_bc_NODE_dof,LL,ivar,fine_node);
                           if (_bc[kdof] != 0) _bc[kdof] = bc_flag[ ivar + _dofmap->_VarOff[LL]];
                        }
                    }
                } // i

           }
   //******************* END ONLY FINE LEVEL *****************

 //******************* ALL LEVELS, ELEM VARS *****************
		 int sum_elems_prev_sd_at_lev = 0;
                 for (uint pr = 0; pr < isubd; pr++) { sum_elems_prev_sd_at_lev += _dofmap->_mesh._off_el[BB][_dofmap->_mesh._NoLevels*pr + Level + 1] - _dofmap->_mesh._off_el[BB][ _dofmap->_mesh._NoLevels*pr + Level]; }

		 int bdry_iel_lev =  iel + sum_elems_prev_sd_at_lev;
		 int vol_iel =  _dofmap->_mesh._el_bdry_to_vol[Level][bdry_iel_lev];

//******************* END ALL LEVELS, ELEM VARS **************

            } // end of element loop
        }  // end subdomain

     }//end nolevels

        std::cout << "\n GenerateBdc: boundary conditions defined by  bc_read" << "\n \n";

    delete []DofOff_Lev_kk;

    return;
}














//here, we could pass a std::VECTOR of the UNKNOWN INTERNAL QUANTITIES,
//not a map because we want things to be ordered.
//this vector must be passed to the EqnNS constructor so that it can be passed
//to the DA.
//so we must do things OUTSIDE the NS constructor, i.e. in the MultiLevelProblemTwo.
//We cannot do INSIDE the NS constructor because it is called AFTER the DA constructor.
//where Velocity is in the first place and pressure is in the second
//and later one can add OTHER Quantities with
//at this point things are done for only two types,
//quadratic and linear
//in general, suppose that you want linear velocity/linear pressure,
//then you do not use the quadratics
//BUT, we decided to have every shape funcs and stuff
//both for QQ and for LL
//So, we'll leave these like that,
//so that linear temperature can get velocity values from QUADRATIC velocity
//but, actually, we should do that only the things for the INVOLVED ORDERS
// (Internally INVOLVED and Externally INVOLVED) are first allocated HERE
// and THEN FILLED in the GenMatRhs
//this would mean doing a LOOP of "new" and correspondingly a LOOP of "delete"
//for the involved quantities
//clearly, it is important that these "new" and delete things are not done INSIDE
//the Gauss loop but OUTSIDE, a priori
//the problem is that in order to know these we must do a loop
//not only over the Internal Quantities but also the EXTERNAL ones
//because they need to be interpolated as well!
//also, once we know which orders are involved, it is not absolutely immediate
// what DERIVATIVE orders need to be actually filled: phi,dphidxezeta, dphidxyz,dphidxyz3D
//While the FE order depends on the Involved QUANTITIES,
//the DERIVATIVE order depends on the Involved OPERATORS.
//Involved Quants or Ops means BOTH Internally AND Externally!
//therefore, we cannot know.
//One Day, we'll do an Equation as a list of Quantities and a list of Operators acting
//on each of them.
//So the rule will be: we FILL all that we ALLOCATE here.
//ok, so we'll do a vector of internal, and a map of external for now...
//in this way we may write in more general form all the things

//===============related to the SINGLE EQUATION
//=========number of variables
//the number of variables of my own equation must be deduced from the INTERNAL QUANTITIES
//now let us consider TWO internal quantities
//the InternalVector actually belongs to the BASE class,
//but n_vars[] belongs to THIS class
//so, we must consider that we have to turn things into QUADRATIC and LINEAR,
//no matter how the Quantities are.
//Suppose that you have FOUR Quantities, some of them are scalar, some vectors,
//some are linear, some are quadratic, we must turn them into
//nvars[QQ] and nvars[LL].
//Then, when trying to pick the DOFS, we only need to know the offsets
//for passing from the FIRST Quantity (some components and some order)
//to the Second Quantity. Alright.
//so now we''ll just do a TWO by TWO. From TWO Quantities to TWO different FE.
//---> the important thing is not to use nvars_q_in and nvars_l_in any longer!

//---> now, we must understand if the QQ and LL variables used later are
//used in the sense "quadratic or linear - PUBLIC DOMAIN"
//or Quadratic for the FIRST Quantity and LINEAR for the SECOND Quantity
//----so, we have to understand where they are called.
//_fe[0] must be called for being the _fe of quadratic,
//not for being the fe of velocity...
//well, actually i think that all the variables are treated here for being
//quadratic or linear variables















//=======================================

// one day this will be   GetElDofs(myvect,iel);
  //the external quantities might go ahead on their own,
//the only thing we need is the Gauss point value
//The thing is simple: I need the value of the external quantity at a gauss point,
//FOR THE ELEMENT I AM IN
//So, both the SEARCH FOR DOFS and the INTERPOLATION AT THE GAUSS POINT
//can be done INDEPENDENTLY BY EACH QUANTITY!
  //for a given element,
   //every EQUATION and EVERY VARIABLE IN THE EQUATION
   //has ITS OWN _node_dof map!!!
//considering how the dof map is filled now, it works only for QUADRATIC and LINEAR FINITE ELEMENTS
//we should REWRITE the _node_dof map FOR EACH FINITE ELEMENT ASSOCIATED TO THE VARIABLES
//for now, we use the connectivity of el_conn which works for BOTH QUADRATIC AND LINEAR.

// TODO THAT THIS FUNCTION is NOT COMPLETE! If you want to pick the PRESSURE from NS,
// you cant do like that!
// you need to know how the dofmap behaves.
// For now, and throughout the code, there is no dofmap class to handle this thing.
//But, the dofmap is always like this: First the quadratic variables, then the linear variables
//also, el_conn is always quadratic and trhis is not a problem, since
//you only have to loop over the vertices which are the FIRST in the connectivity order.
//so, (idim+ var_off)*dof_off_gl
//myvect can be either quadratic or linear
//if myvect == QQ, varoff=0
//if myvect == LL, varoff=nvars_q. So this depends on the DA equation in which I find myself.
// This is ok, because I have an equation pointer, so i know where the variables lies on.
// TODO For now, we cannot pick only one component, but we can pick
//either ALL the QUADRATIC, or ALL the LINEAR components
//for now all our equations behave in this manner, but later we'll change this thing
//this is something to consider.
//What happens if I have both velocity and pressure that are linear,
//and i want to pick only velocity or Only pressure?
//i want to pick only a PART of the linear variables
//we need this distinction because the _node_dof map is divided into the quadratic and the linear part
//this is because we would have TWO DISTINCT PHYSICAL TYPES of LINEAR VARIABLES
//so, if both variables are linear, you have that lin_off = 0, because nodedof starts with linear stuff;
//but then, while Velocity is the First and it starts right at the beginning, Pressure
//starts right after velocity.
//So, for every DA equation, instead of distinguishing between
//the QUADRATIC variables and the LINEAR variables
//We distinguish between a FIRST Quantity and a SECOND Quantity
//(later we will do just a LIST of MORE Ordered Quantities)
//based on them we have the number of variables nvars[FIRST] and nvars[SECOND] and
//so we can turn them into nvars[QQ] and nvars[LL],
//on which all the loops are based
//So, we'll have that for the FIRST Quantity the offset is ZERO,
//for the second quantity the offset is the number of variables of the First Quantity:
//lin_off = nvars[FIRST]
//all the loops like multigrid operators or stuff are based only on _nvars[0] and _nvars[1],
//so they're not aware if the first is velocity and the second is pressure,
//they just dont care, they work BLINDLY for an ARBITRARY QUANTITY of QUADRATIC and LINEAR variables







//**************************************************
//here, we should explicitly pass the THREE QUANTITIES,
// or better a vector of Quantities;
//the structures of nodedof
// and bc_eldofs
// are both ordered
// first by QUADRATIC and
// then by LINEAR variables.
//But, the Vect we provide are ordered
// by the way we provide them...
//so they may be first quadratic then linear then quadratic...
//we should remember the passage from the initial order
//to the QL order...
//Since later you fill the matrix based on the
//QL order, so you explicitly do what the initvectors
//does at the beginning,
//the point is that you must pay attention to
//the order given at the BEGINNING
//So we have a hybrid reasoning:
// by List of Quantities
// and
// by List of Quadratic and Linear Variables.
//here we are picking things from structures
//that are all ordered by QL variables
//but some of those are distinguished by Quantity...




//now, we'll consider that qty zero and qty one are ordered
//this must be done because in the el_dof_indices and bc_eldofs the nodes are ordered, first quadratic then linear,
//or first linear then linear...
//The KeM and FeM are filled as ORDERED, FIRST QTYZERO("velocity") then QTYONE("pressure")
//so the bc_eldofs and el_dof_indices are correspondingly ORDERED
//therefore also the corresponding global structures are ORDERED
//if you want to write FIRST the PRESSURE then the VELOCITY equations you must "REPOSITION" the element matrix.
//if you do linear linear how can you be sure that velocity comes first?
//I must assure that also the GLOBAL bc and _node_dof start from the qty_zero (either quadratic or linear) and then qty_one


//I'll do the version for a vector of Vect
//when you explore the global structures,
//they all have the structure with quadratic and linear
//distinguished, and ordered, first quadratic then linear.
//So we have to explore these with nvars[QQ] and nvars[LL]
//Now, the point is: when we want to explore the xold
//to retrieve the dof values, we have to do the inverse
//path that allowed us to
//the first quadratic goes to the first quantity
//the second quadratic goes to the second quantity
//the third quadratic goes to the third quantity... and so on...
//in that case you have to consider that you may have
//more variables for a single quantity
//we can say:we loop over the quantities,
// in the order we gave
// if they are quadratic we pick the dofs from the quadratic part, in the order we need;
//if they are linear, we pick the dofs from the linear part.
//then, we must keep in mind the current cursor for both the quadratic and the linear parts
//E' possibile fare in modo che tutto sia costruito sulle quantities?
//beh, in qualche modo, quando riempi i pezzi della matrice di elemento,
// ti servira' una distinzione per componenti... anche se sicuramente avere la quantita'
//tutta insieme ti puo' permettere di scrivere in modo compatto
//certi operatori specie su strutture vettoriali
//direi che e' il caso di separare quello che e' suddiviso per QL
// da quello che e' suddiviso per Quantity

//finora l'idea e': inizializzo per quantita' fisiche,
//poi converto tutto subito in QL,
//poi costruisco tutto per QL

// ad esempio, le condizioni al contorno sono messe tutte insieme,
// e anche le condizioni iniziali sono tutte insieme.
// Se io volessi un modo per riempire le condizioni al contorno
// non per QL ma per Quantita',
//dovrei fare il loop sulle quantita' nell'ordine che e' stato dato,
//poi fare il loop sul numero di componenti di ciascuna quantita',
//e in base a quello mettere una incognita oppure un'altra

//TODO: Now, we are MULTILEVEL. Nevertheless, the boundary conditions are imposed
//only on the FINE LEVEL. Now, the point is that the DofObj for
//Now, the boundary conditions are imposed by NODE COORDINATES,
// lopping over ALL THE NODES, also the VOLUME NODES.
// Then, you distinguish the BOUNDARY NODES from the VOlUME nodes
// by the COORDINATES
// In the same way, we will distinguish the VOLUME ELEMENTS from the
// BOUNDARY elements by the COORDINATES of the CENTER.
// the point is that we are talking about the VOLUME ELEMENTS,
// so we do not know how big they are for putting ifs on the coordinates
// Then, we will loop over volume elements;
// if a volume element HAS a FACE ON THE BOUNDARY,
// and if that face satisfies the "LINE CONDITION" on the boundary,
// then SET the corresponding BC.
// The point is that we need the list of faces for each element,
// and then the coordinates of the nodes of each face.
// so, so far, let me avoid all this stuff and put the boundary conditions
// only on the QQ and LL
// The important thing to remember is that the bc_flag
// is a flag that is associated to EACH DOF,
// in particular to EACH ROW.

//We need to understand what this routine receives and what it gives.
//This routine receives DOF CARRIERS, in terms of NODE DOF CARRIERS and ELEM DOF CARRIERS.
//What does it give? it gives DOF POSITIONS, and DOF BC FLAGS.
//Of course for any DOF CARRIER there may be ONE OR MORE DOFS,
//so one or more DOF POSITIONS, and one or more DOF BC FLAGS
//So, given a "dof carrier" we do not know how many dofs are associated to it.
//Therefore, an alternative would be: starting from the dofs, can we retrieve
//what their dof carriers are?
//Well, the point is that for purpose of domain decomposition
// it is natural to loop over elements, and so it is natural
// to say: what are the DOF CARRIERS in this element?
// Well, on the other hand, one might say: what are the DOFS in this element?
//But, since you are in an element, to retrieve the dofs you necessarily start from
//the GEOMETRICAL ENTITIES that are available to you in there.
//So, from "DOF CARRIER" to DOF is a "GEOMETRICAL approach"
//Instead, from DOF to "DOF CARRIER" would be an "ALGEBRAIC approach"...,
//but I dont think that anybody uses that

//So the idea is that I always provide all the possible DOF CARRIERS,
//and then they may be used or not depending on the existence of variables of each FE type








//============================
//glnode may be a multiple of the node number,
//since node->dof is a non-function map
//either because you have one node - one scalar dof, but vector variables
//or because you have one node - many dofs



///this function scales the passed dof vector and takes into account the boundary conditions as well
///a function like that can be useful also for multiplying/dividing by reference values or the like

void BoundaryConditions::Bc_ScaleDofVec(NumericVector* myvec,  double ScaleFac /*, dimension */ ) {
//only works with Level = NoLevels - 1 , because bc is only on the finest level

  //pass the pointer to the array,the dimension of the array, and the scale factor
//remember that these vectors we pick are made by quadratic and linear parts,
//so we actually know what kind of vector these are => the dimension can be taken from _Dim

//if you have to scale, scale  BOTH vector components, otherwise you could have a non div free vector after the scaling,
//because you scaled only one component.
//Even if you're controlling only with one component at the boundary, both components of the lifting function
//will change due to the variation of a single component at the boundary


for (uint i=0; i < _dofmap->_Dim[_dofmap->_mesh._NoLevels-1]; i++) { //loop over all the dofs, both quadratic and linear

  if (_bc[i] == 1 ) {  //if the dofs are not fixed, scale them

    myvec->set( i, (*myvec)(i)*ScaleFac );

  }

}



 return;
}

//add only where boundary conditions are not fixed
void BoundaryConditions::Bc_AddDofVec(NumericVector* vec_in,NumericVector* vec_out ) {

for (uint i=0; i < _dofmap->_Dim[_dofmap->_mesh._NoLevels-1]; i++) {

    if (_bc[i] == 1 ) {

      vec_out->add (i,(*vec_in)(i));

    }

}

return;

}

void BoundaryConditions::Bc_AddScaleDofVec(NumericVector* vec_in,NumericVector* vec_out,const double ScaleFac ) {
//add a vector multiplied by a constant (only where it is not fixed)

for (uint i=0; i < _dofmap->_Dim[_dofmap->_mesh._NoLevels-1]; i++) {

    if (_bc[i] == 1 ) {

      vec_out->add (i,(*vec_in)(i)*ScaleFac);

    }

}

return;

}




//=======================
//the implementation of these boundary conditions is related to the particular Domain
//Here we are picking a Box
//So, you get the domain name from the Domain. If it is not a box, you abort.
//the imposition of the boundary conditions is related to the Equation.
//Clearly, it depends on the domain
//So for different domains we would have different parts here, with if's.
//We cannot associate this function to the Box or the Cylinder because
//it depends on the OPERATORS involved in the EQUATION,
//so it must stay stick to the Equation, which is a bunch of operators
//every application has only one domain, but if you want to use different
//domains in the same equation you have to specify it here...
//also, changing the domain would mean changing the functions in the Physics User Quantities,
//so in general we do not automatically switch the domain so quickly

//So, for every Domain we have a different implementation
// of this function
//The idea is: i have to get the Box from where i put it.
//The point is that i set it as a domain but it is also a Box
//So i have to do a CAST from Domain to Box



// // // // // void BoundaryConditions::elem_bc_read(const double el_xm[],int& surf_id, double value[],int el_flag[]) const {
// // // // // //el_xm[] is the NON-DIMENSIONAL node coordinate // lb,le are NONDIMENSIONALIZED
// // // // //
// // // // // const double bdry_toll = DEFAULT_BDRY_TOLL;
// // // // //
// // // // //
// // // // //
// // // // // Box* box= static_cast<Box*>(_dofmap->_mesh.GetDomain());
// // // // //
// // // // //   std::vector<double>  lb(_dofmap->_mesh.get_dim());
// // // // //   std::vector<double>  le(_dofmap->_mesh.get_dim());
// // // // //   lb[0] = box->_lb[0];//already nondimensionalized
// // // // //   le[0] = box->_le[0];
// // // // //   lb[1] = box->_lb[1];
// // // // //   le[1] = box->_le[1];
// // // // //
// // // // //  if (_dofmap->_mesh.get_dim() == 3)  {
// // // // //   lb[2] = box->_lb[2];
// // // // //   le[2] = box->_le[2];
// // // // //   }
// // // // //
// // // // //   std::vector<double> x_rotshift(_dofmap->_mesh.get_dim());
// // // // //   _dofmap->_mesh._domain->TransformPointToRef(el_xm,&x_rotshift[0]);
// // // // //
// // // // //
// // // // //   if (_dofmap->_mesh.get_dim() == 2)  {
// // // // //
// // // // //   if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
// // // // // surf_id=44;
// // // // //      el_flag[NN]=1;
// // // // //      el_flag[TT]=1;
// // // // //   value[NN]=0.;
// // // // //   value[TT]=0.;/*-4.*/
// // // // //
// // // // //   }
// // // // //
// // // // //
// // // // //  if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0]) -(x_rotshift[0]) < bdry_toll){ //right
// // // // // surf_id=66;
// // // // //      el_flag[NN]=1;
// // // // //      el_flag[TT]=1;
// // // // //        value[NN]=0.;
// // // // //        value[TT]=0.;/*+4.*/
// // // // // }
// // // // //
// // // // //    if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom
// // // // // surf_id=22;
// // // // //
// // // // //       el_flag[NN]=0;    //no normal component
// // // // //       el_flag[TT]=1;    //yes tangential component
// // // // //         value[NN]=0.;
// // // // //         value[TT]=0.;
// // // // // }
// // // // //
// // // // //   if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
// // // // //  surf_id=88;
// // // // //
// // // // //      el_flag[NN]=0;     //no normal component
// // // // //      el_flag[TT]=1;     //yes tangential component
// // // // //        value[NN]=0.;
// // // // //        value[TT]=0.;
// // // // //
// // // // //   }
// // // // //
// // // // //   }  //end dim 2
// // // // //
// // // // //   else if (_dofmap->_mesh.get_dim() == 3)  {
// // // // //
// // // // //
// // // // //  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //left
// // // // // surf_id=44;
// // // // //      el_flag[NN]=1;    //yes normal component
// // // // //      el_flag[TT]=1;    //yes tang component
// // // // //   value[NN]=0.;  //value of the normal
// // // // //   value[1]=0.;  //value of the tangential
// // // // //   value[2]=0.;  //value of the  tangential
// // // // //   value[3]=0.;  //value of the tangential
// // // // //   }
// // // // //
// // // // //  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) -x_rotshift[0] < bdry_toll){ //right
// // // // // surf_id=66;
// // // // //      el_flag[NN]=1;    //yes normal component
// // // // //      el_flag[TT]=1;    //yes tang component
// // // // //   value[NN]=0.;  //value of the normal
// // // // //   value[1]=0.;  //value of the tangential
// // // // //   value[2]=0.;  //value of the  tangential
// // // // //   value[3]=0.;  //value of the tangential
// // // // //
// // // // // }
// // // // //
// // // // //    if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom
// // // // //
// // // // //    surf_id=22;
// // // // //
// // // // //      el_flag[NN]=0;    //no normal component
// // // // //      el_flag[TT]=1;    //yes tang component
// // // // //   value[NN]=0.;  //value of the normal
// // // // //   value[1]=0.;  //value of the tangential
// // // // //   value[2]=0.;  //value of the  tangential
// // // // //   value[3]=0.;  //value of the tangential
// // // // //   }
// // // // //
// // // // //   if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  { //top
// // // // // surf_id=88;
// // // // //
// // // // //      el_flag[NN]=0;    //no normal component
// // // // //      el_flag[TT]=1;    //yes tang component
// // // // //   value[NN]=0.;  //value of the normal
// // // // //   value[1]=0.;  //value of the tangential
// // // // //   value[2]=0.;  //value of the  tangential
// // // // //   value[3]=0.;  //value of the tangential
// // // // //   }
// // // // //
// // // // //  if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //symmetry
// // // // //
// // // // // surf_id=11;
// // // // //
// // // // //      el_flag[NN]=1;    //yes normal (equal to zero)
// // // // //      el_flag[TT]=0;    //no tangential (symmetry)
// // // // //   value[NN]=0.;  //value of the normal
// // // // //   value[1]=0.;  //value of the tangential
// // // // //   value[2]=0.;  //value of the  tangential
// // // // //   value[3]=0.;  //value of the tangential
// // // // //   }
// // // // //   if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
// // // // // surf_id=77;
// // // // //      el_flag[NN]=1;    //yes normal component //it can be zero also i think
// // // // //      el_flag[TT]=0;    //no tang component
// // // // //   value[NN]=0.;  //value of the normal
// // // // //   value[1]=0.;  //value of the tangential
// // // // //   value[2]=0.;  //value of the tangential
// // // // //   value[3]=0.;  //value of the tangential
// // // // //
// // // // //   }
// // // // //
// // // // //   } //end dim 3
// // // // //
// // // // //
// // // // //   return;
// // // // // }





} //end namespace femus
