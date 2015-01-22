// Class includes
#include "IO.hpp"

// std libraries 
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

#include "MeshTwo.hpp"
#include "DofMap.hpp"
#include "SystemTwo.hpp"
#include "FEElemBase.hpp"
#include "NumericVector.hpp"



namespace femus {


namespace IO {
// =============================================================
hid_t read_Dhdf5(hid_t file,const std::string & name,double* data) {

  hid_t  dataset = H5Dopen(file,name.c_str(), H5P_DEFAULT);
  hid_t status=H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t read_Ihdf5(hid_t file,const std::string & name,int* data) {

  hid_t  dataset = H5Dopen(file,name.c_str(), H5P_DEFAULT);
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t read_UIhdf5(hid_t file,const std::string & name,uint* data) {

  hid_t  dataset = H5Dopen(file,name.c_str(), H5P_DEFAULT);
  hid_t status=H5Dread(dataset,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

//=========================================
hid_t print_Dhdf5(hid_t file,const std::string & name, hsize_t dimsf[],double* data) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_DOUBLE,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return status;
}

/// Print int data into dhdf5 file
//TODO can we make a TEMPLATE function that takes either "double" or "int" and uses
//either H5T_NATIVE_DOUBLE or H5T_NATIVE_INT? Everything else is the same
hid_t print_Ihdf5(hid_t file,const std::string & name, hsize_t dimsf[],int* data) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return status;
}

//H5T_NATIVE_UINT
hid_t print_UIhdf5(hid_t file,const std::string & name, hsize_t dimsf[],uint* data) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_UINT,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return status;
}


void PrintXDMFAttribute(std::ofstream& outstream, 
				      std::string hdf5_filename, 
				      std::string hdf5_field,
				      std::string attr_name,
				      std::string attr_type,
				      std::string attr_center,
				      std::string data_type,
				      int data_dim_row,
				      int data_dim_col
 				    ) {

                    outstream << "<Attribute Name=\""<< attr_name << "\" "
                              << " AttributeType=\"" << attr_type << "\" "
			      << " Center=\"" <<  attr_center << "\">\n";
                    outstream << "<DataItem  DataType=\""<< data_type << "\" "
                              << " Precision=\"" << 8 << "\" " 
			      << " Dimensions=\"" << data_dim_row << " " << data_dim_col << "\" "
			      << " Format=\"HDF\">  \n";
                    outstream << hdf5_filename  << ":" << hdf5_field << "\n";
                    outstream << "</DataItem>\n" << "</Attribute>\n";
  
  return;
}


void PrintXDMFTopology(std::ofstream& outfstream,
				     std::string hdf5_file,
				     std::string hdf5_field,
				     std::string top_type,
				     int top_dim,
				     int datadim_n_elems,
				     int datadim_el_nodes
				    ) {
  
    outfstream << "<Topology Type="
               << "\"" <<  top_type  << "\" " 
	       << "Dimensions=\"" << top_dim << "\"> \n";
    outfstream << "<DataStructure DataType= \"Int\""
               << " Dimensions=\" "  <<  datadim_n_elems << " " <<  datadim_el_nodes
	       << "\" Format=\"HDF\">  \n";
    outfstream << hdf5_file << ":" << hdf5_field << " \n";
    outfstream << "</DataStructure> \n" << "</Topology> \n";
   
 return; 
}


void PrintXDMFGeometry(std::ofstream& outfstream,
				     std::string hdf5_file,
				     std::string hdf5_field,
			             std::string coord_lev,
				     std::string geom_type,
				     std::string data_type,
				     int data_dim_one,
				     int data_dim_two) {

    outfstream << "<Geometry Type=\"" << geom_type << "\"> \n";
    for (uint ix=1; ix<4; ix++) {
        outfstream << "<DataStructure DataType=\"" << data_type 
                   << "\" Precision=\"8\" " 
		   << "Dimensions=\"" << data_dim_one << " " << data_dim_two
		   << "\" Format=\"HDF\">  \n";
        outfstream << hdf5_file << ":" << hdf5_field << ix << coord_lev << "\n";
        outfstream << "</DataStructure> \n";
    }
    outfstream << " </Geometry>\n";

  return;
}


// ============================================================
/// This function prints the solution: for quad and linear fem
/// the mesh is assumed to be quadratic, so we must print on quadratic nodes
/// first we print QUADRATIC variables over quadratic nodes, straightforward
/// then, LINEAR variables over quadratic nodes --> we interpolate
/// The NODE numbering of our mesh is so that LINEAR NODES COME FIRST
/// So, to interpolate we do the following:
/// First, we put the values of the LINEAR (coarse) MESH in the FIRST POSITIONS of the sol vector
/// Then, we loop over the elements:
///     we pick the element LINEAR values FROM THE FIRST POSITIONS of the SOL VECTOR itself,
///     we multiply them by the PROLONGATOR,
///     we put the result in the sol vector
/// Therefore the sol values are REWRITTEN every time and OVERWRITTEN because of adjacent elements.
// Ok, allora dobbiamo fare in modo che questa routine stampi Quadratici, Lineari e Costanti
// Ricordiamoci che tutto viene fatto in un mesh "RAFFINATO UNA VOLTA IN PIU' RISPETTO AL LIVELLO PIU' FINE"
// per quanto riguarda gli  elementi, noi abbiamo un valore per ogni elemento,
//e poi dovremo trasferirlo ai suoi figli originati dal raffinamento
// quindi dovremo avere una mappa che per ogni elemento ci da' i suoi FIGLI originati da un raffinamento,
// il tutto ovviamente seguendo gli ordinamenti delle connettivita' di FEMuS.
// Allora, per quanto riguarda gli elementi, viene fatta una mesh_conn_lin che si occupa di mostrare gli 
// elementi del mesh linearizzato
//questa connettivita' linearizzata riguarda soltanto il livello FINE
// finora SIA le VARIABILI QUADRATICHE sia le VARIABILI LINEARI 
//sono stampate sul "MESH FINE QUADRATICO, reso come MESH LINEARIZZATO"
// per le variabili COSTANTI, esse le stampero' sulla CELLA, che ovviamente sara' una CELLA LINEARIZZATA,
// la stessa che uso per il resto, il mio file soluzione .xmf credo che possa utilizzare UNA SOLA TIPOLOGIA.
// Quindi dobbiamo appoggiarci alla stessa topologia "raffinata linearizzata".
// Per fare questo penso che possiamo copiare quello che fa il file case.h5!
    
//The printing of the constant will be very similar to the printing of the PID for the case;
// so we could make a common routine which means "print_cell_property_on_linearized_mesh"

//By the way, observe that with multigrid and with the problems about printing in XDMF format with Paraview
//we are in practice having ONE COARSER LEVEL, for LINEAR COARSE NODES,
// and ONE FINER LEVEL for LINEAR NODES, given by this "auxiliary further mesh refinement"

//TODO Parallel: is this function called only by one proc I hope...
//TODO That information should stay INSIDE a FUNCTION, not outside

//TODO this function works only with Level=NoLevels - 1 !
//facciamola funzionare per ogni livello.
//ora che stampo le connettivita' a tutti i livelli per il mesh LINEARIZZATO, posso fare come voglio.
//Quindi stampo su MSH_CONN di un certo livello
//Ok, dobbiamo distinguere come si esplora la node_dof map e la lunghezza di ogni livello
//
    // ===============================================
//An idea could be: can we make the projection of any FE solution onto this "REFINED LINEARIZED" MESH
// an automatic process, by using the Prol operators?

//Now, the behaviour of this function should be somewhat "parallel" to the construction of the node_dof.

//TODO the group must be REOPENED by EACH EQUATION.. forse Gopen anziche' Gcreate
// I need to have a group that is first created, then closed and reopened without being emptied...
// it should do the same as a File does

//ok, so far we will print the variables with _LEVEL... now we have to fix the wrapper as well

//Ok, now that we make it print things for ALL variables at ALL levels, we need to JUMP on the _node_dof map.
// TODO the NODE DOF NUMBERING is "CONSECUTIVE" ONLY ON THE FINE LEVEL... and only for QUADRATIC variables, 	for (uint i=0;i< n_nodes;i++) 
// which is the same order as the mesh!
//So only on the FINE level you can loop over NODES, in all the other cases 
// you need to loop over ELEMENTS
//For all the other levels we have to JUMP... so we will jump the same way as when we build the node dof!
//Now for all variables, loop over all subdomains, collect all the values, and print them...
// I cannot print them all together

      //TODO TODO TODO here it is more complicated... pos_sol_ivar will sum up to the QUADRATIC NODES,
      //but here we pick the quadratic nodes MORE THAN ONCE, so how can we count them EXACTLY?!!
      //the problem is that when you are not on the fine level you have jumps, so, in order to count the LINEAR nodes,
      //you need a flag for every quadratic node that tells you if it is also linear or not...
      //we should either build a flag field in advance, so that we don't do it now, or maybe later is ok
      //since we don't keep any extra info about each node (we do not have a Node class, or an Elem class),
      //now it's time we have some flags...
      //i guess i should loop over all the quadratic nodes, and let the flag for the linear only
      // ok so when i pick the linear nodes i may set the flag is linear there
      //flag = 1 means: it is linear, otherwise 0

// // //       for (int i = mesh->_off_nd[QQ][off_proc];
// // // 	       i < mesh->_off_nd[QQ][off_proc+Level+1];i++) {
// // //       if (i< mesh->_off_nd[QQ][off_proc]
// // // 	   + mesh->_off_nd[LL][off_proc + Level+1 ]
// // // 	   - mesh->_off_nd[LL][off_proc])
// // //         
// // // 	 }
  //now, in the quadratic nodes, what are the positions of the linear nodes?    
  //if I am not wrong, the first nodes in the COARSE LEVEL in every processor
  //are the LINEAR ONES, that's why the offsets are what they are...
  //in fact, we are looping on the GEomElObjects of each level.

//     int* flag_is_linear = new int[n_nodes]; // we are isolating MESH NODES correspoding to LINEAR DOFS
    //what, but we already know from the node numbering what nodes are linear, because we distinguished them, right?!
//     for (uint i=0;i< n_nodes;i++) flag_is_linear[i] = 0;

//LINEAR VARIABLES
//these are re-dimensionalized //so you dont need to multiply below! //the offset for _node_dof is always quadratic, clearly, because the mesh is quadratic
//        //I am filling two QUADRATIC arrays first by the LINEAR POSITIONS
	  //so pay attention that if you are not setting to ZERO for the linear case, but exactly replacing at the required points

//Always remember that in order to pick the DofValues you have to provide the DofObjects in the correct manner
	  
// QUADRATIC VARIABLES  
// pos_in_mesh_obj gives me the position of the GEomElObject: in fact the quadratic dofs are built on the quadratic GeomEls exactly in this order
//here we are picking the NODES per subd and level, so we are sure don't pass MORE TIMES on the SAME NODE

// PRoblem with the linear variables in write_system_solutions and PrintBc. 
//There is one line that brings to mistake, but TWO different mistakes.
// in PrintBc it seems to be related to HDF5;
// in write_system_solutions it seems to concern PETSC
//so that is the wrong line, if i comment it everything seems to work for any processor.
//with two and three processors it seems to give even different errors...
//when there is an error related to HDF5, the stack has the _start thing...
// with two procs there is an error related to Petsc,
// with three procs there is an error related to HDF5...

//If I only use PrintBc and not write_system_solutions, 
//both with 2 and 3 processors the errors are related to HDF5... this is so absolutely weird...

//Now it seems like I am restricted to that line. That line is responsible for the error.
//Then, what is wrong with that? I would say: I am doing sthg wrong on the array positions.
//Let me put a control on the position. maybe the _el_map is somehow ruined

//Quando hai un segmentation fault, devi concentrarti non tanto sui VALORI ASSEGNATI
// quanto sugli INDICI DEGLI ARRAY !!!

        // to do the interpolation over the fine mesh you loop over the ELEMENTS
        //So of course you will pass through most nodes MORE THAN ONCE
        //after setting correctly the linear nodes, then the quadratic ones are done by projection, without any problem
        //for every element, I take the linear nodes.
        //Then I loop over the quadratic nodes, and the value at each of them is the prolongation
        // of the linear values.
        //So, for every quadratic node, I compute the value and set the value 
        //now, the connectivity is that of the quadratic mesh, but with respect
        //to the FINE node numbering
        //now we have to convert from the fine node numbering to the node numbering at Level!!!
        //do we already have something to do this in the MESH?!?
        //given a quadratic node in FINE NUMBERING, obtained by an element AT LEVEL L,
        //can we obtain the position of that node AT LEVEL L
        //don't we have the connectivities AT LEVEL L,
        //where the numbering goes from ZERO to n_nodes_level?
        //I would say it is exactly the _node_dof FOR THE FIRST QUADRATIC VARIABLE!
        //TODO OF COURSE THAT WOULD IMPLY THAT AT LEAST ONE QUADRATIC VARIABLE is BUILT
        //TODO Also, I have to check that COARSER GEOMETRIES are ASSOCIATED to COARSER TOPOLOGIES!
        // _node_dof in serial is the IDENTITY, BUT NOT IN PARALLEL!!
        //Per passare dal Qnode in FINE NUMBERING al Qnode in SERIAL NUMBERING
        //non basta passarlo alla _node_dof del livello, perche' quando sei in parallelo
        //lui conta prima i dof quadratici, poi quelli lineari, poi quelli costanti, e bla bla bla...
        //e quindi bisogna cambiare il modo di pigliarli!
        //Si' bisogna usare qualcosa di SIMILE alla NODE_DOF, ma NON la node_dof,
        //perche' quella, per dato livello, conta i dof QQ,LL,KK del proc0, poi QQ,LL,KK del proc1, and so on...
        //OK OK OK! Direi che quello che dobbiamo usare e' proprio la node_map, quella che leggiamo dal mesh.h5!!!
        //Quindi _node_dof e' per TUTTI i DOF,
        // mentre _node_map e' solo per i NODI GEOMETRICI!
        //mi sa pero' che la _node_map fa l'opposto di quello che vogliamo noi, cioe'
        //dato il nodo di un certo livello ti restituisce il nodo FINE
        //Noi invece abbiamo il NODO FINE, e vogliamo avere il nodo di un certo livello.
        //Questo si otterrebbe leggendo TUTTA la map e non comprimendola come facciamo...
        //oppure io la ricostruisco rifacendo il loop della node_dof ma solo per UNA variabile quadratica!
        //praticamente, MA NON DEL TUTTO!, la _node_map e' l'inverso della _node_dof!!!
        //TODO ma non c'e' nessun altro posto nel codice in cui devi passare 
        //da QNODI FINI a QNODI di LIVELLO?
        // sembra di no, piu' che altro devi passare da QNODI FINI a DOF di LIVELLO!
        //ora qui, essendo che dobbiamo stampare su griglie di diverso livello,
        //e siccome le CONNETTIVITA' di TUTTI I LIVELLI SONO DATE RISPETTO all'unico FINE NUMBERING
        // (TODO beh, volendo potrei utilizzare le connettivita' "di livello" che stampo nel file msh_conn_lin...
        //Il fatto e' che quelle non sono variabile di classe (in Mesh) e quindi mi limito a calcolare, stampare e distruggere...)
        // ALLORA DEVO FARE IL PASSAGGIO da QNODI FINI a QNODI DI LIVELLO.
        //Questo me lo costruisco io domattina.
        //Siccome e' gia' stato calcolato nel gencase, allora mi conviene evitare di fare questo calcolo.
        //Quando leggo il vettore dal gencase, posso costruire un vettore di PAIRS...
        //ovviamente e' piu' lento, perche' deve cercare un elemento con degli if...
        //allora faccio un array 2x, per cui calcolo io l'indice da estrarre...

// We have to find the LINEAR NODES positions in the QUADRATIC LIST
// the QUADRATIC list at each level is based on 
// the "LINEARIZED" CONNECTIVITIES of that level
//    plus the LIST OF COORDINATES of that level
	    
//Now there is another mistake, in printing the LINEAR CONNECTIVITIES!
//The fact is always that _el_map gives the Qnodes in FINE NUMBERING,
//while I want the Qnode numbering AT EACH LEVEL.
//because at coarser levels I have fewer nodes, so the list of coords is shorter!

//Ok, the printing of the linearized connectivities is WRONG at NON-FINE LEVEL
        
//TODO AAAAAAAAAAAAAAAAA: Well, the thing is this: the connectivities of any level must be expressed 
// in terms of the node numbering OF THAT LEVEL!
//if the connectivities are expressed with the FINE NODE NUMBERING, 
//then we always have to convert to the LEVEL NODE NUMBERING!!!
//I guess the COORDINATES are WRONG...
//Ok, first of all we already have the connectivities at all levels
// of quadratic elements in mesh.h5, in FINE NODE NUMBERING.
//We only have to convert them to LEVEL NODE NUMBERING.

//Ok, I fixed the COORDINATES of the Qnodes at EVERY LEVEL.
//Plus, I have the "LINEARIZED" CONNECTIVITIES at EVERY LEVEL.
//Therefore, it seems like I can print any vector at EVERY LEVEL.
//Now, I still have a problem for the LINEAR variables.
//I guess I'm putting them in the WRONG PLACES before the interpolations.

//allora, quando faccio i loop con _off_nd, quei numeri li' corrispondono alle posizioni nel vettore sol?
//TODO chi da' l'ordine dei nodi del mesh A CIASCUN LIVELLO?
// Le COORDINATE DEL MESH A CIASCUN LIVELLO, quello e' l'ordine di riferimento!

//ok, the position i corresponds to the FINE MESH... then you must translate it for the LEVEL mesh

//For the linear variables we have TWO PROC LOOPS

         //Now I guess I have to pick the position of the linear nodes on the mesh by using the ExtendedLevel on the map...
         //In the following interpolation the loop is on the elements. From the elements you pick the connectivities,
//          which yield you the FINE QUADRATIC NODES, from which you pick the node numbers at level


//  PRINT OF SOL0000
// the first sol, and the case, are printed BEFORE the INITIAL CONDITIONS are set,
// or the initial conditions are set only at the FINE LEVEL?
// I would say the first one, otherwise I would expect different solutions...
// The first sol and the case are supposed to be equal...
//No wait, what i am saying is not true, because the first sol and the case 
// at the FINE LEVEL are printed as they should,
// after the initial conditions,
//but not on COARSER LEVELS!
// Ok now the initial conditions are set only AT THE FINE LEVEL.
// Depending on the kind of multigrid cycle, you may want to set the initial conditions
// ONLY AT THE FINE LEVEL or AT ALL LEVELS...
// For us, let us just do the GenIc function in such a way that all the levels can be treated separately.
// then, if we need it, we call it for ALL LEVELS, or we call it for ONLY THE FINE, or ONLY THE COARSE, or whatever...

//TODO ok, we have to remember one basic principle about our multigrid algorithm:
//the true solution is contained only at the FINE LEVEL!
//at all the other levels we have the DELTAS
//so, the level where to pick the values is the FINE LEVEL.
//Then, of course, we will print at all levels, but the VALUES from x_old are taken from the FINE LEVEL
//no... but wait a second... i am printing at all levels, so that's fine! I wanna print the RESIDUAL for all levels,
//except for the fine level where i print the true solution

// This prints All Variables of One Equation    
void write_system_solutions(const std::string namefile, const MeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn) {

  std::vector<FEElemBase*> fe_in(QL);
  for (int fe=0; fe<QL; fe++)    fe_in[fe] = FEElemBase::build(mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh->_mesh_order)._geomel_id.c_str(),fe);
  
  
    hid_t file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
   
    // ==========================================
    // =========== FOR ALL LEVELS ===============
    // ==========================================
    for (uint Level = 0; Level < mesh->_NoLevels; Level++)  {

      std::ostringstream grname; grname << "LEVEL" << Level;
//      hid_t group_id = H5Gcreate(file_id, grname.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//      hid_t group_id = H5Gopen(file_id, grname.str().c_str(),H5P_DEFAULT);
    
    int NGeomObjOnWhichToPrint[QL];
    NGeomObjOnWhichToPrint[QQ] = mesh->_NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[LL] = mesh->_NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[KK] = mesh->_n_elements_vb_lev[VV][Level]*mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh->_mesh_order).n_se;
    
    const uint n_nodes_lev = mesh->_NoNodesXLev[Level];
    double* sol_on_Qnodes  = new double[n_nodes_lev];  //TODO VALGRIND //this is QUADRATIC because it has to hold  either quadratic or linear variables and print them on a QUADRATIC mesh
    
    // ===================================
    // ========= QUADRATIC ===============
    // ===================================
    for (uint ivar=0; ivar<dofmap->_nvars[QQ]; ivar++)        {
      
      int pos_in_mesh_obj = 0;   
         for (uint isubdom=0; isubdom<mesh->_NoSubdom; isubdom++) {
            uint off_proc=isubdom*mesh->_NoLevels;
     
            for (int fine_node = mesh->_off_nd[QQ][off_proc];
                     fine_node < mesh->_off_nd[QQ][off_proc+Level+1]; fine_node++) {
      
  	int pos_in_sol_vec_lev = dofmap->GetDof(Level,QQ,ivar,fine_node);
	int pos_on_Qnodes_lev  = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ]; 

#ifndef NDEBUG
         if ( pos_on_Qnodes_lev >= (int) n_nodes_lev ) { std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl; abort(); }
#endif
        sol_on_Qnodes[ pos_on_Qnodes_lev/* pos_in_mesh_obj*/ ] = (* eqn->_x_old[Level])(pos_in_sol_vec_lev) * eqn->_refvalue[ ivar + dofmap->_VarOff[QQ] ];
	pos_in_mesh_obj++;
	  }
       }  //end subd
       
#ifndef NDEBUG
	 if (pos_in_mesh_obj != NGeomObjOnWhichToPrint[QQ]) { std::cout << "Wrong counting of quadratic nodes" << std::endl; abort(); }
#endif

     std::ostringstream var_name;  var_name << eqn->_var_names[ ivar + dofmap->_VarOff[QQ] ] << "_" << grname.str(); 	 //         std::string var_name = grname.str() + "/" + _var_names[ivar];
     hsize_t  dimsf[2];  dimsf[0] = NGeomObjOnWhichToPrint[QQ];  dimsf[1] = 1;
     IO::print_Dhdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);   //TODO VALGRIND

     }

    // =================================
    // ========= LINEAR ================
    // =================================
    uint elnds[QL_NODES];
    elnds[QQ] = mesh->GetGeomEl(mesh->get_dim()-1-VV,QQ)._elnds;
    elnds[LL] = mesh->GetGeomEl(mesh->get_dim()-1-VV,LL)._elnds;
    double* elsol_c = new double[elnds[LL]];
    
    for (uint ivar=0; ivar < dofmap->_nvars[LL]; ivar++)        {
      
//               for (uint i=0; i< n_nodes_lev; i++) { sol_on_Qnodes[i] = 0.; }
    
	for (uint isubdom=0; isubdom<mesh->_NoSubdom; isubdom++) {
	     uint off_proc=isubdom*mesh->_NoLevels;
            for (int fine_node = mesh->_off_nd[QQ][off_proc];
                     fine_node < mesh->_off_nd[QQ][off_proc]+
                         mesh->_off_nd[LL][off_proc + Level+1 ]
                       - mesh->_off_nd[LL][off_proc]; fine_node++) {
	      
	    int pos_in_sol_vec_lev = dofmap->GetDof(Level,LL,ivar,fine_node);
 	    int pos_on_Qnodes_lev = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ];

#ifndef NDEBUG
	 if ( pos_in_sol_vec_lev == -1 ) { std::cout << "Not correct DOF number at required level" << std::endl; abort(); }
         if ( pos_on_Qnodes_lev >= (int) n_nodes_lev ) { std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl; abort(); }
#endif

         sol_on_Qnodes[ pos_on_Qnodes_lev ] = (*eqn->_x_old[Level])(pos_in_sol_vec_lev) * eqn->_refvalue[ ivar + dofmap->_VarOff[LL] ];
	 
            }
        }

        //  2bB element interpolation over the fine mesh -----------------------
        // the way you filled linear positions before completely affects what happens next, which is only geometric
        for (uint iproc=0; iproc<mesh->_NoSubdom; iproc++) {
               uint off_proc = iproc*mesh->_NoLevels;
	       int iel_b = mesh->_off_el[VV][off_proc + Level];
	       int iel_e = mesh->_off_el[VV][off_proc + Level + 1];
	       
            for (int iel = 0; iel < (iel_e-iel_b); iel++) {
      
                for (uint in=0; in < elnds[LL]; in++) {
		  int pos_Qnode_fine = mesh->_el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];
		  int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
		  elsol_c[in] = sol_on_Qnodes[ pos_Qnode_lev ];   /**_refvalue[ivar]*/ //Do not multiply here!
		}

                for (uint in=0; in < elnds[QQ]; in++) { //TODO this loop can be done from elnds[LL] instead of from 0
                    double sum=0.;
                    for (uint jn=0; jn<elnds[LL]; jn++) {
                        sum += fe_in[LL]->get_prol(in*elnds[LL]+jn)*elsol_c[jn];
                    }
                    
                    int pos_Qnode_fine = mesh->_el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];       //Qnode in FINE NUMBERING
                    int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];  //Qnode in Level NUMBERING

#ifndef NDEBUG
                    if ( pos_Qnode_lev == -1 ) { std::cout << "Not correct node number at required level" << std::endl; abort(); }
 		    if ( pos_Qnode_lev >= (int) n_nodes_lev ) { std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl; abort(); }
#endif

 		    sol_on_Qnodes[ pos_Qnode_lev ] = sum;
                 }
              }
          } // 2bB end interpolation over the fine mesh --------
        
     std::ostringstream var_name; var_name << eqn->_var_names[ ivar + dofmap->_VarOff[LL] ] << "_" << grname.str();
     hsize_t  dimsf[2]; dimsf[0] = NGeomObjOnWhichToPrint[LL];  dimsf[1] = 1;
     IO::print_Dhdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);
     
    } // ivar linear

      delete []elsol_c;
      delete []sol_on_Qnodes;

     // ===================================
     // ========= CONSTANT ================
     // ===================================
  double *sol_on_cells;   sol_on_cells = new double[ NGeomObjOnWhichToPrint[KK] ];

  for (uint ivar=0; ivar < dofmap->_nvars[KK]; ivar++)        {
      
  int cel=0;
  for (uint iproc=0; iproc<mesh->_NoSubdom; iproc++) {
               uint off_proc = iproc*mesh->_NoLevels;
   
            int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < iproc; pr++) { sum_elems_prev_sd_at_lev += mesh->_off_el[VV][pr*mesh->_NoLevels + Level + 1] - mesh->_off_el[VV][pr*mesh->_NoLevels + Level]; }

	    for (int iel = 0;
              iel <    mesh->_off_el[VV][off_proc + Level+1]
                      - mesh->_off_el[VV][off_proc + Level]; iel++) {
             int elem_lev = iel + sum_elems_prev_sd_at_lev;
	  int dof_pos_lev = dofmap->GetDof(Level,KK,ivar,elem_lev);   
      for (uint is=0; is< mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh->_mesh_order).n_se; is++) {      
	   sol_on_cells[cel*mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh->_mesh_order).n_se + is] = (* eqn->_x_old[Level])(dof_pos_lev) * eqn->_refvalue[ ivar + dofmap->_VarOff[KK] ];
      }
      cel++;
    }
  }
  
  std::ostringstream varname; varname << eqn->_var_names[ ivar + dofmap->_VarOff[KK] ] << "_" << grname.str();         //   std::string varname = grname.str() + "/" + _var_names[_nvars[QQ]+_nvars[LL]+ivar];
  hsize_t dimsf[2]; dimsf[0] = NGeomObjOnWhichToPrint[KK]; dimsf[1] = 1;
  IO::print_Dhdf5(file_id,varname.str(),dimsf,sol_on_cells);   
      
    } //end KK

     delete [] sol_on_cells;
  
//         H5Gclose(group_id);
	
    } //end Level
    
    H5Fclose(file_id);   //TODO VALGRIND

      for (int fe=0; fe<QL; fe++)  {  delete fe_in[fe]; }

    
    return;
}


// ===================================================
/// This function reads the system solution from namefile.h5
//TODO this must be modified in order to take into account KK element dofs
void read_system_solutions(std::string namefile, MeshTwo* mesh, DofMap* dofmap, SystemTwo* eqn) {
//this is done in parallel

  std::cout << "read_system_solutions still has to be written for CONSTANT elements, BEWARE!!! ==============================  " << std::endl;
  
    const uint Level = mesh->_NoLevels-1;
    
    const uint mesh_ord = (int) mesh->GetRuntimeMap().get("mesh_ord");
    const uint offset   =       mesh->_NoNodesXLev[mesh->_NoLevels-1];

    // file to read
    double *sol=new double[offset]; // temporary vector
    hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

    // reading loop over system varables
    for (uint ivar=0;ivar< dofmap->_nvars[LL]+dofmap->_nvars[QQ]; ivar++) {
        uint el_nds = mesh->GetGeomEl(mesh->get_dim()-1-VV,QQ)._elnds;
        if (ivar >= dofmap->_nvars[QQ]) el_nds = mesh->GetGeomEl(mesh->get_dim()-1-VV,LL)._elnds; // quad and linear
        // reading ivar param
       std::ostringstream grname; grname << eqn->_var_names[ivar] << "_" << "LEVEL" << Level;
        IO::read_Dhdf5(file_id,grname.str(),sol);
        double Irefval = 1./eqn->_refvalue[ivar]; // units

        // storing  ivar variables (in parallell)
        for (int iel=0;iel <  mesh->_off_el[0][mesh->_iproc*mesh->_NoLevels+mesh->_NoLevels]
                -mesh->_off_el[0][mesh->_iproc*mesh->_NoLevels+mesh->_NoLevels-1]; iel++) {
            uint elem_gidx=(iel+mesh->_off_el[0][mesh->_iproc*mesh->_NoLevels+mesh->_NoLevels-1])*mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh_ord)._elnds;
            for (uint i=0; i<el_nds; i++) { // linear and quad
                int k=mesh->_el_map[0][elem_gidx+i];   // the global node
                eqn->_x[mesh->_NoLevels-1]->set(dofmap->GetDof(mesh->_NoLevels-1,QQ,ivar,k), sol[k]*Irefval); // set the field
            }
        }
    }

    eqn->_x[mesh->_NoLevels-1]->localize(* eqn->_x_old[mesh->_NoLevels-1]);
    // clean
    H5Fclose(file_id);
    delete []sol;
    
    return;
}



// =====================================================================
/// This function  defines the boundary conditions for  DA systems:
// how does this function behave when e.g. both velocity and pressure are linear?
//this function runs over ALL the domain, so it's not only on a portion.
//ok we must think of this function in terms of the quantities
//we must also keep in mind that the DofMap is built by following a GIVEN ORDER;
// so, if you have 2 quadratic and 1 linear, you choose to put first the quadratic then the linear
// here, we loop over BOUNDARY ELEMENTS
//for every boundary element, loop over its nodes
// pick the coordinate of that node
//for every node, pick all its degrees of freedom from the dofmap
//at THAT COORDINATE, the DoF loop corresponds to a Quantity loop



// ================================================================
/// This function prints the boundary conditions: for quad and linear fem
//ok, this routine prints FLAGS, not VALUES.
//therefore, I'm not interested in INTERPOLATING FLAGS.
//Flags are either one or zero, no intermediate values are requested
//if all the nodes are 0, then it means do the pressure integral ---> put a zero
//if at least one node is 1 ---> put a 1. (anyway, the bc_p are not used in "volume" manner)
//the problem is always that we are printing LINEAR VARIABLES on a QUADRATIC mesh,
// so actually the interpolated values are there just for printing...
//

//TODO I must do in such a way as to print bc just like x_old,
// it is supposed to be exactly the same routine

//TODO I should print these bc flags, as well as the solution vectors, for ALL LEVELS!
// in this way I could see what happens at every level
// I already have the connectivities at all levels, thanks to gencase
// The only problem is that in 3D I do not see the connectivities, so I should convert all the levels to LINEARIZED REFINED

//also, apart from very little things, the routine is very similar to printing any field defined on each level
// And it's an integer and when you do the average you need to use the ceil() function...

//ok, when you print a field on a coarser grid, you need to specify a coarser topology,
// but also a coarser number of points...
// but, when I print only BOUNDARY MESH or VOLUME MESH at all levels,
// I put the FINE numbers and it's ok...
// so how does it pick the coordinates?


// A = NUMBER SEEN in PARAVIEW
// B = NUMBER in HDF5

//Ok, when i load coarser meshes, the number of cells changes but not the number of points...
// but the thing is: it seems that the NODES are ORDERED in SUCH A WAY THAT 
// you have first COARSER then 

//AAA: the equation A=B only holds with the WHOLE MESH!!!
// When you add some FILTER in paraview, for instance if you do a CLIP,
// Then the numbers DO NOT CORRESPOND ANYMORE!!!

//Ok, if i have to print a data attribute on a grid, i need to provide a SMALLER NUMBER of COORDINATES,
// or maybe make a LARGER DATA VECTOR...
//Due to the way we order the nodes, I guess we can say in the XDMF that the array CAN BE CUT! Let's see...
//You cannot put in the XDMF a SMALLER NUMBER than the dimension of the HDF5 ARRAY,
//otherwise it gives SEGMENTATION FAULT!
//So I need to have HDF5 fields for the coordinates of EACH LEVEL!
//I dont wanna put useless values, so I'll just print the coordinates for each level...

//Now there is a problem with the TIME printing... I am including all the solution files but it is picking up
// the COARSE LEVEL instead of the FINEST ONE

//TODO the problem is that now every solution file contains DIFFERENT GRIDS,
// and when I load the TIME file then it loads THE FIRST ONE IT MEETS in the SOL FILES!!!

// AAA, ecco l'errore! la connettivita' che percorriamo per interpolare i valori lineari 
// e' quella FINE, ma noi ora dobbiamo prendere quella DI CIASCUN LIVELLO SEPARATAMENTE!


void write_system_solutions_bc(const std::string namefile, const MeshTwo* mesh, const DofMap* dofmap, const SystemTwo* eqn, const int* bc, int** bc_fe_kk ) {
  
  std::vector<FEElemBase*> fe_in(QL);
  for (int fe=0; fe<QL; fe++)    fe_in[fe] = FEElemBase::build(mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh->_mesh_order)._geomel_id.c_str(),fe);

    hid_t file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

    std::string  bdry_suffix = DEFAULT_BDRY_SUFFIX;

    const uint Lev_pick_bc_NODE_dof = mesh->_NoLevels-1;  //we use the FINE Level as reference
    
    // ==========================================
    // =========== FOR ALL LEVELS ===============
    // ==========================================
    for (uint Level = 0; Level < mesh->_NoLevels; Level++)  {
    
      std::ostringstream grname; grname << "LEVEL" << Level;
  
   int NGeomObjOnWhichToPrint[QL];
    NGeomObjOnWhichToPrint[QQ] = mesh->_NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[LL] = mesh->_NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[KK] = mesh->_n_elements_vb_lev[VV][Level]*mesh->GetGeomEl(mesh->get_dim()-1-VV,QQ).n_se;
  
    const uint n_nodes_lev = mesh->_NoNodesXLev[Level];
    int* sol_on_Qnodes = new int[n_nodes_lev];  //this vector will contain the values of ONE variable on ALL the QUADRATIC nodes
    //for the quadratic variables it'll be just a copy, for the linear also interpolation
    
    // ===================================
    // ========= QUADRATIC ================
    // ===================================
    for (uint ivar=0; ivar < eqn->_dofmap._nvars[QQ]; ivar++)        {
      
      for (uint isubdom=0; isubdom<mesh->_NoSubdom; isubdom++) {
            uint off_proc=isubdom*mesh->_NoLevels;
     
          for (int fine_node = mesh->_off_nd[QQ][off_proc];
	           fine_node < mesh->_off_nd[QQ][off_proc+Level+1]; fine_node ++) {

	int pos_in_sol_vec_lev = eqn->_dofmap.GetDof(Lev_pick_bc_NODE_dof,QQ,ivar,fine_node);
	int pos_on_Qnodes_lev  = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ]; 

	sol_on_Qnodes[ pos_on_Qnodes_lev ] = bc[pos_in_sol_vec_lev];
          }
       }  //end subd

       std::ostringstream var_name; var_name << eqn->_var_names[ ivar + eqn->_dofmap._VarOff[QQ] ] << "_" << grname.str() << bdry_suffix;
       hsize_t dimsf[2];  dimsf[0] = NGeomObjOnWhichToPrint[QQ];  dimsf[1] = 1;
       IO::print_Ihdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);
    }

    // ===================================
    // ========= LINEAR ==================
    // ===================================
    uint elnds[QL_NODES];
    elnds[QQ] =mesh->GetGeomEl(mesh->get_dim()-1-VV,QQ)._elnds;
    elnds[LL] =mesh->GetGeomEl(mesh->get_dim()-1-VV,LL)._elnds;
    double *elsol_c = new double[elnds[LL]];

    for (uint ivar=0; ivar < eqn->_dofmap._nvars[LL]; ivar++)   {

	for (uint isubdom=0; isubdom<mesh->_NoSubdom; isubdom++) {
	              uint off_proc=isubdom*mesh->_NoLevels;

            for (int fine_node =   mesh->_off_nd[QQ][off_proc];
                     fine_node <   mesh->_off_nd[QQ][off_proc]
                    + mesh->_off_nd[LL][off_proc+Level+1]
                    - mesh->_off_nd[LL][off_proc]; fine_node++) {
	      
            int pos_in_sol_vec_lev = eqn->_dofmap.GetDof(Lev_pick_bc_NODE_dof,LL,ivar,fine_node); 
 	    int pos_on_Qnodes_lev  = mesh->_Qnode_fine_Qnode_lev[Level][ fine_node ];

	    sol_on_Qnodes[ pos_on_Qnodes_lev ]=  bc[pos_in_sol_vec_lev]; 
          }
        }
 
        //  2bB element interpolation over the fine mesh 
        for (uint iproc=0; iproc<mesh->_NoSubdom; iproc++) {
               uint off_proc = iproc*mesh->_NoLevels;
	       int iel_b = mesh->_off_el[VV][off_proc + Level];
	       int iel_e = mesh->_off_el[VV][off_proc + Level + 1];
            for (int iel = 0; iel < (iel_e-iel_b); iel++) {

                for (uint in=0; in < elnds[LL]; in++)  {
 		  int pos_Qnode_fine = mesh->_el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];
		  int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
        	  elsol_c[in]= sol_on_Qnodes[ pos_Qnode_lev ];
		}
		
		  for (uint in=0; in < elnds[QQ]; in++) { // mid-points
                    double sum=0;
                    for (uint jn=0; jn<elnds[LL]; jn++) {
                        sum += fe_in[LL]->get_prol(in*elnds[LL]+jn)*elsol_c[jn];
                    }
                    
                    int pos_Qnode_fine = mesh->_el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];
                    int pos_Qnode_lev  = mesh->_Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
              
                    sol_on_Qnodes[ pos_Qnode_lev ] = ceil(sum);   //the ceiling, because you're putting double over int!
                }
            }
        } // 2bB end interpolation over the fine mesh

       std::ostringstream var_name; var_name << eqn->_var_names[ ivar + eqn->_dofmap._VarOff[LL] ] << "_" << grname.str() << bdry_suffix;
       hsize_t  dimsf[2];  dimsf[0] = NGeomObjOnWhichToPrint[LL];  dimsf[1] = 1;
       IO::print_Ihdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);
    } // ivar
    
    delete [] elsol_c;
    delete [] sol_on_Qnodes;

     // ===================================
     // ========= CONSTANT ================
     // ===================================
     for (uint ivar=0; ivar < eqn->_dofmap._nvars[KK]; ivar++)        {
      
  int *sol_on_cells;   sol_on_cells = new int[ NGeomObjOnWhichToPrint[KK] ];
  
  int cel=0;
  for (uint iproc=0; iproc<mesh->_NoSubdom; iproc++) {
               uint off_proc = iproc*mesh->_NoLevels;
	       
            int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < iproc; pr++) { sum_elems_prev_sd_at_lev += mesh->_off_el[VV][pr*mesh->_NoLevels + Level + 1] - mesh->_off_el[VV][pr*mesh->_NoLevels + Level]; }

	    for (int iel = 0;
              iel <    mesh->_off_el[VV][off_proc + Level+1]
                     - mesh->_off_el[VV][off_proc + Level]; iel++) {
      for (uint is=0; is< mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh->_mesh_order).n_se; is++) {      
	sol_on_cells[cel* mesh->GetGeomEl(mesh->get_dim()-1-VV,mesh->_mesh_order).n_se + is] = bc_fe_kk[Level][iel + sum_elems_prev_sd_at_lev + ivar*mesh->_n_elements_vb_lev[VV][Level]]; //this depends on level!
      }
      cel++;
    }
  }
  
  std::ostringstream var_name; var_name << eqn->_var_names[ ivar + eqn->_dofmap._VarOff[KK] ] << "_" << grname.str() << bdry_suffix;
  hsize_t dimsf[2]; dimsf[0] = NGeomObjOnWhichToPrint[KK]; dimsf[1] = 1;
  IO::print_Ihdf5(file_id,var_name.str(),dimsf,sol_on_cells);   
      
      delete [] sol_on_cells;
      
      } //end KK   
    
    } //end Level
    

     H5Fclose(file_id);

      for (int fe=0; fe<QL; fe++)  {  delete fe_in[fe]; }
      
    return;
}






} //end namespace IO


} //end namespace femus


