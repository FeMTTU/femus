// tests the inclusion algorithm for all 2D and 3D elements


#include "FemusInit.hpp"
#include "Marker.hpp"
#include "MultiLevelMesh.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "MyVector.hpp"

using namespace femus;

void PrintStreamline(const std::string output_path, const std::vector< std::vector<double> > &xn);

// double InitalValueU(const std::vector < double >& x) {
//   return 0.5;
// }
// 
// double InitalValueV(const std::vector < double >& x) {
//   return 0.5;
// }
// 
// double InitalValueW(const std::vector < double >& x) {
//   return 0.5;
// }


double InitalValueU(const std::vector < double >& x) {
  return -x[1];
}

double InitalValueV(const std::vector < double >& x) {
  return x[0];
}

double InitalValueW(const std::vector < double >& x) {
  return 0.;
}


// double InitalValueU(const std::vector < double >& x) {
//   return (-x[1]+x[2])/sqrt(3);
// }
// 
// double InitalValueV(const std::vector < double >& x) {
//   return (x[0]-x[2])/sqrt(3);
// }
// 
// double InitalValueW(const std::vector < double >& x) {
//   return (x[1]-x[0])/sqrt(3);
// }



bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if(elemgroupnumber == 6 && level < 4) refine = 1;
  if(elemgroupnumber == 7 && level < 5) refine = 1;
  if(elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  std::vector < double > x(3, 0); // marker
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  std::vector < std::string > variablesToBePrinted;

  /* element types
  0 = HEX
  1 = TET
  2 = WEDGE
  3 = QUAD
  4 = TRI
   */

  unsigned solType = 2;

  std::cout << " --------------------------------------------------     TEST     --------------------------------------------------" << std::endl;

  //mlMsh.ReadCoarseMesh("./input/prism3D.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/tri2.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cubeHex.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/cubeTet.neu", "seventh", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  unsigned dim = mlMsh.GetDimension();



  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  mlSol.Initialize("U" , InitalValueU);
  mlSol.Initialize("V" , InitalValueV);
  if(dim == 3) mlSol.Initialize("W", InitalValueW);

//   //Test 1 (QUAD):
//
  //NOTE tests ran with 4 procs
//   x[0] = -0.46875; //the marker is in element 191 (proc 3 of 4)
//   x[1] = -0.5;
//   x[2] = 0.;
  
    x[0] = 0.125; 
    x[1] = 0.125;
    x[2] = -0.25;



// //Test 1 (TET):  element 20
//     //NOTE Tests ran with 2 procs
//       x[0] = -0.5;
//       x[1] = 0.;
//       x[2] = 0.;

  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
  Marker a1Quad(x, VOLUME, mlMsh.GetLevel(0), solType, true);
  //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
  std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
  std::cout << " The marker type is " <<  a1Quad.GetMarkerType() << std::endl;

  double T = 2 * acos(-1.);
  unsigned n  = 100;


  std::vector < std::vector < double > > xn(n + 1);
  for(unsigned k = 0; k < n; k++) {
    a1Quad.GetMarkerCoordinates(xn[k]);
    a1Quad.Advection(mlSol.GetLevel(0), 2, T / n);
  }
  a1Quad.GetMarkerCoordinates(xn[n]);
  for(unsigned i = 0;  i < xn.size(); i++) {
    for(unsigned d = 0; d < xn[i].size(); d++) {
      std::cout << xn[i][d] << " ";
    }
    std::cout << std::endl;
  }








//   std::vector < MyVector < double > > xn(dim);
//   for(unsigned k = 0; k < n; k++) {
//     a1Quad.GetMarkerCoordinates(xn);
//     a1Quad.Advection(mlSol.GetLevel(0), 2, T / n);
//   }
//   a1Quad.GetMarkerCoordinates(xn);
//   for(unsigned d = 0; d < dim; d++) {
//     xn[d].stack();
//     std::cout << xn[d] << std::endl;
//   }






  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);


  PrintStreamline(DEFAULT_OUTPUTDIR, xn);


  return 0;
}

#include "Files.hpp"
#include <b64/b64.h>

void PrintStreamline(const std::string output_path, const std::vector< std::vector<double> > &xn) {

  // *********** open vtu files *************
  std::ofstream fout;

  std::string dirnamePVTK = "./";
  Files files;
  files.CheckDir(output_path, dirnamePVTK);

  std::string filename_prefix = "streamline";

  std::ostringstream filename;
  filename << output_path << "./" << filename_prefix << ".vtu";

  fout.open(filename.str().c_str());
  if(!fout.is_open()) {
    std::cout << std::endl << " The output file " << filename.str() << " cannot be opened.\n";
    abort();
  }

  unsigned nvt = xn.size();
  unsigned nel = nvt - 1;
  
  const unsigned dim_array_coord [] = { nvt * 3 * sizeof(float) };
  const unsigned dim_array_conn[]   = { nel * 2 * sizeof( int ) };
  const unsigned dim_array_off []   = { nel * sizeof( int ) };
  const unsigned dim_array_type []  = { nel * sizeof( short unsigned ) };

  unsigned buffer_size = ( dim_array_coord[0] > dim_array_conn[0] ) ? dim_array_coord[0] : dim_array_conn[0];
  void* buffer_void = new char [buffer_size];
  char* buffer_char = static_cast <char*>(buffer_void);

  size_t cch;
  cch = b64::b64_encode(&buffer_char[0], buffer_size , NULL, 0);
  vector <char> enc;
  enc.resize(cch);
  char* pt_char;

  // *********** write vtu header ************
  fout << "<?xml version=\"1.0\"?>" << std::endl;
  fout << "<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  fout << "  <UnstructuredGrid>" << std::endl;


  fout  << "    <Piece NumberOfPoints= \"" << nvt 
 	<< "\" NumberOfCells= \"" << nel 
	<< "\" >" << std::endl;

  //-----------------------------------------------------------------------------------------------
  // print coordinates *********************************************Solu*******************************************
  fout  << "      <Points>" << std::endl;
  fout  << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;

  // point pointer to common mamory area buffer of void type;
  float* var_coord = static_cast< float* >(buffer_void);

  for(unsigned i = 0; i < nvt; i++) {
    for(unsigned d = 0; d < 3; d++) {
      var_coord[i * 3 + d] = (xn[i].size() > d) ? xn[i][d] : 0.;
    }
  }

  cch = b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), NULL, 0);
  b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), &enc[0], cch);
  pt_char = &enc[0];
  for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;

  //print coordinates array
  cch = b64::b64_encode(&var_coord[0], dim_array_coord[0] , NULL, 0);
  b64::b64_encode(&var_coord[0], dim_array_coord[0], &enc[0], cch);
  pt_char = &enc[0];
  for(unsigned i = 0; i < cch; i++, pt_char++) fout << *pt_char;
  fout << std::endl;

  fout  << "        </DataArray>" << std::endl;
  fout  << "      </Points>" << std::endl;
  //-----------------------------------------------------------------------------------------------

    // Printing of element connectivity - offset - format type  *
    fout  << "      <Cells>" << std::endl;
    //-----------------------------------------------------------------------------------------------
    //print connectivity
    fout  << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">" << std::endl;
  
    // point pointer to common mamory area buffer of void type;
    int* var_conn = static_cast <int*>( buffer_void );
    unsigned icount = 0;
    for(unsigned iel = 0; iel < nel; iel++ ) {
      for( unsigned j = 0; j < 2; j++ ) {
        var_conn[icount] = iel + j;;
        icount++;
      }
    }

    //print connectivity dimension
    cch = b64::b64_encode( &dim_array_conn[0], sizeof( dim_array_conn ), NULL, 0 );
    b64::b64_encode( &dim_array_conn[0], sizeof( dim_array_conn ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print connectivity array
    cch = b64::b64_encode( &var_conn[0], dim_array_conn[0] , NULL, 0 );
    b64::b64_encode( &var_conn[0], dim_array_conn[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    //------------------------------------------------------------------------------------------------
    fout  << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">" << std::endl;
     // point pointer to common memory area buffer of void type;
    int* var_off = static_cast <int*>( buffer_void );
    // print offset array
    for( int iel = 0; iel < nel; iel++ ) {
      var_off[iel] = (iel+1) * 2;
    }

    //print offset dimension
    cch = b64::b64_encode( &dim_array_off[0], sizeof( dim_array_off ), NULL, 0 );
    b64::b64_encode( &dim_array_off[0], sizeof( dim_array_off ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print offset array
    cch = b64::b64_encode( &var_off[0], dim_array_off[0] , NULL, 0 );
    b64::b64_encode( &var_off[0], dim_array_off[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    fout  << std::endl;

    fout  << "        </DataArray>" << std::endl;

  
 
    //--------------------------------------------------------------------------------------------------

    //Element format type : 23:Serendipity(8-nodes)  28:Quad9-Biquadratic
    fout  << "        <DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << std::endl;
    
    // point pointer to common mamory area buffer of void type;
    unsigned short* var_type = static_cast <unsigned short*>( buffer_void );
    
    for( unsigned iel = 0; iel < nel; iel++ ) {
      var_type[iel] = 3;
    }

    //print element format dimension
    cch = b64::b64_encode( &dim_array_type[0], sizeof( dim_array_type ), NULL, 0 );
    b64::b64_encode( &dim_array_type[0], sizeof( dim_array_type ), &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    //print element format array
    cch = b64::b64_encode( &var_type[0], dim_array_type[0] , NULL, 0 );
    b64::b64_encode( &var_type[0], dim_array_type[0], &enc[0], cch );
    pt_char = &enc[0];
    for( unsigned i = 0; i < cch; i++, pt_char++ ) fout << *pt_char;

    fout  << std::endl;
    fout  << "        </DataArray>" << std::endl;
   

    
    //----------------------------------------------------------------------------------------------------
//
    fout  << "      </Cells>" << std::endl;
  
  
  
  //-----------------------------------------------------------------------------------------------
  // Printing of element connectivity - offset - format type  *
//   fout  << "      <Cells>" << std::endl;
//   fout  << "      </Cells>" << std::endl;

  fout << "    </Piece>" << std::endl;
  fout << "  </UnstructuredGrid>" << std::endl;
  fout << "</VTKFile>" << std::endl;
  fout.close();

  delete [] var_coord;
  //--------------------------------------------------------------------------------------------------------
  return;
}

