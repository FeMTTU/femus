# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "obj_io.hpp"

//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

void obj_face_node_print ( int face_num, int order_max, int face_order[],
  int face_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    OBJ_FACE_NODE_PRINT prints the node indices for each face.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int ORDER_MAX, the maximum number of vertices
//    per face.
//
//    Input, int FACE_ORDER[FACE_NUM], the number of vertices
//    per face.
//
//    Input, int FACE_NODE[ORDER_MAX*FACE_NUM], the nodes that
//    make up each face.
//
{
  int face;
  int i;
  int order;

  cout << "\n";
  cout << "    Face   Order      Nodes\n";
  cout << "\n";

  for ( face = 0; face < face_num; face++ )
  {
    order = face_order[face];
    cout << "  " << setw(6) << face
         << "  " << setw(6) << order;
    for ( i = 0; i < order; i++ )
    {
      cout << "  " << setw(6) << face_node[i+face*order_max];
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void obj_normal_vector_print ( int normal_num, double normal_vector[] )

//****************************************************************************80
//
//  Purpose:
//
//    OBJ_NORMAL_VECTOR_PRINT prints the normal vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NORMAL_NUM, the number of normal vectors.
//
//    Input, double NORMAL_VECTOR[3*NORMAL_NUM], the normal vectors.
//
{
  int face;
  int i;
  int normal;

  cout << "\n";
  cout << "  Normal Vectors:\n";
  cout << "\n";

  for ( normal = 0; normal < normal_num; normal++ )
  {
    cout << "  " << setw(6) << normal;
    for ( i = 0; i < 3; i++ )
    {
      cout << "  " << setw(14) << normal_vector[i+normal*3];
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void obj_node_xyz_print ( int node_num, double node_xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    OBJ_NODE_XYZ_PRINT prints the node coordinates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates
//    of the nodes.
//
{
  int i;
  int node;

  cout << "\n";
  cout << "    Node         Coordinates\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(6) << node;
    for ( i = 0; i < 3; i++ )
    {
       cout << "  " << setw(14) << node_xyz[i+node*3];
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void obj_size ( string filename, int *node_num, int *face_num, int *normal_num,
  int *order_max )

//****************************************************************************80
//
//  Purpose:
//
//    OBJ_SIZE determines sizes of graphics objects in an Alias OBJ file.
//
//  Example:
//
//    #  magnolia.obj
//
//    v -3.269770 -39.572201 0.876128
//    v -3.263720 -39.507999 2.160890
//    ...
//    v 0.000000 -9.988540 0.000000
//
//    vn 1.0 0.0 0.0
//    ...
//    vn 0.0 1.0 0.0
//
//    f 8 9 11 10
//    f 12 13 15 14
//    ...
//    f 788 806 774
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the input file name.
//
//    Output, int *NODE_NUM, the number of points.
//
//    Output, int *FACE_NUM, the number of faces.
//
//    Output, int *NORMAL_NUM, the number of normal vectors.
//
//    Output, int *ORDER_MAX, the maximum face order.
//
{
  ifstream input;
  int n;
  string text;
  int text_num;
//
//  Initialize.
//
  *node_num = 0;
  *face_num = 0;
  *normal_num = 0;
  *order_max = 0;
  text_num = 0;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "OBJ_SIZE - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }
  text_num = 0;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    text_num = text_num + 1;

    if ( s_len_trim ( text ) == 0 )
    {
      continue;
    }
    if ( text[0] == '#' )
    {
      continue;
    }
    if ( text[0] == '$' )
    {
      continue;
    }
//
//  F V1 V2 ... VN
//  Face.
//
    if ( text[0] == 'f' || text[0] == 'F' )
    {
      n = s_word_count ( text );
      *order_max = i4_max ( *order_max, n - 1 );
      *face_num = *face_num + 1;
    }
//
//  VN X Y Z
//  Vertex normals.
//
    else if ( ( text[0] == 'v' || text[0] == 'V' ) &&
              ( text[1] == 'n' || text[1] == 'N' ) )
    {
      *normal_num = *normal_num + 1;
    }
//
//  V X Y Z W
//  Geometric vertex.
//
    else if ( text[0] == 'v' || text[0] == 'V' )
    {
      *node_num = *node_num + 1;
    }

  }
//
//  Close the file.
//
  input.close ( );

  cout << "\n";
  cout << "  Read " << text_num << " lines from \"" << filename << "\".\n";
  return;
}
//****************************************************************************80

void obj_size_print ( string filename, int node_num, int face_num,
  int normal_num, int order_max )

//****************************************************************************80
//
//  Purpose:
//
//    OBJ_SIZE_PRINT prints sizes associated with an OBJ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Input, int NODE_NUM, the number of vertices defined.
//
//    Input, int FACE_NUM, the number of faces defined.
//
//    Input, int NORMAL_NUM, the number of normal
//    vectors defined.
//
//    Input, int ORDER_MAX, the maximum number of vertices
//    per face.
//
{
  cout << "\n";
  cout << "  Object sizes for OBJ file \"" << filename << "\".\n";
  cout << "\n";
  cout << "  Nodes              = " << node_num << "\n";
  cout << "  Faces              = " << face_num << "\n";
  cout << "  Maximum face order = " << order_max << "\n";
  cout << "  Normal vectors     = " << normal_num << "\n";

  return;
}
//****************************************************************************80

void obj_vertex_normal_print ( int order_max, int face_num, int face_order[],
  int vertex_normal[] )

//****************************************************************************80
//
//  Purpose:
//
//    OBJ_VERTEX_NORMAL_PRINT prints the normal vectors indices per vertex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER_MAX, the maximum number of vertices
//    per face.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_ORDER[FACE_NUM], the number of vertices
//    per face.
//
//    Input, int VERTEX_NORMAL[ORDER_MAX*FACE_NUM], the
//    indices of normal vectors per vertex.
//
{
  int face;
  int i;
  int order;

  cout << "\n";
  cout << "  Normal Vector Indices:\n";
  cout << "\n";
  cout << "    Face   Order\n";
  cout << "\n";

  for ( face = 0; face < face_num; face++ )
  {
    order = face_order[face];
    cout << "  " << setw(6) << face
         << "  " << setw(6) << order
         << "  ";
    for ( i = 0; i < order; i++ )
    {
      cout << "  " << setw(6) << vertex_normal[i+face*order_max];
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void obj_write ( string output_filename, int node_num, int face_num,
  int normal_num, int order_max, double node_xyz[], int face_order[],
  int face_node[], double normal_vector[], int vertex_normal[] )

//****************************************************************************80
//
//  Purpose:
//
//    OBJ_WRITE writes graphics information to an Alias OBJ file.
//
//  Discussion:
//
//    If no normal vectors are supplied (NORMAL_NUM <= 0) then
//    a simple format is used for the "F" records.  Otherwise,
//    the "v//vn" format is used.
//
//  Example:
//
//    #  no_normals.obj
//
//    g Group002
//
//    v -3.269770 -39.572201 0.876128
//    v -3.263720 -39.507999 2.160890
//    ...
//    v 0.000000 -9.988540 0.000000
//
//    f 8 9 11 10
//    f 12 13 15 14
//    ...
//    f 788 806 774
//
//    #  normals_supplied.obj
//
//    g Group001
//
//    v -3.269770 -39.572201 0.876128
//    v -3.263720 -39.507999 2.160890
//    ...
//    v 0.000000 -9.988540 0.000000
//
//    vn 0.0 1.0 0.0
//    vn 1.0 0.0 0.0
//    ...
//    vn 0.0 0.0 1.0
//
//    f 8//1 9//2 11//3 10//4
//    f 12//5 13//6 15//7 14//8
//    ...
//    f 788//800 806//803 774//807
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the output file.
//
//    Input, int NODE_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int NORMAL_NUM, the number of normal vectors.
//
//    Input, int ORDER_MAX, the maximum number of vertices
//    per face.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of points.
//
//    Input, int FACE_ORDER[FACE_NUM], the number of vertices
//    per face.
//
//    Input, int FACE_NODE[ORDER_MAX*FACE_NUM], the nodes
//    making faces.
//
//    Input, double NORMAL_VECTOR[3*NORMAL_NUM], normal vectors.
//
//    Input, int VERTEX_NORMAL[ORDER_MAX*FACE_NUM], the
//    indices of normal vectors per vertex.
//
{
  int face;
  int i;
  int j;
  int node;
  int normal;
  ofstream output;
  int text_num;
  int vertex;
  double w;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "OBJ_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file \"" << output_filename << "\".\n";
    exit ( 1 );
  }

  text_num = 0;

  output << "# " << output_filename << "\n";
  output << "# created by obj_io::obj_write.C\n";
  output << "\n";
  output << "g Group001\n";

  text_num = text_num + 4;
//
//  V: vertex coordinates.
//  For some reason, a fourth "coordinate" may be recommended.
//  What is its meaning?
//
  if ( 0 < node_num )
  {
    output << "\n";
    text_num = text_num + 1;
  }

  w = 1.0;
  for ( node = 0; node < node_num; node++ )
  {
    output << "v";
    for ( i = 0; i < 3; i++ )
    {
      output << "  " << node_xyz[i+3*node];
    }
    output << "  " << w << "\n";
    text_num = text_num + 1;
  }
//
//  VN: normal vectors.
//
  if ( 0 < normal_num )
  {
    output << "\n";
    text_num = text_num + 1;

    for ( normal = 0; normal < normal_num; normal++ )
    {
      output << "vn";
      for ( i = 0; i < 3; i++ )
      {
        output << "  " << normal_vector[i+normal*3];
      }
      output << "\n";
      text_num = text_num + 1;
    }
  }
//
//  F: Faces, specified as a list of triples, one triple for each vertex:
//     vertex index/vertex texture index/vertex normal index
//
  if ( 0 < face_num )
  {
    output << "\n";
    text_num = text_num + 1;
  }

  for ( face = 0; face < face_num; face++ )
  {
    output << "f";
    for ( vertex = 0; vertex < face_order[face]; vertex++ )
    {
      output << "  " << face_node[vertex+face*order_max];
      if ( 0 < normal_num )
      {
        output << "//" << vertex_normal[vertex+face*order_max];
      }
    }
    output << "\n";
    text_num = text_num + 1;
  }

  output.close ( );
//
//  Report.
//
  if ( false )
  {
    cout << "\n";
    cout << "OBJ_WRITE:\n";
    cout << "  Wrote " << text_num << " text lines to \""
       << output_filename << "\"\n";
  }

  return;
}
//****************************************************************************80

double *r8vec_cross_product_3d ( double v1[3], double v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], the coordinates of the vectors.
//
//    Output, double R8VEC_CROSS_PRODUCT_3D[3], the cross product vector.
//
{
  double *v3;

  v3 = new double[3];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int s_word_count ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int char_count;
  int i;
  int word_count;

  word_count = 0;
  blank = true;

  char_count = s.length ( );

  for ( i = 0; i < char_count; i++ )
  {
    if ( isspace ( s[i] ) )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_count = word_count + 1;
      blank = false;
    }
  }

  return word_count;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
