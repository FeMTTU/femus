int i4_max ( int i1, int i2 );
void obj_face_node_print ( int face_num, int order_max, int face_order[],
  int face_node[] );
void obj_normal_vector_print ( int normal_num, double normal_vector[] );
void obj_node_xyz_print ( int node_num, double node_xyz[] );
void obj_size ( string filename, int *node_num, int *face_num,
  int *normal_num, int *order_max );
void obj_size_print ( string filename, int node_num, int face_num,
  int normal_num, int order_max );
void obj_vertex_normal_print ( int order_max, int face_num, int face_order[],
  int vertex_normal[] );
void obj_write ( string output_filename, int node_num, int face_num,
  int normal_num, int order_max, double node_xyz[], int face_order[],
  int face_node[], double normal_vector[], int vertex_normal[] );
double *r8vec_cross_product_3d ( double v1[3], double v2[3] );
int s_len_trim ( string s );
int s_word_count ( string s );
void timestamp ( );
