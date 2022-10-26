#ifndef __femus_applications_00_tutorial_hpp__
#define __femus_applications_00_tutorial_hpp__


using namespace femus;


 
const std::vector< Unknown >  systems__provide_list_of_unknowns_lagrangian() {


    std::vector< FEFamily >     feFamily = {LAGRANGE, LAGRANGE, LAGRANGE};
    std::vector< FEOrder >       feOrder = {FIRST, SERENDIPITY, SECOND};
    std::vector< int >        time_order = {0, 0, 0};  //0 = steady, 2 = time-dependent
    std::vector< bool >   is_pde_unknown = {true, true, true};

    assert( feFamily.size() == feOrder.size());
    assert( feFamily.size() == is_pde_unknown.size());
    assert( feFamily.size() == time_order.size());

    std::vector< Unknown >  unknowns(feFamily.size());

    for (unsigned int fe = 0; fe < unknowns.size(); fe++) {

        std::ostringstream unk;
        unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
        unknowns[fe]._name           = unk.str();
        unknowns[fe]._fe_family      = feFamily[fe];
        unknowns[fe]._fe_order       = feOrder[fe];
        unknowns[fe]._time_order     = time_order[fe];
        unknowns[fe]._is_pde_unknown = is_pde_unknown[fe];

    }


    return unknowns;

}



const std::vector< Unknown >  systems__provide_list_of_unknowns_all_fe() {


    std::vector< FEFamily >     feFamily = {LAGRANGE, LAGRANGE, LAGRANGE, DISCONTINUOUS_POLYNOMIAL, DISCONTINUOUS_POLYNOMIAL};
    std::vector< FEOrder >       feOrder = {FIRST, SERENDIPITY, SECOND, ZERO, FIRST};
    std::vector< int >        time_order = {0, 0, 0, 0, 0};  //0 = steady, 2 = time-dependent
    std::vector< bool >   is_pde_unknown = {true, true, true, true, true};

    assert( feFamily.size() == feOrder.size());
    assert( feFamily.size() == is_pde_unknown.size());
    assert( feFamily.size() == time_order.size());

    std::vector< Unknown >  unknowns(feFamily.size());

    for (unsigned int fe = 0; fe < unknowns.size(); fe++) {

        std::ostringstream unk;
        unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
        unknowns[fe]._name           = unk.str();
        unknowns[fe]._fe_family      = feFamily[fe];
        unknowns[fe]._fe_order       = feOrder[fe];
        unknowns[fe]._time_order     = time_order[fe];
        unknowns[fe]._is_pde_unknown = is_pde_unknown[fe];

    }


    return unknowns;

}





namespace  Domain_square_01by01  {
    


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  = pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};


//this solution shows SUPERCONVERGENCE for SERENDIPITY FE, and it is like SUPER PERFECT for BIQUADRATIC FE... it is because of the MESH!
template < class type = double >
class Function_Zero_on_boundary_2 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = (1. - 2. * x[0]) *  x[1] * (1. - x[1]);
        solGrad[1]  = (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
    }



};


//this solution does not have SUPERCONVERGENCE even with the straight mesh
template < class type = double >
class Function_Zero_on_boundary_3 : public Math::Function< type > {

public:

// manufactured Laplacian =============
    type value(const std::vector < type >& x) const {
        
        return    x[0] *  x[0] * (1. - x[0]) * x[1] * (1. - x[1]) ;
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  x[0] * (1. - 2. * x[0]) *  x[1] * (1. - x[1]) +  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
        solGrad[1]  =  x[0] * (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return     x[0] *  (  -2. *  x[1] * (1. - x[1])   ) + 2. *  (1. - 2. * x[0]) *  x[1] * (1. - x[1])   +  x[0] * ( -2. *  x[0] * (1. - x[0])) ;
    }



};



} //end namespace


namespace  Domain_square_m05p05  {

 template < class type = double >
class Function_Zero_on_boundary_4 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
  return cos(pi * x[0]) * cos(pi * x[1]);
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
        solGrad[1]  = -pi * cos(pi * x[0]) * sin(pi * x[1]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};

   
    

std::string quad9_all_mesh_generation_methods(const unsigned method_flag, MultiLevelMesh & ml_mesh) {
    
  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");

    
  double scalingFactor = 1.;

    std::string mesh_name("square");
    
    switch(method_flag) {
        case 0: {
            mesh_name += "_femus";
  // ======= Mesh, coarse gen, I: from function - BEGIN  ========================
  const unsigned int nsub_x = 2;
  const unsigned int nsub_y = 2;
  const unsigned int nsub_z = 0;
  const std::vector<double> xyz_min = {-0.5,-0.5,0.};
  const std::vector<double> xyz_max = { 0.5, 0.5,0.};
  const ElemType geom_elem_type = QUAD9/*TRI6*/;
  ml_mesh.GenerateCoarseBoxMesh(nsub_x, nsub_y, nsub_z, xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], xyz_min[2], xyz_max[2], geom_elem_type, fe_quad_rule.c_str() );
  // ======= Mesh, coarse gen, I: from function - END ========================
 
            break;
        }
        
        case 1: {
            mesh_name += "_salome";
//   // ======= Mesh, coarse gen, III: from Salome - BEGIN  ========================
   ml_mesh.ReadCoarseMesh("./input/square_2x2_centered_at_origin.med", fe_quad_rule.c_str(), scalingFactor);
//   // ======= Mesh, coarse gen, III: from Salome - END ========================
            break;
        }
        case 2: {
            mesh_name += "_gambit";
//   // ======= Mesh, coarse gen, II: from Gambit - BEGIN  ========================
    ml_mesh.ReadCoarseMesh("./input/square_quad.neu", fe_quad_rule.c_str(), scalingFactor);
//     ml_mesh.ReadCoarseMesh("./input/square_tri.neu", fe_quad_rule.c_str(), scalingFactor);
//   //   ml_mesh.ReadCoarseMesh("./input/cube_tet.neu", fe_quad_rule.c_str(), scalingFactor);
  // ======= Mesh, coarse gen, II: from Gambit - END ========================
            break;
        }
        
        default: { abort(); }
    }
    
    
    return mesh_name;
    
}

    
}



#endif
