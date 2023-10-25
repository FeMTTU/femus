#ifndef SQUARE_OR_CUBE_OPT_SYSTEMS_INEQUALITIES_HPP
#define SQUARE_OR_CUBE_OPT_SYSTEMS_INEQUALITIES_HPP




//*********************** Domain and Mesh Dependent - BEGIN *****************************************
namespace femus {




   
namespace ctrl {
    
  namespace square_or_cube {
    
class mixed_state_or_ctrl_inequality {

public:

 ///@@@@@@@@@@@@@@@@@@@todo I believe you need to pass other rows here.........
    // Also, I think you the active flag should be frozen if it was frozen once...!!!
    //I don't think I have to put other mu's for the other constraint equations... there is a mu per constrained variable...
    //just make sure where they go
    //only try to constrain each component, and see how it behaves.
    //So, the active set is correctly found, but mu does not change, and it must! Maybe I didn't put all the pieces from the scalar to the vector routine?!
    //ctrl is not modified correctly...
    //I think the constraint is only acting on the interior nodes but not on the boundary, just check the boundary conditions!!!
 static std::vector<double>  InequalityConstraint(const unsigned n_components_ctrl, const std::vector<double> & dof_obj_coord, const bool upper) {

     const unsigned dim = dof_obj_coord.size();

     std::vector<double> constr_value(n_components_ctrl, 0.);


     double constr_value_upper_0 =  1000000.; // dof_obj_coord[1]*(1. - dof_obj_coord[1]);
     double constr_value_lower_0 = -1000000.; //-3.e-13;
     assert(constr_value_lower_0 < constr_value_upper_0);

     double constr_value_upper_1 =   1000.;
     double constr_value_lower_1 =  -1000.;
     assert(constr_value_lower_1 < constr_value_upper_1);

     double constr_value_upper_2 =  1000.;
     double constr_value_lower_2 = -1000.;
     assert(constr_value_lower_2 < constr_value_upper_2);

    if (upper)   {
                      constr_value[0] = constr_value_upper_0;
                      constr_value[1] = constr_value_upper_1;
       if (dim == 3)  constr_value[2] = constr_value_upper_2;
    }
    else         {
                      constr_value[0] = constr_value_lower_0;
                      constr_value[1] = constr_value_lower_1;
       if (dim == 3)  constr_value[2] = constr_value_lower_2;
    }


  return constr_value;

}


};

//*********************** Domain and Mesh Dependent - END *****************************************

  }



}




}















#endif
