#include "ElemType_template_base.hpp" 


    
   #include "ElemType_template.hpp"
  

namespace femus 

{
    


///@todo I have to put it here because I need to know what the children are... do separate files if needed
// run-time selection
// Since this is a TEMPLATED build function, I cannot put it in the cpp so easily...
// The only possibility is that I do EXPLICIT INSTANTIATION (see below)    
 template <class type, class type_mov>
       elem_type_templ_base<type, type_mov> * elem_type_templ_base<type, type_mov>::build(
                                         const std::string geom_elem, 
                                         const std::string fe_fam,
                                         const std::string order_gauss,
                                         const unsigned space_dimension) {

       
              if  ( geom_elem.compare("hex") == 0)     return  /**(*/ new elem_type_templ<type, type_mov, 3, 3>(geom_elem,  fe_fam, order_gauss) /*)*/;
              else if  (geom_elem.compare("tet") == 0)     return  /**(*/ new elem_type_templ<type, type_mov, 3, 3>(geom_elem,  fe_fam, order_gauss) /*)*/;
              else if  (geom_elem.compare("wedge") == 0)   return  /**(*/ new elem_type_templ<type, type_mov, 3, 3>(geom_elem,  fe_fam, order_gauss) /*)*/;
              else if  (geom_elem.compare("quad") == 0)    return  /**(*/ new elem_type_templ<type, type_mov, 2, 3>(geom_elem,  fe_fam, order_gauss) /*)*/;
              else if  (geom_elem.compare("tri") == 0)    return  /**(*/ new elem_type_templ<type, type_mov, 2, 3>(geom_elem,  fe_fam, order_gauss) /*)*/;
              else if  (geom_elem.compare("line") == 0)   return  /**(*/ new elem_type_templ<type, type_mov, 1, 3>(geom_elem,  fe_fam, order_gauss) /*)*/;
              else {std::cout << "Not implemented" << std::endl; abort(); }
          
          
      }
      
      
      
//  EXPLICIT INSTANTIATION
 template class elem_type_templ_base<double, double>;
 template class elem_type_templ_base<adept::adouble, double>;
 template class elem_type_templ_base<adept::adouble, adept::adouble>;
    

    
    
    
    
}
