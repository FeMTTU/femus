#ifndef __femus_meshGencase_GeomElemBase_hpp__
#define __femus_meshGencase_GeomElemBase_hpp__


#include <string>
#include <vector>
#include <iostream>



namespace femus {



class GeomElemBase  {

public:

// ===  Constructors / Destructor - BEGIN =================
    GeomElemBase();
    
    virtual ~GeomElemBase();

    static  GeomElemBase* build(const std::string geomel_id_in, const uint fe_family);
// ===  Constructors / Destructor - END =================

    virtual unsigned int get_dimension() const { std::cout << "Not implemented FE" << std::endl; abort(); };
    virtual unsigned int n_nodes_linear() const { std::cout << "Not implemented FE" << std::endl; abort(); };
    
    virtual unsigned int n_nodes()       const { std::cout << "Not implemented FE" << std::endl; abort(); };
    
    virtual std::vector<unsigned> get_face(const unsigned f) const { std::cout << "Not implemented FE" << std::endl; abort();  };

    virtual std::string  get_name_med()  const { std::cout << "Not implemented FE" << std::endl; abort(); };
    virtual std::string  get_name_xdmf() const { std::cout << "Not implemented FE" << std::endl; abort(); };
    
// Refinement - BEGIN ===
public:
    
    virtual float get_embedding_matrix(const uint,const uint,const uint) = 0; //[/*NCHILDS*/][/*NNDS*/][/*NNDS*/]
    virtual double get_prol(const uint) = 0;
    
// Refinement - END ===

protected:  
    
   
};


} //end namespace femus



#endif
