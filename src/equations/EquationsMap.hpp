#ifndef __mgequationsmap_h__
#define __mgequationsmap_h__

#include <map>
#include <string>
using namespace std;

#include "Typedefs.hpp"


#include "EqnBase.hpp"
#include "GaussPoints.hpp"


namespace femus {



class Utils;
class Physics;
class MeshTwo;
class FEElemBase;
class elem_type;
class TimeLoop;

class QuantityMap;


class EquationsMap  {

public:
  
    Files&       _files;
    Physics&     _phys;
    QuantityMap& _qtymap;
    MeshTwo&     _mesh;
    std::vector<FEElemBase*> &  _AbstractFE;
    std::vector< std::vector<elem_type*> >  &  _elem_type;
    std::vector<Gauss>       _qrule;

  /// Constructor
    EquationsMap( Files& files_in,
		  Physics& mgphys_in,
		  QuantityMap& qtymap_in,
		  MeshTwo& mgmesh_in,
		  std::vector<FEElemBase*> & absfe_in,
                  std::vector< std::vector<elem_type*> > & elem_type_in,
		  std::vector<Gauss> qrule_in
		);

    
    /// Destructor
  ~EquationsMap(){};
  
  void clean(); ///< Clean all substructures

  // equation get/set
  inline          void  set_eqs(EqnBase* value)            {_equations.insert(make_pair(value->_eqname,value));}
  inline       EqnBase* get_eqs(const string & name)       {return _equations.find(name)->second;}
  inline const EqnBase* get_eqs(const string & name) const {return _equations.find(name)->second;}

  typedef std::map<string, EqnBase*>::iterator iterator;
  typedef std::map<string, EqnBase*>::const_iterator const_iterator;

  inline iterator       begin()       { return _equations.begin();}
  inline iterator         end()       { return _equations.end();}
  inline const_iterator begin() const { return _equations.begin();}
  inline const_iterator   end() const { return _equations.end();}


  void setDofBcOpIc() ;


  // read-print functions -------------------------------------------
  void PrintSol(const uint t_step,const double curr_time) const; ///< Print solution 
  void ReadSol(const uint t_step,double& time_out) const;  ///< Read solution //TODO must be updated
  void PrintCase(const uint t_init) const; ///< Print ic and bc
  
private:
  
 map<string,EqnBase*> _equations;   // system map
    
 void PrintSolXDMF(const uint t_step,const double curr_time) const;
 void PrintSolHDF5(const uint t_flag) const;
 
 void PrintCaseXDMF(const uint t_init) const;
 void PrintCaseHDF5(const uint t_init) const;

 void PrintXDMFTopologyGeometry(std::ofstream& out,const uint Level, const uint vb) const;

 
};


} //end namespace femus



#endif
