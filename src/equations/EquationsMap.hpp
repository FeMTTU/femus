#ifndef __mgequationsmap_h__
#define __mgequationsmap_h__

#include <map>
#include <string>
using namespace std;

#include "Typedefs.hpp"


#include "EqnBase.hpp"


namespace femus {



class Utils;
class Physics;
class Mesh;
class FEElemBase;
class QRule;
class TimeLoop;

class QuantityMap;


class EquationsMap  {

public:
  
    Files&       _files;
    Physics&     _phys;
    QuantityMap& _qtymap;
    Mesh&        _mesh;
    std::vector<FEElemBase*>&  _AbstractFE;
    QRule&       _qrule;
    TimeLoop&    _timeloop;

  /// Constructor
    EquationsMap( Files& files_in,
		  Physics& mgphys_in,
		  QuantityMap& qtymap_in,
		  Mesh& mgmesh_in,
		  std::vector<FEElemBase*>&  absfe_in,
		  QRule& qrule_in,
		  TimeLoop& timeloop_in
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
  void OneTimestepEqnLoop(const double time, const uint delta_t_step_in);

  void TransientSetup();  //initialization of all the equations in the map
  void TransientLoop();  //a standard transient loop in alphabetical order


  // read-print functions -------------------------------------------
  void PrintSol(const uint t_step,const double curr_time) const; ///< Print solution 
  void ReadSol(const uint t_step,double& time_out);  ///< Read solution //TODO must be updated
  
private:
  
 map<string,EqnBase*> _equations;   // system map
    
 void PrintSolXDMF(const uint t_step,const double curr_time) const;
 void PrintSolHDF5(const uint t_flag) const;
 
 void PrintCase(const uint t_init) const; ///< Print ic and bc
 void PrintCaseXDMF(const uint t_init) const;
 void PrintCaseHDF5(const uint t_init) const;

 void PrintXDMFTopologyGeometry(std::ofstream& out,const uint Level, const uint vb) const;

 
};


} //end namespace femus



#endif
