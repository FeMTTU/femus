#ifndef __mgequationsmap_h__
#define __mgequationsmap_h__

#include <map>
#include <string>
using namespace std;

#include "Typedefs_conf.hpp"


#include "EqnBase.hpp"

class Utils;
class Physics;
class Mesh;
class FEElemBase;
class QRule;
class TimeLoop;

class QuantityMap;


class EquationsMap  {

  protected:

    map<string,EqnBase*> _equations;   // system map
public:  
    Utils&       _utils;  // Utils class  pointer
    Physics&     _phys;   // Physics class pointer
    QuantityMap& _qtymap;
    Mesh&        _mesh;   // Mesh  class  pointer
    std::vector<FEElemBase*>&  _AbstractFE;
    QRule&       _qrule;
    TimeLoop&    _timeloop;

  /// Constructor
    EquationsMap( Utils& mgutils_in,
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
  void PrintSol(const uint t_step,const double curr_time); ///< Print solution 
  void PrintCase(const uint t_init); ///< Print ic and bc
  void ReadSol(const uint t_step,double& time_out);  ///< Read solution //TODO must be updated
  
private:
  
 void PrintSolXDMF(const uint t_step,const double curr_time);
 void PrintSolHDF5(const uint t_flag);
 
 void PrintCaseXDMF(const uint t_init);
 void PrintCaseHDF5(const uint t_init);

 void PrintXDMFTopologyGeometry(std::ofstream& out,const uint Level, const uint vb);

 
};

#endif
