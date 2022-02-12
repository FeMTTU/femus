#include "ElemSto.hpp"


namespace femus {




//=====================
  ElemStoBase::ElemStoBase(int nnds_in, uint spacedim_in) {
    
    _nnds =  nnds_in;
    _elnds = new int[nnds_in];
    _spacedim = spacedim_in;
    
  }


//===== may it happen that you dont call the father destructor?
  ElemStoBase::~ElemStoBase() { delete [] _elnds;  }

  
  //========================
  ElemStoBdry::ElemStoBdry(int nnds_in, uint spacedim_in) : ElemStoBase(nnds_in,spacedim_in)  {
 
 _id=  0  ; 
_lev=  0  ;
_subd = 0;
     for (int i=0; i<nnds_in; i++) { _elnds[i]=0;  }

_vol_id=  0  ;
_nside=  0  ;
 
  
  }
  
//==============================
ElemStoBdry::~ElemStoBdry() {  }



//========================
  ElemStoVol::ElemStoVol(int nnds_in, uint spacedim_in) : ElemStoBase(nnds_in,spacedim_in)  {
    
   _id = -1;
  _lev = -1;
  for (int i=0; i<nnds_in; i++) { _elnds[i]=-1;  }

  _subd = -1;
  _par = -1;
  _nch = -1;
    
    uint nchd =  4*(_spacedim-1);
    _elchs = new int[nchd];
      for (uint i=0; i<nchd; i++) { _elchs[i]=-1;  }
      
    
  }

//==============================
ElemStoVol::~ElemStoVol() {  }


//=============================
     NodeSto::NodeSto(int nl1_in,int n_levs_in){
       
          _id = -1  ;
         _lev = nl1_in ;
        _levp = n_levs_in ;
        _subd = -1 ;
         _var = -1 ;

    }
     NodeSto::~NodeSto(){}




} //end namespace femus

