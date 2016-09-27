
#ifndef __femus_mesh_MeshParallelVector_hpp__
#define __femus_mesh_MeshParallelVector_hpp__


namespace femus {

  
  class MeshParalleVector {
    MeshParalleVector(const unsigned &typeSize, const unsigned &size){
      _typeSIZE = typeSize;
      //_gblMemory = new void [typeSize*size];
    }
  private:
    void *_gblMemory;
    void *_localMemory;
    void *_localizedFromJprocMemory;
    unsigned _typeSIZE;
    unsigned _globalSize;
  };
  
  
} //end namespace femus

#endif