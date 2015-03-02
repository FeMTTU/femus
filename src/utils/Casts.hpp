#ifndef __femus_utils_Casts_hpp__
#define __femus_utils_Casts_hpp__

//The appropriate includes...
#include <cassert>

#include <typeinfo>


namespace femus {


using std::bad_cast;

//TODO
template <typename Tnew, typename Told>
inline Tnew libmeshM_cast_ptr (Told* oldvar) {
#ifndef NDEBUG
  Tnew newvar = dynamic_cast<Tnew>(oldvar);
  assert (newvar);
  return newvar;
#else
  return(static_cast<Tnew>(oldvar));
#endif
}


//***********PORCHETTA
template <typename Tnew, typename Told>
inline Tnew libmeshM_cast_ref(Told& oldvar) {
#ifndef NDEBUG
  try {
    Tnew newvar = dynamic_cast<Tnew>(oldvar);
    return newvar;
  } catch (std::bad_cast) {
    assert (false);
  }
#else
  return(static_cast<Tnew>(oldvar));
#endif
}

//**********PORCHETTA

//*********PORCHETTA che non va qui
template <typename ForwardIter, typename T>
void iota (ForwardIter first, ForwardIter last, T value) {
  while (first != last) {
    *first = value++;
    ++first;
  }
}



} //end namespace femus



#endif