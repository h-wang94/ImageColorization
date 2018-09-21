#ifndef TENSOR3D_HH
#define TENSOR3D_HH

#include "defs.h"


class Tensor3d {
public:
  /* Data Members */
  
  int len_x, len_y, len_z , plane, length ;
  real * data ;

  /* Methods */
  
  Tensor3d() { len_x = len_y = len_z = plane = length = 0 ; data = NULL ; }
  Tensor3d(int lx, int ly, int lz) ;
  Tensor3d(const Tensor3d& other) ;

  ~Tensor3d() { if(data) delete[] data ; data = NULL ; }
  void swap(Tensor3d& other) ;
  void empty() { len_x = len_y = plane = length = 0 ; if(data) delete[] data ; data = NULL ; }
  
  Tensor3d& operator=(const Tensor3d& other) ;

  
  real& operator[](int i) {
    
    //if (i>= (len_x*len_y*len_z)){
    //  cout << "OUT OF BOUNDS! , trying to access "<<i<<" elements num"<<len_x*len_y*len_z << endl ; 
    // exit(1) ; 
    //}
    
    return data[i] ;
  }
  real& operator()(int x, int y, int z) {
    
    //if(x >= len_x || y >= len_y || z>=len_z) { 
    //  cout << "OUT OF BOUNDS!, trying to access element ("<<x<<","<<y
    //   <<","<<z<<"), size=("<<len_x<<","<<len_y<<","<<len_z<<")"
    //	   << endl ; 
    //  exit(1) ; 
    //}
    
    return data[x + y * len_x + z * plane] ;
  }
  inline real operator()(double x, double y, double z) ;
  void one_minus() ;

  real cube(double x, double y, double z) ;
    
  void set(int lx, int ly, int lz) ;
  void clear(real val = 0) ;
  
  void restrict_full(Tensor3d& r) ;
  void restrict_half(Tensor3d& r) ;
  void restrict_none(Tensor3d& r) ;
  void prolong_full(Tensor3d& p) ;
  void prolong_none(Tensor3d& p) ;

  double norm_max() ;
  double norm_abs() ;
  double norm_euc() ;

  Tensor3d& operator*=(real val) ;
  Tensor3d& operator+=(Tensor3d& other) ;
  Tensor3d& operator-=(Tensor3d& other) ;
  
 
} ;

/******** Inline definintions ********/



real Tensor3d::operator()(double x, double y, double z) {
  int ix = (int) x ;
  int iy = (int) y ;
  int iz = (int) z ;

  double dx = x - ix ;
  double dy = y - iy ;
  double dz = z - iz ;
  double ddx = 1 - dx ;
  double ddy = 1 - dy ;
  double ddz = 1 - dz ;

  int i = ix + iy * len_x + iz * plane ;

  double vz = data[i] * ddx * ddy + data[i + 1] * dx * ddy +
    data[i + len_x] * ddx * dy + data[i + 1 + len_x] * dx * dy ;

  double vvz = data[i + plane] * ddx * ddy + data[i + 1 + plane] * dx * ddy +
    data[i + len_x + plane] * ddx * dy + data[i + 1 + len_x + plane] * dx * dy ;
  
  return vz * ddz + vvz * dz ;
}


#endif /* TENSOR3D_HH */
