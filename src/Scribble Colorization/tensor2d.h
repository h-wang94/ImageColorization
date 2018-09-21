#ifndef TENSOR2D_HH
#define TENSOR2D_HH

#include "defs.h"

#include <iostream.h>
#include <stdlib.h>
class Tensor2d {
public:
  /* Members */
  
  int len_x, len_y ;
  real * data ;

  /* Methods */

  Tensor2d() { len_x = len_y = 0 ; data = NULL ; }
  Tensor2d(int lx, int ly) ;
  Tensor2d(const Tensor2d& other) ;

  ~Tensor2d() { if(data) delete[] data ; data = NULL ; }

  void point(Tensor2d& other) ;
  void swap(Tensor2d& other) ;
  void empty() { len_x = len_y = 0 ; if(data) delete[] data ; data = NULL ; }
  Tensor2d& operator=(const Tensor2d& other) ;
  
  Tensor2d& comb(Tensor2d& t1, Tensor2d& t2);
  
  real& operator()(int x, int y) {
    //if(x >= len_x || y >= len_y) {
    //   cout << "OUT OF BOUNDS!, trying to access element ("<<x<<","<<y
    //	   <<"), size=("<<len_x<<","<<len_y<<")"
    //	   << endl ; 
    //    exit(1) ; 
    //} 
    return data[x + y * len_x] ;
  }
  
  real& operator[](int i) {
    //if (i>=( len_y* len_x)){
    //  cout << "OUT OF BOUNDS!" << endl ; 
    //  exit(1) ; 
    //}
    return data[i] ;
  }
  inline real operator()(double x, double y) ;

  void set(int lx, int ly) ;
  void set(int lx, int ly, real* datap) ;
  void clear(real val = 0) ;

  
  void save_text(const char * filename) ;
  void read_text(const char * filename, int lx, int ly) ;


  void restrict_full(Tensor2d& r) ;
  void restrict_half(Tensor2d& r) ;
  void restrict_none(Tensor2d& r) ;
  void prolong_full(Tensor2d& p) ;
  void prolong_none(Tensor2d& p) ;

  void Pow(double c) ;
  
  double norm_max() ;
  double norm_abs() ;
  double norm_euc() ;

  real min() ;
  real max() ;
  
  Tensor2d& operator*=(real val) ;
  Tensor2d& operator+=(Tensor2d& other) ;
  Tensor2d& operator-=(Tensor2d& other) ;

  
  
} ;

real Tensor2d::operator()(double x, double y) {
int ix = (int)x ;
int iy = (int)y ;
double dx = x - ix ;
double dy = y - iy ;

return (*this)(ix,iy) * (1-dx) * (1-dy) + (*this)(ix+1,iy) * (dx) *
(1-dy) +(*this)(ix,iy+1) * (1-dx) * (dy) +(*this)(ix+1,iy+1) * (dx) *
(dy)  ;


}

#endif /* TENSOR2D_HH */
