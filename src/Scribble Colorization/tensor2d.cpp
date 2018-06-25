#include <iomanip.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "tensor2d.h"

Tensor2d::Tensor2d(int lx, int ly) {
  data = NULL ;
  set(lx, ly) ;
}

Tensor2d::Tensor2d(const Tensor2d& other) {
  data = NULL ;
  *this = other ;
}




void Tensor2d::swap(Tensor2d& other) {
   int tt = other.len_x ;
   other.len_x = len_x ;
   len_x = tt ;
   
   tt = other.len_y ;
   other.len_y = len_y ;
   len_y = tt ;

   real *t = data ; data = other.data ; other.data = t ;
} 

void Tensor2d::point(Tensor2d& other) {
   if(data)
     delete[] data ;
   
   len_x = other.len_x ;
   len_y = other.len_y ;
   
   data = other.data ;
}


Tensor2d& Tensor2d::operator=(const Tensor2d& other) {
  set(other.len_x, other.len_y) ; 
  int length = len_x * len_y ;
  
  for(int i = 0 ; i < length ; i++)
    data[i] = other.data[i] ;
  
  return *this ;
}


Tensor2d& Tensor2d::comb(Tensor2d& t1, Tensor2d& t2){
  set(t1.len_x,t1.len_y+t2.len_y);
  int length1 = t1.len_x * t1.len_y ;
  int length2 = t2.len_x * t2.len_y ;
  int i;

  for(i = 0 ; i < length1 ; i++)
    data[i] = t1.data[i] ;
  for(i = 0 ; i < length2 ; i++)
    data[i+length1] = t2.data[i] ;
  
  return *this ;
  
}




void Tensor2d::set(int lx, int ly) {
   if(data && lx == len_x && ly == len_y)
     return ;
     
   len_x = lx ; len_y = ly ;

   int length = lx * ly ;
   
  if(data)
    delete[] data ;

  data = new real[length] ;

  if (!data) {
    cerr << "Allocation error of Tensor2d. (size req - " << length << ")" << endl ;
    exit(1) ;
  }
}



void Tensor2d::set(int lx, int ly, real* datap) {
   
     
   len_x = lx ; len_y = ly ;

   int length = lx * ly ;
   
  if(data)
    delete[] data ;

  //data = new real[length] ;
  data=datap;
  
  if (!data) {
    cerr << "Allocation error of Tensor2d. (size req - " << length << ")" << endl ;
    exit(1) ;
  }
}

void Tensor2d::clear(real val ) {
  int length = len_x * len_y ;
  for(int i = 0 ; i < length ; i++)
    data[i] = val ;
}


void Tensor2d::restrict_full(Tensor2d& r) {
  r.set((len_x - 1) / 2 + 1, (len_y - 1) / 2 + 1) ;
   
  int i, j ;
  int x,y;

  for( y = 1 ; y < r.len_y - 1 ; y++) {
    for( x = 1 ; x < r.len_x - 1 ; x++) {

      i = x + y * r.len_x ;
      j = 2 * x + 2 * y * len_x ;

      r.data[i] = 0.25 * data[j] + 0.125 * (data[j+1] + data[j-1] + data[j+len_x] + data[j-len_x]) +
	0.0625 * (data[j+1+len_x] + data[j+1-len_x] + data[j-1+len_x] + data[j-1-len_x]) ;
    }
  }

  for( x = 1 ; x < r.len_x - 1 ; x++) {
    i = x ;
    j = 2 * x ;

    r.data[i] = 0.333 * data[j] + 0.1666 * (data[j+1] + data[j-1] + data[j+len_x]) +
      0.08333 * (data[j+1+len_x] + data[j-1+len_x]) ;

    i += (r.len_y-1) * r.len_x ;
    j += (len_y-1) * len_x ;

    r.data[i] = 0.333 * data[j] + 0.1666 * (data[j+1] + data[j-1] + data[j-len_x]) +
      0.08333 * (data[j+1-len_x] + data[j-1-len_x]) ;
  }
  
  for( y = 1 ; y < r.len_y - 1 ; y++) {
    i = y * r.len_x ;
    j = 2 * y * len_x ;
    
    r.data[i] = 0.333 * data[j] + 0.1666 * (data[j+1] + data[j-len_x] + data[j+len_x]) +
      0.08333 * (data[j+1+len_x] + data[j+1-len_x]) ;

    i += r.len_x - 1 ;
    j += len_x - 1 ;

     r.data[i] = 0.333 * data[j] + 0.1666 * (data[j-1] + data[j-len_x] + data[j+len_x]) +
      0.08333 * (data[j-1+len_x] + data[j-1-len_x]) ;
  }

  r.data[0] = 0.444 * data[0] + 0.222 * (data[1] + data[len_x]) + 0.111 * data[1+len_x] ;

  r.data[r.len_x-1] = 0.444 * data[len_x-1] + 0.222 * (data[len_x-2] + data[2*len_x-1]) + 0.111 * data[2*len_x-2] ;

  r.data[(r.len_y-1)*r.len_x] = 0.444 * data[(len_y-1)*len_x] + 0.222 * (data[(len_y-1)*len_x+1] + data[(len_y-2)*len_x]) +
    0.111 * data[(len_y-2)*len_x+1] ;

  r.data[(r.len_y-1)*r.len_x + r.len_x-1] = 0.444 * data[(len_y-1)*len_x+len_x-1] +
    0.222 * (data[(len_y-1)*len_x+len_x-2] + data[(len_y-1)*len_x-1]) + 0.111 * data[(len_y-1)*len_x-2] ;
}

void Tensor2d::restrict_half(Tensor2d& r) {
  r.set((len_x - 1) / 2 + 1, (len_y - 1) / 2 + 1) ;
  
  int i, j, x,y ;

  for( y = 1 ; y < r.len_y - 1 ; y++) {
    for( x = 1 ; x < r.len_x - 1 ; x++) {

      i = x + y * r.len_x ;
      j = 2 * x + 2 * y * len_x ;

      r.data[i] = 0.3333 * data[j] + 0.16666 * (data[j+1] + data[j-1] + data[j+len_x] + data[j-len_x]) ;
    }
  }

  for( x = 1 ; x < r.len_x - 1 ; x++) {
    i = x ;
    j = 2 * x ;

    r.data[i] = 0.4 * data[j] + 0.2 * (data[j+1] + data[j-1] + data[j+len_x]) ;

    i += (r.len_y-1) * r.len_x ;
    j += (len_y-1) * len_x ;

    r.data[i] = 0.4 * data[j] + 0.2 * (data[j+1] + data[j-1] + data[j-len_x]) ;
  }
  
  for( y = 1 ; y < r.len_y - 1 ; y++) {
    i = y * r.len_x ;
    j = 2 * y * len_x ;
    
    r.data[i] = 0.4 * data[j] + 0.2 * (data[j+1] + data[j-len_x] + data[j+len_x]) ;

    i += r.len_x - 1 ;
    j += len_x - 1 ;

    r.data[i] = 0.4 * data[j] + 0.2 * (data[j-1] + data[j-len_x] + data[j+len_x]) ;
  }

  r.data[0] = 0.5 * data[0] + 0.25 * (data[1] + data[len_x]) ;

  r.data[r.len_x-1] = 0.5 * data[len_x-1] + 0.25 * (data[len_x-2] + data[2*len_x-1]) ;

  r.data[(r.len_y-1)*r.len_x] = 0.5 * data[(len_y-1)*len_x] + 0.25 * (data[(len_y-1)*len_x+1] + data[(len_y-2)*len_x]) ;

  r.data[(r.len_y-1)*r.len_x + r.len_x-1] = 0.5 * data[(len_y-1)*len_x+len_x-1] +
    0.25 * (data[(len_y-1)*len_x+len_x-2] + data[(len_y-1)*len_x-1]) ;
}
  
void Tensor2d::restrict_none(Tensor2d& r) {
  r.set((len_x - 1) / 2 + 1, (len_y - 1) / 2 + 1) ;
  
  for(int y = 0 ; y < r.len_y ; y++) {
    for(int x = 0 ; x < r.len_x ; x++) {
      r.data[x+y*r.len_x] = data[2*x+2*y*len_x] ;
    }
  }
}

void Tensor2d::prolong_full(Tensor2d& p) {
  p.set((len_x - 1) * 2 + 1, (len_y - 1) * 2 + 1) ;

  p.data[0] = data[0] ;
  int x,y;
  for( y = 1 ; y < len_y ; y++) {
    p.data[y*2*p.len_x] = data[y*len_x] ;
    p.data[(y*2-1)*p.len_x] = 0.5 * (data[(y-1)*len_x] + data[y*len_x]) ;
  }
  
  for( x = 1 ; x < len_x ; x++) {
    p.data[x*2] = data[x] ;
    p.data[x*2-1] = 0.5 * (data[x-1] + data[x]) ;
  }

  int i, j ;

  for( y = 1 ; y < len_y ; y++) {
    for( x = 1 ; x < len_x ; x++) {

      i = x + y * len_x ;
      j = 2 * x + 2 * y * p.len_x ;

      p.data[j] = data[i] ;

      p.data[j-1] = 0.5 * (data[i] + data[i-1]) ;
      p.data[j-p.len_x] = 0.5 * (data[i] + data[i-len_x]) ;
      
      p.data[j-1-p.len_x] = 0.25 * (data[i] + data[i-1] + data[i-len_x] + data[i-len_x-1]) ;
    }
  }
}

void Tensor2d::prolong_none(Tensor2d& p) {
  p.set((len_x - 1) * 2 + 1, (len_y - 1) * 2 + 1) ;

  int i, j, y,x ;

  real val ;

  p.data[0] = data[0] ;
 
  for( y = 1 ; y < len_y ; y++) {
    p.data[y*2*p.len_x] = data[y*len_x] ;
    p.data[(y*2-1)*p.len_x] = data[y*len_x] ;
  }
  
  for( x = 1 ; x < len_x ; x++) {
    p.data[x*2] = data[x] ;
    p.data[x*2-1] = data[x] ;
  }
  
  for( y = 1 ; y < len_y ; y++) {
    for( x = 1 ; x < len_x ; x++) {
      
      i = x + y * len_x ;
      j = 2 * x + 2 * y * p.len_x ;

      val = data[i] ;
       
      p.data[j] = val ;
      
      p.data[j-1] = val ;
      p.data[j-p.len_x] = val ;
      
      p.data[j-1-p.len_x] = val ;
    }
  }
}

void Tensor2d::save_text(const char * filename) {
  ofstream file(filename) ;

  file << setiosflags(ios::fixed) ;
  file << setprecision(50) ;

  for(int y = 0 ; y < len_y ; y++) {
    for(int x = 0 ; x < len_x ; x++) {
      file << (*this)(x,y) << " " ;
    }
    file << endl ;
  }

  file.close() ;
}
   


void Tensor2d::read_text(const char * filename, int lx, int ly)
{
   
   
      if(data)
          delete[] data ;
   
      data = new real[lx*ly] ;
   
      len_x = lx ;
      len_y = ly ;
   
        ifstream file(filename) ;
      int i = 0 ;
        for(int y = 0 ; y < len_y ; y++)
     {
	
	
	            for(int x = 0 ; x < len_x ; x++)
	  {
	     
	     
	                        file >> data[i++] ;
	  }
	
	
     }
   
   
   
        file.close() ;
}



                           

double Tensor2d::norm_max() {
  double norm = 0 ;

  int length = len_x * len_y ;
  
  for(int i = 0 ; i < length ; i++)
    if(fabs(data[i]) > norm)
      norm = fabs(data[i]) ;

  return norm ;
}
void Tensor2d::Pow(double c) {
int length = len_x * len_y ;
 
  for(int i = 0 ; i < length ; i++)
     data[i] = pow(data[i],c) ;
}

Tensor2d& Tensor2d::operator*=(real val) {
  int length = len_x * len_y ;
 
  for(int i = 0 ; i < length ; i++)
     data[i] *= val ;
  return *this ;  
}

/*
Tensor2d& Tensor2d::operator+=(Tensor2d& other) {

  int length = len_x * len_y ;

  if(len_x != other.len_x || len_y != other.len_y ) {
    cerr << "Vectors length mismatch! (op+)" << endl ;
    exit(1) ;
  }
  for(int i = 0 ; i < length ; i++)
    data[i] += other.data[i] ;
  return *this ;
}
*/

Tensor2d& Tensor2d::operator+=(Tensor2d& other) {

  int length = len_x * len_y ;

  //if(len_x != other.len_x || len_y != other.len_y ) {
  //  cerr << "Vectors length mismatch! (op+)" << endl ;
  //  exit(1) ;
  //}
  int m_len_x, m_len_y;
  if (len_x <= other.len_x){ m_len_x=len_x;}else{m_len_x=other.len_x;}
  if (len_y <= other.len_y){ m_len_y=len_y;}else{m_len_y=other.len_y;}
  
  for(int y = 0 ; y < m_len_y ; y++){
    for(int x = 0 ; x < m_len_x ;x++){
      data[x + y * len_x]+=other.data[x + y * other.len_x];
    }
  }
  return *this ;
}
/*
Tensor2d& Tensor2d::operator-=(Tensor2d& other) {

  int length = len_x * len_y ;

  if(len_x != other.len_x || len_y != other.len_y ) {
    cerr << "Vectors length mismatch! (op-)" << endl ;
    exit(1) ;
  }
  for(int i = 0 ; i < length ; i++)
    data[i] -= other.data[i] ;

  return *this ;
}
*/
Tensor2d& Tensor2d::operator-=(Tensor2d& other) {

  int length = len_x * len_y ;

  //if(len_x != other.len_x || len_y != other.len_y ) {
  //  cerr << "Vectors length mismatch! (op+)" << endl ;
  //  exit(1) ;
  //}
  int m_len_x, m_len_y;
  if (len_x <= other.len_x){ m_len_x=len_x;}else{m_len_x=other.len_x;}
  if (len_y <= other.len_y){ m_len_y=len_y;}else{m_len_y=other.len_y;}
  
  for(int y = 0 ; y < m_len_y ; y++){
    for(int x = 0 ; x < m_len_x ;x++){
      data[x + y * len_x]-=other.data[x + y * other.len_x];
    }
  }
  return *this ;
}


double Tensor2d::norm_abs() {
  double norm = 0 ;

  int length = len_x * len_y ;
  
  for (int i = 0 ; i < length ; i++)
    norm += fabs(data[i]) ;

  return norm / length ;
}

double Tensor2d::norm_euc() {
  double norm = 0 ;

  int length = len_x * len_y ;
  
  for (int i = 0 ; i < length ; i++)
    norm += data[i] * data[i] ;

  return norm / length ;
}

real Tensor2d::min() 
{
   double min = 1e50 ;

   int length = len_x * len_y ;
   
   for(int i = 0 ; i < length ; i++)
     if(min > data[i])
       min = data[i] ;
   
   return min;
}

real Tensor2d::max() 
{
   double max = -1e50 ;

   int length = len_x * len_y ;
   
   for(int i = 0 ; i < length ; i++)
          if(max < data[i])
              max = data[i] ;
   
      return max;
     
}


