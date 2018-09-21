#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "tensor3d.h"

Tensor3d::Tensor3d(int lx, int ly, int lz) {
  data = NULL ;
  set(lx, ly, lz) ;
}

Tensor3d::Tensor3d(const Tensor3d& other) {
  data = NULL ;
  *this = other ;
}


real Tensor3d::cube(double x, double y, double z) {
  int i,j;
  int ix = (int) x ;
  int iy = (int) y ;
  int iz = (int) z ;

  double dx = x - ix ;
  double dy = y - iy ;
  double dz = z - iz ;

  ix -- ;
  iy -- ;
  iz -- ;

  if(ix < 0 || iy < 0 || iz < 0 || ix+3 >= len_x || iy+3 >= len_y || iz+3 >= len_z)
    return 0 ;

  double A = dz * (dz * (-0.5 * dz + 1) - 0.5) ;
  double B = dz * (dz * (1.5 * dz - 2.5)) + 1 ;
  double C = dz * (dz * (-1.5 * dz + 2) + 0.5) ;
  double D = dz * (dz * (0.5 * dz - 0.5)) ;

  double p[16] ;

  for( i = 0 ; i < 4 ; i++)
    for( j = 0 ; j < 4 ; j++)
      
      p[j+i*4] = operator()(ix + i, iy + j, iz) * A +
	operator()(ix + i, iy + j, iz+1) * B + 
	operator()(ix + i, iy + j, iz+2) * C + 
	operator()(ix + i, iy + j, iz+3) * D ;
  

  A = dy * (dy * (-0.5 * dy + 1) - 0.5) ;
  B = dy * (dy * (1.5 * dy - 2.5)) + 1 ;
  C = dy * (dy * (-1.5 * dy + 2) + 0.5) ;
  D = dy * (dy * (0.5 * dy - 0.5)) ;

  double l[4] ;

  for( i = 0 ; i < 4 ; i++)
    l[i] = p[i*4] * A +  p[i*4+1] * B +  p[i*4+2] * C +  p[i*4+3] * D ;
  
  A = dx * (dx * (-0.5 * dx + 1) - 0.5) ;
  B = dx * (dx * (1.5 * dx - 2.5)) + 1 ;
  C = dx * (dx * (-1.5 * dx + 2) + 0.5) ;
  D = dx * (dx * (0.5 * dx - 0.5)) ;

  return l[0] * A + l[1] * B + l[2] * C + l[3] * D ;
}


void Tensor3d::swap(Tensor3d& other) {
  int tt = other.len_x ;
  
  other.len_x = len_x ;
  len_x = tt ;
  
  tt = other.len_y ;
  other.len_y = len_y ;
  len_y = tt ;

  tt = other.plane ;
  other.plane = plane ;
  plane = tt ;

  tt = other.length ;
  other.length = length ;
  length = tt ;
  
  real *t = data ; data = other.data ; other.data = t ;
}

Tensor3d& Tensor3d::operator=(const Tensor3d& other) {
  len_x = other.len_x ;
  len_y = other.len_y ;
  len_z = other.len_z ;

  length = len_x * len_y * len_z ;
  plane = len_x * len_y ;
  
  if(data)
    delete[] data ;

  data = new real[length] ;

  if (!data) {
    cerr << "Allocation error of Tensor3d. (size req - " << length << ")" << endl ;
    exit(1) ;
  }

  for(int i = 0 ; i < length ; i++)
    data[i] = other.data[i] ;

  return *this ;
}




void Tensor3d::set(int lx, int ly, int lz) {
  if(data && lx == len_x && ly == len_y && lz == len_z)
    return ;
  
  len_x = lx ; len_y = ly ; len_z = lz ;
  
  length = len_x * len_y * len_z ;
  plane = len_x * len_y ;
  
  
  if(data)
    delete[] data ;
  
  data = new real[length] ;
  
  if (!data) {
    cerr << "Allocation error of Tensor3d. (size req - " << length << ")" << endl ;
    exit(1) ;
  }
}

void Tensor3d::clear(real val) {
  for(int i = 0 ; i < length ; i++)
    data[i] = val ;
}

void Tensor3d::restrict_full(Tensor3d& r) {
  r.set((len_x - 1) / 2 + 1, (len_y - 1) / 2 + 1, (len_z - 1) / 2 + 1) ;

  int i, j,x,y,z ;
    
  for( z = 1 ; z < r.len_z - 1 ; z++) {
    for(y = 1 ; y < r.len_y - 1 ; y++) {
      for( x = 1 ; x < r.len_x - 1 ; x++) {
	i = x + y * r.len_x + z * r.plane ;
	j = 2 * x + 2 * y * len_x + 2 * z * plane ;

	r.data[i] = 0.125 * data[j]
	  + 0.0625 * (data[j-1] + data[j+1] + data[j-len_x] + data[j+len_x] + data[j-plane] + data[j+plane])
	  + 0.03125 * (data[j-1-len_x] + data[j+1-len_x] + data[j-1+len_x] + data[j+1+len_x]
		       + data[j-1-plane] + data[j+1-plane] + data[j-1+plane] + data[j+1+plane]
		       + data[j-len_x-plane] + data[j+len_x-plane] + data[j-len_x+plane] + data[j+len_x+plane])
	  + 0.015625 * (data[j-1-len_x-plane] + data[j-1-len_x+plane] + data[j+1-len_x-plane] + data[j+1-len_x+plane] +
			data[j-1+len_x-plane] + data[j-1+len_x+plane] + data[j+1+len_x-plane] + data[j+1+len_x+plane]) ;
      }
    }
  }

  /* Z-Plane (XY) */
  
  for( y = 1 ; y < r.len_y - 1 ; y++) {
    for( x = 1 ; x < r.len_x - 1 ; x++) {
      i = x + y * r.len_x ;
      j = 2 * x + 2 * y * len_x ;

      r.data[i] = 0.1666667 * data[j]
	        + 0.0833333 * (data[j-1] + data[j+1] + data[j-len_x] + data[j+len_x] + data[j+plane])
	        + 0.0416667 * (data[j-1+len_x] + data[j+1+len_x] + data[j-1-len_x] + data[j+1-len_x] +
			       data[j-1+plane] + data[j+1+plane] + data[j+plane-len_x] + data[j+plane+len_x])
   	        + 0.0208333 * (data[j-1-len_x+plane] + data[j+1-len_x+plane] + data[j-1+len_x+plane] + data[j+1+len_x+plane]) ;

      i += (r.len_z - 1) * r.plane ;
      j += (len_z - 1) * plane ;
      
      r.data[i] = 0.1666667 * data[j]
	        + 0.0833333 * (data[j-1] + data[j+1] + data[j-len_x] + data[j+len_x] + data[j-plane])
	        + 0.0416667 * (data[j-1+len_x] + data[j+1+len_x] + data[j-1-len_x] + data[j+1-len_x] +
			       data[j-1-plane] + data[j+1-plane] + data[j-plane-len_x] + data[j-plane+len_x])
   	        + 0.0208333 * (data[j-1-len_x-plane] + data[j+1-len_x-plane] + data[j-1+len_x-plane] + data[j+1+len_x-plane]) ;
    }
  }

  /* Y-Plane (ZX) */
  
  for( z = 1 ; z < r.len_z - 1 ; z++) {
    for( x = 1 ; x < r.len_x - 1 ; x++) {
      i = x + z * r.plane ;
      j = 2 * x + 2 * z * plane ;

      r.data[i] = 0.1666667 * data[j]
	        + 0.0833333 * (data[j-1] + data[j+1] + data[j-plane] + data[j+plane] + data[j+len_x])
	        + 0.0416667 * (data[j-1+plane] + data[j+1+plane] + data[j-1-plane] + data[j+1-plane] +
			       data[j-1+len_x] + data[j+1+len_x] + data[j+len_x-plane] + data[j+len_x+plane])
   	        + 0.0208333 * (data[j-1-plane+len_x] + data[j+1-plane+len_x] + data[j-1+plane+len_x] + data[j+1+plane+len_x]) ;

      i += (r.len_y-1) * r.len_x ;
      j += (len_y-1) * len_x ;
      
      r.data[i] = 0.1666667 * data[j]
	        + 0.0833333 * (data[j-1] + data[j+1] + data[j-plane] + data[j+plane] + data[j-len_x])
	        + 0.0416667 * (data[j-1+plane] + data[j+1+plane] + data[j-1-plane] + data[j+1-plane] +
			       data[j-1-len_x] + data[j+1-len_x] + data[j-len_x-plane] + data[j-len_x+plane])
   	        + 0.0208333 * (data[j-1-plane-len_x] + data[j+1-plane-len_x] + data[j-1+plane-len_x] + data[j+1+plane-len_x]) ;
    }
  }

  /* X-Plane (ZY) */
  
  for( z = 1 ; z < r.len_z - 1 ; z++) {
    for( y = 1 ; y < r.len_y - 1 ; y++) {
      i = y * r.len_x + z * r.plane ;
      j = 2 * y * len_x + 2 * z * plane ;

      r.data[i] = 0.1666667 * data[j]
	        + 0.0833333 * (data[j-len_x] + data[j+len_x] + data[j-plane] + data[j+plane] + data[j+1])
	        + 0.0416667 * (data[j-len_x+plane] + data[j+len_x+plane] + data[j-len_x-plane] + data[j+len_x-plane] +
			       data[j-len_x+1] + data[j+len_x+1] + data[j+1-plane] + data[j+1+plane])
   	        + 0.0208333 * (data[j-len_x-plane+1] + data[j+len_x-plane+1] + data[j-len_x+plane+1] + data[j+len_x+plane+1]) ;

      i += r.len_x - 1 ;
      j += len_x - 1 ;
      
      r.data[i] = 0.1666667 * data[j]
	        + 0.0833333 * (data[j-len_x] + data[j+len_x] + data[j-plane] + data[j+plane] + data[j-1])
	        + 0.0416667 * (data[j-len_x+plane] + data[j+len_x+plane] + data[j-len_x-plane] + data[j+len_x-plane] +
			       data[j-len_x-1] + data[j+len_x-1] + data[j-1-plane] + data[j-1+plane])
   	        + 0.0208333 * (data[j-len_x-plane-1] + data[j+len_x-plane-1] + data[j-len_x+plane-1] + data[j+len_x+plane-1]) ;
    }
  }

  /* X-Axis */
  
  for( x = 1 ; x < r.len_x - 1 ; x++) {
    i = x ;
    j = 2 * x ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+1] + data[j-1] + data[j+len_x] + data[j+plane])
      + 0.055556 * (data[j+len_x+1] + data[j+len_x-1] + data[j+plane+1] + data[j+plane-1] + data[j+plane+len_x])
      + 0.027778 * (data[j+len_x+1+plane] + data[j+len_x-1+plane]) ;

    i = x + (r.len_y-1) * r.len_x ;
    j = 2 * x + (len_y-1) * len_x ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+1] + data[j-1] + data[j-len_x] + data[j+plane])
      + 0.055556 * (data[j-len_x+1] + data[j-len_x-1] + data[j+plane+1] + data[j+plane-1] + data[j+plane-len_x])
      + 0.027778 * (data[j-len_x+1+plane] + data[j-len_x-1+plane]) ;
    
    i = x + (r.len_z-1) * r.plane ;
    j = 2 * x + (len_z-1) * plane ;
    
    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+1] + data[j-1] + data[j+len_x] + data[j-plane])
      + 0.055556 * (data[j+len_x+1] + data[j+len_x-1] + data[j-plane+1] + data[j-plane-1] + data[j-plane+len_x])
      + 0.027778 * (data[j+len_x+1-plane] + data[j+len_x-1-plane]) ;

    i = x + (r.len_z-1) * r.plane + (r.len_y-1) * r.len_x ;
    j = 2 * x + (len_z-1) * plane + (len_y-1) * len_x ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+1] + data[j-1] + data[j-len_x] + data[j-plane])
      + 0.055556 * (data[j-len_x+1] + data[j-len_x-1] + data[j-plane+1] + data[j-plane-1] + data[j-plane-len_x])
      + 0.027778 * (data[j-len_x+1-plane] + data[j-len_x-1-plane]) ;
  }

  /* Y-Axis */

  for( y = 1 ; y < r.len_y - 1 ; y++) {
    
    i = y * r.len_x ;
    j = 2 * y * len_x ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+len_x] + data[j-len_x] + data[j+1] + data[j+plane])
      + 0.055556 * (data[j+1+len_x] + data[j+1-len_x] + data[j+plane+len_x] + data[j+plane-len_x] + data[j+plane+1])
      + 0.027778 * (data[j+1+len_x+plane] + data[j+1-len_x+plane]) ;

    i = y * r.len_x + r.len_x - 1 ;
    j = 2 * y * len_x + len_x - 1 ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+len_x] + data[j-len_x] + data[j-1] + data[j+plane])
      + 0.055556 * (data[j-1+len_x] + data[j-1-len_x] + data[j+plane+len_x] + data[j+plane-len_x] + data[j+plane-1])
      + 0.027778 * (data[j-1+len_x+plane] + data[j-1-len_x+plane]) ;

    i = y * r.len_x + (r.len_y-1) * r.plane ;
    j = 2 * y * len_x + (len_y-1) * plane ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+len_x] + data[j-len_x] + data[j+1] + data[j-plane])
      + 0.055556 * (data[j+1+len_x] + data[j+1-len_x] + data[j-plane+len_x] + data[j-plane-len_x] + data[j-plane+1])
      + 0.027778 * (data[j+1+len_x-plane] + data[j+1-len_x-plane]) ;
    
    i = y * r.len_x + (r.len_y-1) * r.plane + r.len_x - 1 ;
    j = 2 * y * len_x + (len_y-1) * plane + len_x - 1 ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j-len_x] + data[j-len_x] + data[j+1] + data[j-plane])
      + 0.055556 * (data[j+1-len_x] + data[j+1-len_x] + data[j-plane-len_x] + data[j-plane-len_x] + data[j-plane+1])
      + 0.027778 * (data[j+1-len_x-plane] + data[j+1-len_x-plane]) ;
  }

  /* Z-Axis */
     
  for( z = 1 ; z < r.len_z - 1 ; z++) {
    i = z * r.plane ;
    j = 2 * z * plane ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+plane] + data[j-plane] + data[j+len_x] + data[j+1])
      + 0.055556 * (data[j+len_x+plane] + data[j+len_x-plane] + data[j+1+plane] + data[j+1-plane] + data[j+1+len_x])
      + 0.027778 * (data[j+len_x+plane+1] + data[j+len_x-plane+1]) ;

    i = z * r.plane + r.len_x - 1 ;
    j = 2 * z * plane +len_x - 1 ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+plane] + data[j-plane] + data[j+len_x] + data[j-1])
      + 0.055556 * (data[j+len_x+plane] + data[j+len_x-plane] + data[j-1+plane] + data[j-1-plane] + data[j-1+len_x])
      + 0.027778 * (data[j+len_x+plane-1] + data[j+len_x-plane-1]) ;

    i = z * r.plane + (r.len_y-1) * r.len_x ;
    j = 2 * z * plane + (len_y-1) * len_x ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+plane] + data[j-plane] + data[j-len_x] + data[j+1])
      + 0.055556 * (data[j-len_x+plane] + data[j-len_x-plane] + data[j+1+plane] + data[j+1-plane] + data[j+1-len_x])
      + 0.027778 * (data[j-len_x+plane+1] + data[j-len_x-plane+1]) ;

    i = z * r.plane + (r.len_y-1) * r.len_x + r.len_x - 1 ;
    j = 2 * z * plane + (len_y-1) * len_x + len_x - 1 ;

    r.data[i] = 0.222222 * data[j]
      + 0.111111 * (data[j+plane] + data[j-plane] + data[j-len_x] + data[j-1])
      + 0.055556 * (data[j-len_x+plane] + data[j-len_x-plane] + data[j-1+plane] + data[j-1-plane] + data[j-1-len_x])
      + 0.027778 * (data[j-len_x+plane-1] + data[j-len_x-plane-1]) ;
  }

  i = 0 ;
  j = 0 ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j+1] + data[j+len_x] + data[j+plane])
    + 0.074 * (data[j+1+len_x] + data[j+1+plane] + data[j+len_x+plane])
    + 0.037 * data[j+1+len_x+plane] ;

  i = r.len_x - 1 ;
  j = len_x - 1 ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j-1] + data[j+len_x] + data[j+plane])
    + 0.074 * (data[j-1+len_x] + data[j-1+plane] + data[j+len_x+plane])
    + 0.037 * data[j-1+len_x+plane] ;

  i = (r.len_y - 1) * r.len_x ;
  j = (len_y - 1) * len_x  ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j+1] + data[j-len_x] + data[j+plane])
    + 0.074 * (data[j+1-len_x] + data[j+1+plane] + data[j-len_x+plane])
    + 0.037 * data[j+1-len_x+plane] ;
  
  i = (r.len_y - 1) * r.len_x + r.len_x - 1 ; 
  j = (len_y - 1) * len_x + len_x - 1 ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j-1] + data[j-len_x] + data[j+plane])
    + 0.074 * (data[j-1-len_x] + data[j-1+plane] + data[j-len_x+plane])
    + 0.037 * data[j-1-len_x+plane] ;

  i = (r.len_z - 1) * r.plane ;
  j = (len_z - 1 ) * plane ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j+1] + data[j+len_x] + data[j-plane])
    + 0.074 * (data[j+1+len_x] + data[j+1-plane] + data[j+len_x-plane])
    + 0.037 * data[j+1+len_x-plane] ;

  i = r.len_x - 1 + (r.len_z - 1) * r.plane ;
  j = len_x - 1 + (len_z - 1 ) * plane ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j-1] + data[j+len_x] + data[j-plane])
    + 0.074 * (data[j-1+len_x] + data[j-1-plane] + data[j+len_x-plane])
    + 0.037 * data[j-1+len_x-plane] ;

  i = (r.len_y - 1) * r.len_x + (r.len_z - 1) * r.plane ;
  j = (len_y - 1) * len_x + (len_z - 1 ) * plane  ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j+1] + data[j-len_x] + data[j-plane])
    + 0.074 * (data[j+1-len_x] + data[j+1-plane] + data[j-len_x-plane])
    + 0.037 * data[j+1-len_x-plane] ;
  
  i = (r.len_y - 1) * r.len_x + r.len_x - 1 + (r.len_z - 1) * r.plane ; 
  j = (len_y - 1) * len_x + len_x - 1 + (len_z - 1 ) * plane ;
  
  r.data[i] = 0.296 * data[j]
    + 0.148 * (data[j-1] + data[j-len_x] + data[j-plane])
    + 0.074 * (data[j-1-len_x] + data[j-1-plane] + data[j-len_x-plane])
    + 0.037 * data[j-1-len_x-plane] ;


}

void Tensor3d::restrict_half(Tensor3d& r) {
 int x,y, z;
  r.set((len_x - 1) / 2 + 1, (len_y - 1) / 2 + 1, (len_z - 1) / 2 + 1) ;

  int i, j ;
  
  for(z = 1 ; z < r.len_z - 1 ; z++) {
    for( y = 1 ; y < r.len_y - 1 ; y++) {
      for( x = 1 ; x < r.len_x - 1 ; x++) {
	i = x + y * r.len_x + z * r.plane ;
	j = 2 * x + 2 * y * len_x + 2 * z * plane ;

	r.data[i] = 0.25 * data[j]
	  + 0.125 * (data[j-1] + data[j+1] + data[j-len_x] + data[j+len_x] + data[j-plane] + data[j+plane]) ;
      }
    }
  }

  /* Z-Plane (XY) */
  
  for(y = 1 ; y < r.len_y - 1 ; y++) {
    for( x = 1 ; x < r.len_x - 1 ; x++) {
      i = x + y * r.len_x ;
      j = 2 * x + 2 * y * len_x ;

      r.data[i] = 0.2857 * data[j]
	+ 0.14285 * (data[j-1] + data[j+1] + data[j-len_x] + data[j+len_x] + data[j+plane]) ;

      i += (r.len_z - 1) * r.plane ;
      j += (len_z - 1) * plane ;
      
      r.data[i] = 0.2857 * data[j]
	+ 0.14285 * (data[j-1] + data[j+1] + data[j-len_x] + data[j+len_x] + data[j-plane]) ;
    }
  }

  /* Y-Plane (ZX) */
  
  for( z = 1 ; z < r.len_z - 1 ; z++) {
    for( x = 1 ; x < r.len_x - 1 ; x++) {
      i = x + z * r.plane ;
      j = 2 * x + 2 * z * plane ;

      r.data[i] = 0.2857 * data[j]
	+ 0.14285 * (data[j-1] + data[j+1] + data[j-plane] + data[j+plane] + data[j+len_x]) ;

      i += (r.len_y-1) * r.len_x ;
      j += (len_y-1) * len_x ;
      
      r.data[i] = 0.2857 * data[j]
	+ 0.14285 * (data[j-1] + data[j+1] + data[j-plane] + data[j+plane] + data[j-len_x]) ;
    }
  }

  /* X-Plane (ZY) */
  
  for( z = 1 ; z < r.len_z - 1 ; z++) {
    for( y = 1 ; y < r.len_y - 1 ; y++) {
      i = y * r.len_x + z * r.plane ;
      j = 2 * y * len_x + 2 * z * plane ;

      r.data[i] = 0.2857 * data[j]
	+ 0.14285 * (data[j-len_x] + data[j+len_x] + data[j-plane] + data[j+plane] + data[j+1]) ;

      i += r.len_x - 1 ;
      j += len_x - 1 ;
      
      r.data[i] = 0.2857 * data[j]
	+ 0.14285 * (data[j-len_x] + data[j+len_x] + data[j-plane] + data[j+plane] + data[j-1]) ;
    }
  }

  /* X-Axis */
  
  for( x = 1 ; x < r.len_x - 1 ; x++) {
    i = x ;
    j = 2 * x ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+1] + data[j-1] + data[j+len_x] + data[j+plane]) ;

    i = x + (r.len_y-1) * r.len_x ;
    j = 2 * x + (len_y-1) * len_x ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+1] + data[j-1] + data[j-len_x] + data[j+plane]) ;
    
    i = x + (r.len_z-1) * r.plane ;
    j = 2 * x + (len_z-1) * plane ;
    
    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+1] + data[j-1] + data[j+len_x] + data[j-plane]) ;

    i = x + (r.len_z-1) * r.plane + (r.len_y-1) * r.len_x ;
    j = 2 * x + (len_z-1) * plane + (len_y-1) * len_x ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+1] + data[j-1] + data[j-len_x] + data[j-plane]) ;

  }

  /* Y-Axis */

  for( y = 1 ; y < r.len_y - 1 ; y++) {
    
    i = y * r.len_x ;
    j = 2 * y * len_x ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+len_x] + data[j-len_x] + data[j+1] + data[j+plane]) ;

    i = y * r.len_x + r.len_x - 1 ;
    j = 2 * y * len_x + len_x - 1 ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+len_x] + data[j-len_x] + data[j-1] + data[j+plane]) ;

    i = y * r.len_x + (r.len_y-1) * r.plane ;
    j = 2 * y * len_x + (len_y-1) * plane ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+len_x] + data[j-len_x] + data[j+1] + data[j-plane]) ;
    
    i = y * r.len_x + (r.len_y-1) * r.plane + r.len_x - 1 ;
    j = 2 * y * len_x + (len_y-1) * plane + len_x - 1 ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j-len_x] + data[j-len_x] + data[j+1] + data[j-plane]) ;
  }
  
  /* Z-Axis */
     
  for( z = 1 ; z < r.len_z - 1 ; z++) {
    i = z * r.plane ;
    j = 2 * z * plane ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+plane] + data[j-plane] + data[j+len_x] + data[j+1]) ;
  
    i = z * r.plane + r.len_x - 1 ;
    j = 2 * z * plane +len_x - 1 ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+plane] + data[j-plane] + data[j+len_x] + data[j-1]) ;

    i = z * r.plane + (r.len_y-1) * r.len_x ;
    j = 2 * z * plane + (len_y-1) * len_x ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+plane] + data[j-plane] + data[j-len_x] + data[j+1]) ;

    i = z * r.plane + (r.len_y-1) * r.len_x + r.len_x - 1 ;
    j = 2 * z * plane + (len_y-1) * len_x + len_x - 1 ;

    r.data[i] = 0.3333 * data[j]
      + 0.16666 * (data[j+plane] + data[j-plane] + data[j-len_x] + data[j-1]) ;
  }

  i = 0 ;
  j = 0 ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j+1] + data[j+len_x] + data[j+plane]) ;

  i = r.len_x - 1 ;
  j = len_x - 1 ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j-1] + data[j+len_x] + data[j+plane]) ;

  i = (r.len_y - 1) * r.len_x ;
  j = (len_y - 1) * len_x  ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j+1] + data[j-len_x] + data[j+plane]) ;
  
  i = (r.len_y - 1) * r.len_x + r.len_x - 1 ; 
  j = (len_y - 1) * len_x + len_x - 1 ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j-1] + data[j-len_x] + data[j+plane]) ;

  i = (r.len_z - 1) * r.plane ;
  j = (len_z - 1 ) * plane ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j+1] + data[j+len_x] + data[j-plane]) ;

  i = r.len_x - 1 + (r.len_z - 1) * r.plane ;
  j = len_x - 1 + (len_z - 1 ) * plane ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j-1] + data[j+len_x] + data[j-plane]) ;

  i = (r.len_y - 1) * r.len_x + (r.len_z - 1) * r.plane ;
  j = (len_y - 1) * len_x + (len_z - 1 ) * plane  ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j+1] + data[j-len_x] + data[j-plane]) ;
  
  i = (r.len_y - 1) * r.len_x + r.len_x - 1 + (r.len_z - 1) * r.plane ; 
  j = (len_y - 1) * len_x + len_x - 1 + (len_z - 1 ) * plane ;
  
  r.data[i] = 0.4 * data[j]
    + 0.2 * (data[j-1] + data[j-len_x] + data[j-plane]) ;
}

void Tensor3d::restrict_none(Tensor3d& r) {
  r.set((len_x - 1) / 2 + 1, (len_y - 1) / 2 + 1, (len_z - 1) / 2 + 1) ;

  int i, j,x,y,z ;
   
  for( z = 0 ; z < r.len_z ; z++) {
    for(y = 0 ; y < r.len_y ; y++) {
      for( x = 0 ; x < r.len_x ; x++) {
	i = x + y * r.len_x + z * r.plane ;
	j = 2 * x + 2 * y * len_x + 2 * z * plane ;

	r.data[i] = data[j] ;
      }
    }
  }
}

void Tensor3d::prolong_full(Tensor3d& p) {
  p.set((len_x - 1) * 2 + 1, (len_y - 1) * 2 + 1, (len_z - 1) * 2 + 1) ;

  int i, j ,x,y,z;
  
  p.data[0] = data[0] ;
  
  for(x = 1 ; x < len_x ; x++) {
    i = x ;
    j = 2 * x ;

    p.data[j] = data[i] ;
    p.data[j-1] = 0.5 * (data[i] + data[i-1]) ;
  }
  
  for( y = 1 ; y < len_y ; y++) {
    i = y * len_x ;
    j = 2 * y * p.len_x  ;

    p.data[j] = data[i] ;
    p.data[j-p.len_x] = 0.5 * (data[i] + data[i-len_x]) ;
  }

  
  for( z = 1 ; z < len_z ; z++) {
    i = z * plane ;
    j = 2 * z * p.plane  ;

    p.data[j] = data[i] ;
    p.data[j-p.plane] = 0.5 * (data[i] + data[i-plane]) ;
  }

  for( y = 1 ; y < len_y ; y++) {
    for( x = 1 ; x < len_x ; x++) {
      i = x + y * len_x ;
      j = 2 * x + 2 * y * p.len_x  ;
      
      p.data[j] = data[i] ;
      
      p.data[j-1] = 0.5 * (data[i] + data[i-1]) ;
      p.data[j-p.len_x] = 0.5 * (data[i] + data[i-len_x]) ;

      p.data[j-1-p.len_x] = 0.25 * (data[i] + data[i-len_x] + data[i-1] + data[i-1-len_x]) ;
    }
  }
  
  for(z = 1 ; z < len_z ; z++) {
    for( x = 1 ; x < len_x ; x++) {
      i = x + z * plane ;
      j = 2 * x + 2 * z * p.plane  ;
      
      p.data[j] = data[i] ;
      
      p.data[j-1] = 0.5 * (data[i] + data[i-1]) ;
      p.data[j-p.plane] = 0.5 * (data[i] + data[i-plane]) ;

      p.data[j-1-plane] = 0.25 * (data[i] + data[i-plane] + data[i-1] + data[i-1-plane]) ;
    }
  }
  
  for( z = 1 ; z < len_z ; z++) {
    for( y = 1 ; y < len_y ; y++) {
      i = y * len_x + z * plane ;
      j = 2 * y * p.len_x + 2 * z * p.plane  ;
      
      p.data[j] = data[i] ;
      
      p.data[j-p.len_x] = 0.5 * (data[i] + data[i-len_x]) ;
      p.data[j-p.plane] = 0.5 * (data[i] + data[i-plane]) ;

      p.data[j-p.len_x-plane] = 0.25 * (data[i] + data[i-plane] + data[i-len_x] + data[i-len_x-plane]) ;
    }
  }
    
  for( z = 1 ; z < len_z ; z++) {
    for( y = 1 ; y < len_y ; y++) {
      for( x = 1 ; x < len_x ; x++) {
	i = x + y * len_x + z * plane ;
	j = 2 * x + 2 * y * p.len_x + 2 * z * p.plane ;

	p.data[j] = data[i] ;

	p.data[j-1] = 0.5 * (data[i-1] + data[i]) ;
	p.data[j-p.len_x] = 0.5 * (data[i-len_x] + data[i]) ;
	p.data[j-p.plane] = 0.5 * (data[i-plane] + data[i]) ;

	p.data[j-1-p.len_x] = 0.25 * (data[i-1] + data[i] + data[i-len_x] + data[i-len_x-1]) ;
	p.data[j-1-p.plane] = 0.25 * (data[i-1] + data[i] + data[i-plane] + data[i-plane-1]) ;
	p.data[j-p.len_x-p.plane] = 0.25 * (data[i-len_x] + data[i] + data[i-plane] + data[i-plane-len_x]) ;

	p.data[j-1-p.len_x-p.plane] = 0.125 * (data[i-len_x] + data[i] + data[i-plane] + data[i-plane-len_x] +
					       data[i-len_x-1] + data[i-1-plane] + data[i-1] + data[i-plane-len_x-1]) ;
      }
    }
  }
}

void Tensor3d::prolong_none(Tensor3d& p) {
  p.set((len_x - 1) * 2 + 1, (len_y - 1) * 2 + 1, (len_z - 1) * 2 + 1) ;

  int i, j,x,y,z ;
  
  p.data[0] = data[0] ;
  
  for(x = 1 ; x < len_x ; x++) {
    i = x ;
    j = 2 * x ;

    p.data[j-1] = p.data[j] = data[i] ;
  }
  
  for(y = 1 ; y < len_y ; y++) {
    i = y * len_x ;
    j = 2 * y * p.len_x  ;

    p.data[j-p.len_x] = p.data[j] = data[i] ;
  }

  
  for( z = 1 ; z < len_z ; z++) {
    i = z * plane ;
    j = 2 * z * p.plane  ;

    p.data[j-p.plane] = p.data[j] = data[i] ;
  }

  for( y = 1 ; y < len_y ; y++) {
    for( x = 1 ; x < len_x ; x++) {
      i = x + y * len_x ;
      j = 2 * x + 2 * y * p.len_x  ;
      
      p.data[j-1-p.len_x] = p.data[j-p.len_x] = p.data[j-1] = p.data[j] = data[i] ;
    }
  }
  
  for( z = 1 ; z < len_z ; z++) {
    for( x = 1 ; x < len_x ; x++) {
      i = x + z * plane ;
      j = 2 * x + 2 * z * p.plane  ;
      
      p.data[j-1-plane] = p.data[j-p.plane] = p.data[j-1] = p.data[j] = data[i] ;
    }
  }
  
  for( z = 1 ; z < len_z ; z++) {
    for( y = 1 ; y < len_y ; y++) {
      i = y * len_x + z * plane ;
      j = 2 * y * p.len_x + 2 * z * p.plane  ;
      
      p.data[j-p.len_x-plane] = p.data[j-p.plane] = p.data[j-p.len_x] = p.data[j] = data[i] ;
    }
  }
    
  for( z = 1 ; z < len_z ; z++) {
    for( y = 1 ; y < len_y ; y++) {
      for( x = 1 ; x < len_x ; x++) {
	i = x + y * len_x + z * plane ;
	j = 2 * x + 2 * y * p.len_x + 2 * z * p.plane ;
	
	p.data[j-1-p.len_x-p.plane] = p.data[j-p.len_x-p.plane] = p.data[j-1-p.plane] = p.data[j-1-p.len_x] = p.data[j-p.plane] = p.data[j-p.len_x] = p.data[j-1] = p.data[j] = data[i] ;
      }
    }
  }
}

  
double Tensor3d::norm_max() {
  double norm = 0 ;

  for(int i = 0 ; i < length ; i++)
    if(fabs(data[i]) > norm)
      norm = fabs(data[i]) ;

  return norm ;
}

void Tensor3d::one_minus() {
  for(int i = 0 ; i < length ; i++)
    data[i] = 1.0 - data[i] ;
}
 
Tensor3d& Tensor3d::operator*=(real val) {
  for(int i = 0 ; i < length ; i++)
     data[i] *= val ;
  return *this ;
}

Tensor3d& Tensor3d::operator+=(Tensor3d& other) {
  if(length != other.length) {
    cerr << "Vectors length mismatch! (op+)" << endl ;
    exit(1) ;
  }
  
  for(int i = 0 ; i < length ; i++)
    data[i] += other.data[i] ;
  return *this ;  
}

Tensor3d& Tensor3d::operator-=(Tensor3d& other) {
  if(length != other.length) {
    cerr << "Vectors length mismatch! (op-)" << endl ;
    exit(1) ;
  }
  
  for(int i = 0 ; i < length ; i++)
    data[i] -= other.data[i] ;
  return *this ;
}

double Tensor3d::norm_abs() {
  double norm = 0 ;
  
  for (int i = 0 ; i < length ; i++)
    norm += fabs(data[i]) ;

  return norm / length ;
}

double Tensor3d::norm_euc() {
  double norm = 0 ;
  
  for (int i = 0 ; i < length ; i++)
    norm += data[i] * data[i] ;

  //return norm / length ;
  return norm ;
}
