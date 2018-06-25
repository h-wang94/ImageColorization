#include "fmg.h"
#include "tensor3d.h"
#include <fstream.h>
#include <vector>
#include <algorithm>
#ifndef VSMKGRD_HH
#define VSMKGRD_HH


//int min(int a, int b){	
//	return( (a)<(b) ? (a) : (b) );
//}

//int max(int a, int b){	
//	return( (a)>(b) ? (a) : (b) );
//}

using namespace std;

class MG : public TensorField{
public:

  MG() {_nebl=NULL; }
  MG(int nx, int ny,int nz, int depth) {
    _nebl=NULL;
    set(nx, ny,nz, depth);
  }

  ~MG() {
    if (_nebl){
      delete [] _nebl;
    }

      
  }
  
  
  void set(int nx, int ny,int nz,int depth) {
    
    _deg=1;
    _depth=depth;
     
    
  
    _len_x = nx ;
    _len_y = ny ;
    _len_z = nz;
    int i,j;
    
    if (_nebl){
      delete [] _nebl;
    }
    
    _nebl=new Tensore3dVP[ _depth ];
  
    
    _P.resize(_depth);
    _Div.resize(_depth);
    _I.resize(_depth);
    _G.resize(_depth);
   
    for ( i=0; i<_depth; i++){
      _G[i].resize(_deg+1);
    }
    
    
    for(i = 0 ; i < _depth ; i++) {
      _P[i].set((int)(_len_x * pow(2,-i)), (int)(_len_y * pow(2,-i)), _len_z) ;
      _P[i].clear(0) ;
      _I[i].set((int)(_len_x * pow(2,-i)), (int)(_len_y * pow(2,-i)),_len_z) ;
      _I[i].clear(0) ;
      _Div[i].set((int)(_len_x * pow(2,-i)), (int)(_len_y * pow(2,-i)), _len_z) ;
      _Div[i].clear(0) ;
      
     
      for ( j=0; j<=_deg; j++){
	_G[i][j].set((int)(_len_x * pow(2,-i)), (int)(_len_y * pow(2,-i)),_len_z) ;
	_G[i][j].clear(0);
      }
    }

    _perms_num=1;   
    _perms.resize(_perms_num);
    _perms[0].resize(27);  
    _perms_len=27;
    _perms[0][0].resize(0); 
    _perms[0][0].push_back(0);_perms[0][0].push_back(0);_perms[0][0].push_back(0);
    _perms[0][1].resize(0);
    _perms[0][1].push_back(1);_perms[0][1].push_back(0);_perms[0][1].push_back(0);
    _perms[0][2].resize(0);
    _perms[0][2].push_back(0);_perms[0][2].push_back(1);_perms[0][2].push_back(0);
    _perms[0][3].resize(0);
    _perms[0][3].push_back(1);_perms[0][3].push_back(1);_perms[0][3].push_back(0);
    _perms[0][4].resize(0);
    _perms[0][4].push_back(0);_perms[0][4].push_back(-1);_perms[0][4].push_back(0);
    _perms[0][5].resize(0);
    _perms[0][5].push_back(-1);_perms[0][5].push_back(0);_perms[0][5].push_back(0);
    _perms[0][6].resize(0);
    _perms[0][6].push_back(-1);_perms[0][6].push_back(-1);_perms[0][6].push_back(0);
    _perms[0][7].resize(0);
    _perms[0][7].push_back(1);_perms[0][7].push_back(-1);_perms[0][7].push_back(0);
    _perms[0][8].resize(0);
    _perms[0][8].push_back(-1);_perms[0][8].push_back(1);_perms[0][8].push_back(0);
    _perms[0][9].resize(0); 
    _perms[0][9].push_back(0);_perms[0][9].push_back(0);_perms[0][9].push_back(1);
    _perms[0][10].resize(0);
    _perms[0][10].push_back(1);_perms[0][10].push_back(0);_perms[0][10].push_back(1);
    _perms[0][11].resize(0);
    _perms[0][11].push_back(0);_perms[0][11].push_back(1);_perms[0][11].push_back(1);
    _perms[0][12].resize(0);
    _perms[0][12].push_back(1);_perms[0][12].push_back(1);_perms[0][12].push_back(1);
    _perms[0][13].resize(0);
    _perms[0][13].push_back(0);_perms[0][13].push_back(-1);_perms[0][13].push_back(1);
    _perms[0][14].resize(0);
    _perms[0][14].push_back(-1);_perms[0][14].push_back(0);_perms[0][14].push_back(1);
    _perms[0][15].resize(0);
    _perms[0][15].push_back(-1);_perms[0][15].push_back(-1);_perms[0][15].push_back(1);
    _perms[0][16].resize(0);
    _perms[0][16].push_back(1);_perms[0][16].push_back(-1);_perms[0][16].push_back(1);
    _perms[0][17].resize(0);
    _perms[0][17].push_back(-1);_perms[0][17].push_back(1);_perms[0][17].push_back(1);
    
    _perms[0][18].resize(0); 
    _perms[0][18].push_back(0);_perms[0][18].push_back(0);_perms[0][18].push_back(-1);
    _perms[0][19].resize(0);
    _perms[0][19].push_back(1);_perms[0][19].push_back(0);_perms[0][19].push_back(-1);
    _perms[0][20].resize(0);
    _perms[0][20].push_back(0);_perms[0][20].push_back(1);_perms[0][20].push_back(-1);
    _perms[0][21].resize(0);
    _perms[0][21].push_back(1);_perms[0][21].push_back(1);_perms[0][21].push_back(-1);
    _perms[0][22].resize(0);
    _perms[0][22].push_back(0);_perms[0][22].push_back(-1);_perms[0][22].push_back(-1);
    _perms[0][23].resize(0);
    _perms[0][23].push_back(-1);_perms[0][23].push_back(0);_perms[0][23].push_back(-1);
    _perms[0][24].resize(0);
    _perms[0][24].push_back(-1);_perms[0][24].push_back(-1);_perms[0][24].push_back(-1);
    _perms[0][25].resize(0);
    _perms[0][25].push_back(1);_perms[0][25].push_back(-1);_perms[0][25].push_back(-1);
    _perms[0][26].resize(0);
    _perms[0][26].push_back(-1);_perms[0][26].push_back(1);_perms[0][26].push_back(-1);
    
 
    
   
  }
  
  Tensor3d& Div(int level=0) { return _Div[level]; } ;
  void setDepth(int d){ _depth=d;};
  
  void setI(Tensor3d& II) {
    _I[0] = II ;

    int nx = _len_x ;
    int ny = _len_y;
    int nz = _len_z;
    
    for(int d = 1 ; d < _depth ; d++) {
    
      nx /= 2 ;
      ny /= 2 ;
      for(int z = 0 ; z < nz ; z++)
	for(int y = 0 ; y < ny ; y++)
	  for(int x = 0 ; x < nx ; x++) {
	  
	    if(_I[d-1](2*x,2*y,z) +  _I[d-1](2*x+1,2*y,z) + _I[d-1](2*x,2*y+1,z) + _I[d-1](2*x+1,2*y+1,z) < 4)
	      _I[d](x,y,z) = 0 ;
	    else
	      _I[d](x,y,z) = 1 ;
	  }
    }
  }


  void setG(Tensor3d& GG) {
    _G[0][1] = GG ;
    _G[0][0].clear(1);

    int nx = _len_x ;
    int ny = _len_y ;
    int nz = _len_z;
    
    for(int lvl = 1 ; lvl < _depth ; lvl++) {
      _G[lvl][0].clear(1);
       for (int d=1; d<2; d++){
	 restrict(_G[lvl-1][d],_G[lvl][d]);
       }
    }
    
  }

  void setFlow(Tensor3d& dx,Tensor3d& dy,Tensor3d& idx,Tensor3d& idy);
  
  void smooth(int level)  ;
  void set_init_guess(void)  ;
  
  void calc_next_level_residual(int level) ;
  void zero_next_level(int level) ;
  void add_prolonged_prev_level(int level) ;
  void advance_in_time() ;
  double residual() { return  _residual ; }

  Tensor3d& P(int level=0) { return _P[level] ; } ;
  
  void prolong(Tensor3d& P, Tensor3d& pP) ;
  void restrict(Tensor3d& P, Tensor3d& rP) ;
  
private:

  double _residual;
  
  vector< vector< Tensor3d> > _G;
  
  vector<Tensor3d> _Div;
  vector<Tensor3d> _I;
  vector<Tensor3d> _P;
  
  double _dx, _dy , _dz;
  int _perms_num,_perms_len, _deg;
  
  vector< vector < vector< int > > > _perms;
  
  
  int _len_z;
  
  typedef double* index;
  struct Tensore3dVI{
    void set(int len_x,int len_y,int len_z){
      VV.clear();
      VV.resize( len_x* len_y*len_z);
      _len_x=len_x;
      _len_y=len_y;
      _len_z=len_z;
    }
  
    vector<int>& operator()(int x,int y,int z){return VV[x+y*_len_x+z*_len_x*_len_y]; }
    vector<int>& operator[](int x){return VV[x]; }
    
    vector< vector<int> > VV;
    int _len_x;
    int _len_y;
    int _len_z;
  };

  struct Tensore3dI{
    void set(int len_x,int len_y,int len_z){
      VV.clear();
      VV.resize( len_x* len_y*len_z);
      _len_x=len_x;
      _len_y=len_y;
      _len_z=len_z;
    }
  
    int& operator()(int x,int y,int z){return VV[x+y*_len_x+z*_len_x*_len_y]; }
    vector<int>  VV;
    int _len_x;
    int _len_y;
    int _len_z;
  };
  
  
  struct Tensore3dVD{
    void set(int len_x,int len_y,int len_z){
      VV.clear();
      VV.resize( len_x* len_y*len_z);
      _len_x=len_x;
      _len_y=len_y;
      _len_z=len_z;
    }
  
    vector<float>& operator()(int x,int y,int z){return VV[x+y*_len_x+z*_len_x*_len_y]; }
    vector< vector<float> > VV;
    int _len_x;
    int _len_y;
    int _len_z;
  };

   
  struct Tensore3dVP{
    void set(int len_x,int len_y,int len_z){
      VV.clear();
      VV.resize( len_x* len_y*len_z);
      _len_x=len_x;
      _len_y=len_y;
      _len_z=len_z;
    }
  
    vector<pair<index,float> >& operator()(int x,int y,int z){return VV[x+y*_len_x+z*_len_x*_len_y]; }

    vector<pair<index,float> >& operator[](int x){return VV[x]; }
    
    vector< vector< pair<index,float> > > VV;
    int _len_x;
    int _len_y;
    int _len_z;
  };


  struct Tensore3dMD{
    Tensore3dMD(){
      VV=NULL;
    }
    ~Tensore3dMD(){
      if (VV!=0){
	delete [] VV;
      }
      VV=NULL;
    }
    void set(int len_x,int len_y,int len_z){
      if (VV!=0){
	delete [] VV;
      }
      VV = new  vector< vector< double> >[len_x*len_y*len_z];
      _len_x=len_x;
      _len_y=len_y;
      _len_z=len_z;
    }
  
    vector< vector<double> >& operator()(int x,int y,int z){return VV[x+y*_len_x+z*_len_x*_len_y]; }
    vector< vector<double> >  *VV;
    int _len_x;
    int _len_y;
    int _len_z;
  };
  
  
  Tensore3dVP * _nebl;
} ;

#endif 
