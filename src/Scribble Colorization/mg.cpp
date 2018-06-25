#include "mg.h"
#include <iostream.h>

 

void MG::smooth(int level) {
  Tensor3d& P = _P[level] ;
  Tensor3d& D = _Div[level];
  
  
  int len_x = (int)(_len_x * pow(2,-level)) ;
  int len_y = (int)(_len_y * pow(2,-level)) ;
  
  
 
  double v,av;
  int nanind=-1;
  vector<pair<index,float> >::iterator cnebl_itr;
  
  bool first_nan=1;
  
  for(int z = 0 ; z< _len_z ; z++){
    for(int y = 0 ; y < len_y ; y++){
      for(int x = 0 ; x < len_x ; x++) {
	
	if ((level==0)&(!_I[level](x,y,z))){
	  P(x,y,z)=D(x,y,z);
	  continue;
	}
	if ((level>0)&(!_I[level](x,y,z))){
	  continue;
	}
	vector<pair<index,float> >& cnebl=_nebl[level](x,y,z);
	  
	v=0;
        av=cnebl.begin()->second;

        cnebl_itr=cnebl.begin(); cnebl_itr++;
	for (; cnebl_itr!=cnebl.end();cnebl_itr++){
	  v+=cnebl_itr->second*(*cnebl_itr->first);
	}
	P(x,y,z)=(D(x,y,z)-v)/(av);
	
      }
    } 
  }
}

void MG::set_init_guess(void) {
}
  
void MG::calc_next_level_residual(int level) {
  Tensor3d& D = _Div[level] ;
  int len_x = (int)(_len_x * pow(2,-level)) ;
  int len_y = (int)(_len_y * pow(2,-level)) ;

  Tensor3d R(len_x, len_y,_len_z) ;
  vector<pair<index,float> >::iterator cnebl_itr;
  
  double v ;
  
  if (level==0){
    _residual=0;
  }
 
  R.clear(0);
  for(int z = 0 ; z < _len_z ; z++){
    for(int y = 0 ; y < len_y ; y++){
      for(int x = 0 ; x < len_x ; x++) {
	if(_I[level](x,y,z) < 1) {
	  R(x,y,z) = 0 ;
	  continue;
	}
	v=0;
	vector<pair<index,float> >& cnebl=_nebl[level](x,y,z);
	cnebl_itr=cnebl.begin();
	for (; cnebl_itr!=cnebl.end();cnebl_itr++){
	  v+=cnebl_itr->second*(*cnebl_itr->first);
	}
	R(x,y,z)=D(x,y,z)-v;
      }
    }
  }
  if (level==0){
    _residual =_residual+ R.norm_euc() ;
  }
  if ( (level+1)<_depth){
    restrict(_P[level],_P[level+1]) ;
  }
  
  int nx = len_x / 2 ;
  int ny = len_y / 2 ;
  
  if(level == 0) {
    //  cout <<"level "<<level<<" residual "<< _residual << endl ;
  }
  
}

void MG::zero_next_level(int level) {
}


void MG::add_prolonged_prev_level(int level) {
  Tensor3d pP,tpP; 

  restrict(_P[level],tpP);
  tpP-=_P[level+1];
  prolong(tpP, pP);  
  _P[level] -= pP ;
}


void MG::advance_in_time() {
}


void MG::prolong(Tensor3d& P, Tensor3d& pP) {
  {
    int nx = P.len_x ;
    int ny = P.len_y ;
    int nz = P.len_z ;
    
    int px = nx ;
    int py = ny ;
    int pz = nz ;
    
    int ppx = nx * 2 ;
    int ppy = ny * 2 ;
    int ppz = nz ;
    
    pP.set(ppx, ppy ,ppz) ;

    /* prolonging P */
    for (int z=0; z<nz; z++){
      
      int xxyy = 1 + ppx +z*ppx*ppy;
      int xxpyy = xxyy + 1 ;
      int xxyyp = xxyy + ppx ;
      int xxpyyp = xxyy + 1 + ppx ;
    
      int xy = 0 +z*px*py;
      int xpy = 1 +z*px*py;
      int xyp = px +z*px*py;
      int xpyp = 1 + px +z*px*py ;
      int x,y;
      for( y = 1 ; y < py ; y++) {
	for(x = 1 ; x < px ; x++) {
         
	  pP[xxyy] = 0.5 * P[xy] + 0.25 * (P[xpy] + P[xyp]) ;
	  pP[xxpyy] = 0.5 * P[xpy] + 0.25 * (P[xy] + P[xpyp]) ;
	  pP[xxyyp] = 0.5 * P[xyp] + 0.25 * (P[xpyp] + P[xy]) ;
	  pP[xxpyyp] = 0.5 * P[xpyp] + 0.25 * (P[xpy] + P[xyp]) ;
	  xy ++ ; xpy ++ ; xyp ++ ; xpyp ++ ;
	  xxyy += 2 ; xxpyy += 2 ; xxyyp += 2 ; xxpyyp += 2 ; 
	}
	xy ++ ; xpy ++ ; xyp ++ ; xpyp ++ ;
	xxyy += ppx + 2 ; xxpyy += ppx + 2 ; xxyyp += ppx + 2 ; xxpyyp += ppx + 2 ;
      }

      for( x = 1 ; x < px ; x++) {
	pP(2*x-1,0,z) = P(x-1,0,z) * 0.75 +  P(x,0,z) * 0.25 ;
	pP(2*x-1,ppy-1,z) = P(x-1,py-1,z) * 0.75 +  P(x,py-1,z) * 0.25 ;

	pP(2*x,0,z) = P(x-1,0,z) * 0.25 +  P(x,0,z) * 0.75 ;
	pP(2*x,ppy-1,z) = P(x-1,py-1,z) * 0.25 +  P(x,py-1,z) * 0.75 ;

      }
      for( y = 1 ; y < py ; y++) {
	pP(0,2*y-1,z) = P(0,y-1,z) * 0.75 +  P(0,y,z) * 0.25 ;
	pP(ppx-1,2*y-1,z) = P(px-1,y-1,z) * 0.75 +  P(px-1,y,z) * 0.25 ;

	pP(0,2*y,z) = P(0,y-1,z) * 0.25 + P(0,y,z) * 0.75 ;
	pP(ppx-1,2*y,z) = P(px-1,y-1,z) * 0.25 +  P(px-1,y,z) * 0.75 ;
      }
     
      pP(0,0,z) = P(0,0,z) ;
      pP(ppx-1,0,z) = P(px-1,0,z) ;
      pP(0,ppy-1,z) = P(0,py-1,z) ;
      pP(ppx-1,ppy-1,z) = P(px-1,py-1,z) ;
    }
  }

}



void MG::restrict(Tensor3d& P, Tensor3d& rP) {
 
  int nx = P.len_x ;
  int ny = P.len_y ;
  int nz = P.len_z;
  
  int px = nx / 2 ;
  int py = ny / 2 ;
  int pz = nz;
  rP.set(px, py, pz) ;
  
  int ppx = nx ;

  for (int z=0; z<nz; z++)
  {
    int xy = z*px*py;
    int xxyy = z*nx*ny ;
    int xxpyy = 1 + z*nx*ny;
    int xxyyp = 1 * ppx + z*nx*ny ;
    int xxpyyp = 1 + 1 * ppx + z*nx*ny ;
    
    for(int y = 0 ; y < py ; y++) {
      for(int x = 0 ; x < px ; x++) {
	rP[xy] = (P[xxyy] + P[xxpyy] + P[xxyyp] + P[xxpyyp]) * 0.25 ;
	xy ++ ;
	xxyy += 2 ; xxpyy += 2 ; xxyyp += 2 ; xxpyyp += 2 ;
      }   
      xxyy += ppx  ; xxpyy += ppx ; xxyyp += ppx ; xxpyyp += ppx ;
    }
  }
}



void MG::setFlow(Tensor3d& dx,Tensor3d& dy,Tensor3d& idx,Tensor3d& idy){
  int x0, y0, z0, cx0, cy0, x,y,z,level;
  int len_x=2*_len_x;
  int len_y=2*_len_y;
  double tw;

  Tensor3d pdx,pdy,pidx,pidy,td;
  pdx=dx;
  pdy=dy;
  pidx=idx;
  pidy=idy;
  vector< vector< Tensore3dVI> > neb;
  neb.resize(_depth);
  
  for ( level=0; level<_depth; level++){
    len_x=len_x/2;
    len_y=len_y/2;
    neb[level].resize(_perms_num);
    
    for (int p=0; p<_perms_num;p++){
      neb[level][p].set(len_x,len_y,_len_z);
      
      for ( z=0; z<_len_z; z++){ 
	for ( y=0; y<len_y; y++){
	  for ( x=0; x<len_x; x++){
	    neb[level][p](x,y,z).reserve( _perms_len);
	  }
	}
      }
    
      for ( z=0; z<_len_z; z++){ 
	for ( y=0; y<len_y; y++){
	  for ( x=0; x<len_x; x++){
       
	    for (int q=0; q<_perms_len; q++){
	      x0=x+_perms[p][q][0];
	      y0=y+_perms[p][q][1];
	      z0=z+_perms[p][q][2];
	      
	      if ((x0<0)|(y0<0)|(z0<0)|(x0>=len_x)|(y0>=len_y)|(z0>=_len_z)){
		continue;
	      }
	      if (z0<z){
		 cx0=int(x0-pdx(x0,y0,z0)+0.5);
		 cy0=int(y0-pdy(x0,y0,z0)+0.5);
	      }
	      if (z0>z){
		 cx0=int(x0-pidx(x0,y0,z)+0.5);
		 cy0=int(y0-pidy(x0,y0,z)+0.5);
	      }
	      if (z0==z){
		cx0=x0; cy0=y0;
	      }
	      if ((cx0<0)|(cy0<0)|(z0<0)|(cx0>=len_x)|(cy0>=len_y)|(z0>=_len_z)){
		continue;
	      } 
	      
	      neb[level][p](x,y,z).push_back(cx0+cy0*len_x+z0*len_x*len_y);
	   
	    }  
	  }
	}
      }
    }


    restrict(pdx,td);
    pdx=td;
    restrict(pdy,td);
    pdy=td;
    restrict(pidx,td);
    pidx=td;
    restrict(pidy,td);
    pidy=td;
  }
  



  vector<int> tlist, tneb;
  vector<int>::iterator tneb_itr,tlist_itr,end_itr;
  
  vector<pair<index,float> > cnebl;
  vector<pair<index,float> > *pcnebl;
  pair<index,float> tpair;
  int prev_ind;
  len_x=2*_len_x; len_y=2*_len_y;
  
  for ( level=0; level<_depth; level++){
    len_x=len_x/2;
    len_y=len_y/2;
    _nebl[level].set(len_x,len_y,_len_z);
 
    for (int z=0; z<_len_z; z++){ 
      for (int y=0; y<len_y; y++){
	for (int x=0; x<len_x; x++){
	  
          tlist.resize(0);
	  for (int p=0; p<_perms_num;p++){
	    tneb=neb[level][p](x,y,z);
	    for (tneb_itr=tneb.begin(); tneb_itr!=tneb.end();tneb_itr++){
	      tlist.push_back(*tneb_itr); 
	    }
	    
	  }
	  
          sort(tlist.begin(),tlist.end());
	
 
	  if ( _len_z==1){
	    _nebl[level](x,y,z).reserve(9);
	  }else{
	    _nebl[level](x,y,z).reserve(27);
	  }
	  
 
	  pcnebl=&_nebl[level](x,y,z);
	  
	  prev_ind=-1;
	  for(tlist_itr=tlist.begin(); tlist_itr!=tlist.end(); tlist_itr++){
	    if (prev_ind< (*tlist_itr)){
              tpair.first=&(_P[level][*tlist_itr]);
	      tpair.second=0;		   
	      pcnebl->push_back(tpair);
	      prev_ind= (*tlist_itr);   
	    }
	  }
	 
	  
	}

      }
    }
    
  }

  
  
   
  double var, mean, g_val, sum_w, t_val;
  vector<double> g_vec,w_vec,b_g_vec;
  vector<index> p_vec;
  vector<pair<index, pair<double,double> > > t_vec;
  vector<pair<index, pair<double,double> > >::iterator t_vec_itr;
  pair<index, pair<double,double> > t_data;
  vector<pair<index,float> >::iterator cnebl_itr;
  len_x=2*_len_x; len_y=2*_len_y;
  bool first_nan=1;
  
  
  for ( level=0; level<_depth; level++){
    len_x=len_x/2;
    len_y=len_y/2;
    
    for (int z=0; z<_len_z; z++){ 
      for (int y=0; y<len_y; y++){
	for (int x=0; x<len_x; x++){
	 
	  for (int p=0; p<_perms_num;p++){
	   
	    tneb=neb[level][p](x,y,z);
	   
	    sort(tneb.begin(),tneb.end());
	    end_itr=unique(tneb.begin(),tneb.end());

	   
	    
	    double csig;
	    int i;
            var=0; t_val=_G[level][1](x,y,z);
	    
	   
	    var=0; mean=0; sum_w=0;
	    g_vec.resize(0);
	    tneb_itr=tneb.begin();
	    for (i=0; tneb_itr!=end_itr;tneb_itr++, i++){
	      g_val=_G[level][1][*tneb_itr];
	      g_vec.push_back((g_val-t_val)*(g_val-t_val)); 
	      sum_w+=1;
	      var+=g_val*g_val;
	      mean+=g_val;
	    }

	    if (g_vec.size()<2){
	      continue;
	    }
	    
	    var/=sum_w;
	    mean/=sum_w;
	    var-=mean*mean;
  
	    sort(g_vec.begin(),g_vec.end());
	    csig=var*0.6;
	    
	    if (csig<(-g_vec[1]/log(0.01))){
	      csig=-g_vec[1]/log(0.01);
	    }
	    if (csig<0.000002){
	      csig=0.000002;
	    }
            sum_w=0;
	    for(i=0; i<g_vec.size(); i++){
	      sum_w+=exp(-g_vec[i]/csig);
	    }
	 
	    t_vec.resize(0);
	    tneb_itr=tneb.begin();   
	    
	    for (i=0; tneb_itr!=end_itr;tneb_itr++, i++){	
	      g_val=_G[level][1][*tneb_itr];
	       
	      tw=exp(-(g_val-t_val)*(g_val-t_val)/csig)/sum_w;
	      
	      if(  (&_P[level][tneb[i]])==(&_P[level](x,y,z))){
		tw-=1;
	      }
	      t_data.first=&_P[level][tneb[i]];
	      t_data.second.first=g_val;
	      t_data.second.second=tw;
	      t_vec.push_back(t_data); 
	
	    }
  
            sort(t_vec.begin(),t_vec.end());
	    
	    pcnebl=&_nebl[level](x,y,z);
	    cnebl_itr=pcnebl->begin();
	  
	    for (t_vec_itr=t_vec.begin();t_vec_itr!=t_vec.end();t_vec_itr++){
	      for(;cnebl_itr->first<t_vec_itr->first;++cnebl_itr);	
	
	      cnebl_itr->second+=(t_vec_itr->second.second);
	      cnebl_itr++; 
	    }
	     
	  }
	   
	  
	}
      }
    }
  }
  
  index cp;
  vector<pair<index,float> >::iterator pnebl_itr;
  len_x=2*_len_x; len_y=2*_len_y;
  for (level=0; level<_depth; level++){
    len_x=len_x/2;
    len_y=len_y/2;
    
    for (int z=0; z<_len_z; z++){ 
      for (int y=0; y<len_y; y++){
	for (int x=0; x<len_x; x++){
	  pcnebl=&_nebl[level](x,y,z);
	  cp=&_P[level](x,y,z);
	 
	  for (cnebl_itr=pcnebl->begin(); cnebl_itr->first!=cp; cnebl_itr++);
	  tpair=*cnebl_itr;
	  pnebl_itr=cnebl_itr;
	   
	  for (; cnebl_itr!=pcnebl->begin();--cnebl_itr){
	    --pnebl_itr;
	    *cnebl_itr=*pnebl_itr;
	  }
	   
	  *(pcnebl->begin())=tpair;
	}

      }
    }
  }  
}


