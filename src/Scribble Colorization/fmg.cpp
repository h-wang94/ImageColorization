#include "fmg.h"
#include <iostream.h>


void Field_MGN(TensorField *tf, int iters, int ncycle) {

  int depth = tf->get_depth() ;
  double prev_residual=100000;
  
  tf->set_init_guess() ;
 
  int up_to_level = depth - 1 ;
  int nc = 0 ;
  int level,i;
  
  for(nc = 0 ; nc < ncycle ; nc++) {
    
    /* Going down */
    for(level = 0 ; level < up_to_level ; level ++) {
      for( i = 0 ; i < iters ; i++) {
	tf->smooth(level) ; // error, residual 
      }

      tf->calc_next_level_residual(level) ; // calc next level resudual by restricting the current residual & set init guess is 0)
      
    }

    /* End of goind down */
    for( i = 0 ; i < iters * 5 ; i++) {
      tf->smooth(up_to_level) ; // error, residual 
    }
    
    /* Going up*/
    for(level = up_to_level-1 ; level >= 0 ; level --) {

      tf->add_prolonged_prev_level(level) ;
	
      for( i = 0 ; i < iters ; i++) {
	tf->smooth(level) ; // error, residual 
      }
    }
   
    if ((prev_residual<tf->residual())) {     
      break ;
    }
    prev_residual=tf->residual();
  
    
  }
}

    
