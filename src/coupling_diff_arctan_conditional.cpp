/** @file coupling_diff_arctan.cpp
  @brief A new Coupling class, where nu_ab follows a smooth function.

  The smooth variation of nu_ab is defined by the difference of
  two arctangent functions.


  @author Farah Deeba,
*/

// Main module header
#include "coupling_diff_arctan_conditional.h" // CouplingDiffArctan;

// Other neurofield  headers
#include "configf.h"    // Configf;
#include "coupling.h"   // Coupling;
#include "population.h" // Population;
#include "propagator.h" // Propagator;

// C++ standard library headers
#include <algorithm> // std::min_element; std::max_element;
#include <cmath>     // std::atan;
#include <iostream>  // std::cerr; std::endl;
using std::atan;
using std::cerr;
using std::endl;
using std::cout;
void CouplingDiffArctanConditional::init( Configf& configf ) {
  //read initial ramp conditions from config file
  configf.param("nu_min",nu_min);
  configf.param("nu_max",nu_max);
  configf.param("delta",delta);
  configf.param("t_half_up",t_half_up);
  configf.param("t_half_down",t_half_down);
  configf.param("break_point",break_point);
  
  if(t_half_up > t_half_down) {
    cerr<<"t_half_up must be less than t_half_down" <<endl;
  } else if(t_half_up == t_half_down) {
    cerr<<"t_half_up must be less than t_half_down" <<endl;
    exit (EXIT_FAILURE);
  }

  //size of time vector:
  time_int = time_tot/deltat;

  //pre-compute the profile
  for( int i=0; i<time_tot; i++ ) {
    diff_atan = (atan((i-t_half_up)/delta)-atan((i-t_half_down)/delta));
    deltanu.push_back(diff_atan);
    if(i==break_point){
        nu_break=diff_atan;
        
    }
  }
  //min and max value -  standardise values between 0 and 1
  
  diff_atan_min=*std::min_element(deltanu.begin(), deltanu.end());
  diff_atan_max=*std::max_element(deltanu.begin(), deltanu.end());
 
  
  nu.clear();
  nu.resize(nodes, nu_min);
  pos = (nu_min>0)?1:-1;
  
  
  // value of nu at break_point
  for( size_type i=0; i<nodes; i++ ) {
      nu_break_norm= nu_min + (nu_max-nu_min)*((nu_break-diff_atan_min)/(diff_atan_max-diff_atan_min));
      temp.push_back(nu_break_norm);
 }
  cout << "Value of nu at break_point is : " << nu_break_norm << endl;    
  
  for( size_type i=0; i<nodes; i++) {
    P[i] = nu[i]*prepropag.phiinit(configf);
  }

  //initial time
  time = 0;
}

void CouplingDiffArctanConditional::step(void) {
  // Return the right value at each time point
  time += deltat;

  for( size_type i=0; i<nodes; i++ ) {
      if ( time<break_point){
    nu[i]= nu_min + (nu_max-nu_min)*(((atan((time-t_half_up)/delta)-atan((time-t_half_down)/delta))-diff_atan_min)/(diff_atan_max-diff_atan_min));
      }
      else
      {nu[i]=temp[i];}      
  }
  Coupling::step();
}


CouplingDiffArctanConditional::CouplingDiffArctanConditional( size_type nodes, double deltat, size_type index,
                                        const Propagator& prepropag, const Population& postpop,
                                        double tempf )
  : Coupling(nodes,deltat,index,prepropag,postpop) {
  // total simulation time is stored in tempf (defined in solver.cpp line 73)
  time_tot = tempf;
}

CouplingDiffArctanConditional::~CouplingDiffArctanConditional(void) {
}
