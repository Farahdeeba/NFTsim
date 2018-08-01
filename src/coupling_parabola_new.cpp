/** @file coupling_parabola_new.cpp
  @brief A new Coupling class, where nu_ab follows a parabolic variation function.

  The smooth variation of nu_ab is defined by nu_theta+4*(nu_max-nu_theta)*(t-t_1)(t_2-t)/(t_2-t_1)^2.


  @author Farah Deeba,
*/

// Main module header
#include "coupling_parabola_new.h" // CouplingDiffArctan;

// Other neurofield  headers
#include "configf.h"    // Configf;
#include "coupling.h"   // Coupling;
#include "population.h" // Population;
#include "propagator.h" // Propagator;

// C++ standard library headers
#include <algorithm> // std::min_element; std::max_element;
#include <cmath>     // std::atan;
#include <iostream>  // std::cerr; std::endl;
#include <math.h>    // if ((x>=a) && (x<=b))
#include <stdlib.h>   //abs (x)

using std::atan;
using std::cerr;
using std::endl;
using std::cout;
void CouplingParabolaNew::init( Configf& configf ) {
  //read initial ramp conditions from config file
  
  configf.param("nu_thr",nu_thr);
  configf.param("nu_max",nu_max);

  configf.param("t_thr_up",t_thr_up);
  configf.param("t_thr_down",t_thr_down);
  
   
  nu.clear();
  nu.resize(nodes, nu_thr);
  pos = (nu_thr>0)?1:-1;
  
  
  // value of nu at rest
  for( size_type i=0; i<nodes; i++ ) {
     temp_i.push_back(nu_thr); //input the value of nu_min
 }
  
  
  
  for( size_type i=0; i<nodes; i++) {
    P[i] = nu[i]*prepropag.phiinit(configf);
  }

  
  
  //initial time
  time = 0;
}

void CouplingParabolaNew::step() {
  // Return the right value at each time point
  time += deltat;

  for( size_type i=0; i<nodes; i++ ) {
      
      if ((time>=t_thr_up)&& (time<=t_thr_down)){
       // it we want start from nu_0
         //nu[i]= nu_min+4*(nu_max-nu_thr)*(((t_start-t_thr_up)*(t_thr_down-t_end))/pow(t_thr_down-t_thr_up,2)) + 4*(nu_max-nu_thr)*(((time-t_thr_up)*(t_thr_down-time))/pow(t_thr_down-t_thr_up,2));
         nu[i]= nu_thr+ 4*(nu_max-nu_thr)*(((time-t_thr_up)*(t_thr_down-time))/pow(t_thr_down-t_thr_up,2)); 
         
      }
      
      else
      {nu[i]=temp_i[i];} 
  }
  
  
    
  Coupling::step();
}


CouplingParabolaNew::CouplingParabolaNew( size_type nodes, double deltat, size_type index,
                                        const Propagator& prepropag, const Population& postpop,
                                        double tempf )
  : Coupling(nodes,deltat,index,prepropag,postpop) {
  // total simulation time is stored in tempf (defined in solver.cpp line 73)
  time_tot = tempf;
}

CouplingParabolaNew::~CouplingParabolaNew() {
}
