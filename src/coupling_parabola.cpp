/** @file coupling_parabola.cpp
  @brief A new Coupling class, where nu_ab follows a smooth function.

  The smooth variation of nu_ab is defined by the difference of
  two arctangent functions.


  @author Farah Deeba,
*/

// Main module header
#include "coupling_parabola.h" // CouplingDiffArctan;

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
void CouplingParabola::init( Configf& configf ) {
  //read initial ramp conditions from config file
  configf.param("nu_min",nu_min);
  configf.param("nu_max",nu_max);
  configf.param("t_start",t_start);
  configf.param("t_end",t_end);
  //configf.param("break_point",break_point);


  //size of time and nu vector:
  //time_int = t_end-t_start;     //duration of the variation
  
  time_int = time_tot/deltat;
  t_mid = (t_end+t_start)/2;    // mid value of the variation
  nu_start = nu_min;            // Rest value of nu_se 
  nu_mid   = nu_max;             // Maximum value at the max
  nu_end = nu_min;               // Rest value of nu_se at the end
  
  // buid the Y column, i.e. nu_ary as a function of Y=x^2+x+1
  Arry[0][0]=pow(t_start,2);
  Arry[0][1]=t_start;
  Arry[0][2]=1;
  
  Arry[1][0]=pow(t_mid,2);
  Arry[1][1]=t_mid;
  Arry[1][2]=1;
  
  Arry[2][0]=pow(t_end,2);
  Arry[2][1]=t_end;
  Arry[2][2]=1;
  
  // Find the determinent of Arry
  Arry_det= Arry[0][0]*(Arry[1][1]*Arry[2][2]-Arry[1][2]*Arry[2][1])-Arry[0][1]*(Arry[1][0]*Arry[2][2]-Arry[1][2]*Arry[2][0])+Arry[0][2]*(Arry[1][0]*Arry[2][1]-Arry[1][1]*Arry[2][0]);
  Arry_d=abs(Arry_det);  //absolute value of Arry_det
  
  //cout << "Value of det is :"<<Arry_det<< endl;
  
  // Find the adjoint matrix of Arry in transpose order
  
  Arry_adj[0][0]=Arry[1][1]*Arry[2][2]-Arry[1][2]*Arry[2][1];
  Arry_adj[1][0]=Arry[1][2]*Arry[2][0]-Arry[1][0]*Arry[2][2];
  Arry_adj[2][0]=Arry[1][0]*Arry[2][1]-Arry[1][1]*Arry[2][0];
    
  Arry_adj[0][1]=Arry[0][2]*Arry[2][1]-Arry[0][1]*Arry[2][2];
  Arry_adj[1][1]=Arry[0][0]*Arry[2][2]-Arry[0][2]*Arry[2][0];
  Arry_adj[2][1]=Arry[0][1]*Arry[2][0]-Arry[2][1]*Arry[0][0];
          
  Arry_adj[0][2]=Arry[0][1]*Arry[1][2]-Arry[0][2]*Arry[1][1];
  Arry_adj[1][2]=Arry[0][2]*Arry[1][0]-Arry[0][0]*Arry[1][2];
  Arry_adj[2][2]=Arry[0][0]*Arry[1][1]-Arry[0][1]*Arry[1][0];
  
  
  // Find the inverse matrix
  
  Arry_inv[0][0]=Arry_adj[0][0]/Arry_det; 
  Arry_inv[1][0]=Arry_adj[1][0]/Arry_det;
  Arry_inv[2][0]=Arry_adj[2][0]/Arry_det;

  Arry_inv[0][1]=Arry_adj[0][1]/Arry_det;
  Arry_inv[1][1]=Arry_adj[1][1]/Arry_det;
  Arry_inv[2][1]=Arry_adj[2][1]/Arry_det;

  Arry_inv[0][2]=Arry_adj[0][2]/Arry_det;
  Arry_inv[1][2]=Arry_adj[1][2]/Arry_det;
  Arry_inv[2][2]=Arry_adj[2][2]/Arry_det;
 
  
  // Find the time matrix
  time_Arry[0]= Arry_inv[0][0]*nu_start+Arry_inv[0][1]*nu_mid+Arry_inv[0][2]*nu_end;
  time_Arry[1]= Arry_inv[1][0]*nu_start+Arry_inv[1][1]*nu_mid+Arry_inv[1][2]*nu_end;
  time_Arry[2]= Arry_inv[2][0]*nu_start+Arry_inv[2][1]*nu_mid+Arry_inv[2][2]*nu_end;
  

  
  
 //pre-compute the profile
  
  for( int i=t_start; i<t_end; i++ ) {
      para_nu = time_Arry[0]*pow(i,2)+ time_Arry[1]*(i);
    deltanu.push_back(para_nu);
      }
  
  //min and max value -  standardise values between 0 and 1
  
  para_nu_min=*std::min_element(deltanu.begin(), deltanu.end());
  para_nu_max=*std::max_element(deltanu.begin(), deltanu.end());
 
  nu.clear();
  nu.resize(nodes, nu_min);
  pos = (nu_min>0)?1:-1;
  
  
  // value of nu at rest
  for( size_type i=0; i<nodes; i++ ) {
     temp_i.push_back(nu_min);
     temp_f.push_back(nu_end);
 }
  
  
  
  for( size_type i=0; i<nodes; i++) {
    P[i] = nu[i]*prepropag.phiinit(configf);
  }

  
  
  //initial time
  time = 0;
}

void CouplingParabola::step(void) {
  // Return the right value at each time point
  time += deltat;

  for( size_type i=0; i<nodes; i++ ) {
      
      if ((time>=t_start)&& (time<=t_end)){
       nu[i]= nu_min + (nu_max-nu_min)*((((time_Arry[0]*pow(time,2))+(time_Arry[1]*time))-para_nu_min)/(para_nu_max-para_nu_min)); 
      }
      
      else
      {nu[i]=temp_i[i];} 
  }
  
  
    
  Coupling::step();
}


CouplingParabola::CouplingParabola( size_type nodes, double deltat, size_type index,
                                        const Propagator& prepropag, const Population& postpop,
                                        double tempf )
  : Coupling(nodes,deltat,index,prepropag,postpop) {
  // total simulation time is stored in tempf (defined in solver.cpp line 73)
  time_tot = tempf;
}

CouplingParabola::~CouplingParabola(void) {
}
