/** @file coupling_parabola_continuous.cpp
  @brief A new Coupling class, where nu_ab follows a parabolic function and repeat 4 times.

  The smooth variation of nu_ab is govern by the equation for the parabola: x^2+x+1.


  @author Farah Deeba,
*///26/06/2018

// Main module header
#include "coupling_parabola_continuous.h" // CouplingDiffArctan;

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
void CouplingParabolaContinuous::init( Configf& configf ) {
  //read initial ramp conditions from config file
  configf.param("nu_min",nu_min);
  configf.param("nu_max",nu_max);
  configf.param("time_start",time_start);
  configf.param("time_end",time_change);
  configf.param("amp_inc",amp_inc);          // increament per cycle
  configf.param("time_laps",time_laps);      //time gap between the cycle
  


  //size of time and nu vector:
      
  
  time_int = time_tot/deltat;       //duration of the variation
  time_mid = (time_change+time_start)/2;    // mid value of the variation

  nu_start = nu_min;            // Rest value of nu_se  
  nu_mid   = nu_max;             // Maximum value at the max
  nu_end = nu_min;               // Rest value of nu_se at the end
  


  // Compute the subsequent amplitudes
  nu_mid_array[1]= nu_max+amp_inc;                         // Maximum value of 2nd cycle
  nu_mid_array[2]= nu_max+amp_inc+amp_inc;                 // Maximum value of 3rd cycle 
  nu_mid_array[3]= nu_max+amp_inc+amp_inc+amp_inc;         // Maximum value of 4th cycle
  nu_mid_array[4]= nu_max+amp_inc+amp_inc+amp_inc+amp_inc; // Maximum value of 5th cycle

  // Initialize variables to compute time periods of the subsequent cycles
  
  int cycle=5;                   //Number of cycle. 
  int time_st_1=time_start;
  int time_ch_1=time_change;
  
  for (int j=0; j<cycle; j++){ 
          time_start_array[j]=time_st_1;
          time_change_array[j]=time_ch_1;        
          time_st_1=time_st_1+time_laps;
          time_ch_1=time_ch_1+time_laps;           
      
  }

  time_mid_2 = (time_start_array[1]+time_change_array[1])/2;    // mid time of 2nd cycle
  time_mid_3 = (time_start_array[2]+time_change_array[2])/2;    // mid time of 3rd cycle
  time_mid_4 = (time_start_array[3]+time_change_array[3])/2;    // mid time of 4th cycle
  time_mid_5 = (time_start_array[4]+time_change_array[4])/2;    // mid time of 5th cycle
  
  
  
// buid the Y column, i.e. nu_ary as a function of Y=x^2+x+1 for 1st cycle
  Arry[0][0]=pow(time_start,2);
  Arry[0][1]=time_start;
  Arry[0][2]=1;
  
  Arry[1][0]=pow(time_mid,2);
  Arry[1][1]=time_mid;
  Arry[1][2]=1;
  
  Arry[2][0]=pow(time_change,2);
  Arry[2][1]=time_change;
  Arry[2][2]=1;
  
  // Find the determinent of Arry 
  Arry_det= Arry[0][0]*(Arry[1][1]*Arry[2][2]-Arry[1][2]*Arry[2][1])-Arry[0][1]*(Arry[1][0]*Arry[2][2]-Arry[1][2]*Arry[2][0])+Arry[0][2]*(Arry[1][0]*Arry[2][1]-Arry[1][1]*Arry[2][0]);
  Arry_d=abs(Arry_det);  //absolute value of Arry_det
  
  
  // Find the adjoint matrix of Arry in transpose order for 1st cycle
  
  Arry_adj[0][0]=Arry[1][1]*Arry[2][2]-Arry[1][2]*Arry[2][1];
  Arry_adj[1][0]=Arry[1][2]*Arry[2][0]-Arry[1][0]*Arry[2][2];
  Arry_adj[2][0]=Arry[1][0]*Arry[2][1]-Arry[1][1]*Arry[2][0];
    
  Arry_adj[0][1]=Arry[0][2]*Arry[2][1]-Arry[0][1]*Arry[2][2];
  Arry_adj[1][1]=Arry[0][0]*Arry[2][2]-Arry[0][2]*Arry[2][0];
  Arry_adj[2][1]=Arry[0][1]*Arry[2][0]-Arry[2][1]*Arry[0][0];
          
  Arry_adj[0][2]=Arry[0][1]*Arry[1][2]-Arry[0][2]*Arry[1][1];
  Arry_adj[1][2]=Arry[0][2]*Arry[1][0]-Arry[0][0]*Arry[1][2];
  Arry_adj[2][2]=Arry[0][0]*Arry[1][1]-Arry[0][1]*Arry[1][0];
  
  
  // Find the inverse matrix for 1st cycle
  
  Arry_inv[0][0]=Arry_adj[0][0]/Arry_det; 
  Arry_inv[1][0]=Arry_adj[1][0]/Arry_det;
  Arry_inv[2][0]=Arry_adj[2][0]/Arry_det;

  Arry_inv[0][1]=Arry_adj[0][1]/Arry_det;
  Arry_inv[1][1]=Arry_adj[1][1]/Arry_det;
  Arry_inv[2][1]=Arry_adj[2][1]/Arry_det;

  Arry_inv[0][2]=Arry_adj[0][2]/Arry_det;
  Arry_inv[1][2]=Arry_adj[1][2]/Arry_det;
  Arry_inv[2][2]=Arry_adj[2][2]/Arry_det;
 
  
  // Find the time matrix for 1st cycle
  
  time_Arry[0]= Arry_inv[0][0]*nu_start+Arry_inv[0][1]*nu_max+Arry_inv[0][2]*nu_end;
  time_Arry[1]= Arry_inv[1][0]*nu_start+Arry_inv[1][1]*nu_max+Arry_inv[1][2]*nu_end;
  time_Arry[2]= Arry_inv[2][0]*nu_start+Arry_inv[2][1]*nu_max+Arry_inv[2][2]*nu_end;
  
  
  //Profile for second Cycle
  // buid the Y column, i.e. nu_ary as a function of Y=x^2+x+1 for 2nd cycle
  Arry_2[0][0]=pow(time_start_array[1],2);
  Arry_2[0][1]=time_start_array[1];
  Arry_2[0][2]=1;
  
  Arry_2[1][0]=pow(time_mid_2,2);
  Arry_2[1][1]=time_mid_2;
  Arry_2[1][2]=1;
  
  Arry_2[2][0]=pow(time_change_array[1],2);
  Arry_2[2][1]=time_change_array[1];
  Arry_2[2][2]=1;
  
  // Find the determinent of Arryfor 2nd cycle
  
  Arry_det_2= Arry_2[0][0]*(Arry_2[1][1]*Arry_2[2][2]-Arry_2[1][2]*Arry_2[2][1])-Arry_2[0][1]*(Arry_2[1][0]*Arry_2[2][2]-Arry_2[1][2]*Arry_2[2][0])+Arry_2[0][2]*(Arry_2[1][0]*Arry_2[2][1]-Arry_2[1][1]*Arry_2[2][0]);
  Arry_d_2=abs(Arry_det_2);  //absolute value of Arry_det
 
  
  // Find the adjoint matrix of Arry in transpose order for 2nd cycle
  
  Arry_adj_2[0][0]=Arry_2[1][1]*Arry_2[2][2]-Arry_2[1][2]*Arry_2[2][1];
  Arry_adj_2[1][0]=Arry_2[1][2]*Arry_2[2][0]-Arry_2[1][0]*Arry_2[2][2];
  Arry_adj_2[2][0]=Arry_2[1][0]*Arry_2[2][1]-Arry_2[1][1]*Arry_2[2][0];
    
  Arry_adj_2[0][1]=Arry_2[0][2]*Arry_2[2][1]-Arry_2[0][1]*Arry_2[2][2];
  Arry_adj_2[1][1]=Arry_2[0][0]*Arry_2[2][2]-Arry_2[0][2]*Arry_2[2][0];
  Arry_adj_2[2][1]=Arry_2[0][1]*Arry_2[2][0]-Arry_2[2][1]*Arry_2[0][0];
          
  Arry_adj_2[0][2]=Arry_2[0][1]*Arry_2[1][2]-Arry_2[0][2]*Arry_2[1][1];
  Arry_adj_2[1][2]=Arry_2[0][2]*Arry_2[1][0]-Arry_2[0][0]*Arry_2[1][2];
  Arry_adj_2[2][2]=Arry_2[0][0]*Arry_2[1][1]-Arry_2[0][1]*Arry_2[1][0];
  
  
  // Find the inverse matrix for 2nd cycle
  
  Arry_inv_2[0][0]=Arry_adj_2[0][0]/Arry_det_2; 
  Arry_inv_2[1][0]=Arry_adj_2[1][0]/Arry_det_2;
  Arry_inv_2[2][0]=Arry_adj_2[2][0]/Arry_det_2;

  Arry_inv_2[0][1]=Arry_adj_2[0][1]/Arry_det_2;
  Arry_inv_2[1][1]=Arry_adj_2[1][1]/Arry_det_2;
  Arry_inv_2[2][1]=Arry_adj_2[2][1]/Arry_det_2;

  Arry_inv_2[0][2]=Arry_adj_2[0][2]/Arry_det_2;
  Arry_inv_2[1][2]=Arry_adj_2[1][2]/Arry_det_2;
  Arry_inv_2[2][2]=Arry_adj_2[2][2]/Arry_det_2;
 
  
  // Find the time matrix for 2nd cycle
  
  time_Arry_2[0]= Arry_inv_2[0][0]*nu_start+Arry_inv_2[0][1]*nu_mid_array[1]+Arry_inv_2[0][2]*nu_end;
  time_Arry_2[1]= Arry_inv_2[1][0]*nu_start+Arry_inv_2[1][1]*nu_mid_array[1]+Arry_inv_2[1][2]*nu_end;
  time_Arry_2[2]= Arry_inv_2[2][0]*nu_start+Arry_inv_2[2][1]*nu_mid_array[1]+Arry_inv_2[2][2]*nu_end;
  
  
  //Profile for third slope
  // buid the Y column, i.e. nu_ary as a function of Y=x^2+x+1 for 3rd cycle
  Arry_3[0][0]=pow(time_start_array[2],2);
  Arry_3[0][1]=time_start_array[2];
  Arry_3[0][2]=1;
  
  Arry_3[1][0]=pow(time_mid_3,2);
  Arry_3[1][1]=time_mid_3;
  Arry_3[1][2]=1;
  
  Arry_3[2][0]=pow(time_change_array[2],2);
  Arry_3[2][1]=time_change_array[2];
  Arry_3[2][2]=1;
  
  
  // Find the determinent of Arry  for 3rd cycle
  
  Arry_det_3= Arry_3[0][0]*(Arry_3[1][1]*Arry_3[2][2]-Arry_3[1][2]*Arry_3[2][1])-Arry_3[0][1]*(Arry_3[1][0]*Arry_3[2][2]-Arry_3[1][2]*Arry_3[2][0])+Arry_3[0][2]*(Arry_3[1][0]*Arry_3[2][1]-Arry_3[1][1]*Arry_3[2][0]);
  Arry_d_3=abs(Arry_det_3);  //absolute value of Arry_det
  
  
  // Find the adjoint matrix of Arry in transpose order for 3rd cycle
  
  Arry_adj_3[0][0]=Arry_3[1][1]*Arry_3[2][2]-Arry_3[1][2]*Arry_3[2][1];
  Arry_adj_3[1][0]=Arry_3[1][2]*Arry_3[2][0]-Arry_3[1][0]*Arry_3[2][2];
  Arry_adj_3[2][0]=Arry_3[1][0]*Arry_3[2][1]-Arry_3[1][1]*Arry_3[2][0];
    
  Arry_adj_3[0][1]=Arry_3[0][2]*Arry_3[2][1]-Arry_3[0][1]*Arry_3[2][2];
  Arry_adj_3[1][1]=Arry_3[0][0]*Arry_3[2][2]-Arry_3[0][2]*Arry_3[2][0];
  Arry_adj_3[2][1]=Arry_3[0][1]*Arry_3[2][0]-Arry_3[2][1]*Arry_3[0][0];
          
  Arry_adj_3[0][2]=Arry_3[0][1]*Arry_3[1][2]-Arry_3[0][2]*Arry_3[1][1];
  Arry_adj_3[1][2]=Arry_3[0][2]*Arry_3[1][0]-Arry_3[0][0]*Arry_3[1][2];
  Arry_adj_3[2][2]=Arry_3[0][0]*Arry_3[1][1]-Arry_3[0][1]*Arry_3[1][0];
  
  
  // Find the inverse matrix for 3rd cycle
  
  Arry_inv_3[0][0]=Arry_adj_3[0][0]/Arry_det_3; 
  Arry_inv_3[1][0]=Arry_adj_3[1][0]/Arry_det_3;
  Arry_inv_3[2][0]=Arry_adj_3[2][0]/Arry_det_3;

  Arry_inv_3[0][1]=Arry_adj_3[0][1]/Arry_det_3;
  Arry_inv_3[1][1]=Arry_adj_3[1][1]/Arry_det_3;
  Arry_inv_3[2][1]=Arry_adj_3[2][1]/Arry_det_3;

  Arry_inv_3[0][2]=Arry_adj_3[0][2]/Arry_det_3;
  Arry_inv_3[1][2]=Arry_adj_3[1][2]/Arry_det_3;
  Arry_inv_3[2][2]=Arry_adj_3[2][2]/Arry_det_3;
 
  
  // Find the time matrix for 3rd cycle
  
  time_Arry_3[0]= Arry_inv_3[0][0]*nu_start+Arry_inv_3[0][1]*nu_mid_array[2]+Arry_inv_3[0][2]*nu_end;
  time_Arry_3[1]= Arry_inv_3[1][0]*nu_start+Arry_inv_3[1][1]*nu_mid_array[2]+Arry_inv_3[1][2]*nu_end;
  time_Arry_3[2]= Arry_inv_3[2][0]*nu_start+Arry_inv_3[2][1]*nu_mid_array[2]+Arry_inv_3[2][2]*nu_end;
  
  
  //Profile for fourth cycle
  // buid the Y column, i.e. nu_ary as a function of Y=x^2+x+1 for 4th cycle
  Arry_4[0][0]=pow(time_start_array[3],2);
  Arry_4[0][1]=time_start_array[3];
  Arry_4[0][2]=1;
  
  Arry_4[1][0]=pow(time_mid_4,2);
  Arry_4[1][1]=time_mid_4;
  Arry_4[1][2]=1;
  
  Arry_4[2][0]=pow(time_change_array[3],2);
  Arry_4[2][1]=time_change_array[3];
  Arry_4[2][2]=1;
  
  // Find the determinent of Arry for 4th cycle
  
  Arry_det_4= Arry_4[0][0]*(Arry_4[1][1]*Arry_4[2][2]-Arry_4[1][2]*Arry_4[2][1])-Arry_4[0][1]*(Arry_4[1][0]*Arry_4[2][2]-Arry_4[1][2]*Arry_4[2][0])+Arry_4[0][2]*(Arry_4[1][0]*Arry_4[2][1]-Arry_4[1][1]*Arry_4[2][0]);
  Arry_d_4=abs(Arry_det_4);  //absolute value of Arry_det
  
  
  // Find the adjoint matrix of Arry in transpose order for 4th cycle
  
  Arry_adj_4[0][0]=Arry_4[1][1]*Arry_4[2][2]-Arry_4[1][2]*Arry_4[2][1];
  Arry_adj_4[1][0]=Arry_4[1][2]*Arry_4[2][0]-Arry_4[1][0]*Arry_4[2][2];
  Arry_adj_4[2][0]=Arry_4[1][0]*Arry_4[2][1]-Arry_4[1][1]*Arry_4[2][0];
    
  Arry_adj_4[0][1]=Arry_4[0][2]*Arry_4[2][1]-Arry_4[0][1]*Arry_4[2][2];
  Arry_adj_4[1][1]=Arry_4[0][0]*Arry_4[2][2]-Arry_4[0][2]*Arry_4[2][0];
  Arry_adj_4[2][1]=Arry_4[0][1]*Arry_4[2][0]-Arry_4[2][1]*Arry_4[0][0];
          
  Arry_adj_4[0][2]=Arry_4[0][1]*Arry_4[1][2]-Arry_4[0][2]*Arry_4[1][1];
  Arry_adj_4[1][2]=Arry_4[0][2]*Arry_4[1][0]-Arry_4[0][0]*Arry_4[1][2];
  Arry_adj_4[2][2]=Arry_4[0][0]*Arry_4[1][1]-Arry_4[0][1]*Arry_4[1][0];
  
  
  // Find the inverse matrix for 4th cycle
  
  Arry_inv_4[0][0]=Arry_adj_4[0][0]/Arry_det_4; 
  Arry_inv_4[1][0]=Arry_adj_4[1][0]/Arry_det_4;
  Arry_inv_4[2][0]=Arry_adj_4[2][0]/Arry_det_4;

  Arry_inv_4[0][1]=Arry_adj_4[0][1]/Arry_det_4;
  Arry_inv_4[1][1]=Arry_adj_4[1][1]/Arry_det_4;
  Arry_inv_4[2][1]=Arry_adj_4[2][1]/Arry_det_4;

  Arry_inv_4[0][2]=Arry_adj_4[0][2]/Arry_det_4;
  Arry_inv_4[1][2]=Arry_adj_4[1][2]/Arry_det_4;
  Arry_inv_4[2][2]=Arry_adj_4[2][2]/Arry_det_4;
 
  
  // Find the time matrix for 4th cycle
  
  time_Arry_4[0]= Arry_inv_4[0][0]*nu_start+Arry_inv_4[0][1]*nu_mid_array[3]+Arry_inv_4[0][2]*nu_end;
  time_Arry_4[1]= Arry_inv_4[1][0]*nu_start+Arry_inv_4[1][1]*nu_mid_array[3]+Arry_inv_4[1][2]*nu_end;
  time_Arry_4[2]= Arry_inv_4[2][0]*nu_start+Arry_inv_4[2][1]*nu_mid_array[3]+Arry_inv_4[2][2]*nu_end;

  
  //Profile for fifth cycle
  // buid the Y column, i.e. nu_ary as a function of Y=x^2+x+1 for 4th cycle
  Arry_5[0][0]=pow(time_start_array[4],2);
  Arry_5[0][1]=time_start_array[4];
  Arry_5[0][2]=1;
  
  Arry_5[1][0]=pow(time_mid_5,2);
  Arry_5[1][1]=time_mid_5;
  Arry_5[1][2]=1;
  
  Arry_5[2][0]=pow(time_change_array[4],2);
  Arry_5[2][1]=time_change_array[4];
  Arry_5[2][2]=1;
  
  // Find the determinent of Arry for 4th cycle
  
  Arry_det_5= Arry_5[0][0]*(Arry_5[1][1]*Arry_5[2][2]-Arry_5[1][2]*Arry_5[2][1])-Arry_5[0][1]*(Arry_5[1][0]*Arry_5[2][2]-Arry_5[1][2]*Arry_5[2][0])+Arry_5[0][2]*(Arry_5[1][0]*Arry_5[2][1]-Arry_5[1][1]*Arry_5[2][0]);
  Arry_d_5=abs(Arry_det_5);  //absolute value of Arry_det
  
  
  // Find the adjoint matrix of Arry in transpose order for 4th cycle
  
  Arry_adj_5[0][0]=Arry_5[1][1]*Arry_5[2][2]-Arry_5[1][2]*Arry_5[2][1];
  Arry_adj_5[1][0]=Arry_5[1][2]*Arry_5[2][0]-Arry_5[1][0]*Arry_5[2][2];
  Arry_adj_5[2][0]=Arry_5[1][0]*Arry_5[2][1]-Arry_5[1][1]*Arry_5[2][0];
    
  Arry_adj_5[0][1]=Arry_5[0][2]*Arry_5[2][1]-Arry_5[0][1]*Arry_5[2][2];
  Arry_adj_5[1][1]=Arry_5[0][0]*Arry_5[2][2]-Arry_5[0][2]*Arry_5[2][0];
  Arry_adj_5[2][1]=Arry_5[0][1]*Arry_5[2][0]-Arry_5[2][1]*Arry_5[0][0];
          
  Arry_adj_5[0][2]=Arry_5[0][1]*Arry_5[1][2]-Arry_5[0][2]*Arry_5[1][1];
  Arry_adj_5[1][2]=Arry_5[0][2]*Arry_5[1][0]-Arry_5[0][0]*Arry_5[1][2];
  Arry_adj_5[2][2]=Arry_5[0][0]*Arry_5[1][1]-Arry_5[0][1]*Arry_5[1][0];
  
  
  // Find the inverse matrix for 4th cycle
  
  Arry_inv_5[0][0]=Arry_adj_5[0][0]/Arry_det_5; 
  Arry_inv_5[1][0]=Arry_adj_5[1][0]/Arry_det_5;
  Arry_inv_5[2][0]=Arry_adj_5[2][0]/Arry_det_5;

  Arry_inv_5[0][1]=Arry_adj_5[0][1]/Arry_det_5;
  Arry_inv_5[1][1]=Arry_adj_5[1][1]/Arry_det_5;
  Arry_inv_5[2][1]=Arry_adj_5[2][1]/Arry_det_5;

  Arry_inv_5[0][2]=Arry_adj_5[0][2]/Arry_det_5;
  Arry_inv_5[1][2]=Arry_adj_5[1][2]/Arry_det_5;
  Arry_inv_5[2][2]=Arry_adj_5[2][2]/Arry_det_5;
 
  
  // Find the time matrix for 4th cycle
  
  time_Arry_5[0]= Arry_inv_5[0][0]*nu_start+Arry_inv_5[0][1]*nu_mid_array[3]+Arry_inv_5[0][2]*nu_end;
  time_Arry_5[1]= Arry_inv_5[1][0]*nu_start+Arry_inv_5[1][1]*nu_mid_array[3]+Arry_inv_5[1][2]*nu_end;
  time_Arry_5[2]= Arry_inv_5[2][0]*nu_start+Arry_inv_5[2][1]*nu_mid_array[3]+Arry_inv_5[2][2]*nu_end;

 //pre-compute the profile for first cycle
  
  for( int i=time_start; i<time_change; i++ ) {
      para_nu = time_Arry[0]*pow(i,2)+ time_Arry[1]*(i);
    deltanu.push_back(para_nu);
      }
  
  //min and max value -  standardise values between 0 and 1
  
  para_nu_min=*std::min_element(deltanu.begin(), deltanu.end());
  para_nu_max=*std::max_element(deltanu.begin(), deltanu.end());
  
  para_nu_min_array[0]=para_nu_min;
 
  para_nu_max_array[0]=para_nu_max;
  
  
  //pre-compute the profile for second cycle
  
  for( int i=time_start_array[1]; i<time_change_array[1]; i++ ) {
      para_nu_array[1] = time_Arry_2[0]*pow(i,2)+ time_Arry_2[1]*(i);
    deltanu_array[1].push_back(para_nu_array[1]);
      }
  
  //min and max value -  standardise values between 0 and 1
  
  para_nu_min_array[1]=*std::min_element(deltanu_array[1].begin(), deltanu_array[1].end());
  para_nu_max_array[1]=*std::max_element(deltanu_array[1].begin(), deltanu_array[1].end());
  
  
  //pre-compute the profile for third cycle
  
  for( int i=time_start_array[2]; i<time_change_array[2]; i++ ) {
      para_nu_array[2] = time_Arry_3[0]*pow(i,2)+ time_Arry_3[1]*(i);
    deltanu_array[2].push_back(para_nu_array[2]);
      }
  
  //min and max value -  standardise values between 0 and 1
  
  para_nu_min_array[2]=*std::min_element(deltanu_array[2].begin(), deltanu_array[2].end());
  para_nu_max_array[2]=*std::max_element(deltanu_array[2].begin(), deltanu_array[2].end());
  
  
  //pre-compute the profile for fourth cycle
  
  for( int i=time_start_array[3]; i<time_change_array[3]; i++ ) {
      para_nu_array[3] = time_Arry_4[0]*pow(i,2)+ time_Arry_4[1]*(i);
    deltanu_array[3].push_back(para_nu_array[3]);
      }
  
  //min and max value -  standardise values between 0 and 1
  
  para_nu_min_array[3]=*std::min_element(deltanu_array[3].begin(), deltanu_array[3].end());
  para_nu_max_array[3]=*std::max_element(deltanu_array[3].begin(), deltanu_array[3].end());
  
  
//pre-compute the profile for fifth cycle
  
  for( int i=time_start_array[4]; i<time_change_array[4]; i++ ) {
      para_nu_array[4] = time_Arry_5[0]*pow(i,2)+ time_Arry_5[1]*(i);
    deltanu_array[4].push_back(para_nu_array[4]);
      }
  
  //min and max value -  standardise values between 0 and 1
  
  para_nu_min_array[4]=*std::min_element(deltanu_array[4].begin(), deltanu_array[4].end());
  para_nu_max_array[4]=*std::max_element(deltanu_array[4].begin(), deltanu_array[4].end());
  
  // clear previous values
  nu.clear();
  nu.resize(nodes, nu_min);
  pos = (nu_min>0)?1:-1;
  
  
  // value of nu at rest
  for( size_type i=0; i<nodes; i++ ) {
     temp_i.push_back(nu_min);
     temp_f.push_back(nu_end);
 }
  
  
  
  
// Compute the time of subsequent cycles
  for( size_type i=0; i<nodes; i++) {
    P[i] = nu[i]*prepropag.phiinit(configf);
  } 
  
  
  //initial time
  time = 0;
}

void CouplingParabolaContinuous::step(void) {
  // Return the right value at each time point
  time += deltat;

  for( size_type i=0; i<nodes; i++ ) {
      if ((time>=time_start)&& (time<=time_change)){
      nu[i]= nu_min + (nu_max-nu_min)*((((time_Arry[0]*pow(time,2))+(time_Arry[1]*time))-para_nu_min)/(para_nu_max-para_nu_min)); 
      }
      
      else if ((time>=time_start_array[1])&& (time<=time_change_array[1])){
       nu[i]= nu_min + (nu_mid_array[1]-nu_min)*((((time_Arry_2[0]*pow(time,2))+(time_Arry_2[1]*time))-para_nu_min_array[1])/(para_nu_max_array[1]-para_nu_min_array[1])); 
      }

      
      else if ((time>=time_start_array[2])&& (time<=time_change_array[2])){
       nu[i]= nu_min + (nu_mid_array[2]-nu_min)*((((time_Arry_3[0]*pow(time,2))+(time_Arry_3[1]*time))-para_nu_min_array[2])/(para_nu_max_array[2]-para_nu_min_array[2])); 
      }

      
      
      else if ((time>=time_start_array[3])&& (time<=time_change_array[3])){
       nu[i]= nu_min + (nu_mid_array[3]-nu_min)*((((time_Arry_4[0]*pow(time,2))+(time_Arry_4[1]*time))-para_nu_min_array[3])/(para_nu_max_array[3]-para_nu_min_array[3])); 
      }
      
      
      else if ((time>=time_start_array[4])&& (time<=time_change_array[4])){
       nu[i]= nu_min + (nu_mid_array[4]-nu_min)*((((time_Arry_5[0]*pow(time,2))+(time_Arry_5[1]*time))-para_nu_min_array[4])/(para_nu_max_array[4]-para_nu_min_array[4])); 
      }

      
      else
      {nu[i]=temp_i[i];} 
      
      
      
  }
  
  
    
  Coupling::step();
}


CouplingParabolaContinuous::CouplingParabolaContinuous( size_type nodes, double deltat, size_type index,
                                        const Propagator& prepropag, const Population& postpop,
                                        double tempf )
  : Coupling(nodes,deltat,index,prepropag,postpop) {
  // total simulation time is stored in tempf (defined in solver.cpp line 73)
  time_tot = tempf;
}

CouplingParabolaContinuous::~CouplingParabolaContinuous(void) {
}
