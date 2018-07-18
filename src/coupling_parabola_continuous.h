/** @file coupling_parabola.h
  @brief A new Coupling class, where the nu_ab vary like a parabola or half cycle.

  The smooth variation of nu_ab is defined by a parabolic function y=ax^2+bx+c.


  @author Farah Deeba, 
*/

#ifndef COUPLING_PARABOLA_CONTINUOUS_H
#define COUPLING_PARABOLA_CONTINUOUS_H

// Other neurofield headers
#include "configf.h"    // Configf;
#include "coupling.h"   // Coupling
#include "propagator.h" // Propagator;
#include "population.h" // Population;


// C++ standard library headers
#include <vector> // std::vector;

class CouplingParabolaContinuous : public Coupling {
  CouplingParabolaContinuous();
  CouplingParabolaContinuous(CouplingParabolaContinuous&);
 protected:
  double time_start, time_mid, time_end, time_mid_2, time_mid_3, time_mid_4, time_mid_5;
  double time_change;
  double nu_start, nu_mid, nu_end, nu_mid_2, nu_mid_3, nu_mid_4, nu_mid_5;
  double Arry[3][3];
  double nu_min = 0.0, nu_max = 0.0;
  
  double Arry_adj[3][3];
  double Arry_det;
  double Arry_d;
  double Arry_inv[3][3];
  double time_Arry[3];  
  
  double Arry_2[3][3];
  double Arry_adj_2[3][3];
  double Arry_det_2;
  double Arry_d_2;
  double Arry_inv_2[3][3];
  double time_Arry_2[3];
  
  double Arry_3[3][3];
  double Arry_adj_3[3][3];
  double Arry_det_3;
  double Arry_d_3;
  double Arry_inv_3[3][3];
  double time_Arry_3[3];  
  
  double Arry_4[3][3];
  double Arry_adj_4[3][3];
  double Arry_det_4;
  double Arry_d_4;
  double Arry_inv_4[3][3];
  double time_Arry_4[3];  
  
  double Arry_5[3][3];
  double Arry_adj_5[3][3];
  double Arry_det_5;
  double Arry_d_5;
  double Arry_inv_5[3][3];
  double time_Arry_5[3];  
  
 
  double para_nu;
  double para_nu_array[5];
  double time_tot;
  double para_nu_min = 0.0,  para_nu_max = 0.0;

  double para_nu_2;
  double para_nu_min_2 = 0.0,  para_nu_max_2 = 0.0;
  
  double para_nu_min_array[5],  para_nu_max_array[5];
  //double para_nu_neg_array[5];
  //double para_nu_min_neg_array[5],  para_nu_max_neg_array[5];
  double time = 0.0;
  double time_int = 0.0;
  double amp_inc;
  double time_laps;
  double time_start_array[5];
  double time_change_array[5];
  double nu_mid_array[5];
  std::vector<double> temp_i, temp_f, temp_nf;
  std::vector<double> deltanu;   
  std::vector<double> deltanu_array[5];
  //std::vector<double> deltanu_2;  
     
public:
  void init( Configf& configf );
  void step(void);
  void find(void);
  CouplingParabolaContinuous( size_type nodes, double deltat, size_type index,
                    const Propagator& prepropag, const Population& postpop, double tempf );
  virtual ~CouplingParabolaContinuous(void);
};

#endif //COUPLING_PARABOLA_H

