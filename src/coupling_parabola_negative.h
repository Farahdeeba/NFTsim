/** @file coupling_parabola.h
  @brief A new Coupling class, where the nu_ab vary like a parabola or half cycle.

  The smooth variation of nu_ab is defined by a parabolic function y=ax^2+bx+c.


  @author Farah Deeba, 
*/

#ifndef COUPLING_PARABOLA_NEGATIVE_H
#define COUPLING_PARABOLA_NEGATIVE_H

// Other neurofield headers
#include "configf.h"    // Configf;
#include "coupling.h"   // Coupling
#include "propagator.h" // Propagator;
#include "population.h" // Population;


// C++ standard library headers
#include <vector> // std::vector;

class CouplingParabolaNegative : public Coupling {
  CouplingParabolaNegative();
  CouplingParabolaNegative(CouplingParabolaNegative&);
 protected:
  double t_start, t_mid, t_end, t_mid_neg;
  double t_change;
  double nu_start, nu_mid, nu_end, nu_mid_neg;
  double Arry[3][3];
  double nu_min = 0.0, nu_max = 0.0, nu_max_neg=0.0;
  double Arry_adj[3][3];
  double Arry_det;
  double Arry_d;
  double Arry_inv[3][3];
  double time_Arry[3];  
  double Arry_neg[3][3];
  double Arry_adj_neg[3][3];
  double Arry_det_neg;
  double Arry_d_neg;
  double Arry_inv_neg[3][3];
  double time_Arry_neg[3];  
  double para_nu;
  double time_tot;
  double para_nu_min = 0.0,  para_nu_max = 0.0;
  double para_nu_neg;
  double para_nu_min_neg = 0.0,  para_nu_max_neg = 0.0;
  double time = 0.0;
  double time_int = 0.0;
  std::vector<double> temp_i, temp_f, temp_nf;
  std::vector<double> deltanu;   
  std::vector<double> deltanu_neg;  
     
public:
  void init( Configf& configf );
  void step(void);
  void find(void);
  CouplingParabolaNegative( size_type nodes, double deltat, size_type index,
                    const Propagator& prepropag, const Population& postpop, double tempf );
  virtual ~CouplingParabolaNegative(void);
};

#endif //COUPLING_PARABOLA_H

