/** @file coupling_parabola.h
  @brief A new Coupling class, where the nu_ab vary like a parabola or half cycle.

  The smooth variation of nu_ab is defined by a parabolic function y=ax^2+bx+c.


  @author Farah Deeba, 
*/

#ifndef COUPLING_PARABOLA_H
#define COUPLING_PARABOLA_H

// Other neurofield headers
#include "configf.h"    // Configf;
#include "coupling.h"   // Coupling
#include "propagator.h" // Propagator;
#include "population.h" // Population;


// C++ standard library headers
#include <vector> // std::vector;

class CouplingParabola : public Coupling {
  //CouplingParabola();
  //CouplingParabola(CouplingParabola&);
 protected:
  double t_start, t_mid, t_end;
  double nu_start, nu_mid, nu_end;
  double Arry[3][3];
  double nu_min = 0.0, nu_max = 0.0;
  double Arry_adj[3][3];
  double Arry_det;
  double Arry_d;
  double Arry_inv[3][3];
  double time_Arry[3];  
  double para_nu;
  double time_tot;
  double para_nu_min = 0.0,  para_nu_max = 0.0;
  double time = 0.0;
  double time_int = 0.0;
  std::vector<double> temp_i, temp_f;
  std::vector<double> deltanu;   
     
     
//public:
    
  //void init( Configf& configf );
  //void step(void);
 // void find(void);
  //CouplingParabola( size_type nodes, double deltat, size_type index,
                 //   const Propagator& prepropag, const Population& postpop, double tempf );
  //virtual ~CouplingParabola(void);
//};

 public:
 CouplingParabola() = delete;                          // No default constructor allowed.
  CouplingParabola(const CouplingParabola&) = delete; // No copy constructor allowed.

  void init( Configf& configf ) override;
  void step() override;
  void find();
  CouplingParabola( size_type nodes, double deltat, size_type index,
                    const Propagator& prepropag, const Population& postpop, double tempf );
  virtual ~CouplingParabola() override;
};

#endif //COUPLING_PARABOLA_H

