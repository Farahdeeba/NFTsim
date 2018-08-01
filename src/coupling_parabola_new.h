/** @file coupling_parabola.h
  @brief A new Coupling class, where the nu_ab vary like a parabola or half cycle.

  The smooth variation of nu_ab is defined by a parabolic function y=ax^2+bx+c.


  @author Farah Deeba, 
*/

#ifndef COUPLING_PARABOLA_NEW_H
#define COUPLING_PARABOLA_NEW_H

// Other neurofield headers
#include "configf.h"    // Configf;
#include "coupling.h"   // Coupling
#include "propagator.h" // Propagator;
#include "population.h" // Population;


// C++ standard library headers
#include <vector> // std::vector;

class CouplingParabolaNew : public Coupling {
 protected:
  double t_thr_up, t_thr_down;
  double nu_thr = 0.0, nu_max = 0.0;
  double time_tot;
  double time = 0.0;
  double time_int = 0.0;
  std::vector<double> temp_i;
  std::vector<double> deltanu;   
     
     

 public:
 CouplingParabolaNew() = delete;                          // No default constructor allowed.
  CouplingParabolaNew(const CouplingParabolaNew&) = delete; // No copy constructor allowed.

  void init( Configf& configf ) override;
  void step() override;
  void find();
  CouplingParabolaNew( size_type nodes, double deltat, size_type index,
                    const Propagator& prepropag, const Population& postpop, double tempf );
  virtual ~CouplingParabolaNew() override;
};

#endif //COUPLING_PARABOLA_H

