#ifndef PROPAG_H
#define PROPAG_H

#include"population.h"

using std::string;

class Propag : public NF
{
  Propag(); // no default constructor
  Propag(Propag&); // no copy constructor
protected:
  void init( Configf& configf );
  //void restart( Restartf& restartf );
  //void dump( Dumpf& dumpf ) const;

  Population& prepop;
  Population& postpop;
  Tau tau; // tau_ab
  int longside;
  vector<double> p; // phi_ab
public: 
  Propag( int nodes, double deltat, int index, Population& prepop,
      Population& postpop, int longside, string topology );
  virtual ~Propag(void);
  virtual void step(void); 
  double phiinit( Configf& configf ) const;
  virtual const vector<double>& phi(void) const;
  inline double operator[]( int node ) const;
  virtual void output( Output& output ) const;
};

double Propag::operator[]( int node ) const
{
  return p[node];
}

#endif