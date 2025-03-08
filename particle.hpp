#ifndef _PARTICLE_
#define _PARTICLE_

#include "./pvector.hpp"
#include "params.hpp"


using ntype = double;

class particle
{
  pvector<ntype,3> rold; // to store particle's position 
protected:
  ntype vcut;
public:
  ntype sigma, epsilon, rc, n, m;
  pvector<ntype,3> r, v, f; // particle's position, velocity and force

  //basic constructor
  particle();

  //constructor using params
  explicit particle(const simpars& params);

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  static void set_vcut();

  void set_sigma(ntype sig);

  void set_epsilon(ntype eps);

  void set_rcut(ntype rcut);

  //methods for MC

  // saving the MC trial move
  void tra_move(const pvector<ntype,3>& delr);

  // storing the position if the MC move is accepted
  void store();

  //restore the old position if the MC move is rejected
  void restore();

  // methods for MD

  // implementation of action of operator exp(iLp*dt)
  void expiLp(ntype dt);

  // implementation of action of operator exp(iLq*dt)
  void expiLq(ntype dt);
};


//Lennard-Jones interacting particles
class particleLJ: public particle
{
public:
  ntype vij(const particleLJ& P, const pvector<ntype,3>& L) const;

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut();

  //calculate the force
  pvector<ntype,3> fij(const particleLJ& P, const pvector<ntype,3>& L, ntype &vij, ntype& vijs, ntype &wij) const;

};


#endif 
