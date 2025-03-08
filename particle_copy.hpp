#ifndef _PARTICLE_
#define _PARTICLE_


#include "./pvector.hpp"

using ntype = double;
class particle
{
  pvector<ntype,3> rold; // to store particle's position 
protected:
  ntype vcut;
public:
  ntype sigma, epsilon, rc, n, m;
  pvector<ntype,3> r, v, f; // particle's position, velocity and force

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  void set_vcut();

  // saving the MC trial move
  void tra_move(const pvector<ntype,3>& delr)
    {
      r+=delr; 
    }

  // storing the position if the MC move is accepted
  void store()
    {
      rold = r;
    }

  //restore the old position if the MC move is rejected
  void restore()
    {
      r = rold;
    }

  void set_sigma(const ntype sig)
    {
      sigma = sig;
    }
  void set_epsilon(const ntype eps)
    {
      epsilon = eps;
    }
  void set_rcut(const ntype rcut)
    {
      rc = rcut;
    }

  // methods for MD

  // implementation of action of operator exp(iLp*dt)
  void expiLp(const ntype dt)
  {
    v += f/m*dt;
  }

  // implementation of action of operator exp(iLq*dt)
  void expiLq(const ntype dt)
  {
    r += v*dt;
  }

  //basic constructor
  particle()
    {
      sigma=1.0;
      epsilon=1.0;
      n = 12;
      rc=2.5;
      vcut = 0.0;
    }

  //constructor using params
  particle(const simpars& params)
  {
    sigma = params.sigma;
    epsilon = params.epsilon;
    rc = params.rc;
    n = params.n;
    vcut = 0.0;
  }
};


//Lennard-Jones interacting particles
class particleLJ: public particle
{
public:
  ntype vij(const particleLJ& P, const pvector<ntype,3>& L) const
  {
      ntype ene=0.0;

      ntype dist = (r - P.r - L.mulcw(rint((r - P.r) / L))).norm(); //calculate the distance between particles (minimum image convention)

      if (dist<=rc) //if distance less than cutoff
      {
        ene = 4.0*epsilon*( pow(sigma/dist, 12) - pow(sigma/dist, 6) );
        //ene -= vcut; //to make a continuous potential
      }

      return ene;
    }

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut()
    {
      const ntype rijsq = rc * rc;
      const ntype srij2 = sigma * sigma / rijsq;
      const ntype srij6 = srij2 * srij2 * srij2;
      const ntype srij12 = srij6 * srij6;
      vcut = 4.0*epsilon*(srij12 - srij6);
    }


  //calculate the force
  pvector<ntype,3> fij(const particleLJ P, const pvector<ntype,3>& L, ntype &vij, ntype& vijs, ntype &wij) const
  {
    pvector<ntype,3> fijv = {0, 0, 0}, rij = r - P.r;

    // MINIMUM IMAGE CONVENTION
    rij = rij - L.mulcw(rint(rij.divcw(L))); // Dr - L*rint(Dr/L)

    // potential and force calculation
    const ntype rijsq = rij * rij;


    if (rijsq < rc*rc) // interaction potential cut-off
    {
      const ntype srij2 = sigma * sigma / rijsq;
      const ntype srij6 = srij2 * srij2 * srij2;
      const ntype srij12 = srij6 * srij6;

      vij = 4*epsilon*(srij12 - srij6);
      vijs = vij - vcut;
      fijv = 24 * epsilon * (2 * srij12 - srij6) / rijsq * rij;
    }else
    {
      vij = 0;
      vijs = 0;
    }
    wij = fijv*rij;

    return fijv;
  }

};




//Soft Sphere interacting particles
class particleSS: public particle
{
  public:
  ntype vij(const particleSS& P, const pvector<ntype,3>& L) const
  {
      ntype ene = 0.0;

      ntype dist = (r - P.r - L.mulcw(rint((r - P.r) / L))).norm(); //calculate the distance between particles (minimum image convention)

      if (dist<=rc) //if distance less than cutoff
      {
        ene = epsilon*pow( sigma/dist, n );
        //ene -= vcut; //to make a continuous potential
      }
      return ene;
    }

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut()
    {
      vcut = epsilon*pow( sigma/rc, n );
    }


  pvector<ntype,3> fij(const particleSS& P, const pvector<ntype,3>& L, ntype &vij, ntype& vijs, ntype &wij) const
  {
  pvector<ntype,3> fijv = {0, 0, 0}, rij = r - P.r;

  // MINIMUM IMAGE CONVENTION
  rij = rij - L.mulcw(rint(rij.divcw(L))); // Dr - L*rint(Dr/L)

  // potential and force calculation
  const ntype rijsq = rij * rij;

  if (rijsq < rc*rc) // interaction potential cut-off
  {
    vij = epsilon*pow(sigma/rij.norm(), n);
    vijs = vij - vcut;
    fijv = n * vij / rijsq * rij;
  }else
  {
    vij = 0;
    vijs = 0;
  }
  wij = fijv*rij;

  return fijv;
}

};
#endif 
