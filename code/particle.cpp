# include "./particle.hpp"


// vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
void particle::set_vcut()
{
}

// saving the MC trial move
void particle::tra_move(const pvector<ntype,3>& delr)
  {
    r+=delr;
  }

// storing the position if the MC move is accepted
void particle::store()
  {
    rold = r;
  }

//restore the old position if the MC move is rejected
void particle::restore()
  {
    r = rold;
  }

void particle::set_sigma(const ntype sig)
  {
    sigma = sig;
  }

void particle::set_epsilon(const ntype eps)
  {
    epsilon = eps;
  }

void particle::set_rcut(const ntype rcut)
  {
    rc = rcut;
  }

  // methods for MD

// implementation of action of operator exp(iLp*dt)
void particle::expiLp(const ntype dt)
{
  v += f/m*dt;
}

// implementation of action of operator exp(iLq*dt)
void particle::expiLq(const ntype dt)
{
  r += v*dt;
}

//basic constructor
particle::particle()
  {
    sigma=1.0;
    epsilon=1.0;
    n = 12;
    m=1.0;
    rc=2.5;
    vcut = 0.0;
  }

//constructor using params
particle::particle(const simpars& params)
{
  sigma = params.sigma;
  epsilon = params.epsilon;
  rc = params.rc;
  n = params.n;
  vcut = 0.0;
}


//calculate the potential
ntype particleLJ::vij(const particleLJ& P, const pvector<ntype,3>& L) const
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
void particleLJ::set_vcut()
  {
    const ntype rijsq = rc * rc;
    const ntype srij2 = sigma * sigma / rijsq;
    const ntype srij6 = srij2 * srij2 * srij2;
    const ntype srij12 = srij6 * srij6;
    vcut = 4.0*epsilon*(srij12 - srij6);
  }


//calculate the force
pvector<ntype,3> particleLJ::fij(const particleLJ& P, const pvector<ntype,3>& L, ntype &vij, ntype& vijs, ntype &wij) const
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




