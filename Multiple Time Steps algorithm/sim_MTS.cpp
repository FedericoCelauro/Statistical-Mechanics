#include "sim_MTS.hpp"
#include <fstream>
#include <iomanip> // for setprecision()
#include <sys/stat.h>


using ntype=double;
using simp = simpars;

//generic simulation class


// calculate energy of particle i
// if opt=1 calculate energies only for i < j
ntype sim::calcenergyi(int i, const int opt=1) const
{
    ntype enei=0.0;
    for (int j = 0; j < pars.Np; j++) // pars.Np total number of particles
      {
        if (opt==1 && i >= j)
          continue;
        if (i==j)
          continue;
        enei += parts[i].vij(parts[j], pars.L);
      }
    return enei;
  }


//calculate total energy
ntype sim::totenergy() const
{
    ntype ene=0.0;
    for (int i=0; i < pars.Np; i++)
      {
        ene+=calcenergyi(i, 1);
      }
    return ene;
  }


void sim::pbc(int i) //implements periodic boundary conditions
  {
    auto Dr = parts[i].r;
    Dr = pars.L.mulcw(rint(Dr.divcw(pars.L))); // L*rint(Dr/L)
    parts[i].r = parts[i].r - Dr;
  }


// saves all particles positions
void sim::save_mgl_snapshot(const long int t)
  {
     std::fstream f;

     struct stat buffer{};
     if (stat("mgl snapshots", &buffer) != 0)
        if(mkdir("mgl snapshots", 0777) == -1)
        {
          std::cerr << "Failed to create folder 'mgl snapshots'" << std::endl;
          exit(1);
        }

    const std::string s = "mgl snapshots/cnf-" + std::to_string(t) + ".mgl";

     f.open(s, std::ios::out|std::ios::trunc);
     for (int i=0; i < pars.Np; i++)
       {
         f << parts[i].r(0) << " " << parts[i].r(1) << " " <<
           parts[i].r(2) << " @ " << pars.sigma*0.5 << "\n";
       }
     f.close();
  }


//read parameters from file
void sim::read_params(const std::string& filename)
{
  const simpars tempParams(filename);

  pars = tempParams;
  std::cout << pars << std::endl;
}


//get n_iter
int sim::get_n_iter() const
{
  return pars.n_iter;
}


// prepares initial configuration
void sim::prepare_initial_conf()
  {
    // initial configuration as a SC lattice
    int cc=0;

    //in a SC lattice each lattice point holds a particle
    pars.Np = pars.nx*pars.ny*pars.nz;
    parts.resize(pars.Np);

    ntype vcell = pow(pars.sigma,3.0); //initial volume that holds a single particle
    ntype rhomax = 1.0/vcell; //maximum possible density (if Soft Sphere or LN we won't have interactions)
    ntype sf = cbrt(rhomax / pars.rho); //scaling factor

    pars.L ={static_cast<ntype>(pars.nx), static_cast<ntype>(pars.ny), static_cast<ntype>(pars.nz)}; //the box dimensions in units of nx, ny, nzz
    ntype clen = sf*pars.sigma; //scaled distance between the particles
    pars.L *= clen; // scaled length of the box

    //particle initialization
    for (int ix = 0; ix < pars.nx; ix++)
      for (int iy = 0; iy < pars.ny; iy++)
        for (int iz = 0; iz < pars.nz; iz++)
          {
            parts[cc].r = {ix*clen, iy*clen, iz*clen};
            parts[cc].r -= pars.L*0.5;
            parts[cc].set_sigma(pars.sigma);
            parts[cc].set_epsilon(pars.epsilon);
            parts[cc].set_rcut(pars.rc);
            cc++;
          }
    // ...or BCC or FCC lattice
  }



//initializing random number generator
void sim::init_rng() const
{
    if (pars.seed < 0)
      rng.rseed();
    else
      rng.seed(pars.seed);
  }

void sim::run()
  {
    // intentionally void
  };


//start measurements
void mdsim::init_measures() const
{
  // open files in writing mode to truncate them to 0
  std::fstream f;

  f.open("totenergy"+std::to_string(i_print)+".txt", std::ios::out|std::ios::trunc);
  if(!f.is_open())
  {
    std::cerr << "Error: can't open file 'totenergy"<<std::to_string(i_print)<<".txt'" << std::endl;
    exit(1);
  }
  f.close();
}


// calculate total kinetic energy
ntype mdsim::calcK() const
{
    ntype K=0.0;
    for (int i=0; i < pars.Np; i++)
      {
        K += 0.5*parts[i].m*parts[i].v*parts[i].v;
      }
    return K;
  }



//save measurements
void mdsim::save_measures(const long int t) const
{
  std::fstream f;
  f.open("totenergy"+std::to_string(i_print)+".txt", std::ios::out|std::ios::app);

  const ntype K=calcK();
  if (f.is_open())
  {
    // K+Us is the total conserved energy, where Us is the shifted potential energy
    f << t << " " << std::setprecision(15) << K+Us << "\n";

  }else
  {
    std::cerr << "Error: can't open file 'totenergy.txt'" << std::endl;
    exit(1);
  }
  f.sync();
  f.close();
}


// calculate force acting on each particle
void mdsim::calc_forces( int i_start,  int i_end, ntype& Us, ntype& U)
  {
    ntype energy, energys, w;
    pvec3d force;

    //reset the energy
    U = 0.0;
    Us = 0.0;

    //reset the total force to zero
    for (int i=0; i < pars.Np; i++)
      parts[i].f = {0, 0, 0};


    // loop on heavier particles (only works if i_end > N1)
    for (int i=pars.N1; i < i_end; i++)
    {
      //force contribute of heavier particles on heavier particle i
      for (int j=i+1; j < pars.Np; j++)
      {
        force = parts[i].fij(parts[j], pars.L, energy, energys, w);
        parts[i].f += force;
        parts[j].f -= force;
        U += energy;
        Us += energys;
      }
      //saving force of heavier particles on heavier particles
      parts[i].f_id = parts[i].f;
      U_id = Us;
    }


    // loop on lighter particles (only works if i_start < N1)
    for (int i=i_start; i < pars.N1; i++)
    {
      //force contribute of lighter particles on lighter particle i
      for (int j=i+1; j < pars.N1; j++)
      {
        force = parts[i].fij(parts[j], pars.L, energy, energys, w);
        parts[i].f += force;
        parts[j].f -= force;
        U += energy;
        Us += energys;
      }
      //saving force of lighter particles on lighter particles
      parts[i].f_id = parts[i].f;
    }


  //force contribute of heavier particles on lighter particles
  for (int i=0; i < pars.N1; i++)
  {
    for (int j=pars.N1; j < pars.Np; j++)
    {
      force = parts[i].fij(parts[j], pars.L, energy, energys, w);
      parts[i].f += force;
      parts[j].f -= force;
      U += energy;
      Us += energys;
    }
  }

    //now we can add to the total force of the "other particles" the force term between them that
    //was not updated during the above loop (because these particles don't move, it doesn't change)
    for (int i=pars.N1 - i_start; i < pars.Np - i_end + pars.N1; i++)
      parts[i].f += parts[i].f_id;
  }


// Generate gaussian number by Box Muller algorithm
ntype mdsim::gauss()
  {
    // gaussian of variance 1 and 0 mean
    double x1, x2;
    do
      {
        x1 = rng.ranf();
        x2 = rng.ranf();
      }
    while (x1==0 || x2==0);
    return cos(2*M_PI*x2)*sqrt(-2.0*log(x1));
  }

void mdsim::prepare_initial_conf()
    {
      bc::prepare_initial_conf(); // place particles on a SC lattice
      pvec3d totmomentum = {0, 0, 0};
      ntype totM = 0;

      if(pars.mass2 < pars.mass1) //so that mass1 always less than mass2
      {
        ntype tempMass = pars.mass1;
        pars.mass1 = pars.mass2;
        pars.mass2 = tempMass;
      }

      pars.N1 = pars.Nratio * pars.Np;


      for(int i=0; i < pars.N1; i++) //first N1 particles with mass1
      {
        parts[i].m = pars.mass1;

        //initialize particle's velocity from Maxwell-Boltzmann distribution
        parts[i].v(0) = sqrt(pars.T / parts[i].m)*gauss();
        parts[i].v(1) = sqrt(pars.T / parts[i].m)*gauss();
        parts[i].v(2) = sqrt(pars.T / parts[i].m)*gauss();

        totmomentum += parts[i].v * parts[i].m;
        totM += parts[i].m;

        parts[i].set_vcut(); //set vcut
      }

      for(int i=pars.N1; i < pars.Np; i++) //N2 = Np - N1 particles with mass2
      {
        parts[i].m = pars.mass2;

        //initialize particle's velocity from Maxwell-Boltzmann distribution
        parts[i].v(0) = sqrt(pars.T / parts[i].m)*gauss();
        parts[i].v(1) = sqrt(pars.T / parts[i].m)*gauss();
        parts[i].v(2) = sqrt(pars.T / parts[i].m)*gauss();

        totmomentum += parts[i].v * parts[i].m;
        totM += parts[i].m;

        parts[i].set_vcut(); //set vcut
      }

      //set total momentum to zero
      totmomentum /= totM; //total velocity
      for(int i=0; i < pars.Np; i++)
        parts[i].v -= totmomentum;
    }

  void mdsim::Measures(long int t)
  {
    if (pars.outstps > 0 && t % pars.outstps== 0)
    {
      const ntype K=calcK();

      std::cout << std::setprecision(15) << "[step #" << t << "] \nTotal energy=" << K+Us << "\n";
      std::cout << std::setprecision(8) << "Excess energy=" << Us << "\n";
      std::cout << std::setprecision(8) << "Kinetic energy=" << K << " T="
        << 2*K/(3*(pars.Np-3)) << "\n";    }

    if (t >= pars.eqstps && pars.savemeasure > 0 && (t % pars.savemeasure == 0))
    {
      save_measures(t);
    }


    if (t > 0 && pars.save_mgl_snapshot > 0 && t % pars.save_mgl_snapshot == 0)
      bc::save_mgl_snapshot(t);

  }


  void mdsim::VelocityVerlet(const int i_start, const int i_end, ntype dt, ntype &Us, ntype &U)
  {
    for(int i=i_start; i < i_end; i++)
    {
      parts[i].expiLp(dt*0.5); //apply exp(iLp*dt/2)
      parts[i].expiLq(dt); //apply exp(iLq*dt)
      pbc(i);
    }
    calc_forces(i_start, i_end, Us, U);  //calculate the forces
    for(int i=i_start; i < i_end; i++)
      parts[i].expiLp(dt*0.5); //apply exp(iLp*dt/2)
  }



  void mdsim::MultipleTimeStep(ntype &Us, ntype &U)
  {
    for (long int t=0; t < pars.totsteps; t++)
    {
      for(int i=0; i < pars.M/2; i++)
        VelocityVerlet(0, pars.N1, pars.deltat,Us,U);

      VelocityVerlet(pars.N1, pars.Np, pars.dt,Us,U);
      Us = 0; //gives the right excess energy if M=1

      for(int i=0; i < pars.M/2; i++)
        VelocityVerlet(0, pars.N1, pars.deltat,Us,U);

      Us += U_id;
      Measures(t);
    }
  }


  void mdsim::NormalTimeStepIntegration(ntype &Us, ntype &U)
  {
    for (long int t=0; t < pars.totsteps; t++)
    {
      for(int j=0; j < pars.M; j++)
        VelocityVerlet(0, pars.Np, pars.deltat,Us,U);

      Measures(t);
    }
  }


  void mdsim::run(int i)
    {
      calc_forces(0, pars.Np, Us, U);
      std::cout << std::setprecision(15) << "\nTotal energy=" << calcK()+Us << "\n";
      std::cout << "Initial Excess energy=" << Us << "\n";
      std::cout << "Initial Kinetic energy=" << calcK() << "\n";

      i_print = i;

      std::fstream fo;
      init_measures();

      //deltat calculation
      pars.M = sqrt(pars.mass2/pars.mass1); //mass1 < mass2 always
      pars.deltat = pars.dt/pars.M;
      std::cout << "deltat = " << pars.deltat << ", M = " << pars.M <<  ", N1 = " << pars.N1 << "\n" << std::endl;


      if(pars.algotype == 0) //normal Velocity Verlet
        NormalTimeStepIntegration(Us, U);

      else if (pars.algotype == 1) // multiple time step algorithm
      {
        MultipleTimeStep(Us, U);

      }else
        std::cerr << "Wrong value for algotype variable.\n"
                     "This can take the following values: 0 for normal Velocity Verlet, 1 for multiple time step algorithm."<< std::endl;
    }