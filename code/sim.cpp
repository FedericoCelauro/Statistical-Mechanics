#include "sim.hpp"

#include <string>
#include <fstream>
#include <iomanip> // for setprecision()
#include <sys/stat.h>


using ntype=double;
using simp = simpars;


//generic simulation class

// calculate energy of particle i
// if opt=1 calculate energies only for i < j
ntype sim::calcenergyi(const int i, const int opt) const
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

//implements periodic boundary conditions
void sim::pbc(const int i)
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

     /*boost::filesystem::path folderPath = "mgl snapshots";

     if(!exists(folderPath))
       if (!create_directory(folderPath))
         std::cerr << "Failed to create folder: " << folderPath << std::endl;
*/
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

    //particle initialization (Simple Cubic)
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
  }



//MC simulation

//  trial move
void mcsim::alpha(const int i)
  {
    //move the particle using as max displacement pars.deltra (uniform distribution)
    pvector<ntype, 3> delr = {
      pars.deltra * 2.0 * (rng.ranf() - 0.5), pars.deltra * 2.0 * (rng.ranf() - 0.5),
      pars.deltra * 2.0 * (rng.ranf() - 0.5)
    };
    parts[i].tra_move(delr);
    pbc(i); //apply pbc
  }


// acceptance move
void mcsim::acc(const int i, const ntype eno)
  {
    // assume kB=1 (Boltzmann constants) so that beta=1/T

    const ntype delu= calcenergyi(i)-eno; //calculate the energy difference

    //reject the move if deltaU > 0 and >= exp(-deltaU/T)
    if (delu > 0.0 && rng.ranf() >= exp(-delu/pars.T)) //Metropolis algorithm
      {
        // reject move
        tra_rej++;
        parts[i].restore();
      }
  }


// NTV move
void mcsim::move_NTV(const int i)
  {
    // 1) store the energy of the i-th particle before the move
    const ntype eno = calcenergyi(i);
    // 2) store the position of the i-th particle before the move
    parts[i].store();
    // 3) trial move
    alpha(i);
    // 4) acceptance move for the i-th particle
    acc(i, eno);
  }


// adjust the deltra parameter looking at the acceptance rates
void mcsim::calc_acceptance_and_adjust()
  {
    if (tot_tra > 0)
      {
        // acceptance rate
        const ntype r_tra = static_cast<ntype>(tot_tra - tra_rej) / tot_tra;

        std::cout << "\nrate tra: " << r_tra << " deltra = " << pars.deltra << "\n";
        if (r_tra > 0.5)
          {
            pars.deltra *= 1.1;
          }
        else if (r_tra < 0.4)
          {
            pars.deltra /= 1.1;
          }
        tot_tra=tra_rej = 0;
      }

    // adjust maximum volume "displacement" in alpha_box
    // so that acceptance rate of box move is around 0.5
    if (tot_vol > 0) // if we don't make volume moves, tot_vol = 0
    {

      const ntype r_vol = static_cast<ntype>(tot_vol - vol_rej) / tot_vol;

      std::cout << "rate vol: " << r_vol << " vmax=" << pars.vmax << "\n";
      if (r_vol > 0.5)
      {
        pars.vmax *= 1.1;
      }
      else if (r_vol < 0.4)
      {
        pars.vmax /= 1.1;
      }
      tot_vol=vol_rej=0.0;
    }
  }


//start measurements
void mcsim::init_measures() const
{
    // open files in writing mode to truncate them to 0
    std::fstream f;

    f.open("energy.txt", std::ios::out|std::ios::trunc);
    if(!f.is_open())
    {
      std::cerr << "Error: can't open file 'energy.txt'" << std::endl;
      exit(1);
    }
    f.close();

    if (pars.simtype==1) // simtype == 1 means an NPT simluations
    {
      // clear file used to store densities in NPT simulation
      f.open("density.txt", std::ios::out|std::ios::trunc);
      if(!f.is_open())
      {
        std::cerr << "Error: can't open file 'density.txt'" << std::endl;
        exit(1);
      }
      f.close();
    }
  }


//save measurements
void mcsim::save_measures(const long int t)
  {
    std::fstream f;
    f.open("energy.txt", std::ios::out|std::ios::app);

    const ntype Vtail = 0;//8.0*M_PI/3*pars.epsilon*pow(pars.sigma,3)*pars.Np/(pars.L(0)*pars.L(1)*pars.L(2))*(pow(pars.sigma/pars.rc, 9)/3 - pow(pars.sigma/pars.rc, 3));


    if (f.is_open())
    {
      // save potential energy per particle
      f << t << " " << (totenergy() + Vtail)/pars.Np << "\n";
    }else
    {
      std::cerr << "Error: can't open file 'energy.txt'" << std::endl;
      exit(1);
    }
    f.close();

    if (pars.simtype==1) // 1 means NPT simulation, save density in this case
    {
      f.open("density.txt", std::ios::out|std::ios::app);
      if (f.is_open())
      {
        // save potential energy per particle
        f << t << " " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
      }else
      {
        std::cerr << "Error: can't open file 'density.txt'" << std::endl;
        exit(1);
      }
      f.close();
    }
  }


// restore all particle positions
void mcsim::restore_all_pars()
{
  for (int i=0; i < pars.Np; i++)
    parts[i].restore();
}


// store all particle position
void mcsim::store_all_pars()
{
  for (int i=0; i < pars.Np; i++)
    parts[i].store();
}


// trial box move
void mcsim::alpha_box(ntype& DG, ntype& fact)
{

  //calculate old volume and old internal energy
  const ntype oldV =  pars.L(0)*pars.L(1)*pars.L(2);
  const ntype oldE = totenergy();

  // choose randomly a new volume
  const ntype newV = oldV +  pars.vmax * 2.0 * (rng.ranf() - 0.5);
  fact = cbrt(newV/oldV);

  // scale all particle positions and box size by factor "fact"
  pars.L *= fact;

  for (int i=0; i < pars.Np; i++)
    parts[i].r *= fact;


  // calculate \DeltaG (see pdf of lectures)
  DG = pars.P*(newV - oldV) + (totenergy() - oldE) - pars.Np*pars.T*log(newV/oldV);
}

// accept or reject box move
void mcsim::acc_box(const ntype DG, const ntype fact)
{
  //reject the move if DeltaG > 0 and >= exp(-DeltaG/T)
  if (DG > 0.0 && rng.ranf() >= exp(-DG/pars.T)) //Metropolis algorithm
  {
    // reject move
    vol_rej++;
    restore_all_pars();
    pars.L /= fact;
  }
}

//box move of isobaric enseble
void mcsim::move_box()
   {
  ntype DG, fact; // DG is \Delta G as discussed during lectures (see pdf of lectures)
  store_all_pars(); // store all particle positions before attempting a box move
  alpha_box(DG, fact);
  acc_box(DG, fact);
   }

//run the simulation
void mcsim::run()
  {
    tot_tra = tot_vol = tra_rej = vol_rej = 0;

    //open files
    init_measures();

    //loop on total number of steps
    for (long int t = 0; t < pars.totsteps; t++)
      {

      // attempt to move all particles (each MC move is Np possible particle moves)
        for (int i = 0; i < pars.Np; i++)
          {
            // choose between box move and particle move
            if (pars.simtype==1 && rng.ranf() < 1.0/(pars.Np+1))
            {
              move_box();
              tot_vol++;
            }
            else
            {
              const int ip = rng.ranf()*pars.Np;
              move_NTV(ip);
              tot_tra++;
            }
          }

        //adjust the deltra parameter
        if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 &&
            t < pars.maxadjstps)
        {
          calc_acceptance_and_adjust();
        }

        // saving measures
        if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0 && t > pars.maxadjstps && t > pars.eqstps)
          {
            save_measures(t);
          }

        // print outputs
        if (t > 0 && t % pars.outstps == 0)
          {
            std::cout << "Step #" << t << "\n";
            // per confrontarsi con il Johnson si deve calcolare l'energia interna ovvero sommare il contributo
            // di energia cinetica
            std::cout << "total energy per particle is " << totenergy()/pars.Np << "\n";
            std::cout << "box size: " << pars.L(0) << " " << pars.L(1) << " " << pars.L(2) << "\n" <<std::endl;
          }

        // save mgl snapshot
        if (t > 0 && pars.save_mgl_snapshot > 0 &&
            t % pars.save_mgl_snapshot == 0 && t > pars.maxadjstps)
          {
            save_mgl_snapshot(t);
          }
      }
  }

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
void mdsim::calc_forces(ntype& Us, ntype& U)
  {
    ntype energy, energys, w;
    pvec3d force;

    //reset the energy
    U = 0.0;
    Us = 0.0;


    //reset the total force to zero
    for (int i=0; i < pars.Np; i++)
      parts[i].f = {0, 0, 0};

    // loop on particles
    for (int i=0; i < pars.Np; i++)
    {
      for (int j=i+1; j < pars.Np; j++)
      {
        force = parts[i].fij(parts[j], pars.L, energy, energys, w);
        parts[i].f += force;
        parts[j].f -= force;
        U += energy;
        Us += energys;
      }
    }
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


    for(int i=0; i < pars.Np; i++) //first N1 particles with mass1
    {
      parts[i].m = pars.mass;

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

void mdsim::Measures(const long int t)
{
  if (pars.outstps > 0 && t % pars.outstps== 0)
  {
    const ntype K=calcK();

    std::cout << std::setprecision(15) << "[step #" << t << "] \nTotal energy=" << K+Us << "\n";
    std::cout << std::setprecision(8) << "Excess energy=" << Us << "\n";
    std::cout << std::setprecision(8) << "Kinetic energy=" << K << "\n";
  }

  if (t >= pars.eqstps && pars.savemeasure > 0 && (t % pars.savemeasure == 0))
    save_measures(t);

  if (t > 0 && pars.save_mgl_snapshot > 0 && t % pars.save_mgl_snapshot == 0)
    bc::save_mgl_snapshot(t);
}


void mdsim::Measures_LF(const long int t, const ntype dt)
{

  //adjust velocity to the right step
  //(works best if savemeasure = 1)
  for(int i=0; i < pars.Np; i++)
    parts[i].expiLp(dt*0.5);

  if (pars.outstps > 0 && t % pars.outstps== 0)
  {
    const ntype K=calcK();

    std::cout << std::setprecision(15) << "[step #" << t << "] \nTotal energy=" << K+Us << "\n";
    std::cout << std::setprecision(8) << "Excess energy=" << Us << "\n";
    std::cout << std::setprecision(8) << "Kinetic energy=" << K << "\n";
  }

  if (t >= pars.eqstps && pars.savemeasure > 0 && (t % pars.savemeasure == 0))
    save_measures(t);

  if (t > 0 && pars.save_mgl_snapshot > 0 && t % pars.save_mgl_snapshot == 0)
    save_mgl_snapshot(t);

  //adjust velocity to the right step
  for(int i=0; i < pars.Np; i++)
    parts[i].expiLp(-dt*0.5);
}

void mdsim::VelocityVerlet(const ntype dt, ntype &Us, ntype &U)
{
  for(int i=0; i < pars.Np; i++)
  {
    parts[i].expiLp(dt*0.5); //apply exp(iLp*dt/2)
    parts[i].expiLq(dt); //apply exp(iLq*dt)
    pbc(i);
  }
  calc_forces(Us, U);  //calculate the forces
  for(int i=0; i < pars.Np; i++)
    parts[i].expiLp(dt*0.5); //apply exp(iLp*dt/2)
}


void mdsim::PositionVerlet(const ntype dt, ntype &Us, ntype &U)
{
  for(int i=0; i < pars.Np; i++)
  {
    parts[i].expiLq(dt*0.5); //apply exp(iLq*dt/2)
    pbc(i);
  }
  calc_forces(Us, U);
  for(int i=0; i < pars.Np; i++)
  {
    parts[i].expiLp(dt); //apply exp(iLp*dtù)
    parts[i].expiLq(dt*0.5); //apply exp(iLq*dt/2)
  }
  calc_forces(Us, U);

}


void mdsim::LeapFrog(const ntype dt, ntype &Us, ntype &U)
{
  for (int i=0; i < pars.Np; i++)
  {
    parts[i].expiLp(dt);
    parts[i].expiLq(dt);
    pbc(i);
  }
  calc_forces(Us, U);
}



void mdsim::run(const int i_print_var)
  {
    calc_forces(Us, U);
    std::cout << std::setprecision(15) << "\nTotal energy=" << calcK()+Us << "\n";
    std::cout << "Initial Excess energy=" << Us << "\n";
    std::cout << "Initial Kinetic energy=" << calcK() << "\n";

    i_print = i_print_var;

    std::fstream fo;
    init_measures();


    if(pars.algotype == 0) //Velocity Verlet
    {
      for (long int t=0; t < pars.totsteps; t++)
      {
        VelocityVerlet(pars.dt,Us,U);
        Measures(t);
      }
    }else if (pars.algotype == 1) // position Verlet
    {
      for (long int t=0; t < pars.totsteps; t++)
      {
        PositionVerlet(pars.dt,Us,U);
        Measures(t);
      }
    }else if(pars.algotype == 2) //Leap-Frog
      {

      //starting conditions v(t-dt/2), r(t)
      for (int j=0; j < pars.Np; j++)
        parts[j].expiLp(-pars.dt*0.5);

      for (long int t=0; t < pars.totsteps; t++)
      {
        LeapFrog(pars.dt,Us,U);
        Measures_LF(t, pars.dt);
      }
    }else
      std::cerr << "Wrong value for algotype variable.\n"
                   "This can take the following values: 0 for Velocity Verlet, 1 for Position Verlet and 2 for Leap-Frog."<< std::endl;
  }

