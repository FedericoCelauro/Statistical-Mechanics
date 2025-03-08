#ifndef _SIMCLASS_
#define _SIMCLASS_

#include <vector>
#include "./params.hpp"
#include "./particle.hpp"


using ntype=double;


//generic simulation class
class sim
{
  using simp = simpars;


protected:
  simp pars;  // parameters of the simulation
  std::vector<particleLJ> parts; //parts contains all particles in the simulation

  // calculate energy of particle i
  // if opt=1 calculate energies only for i < j
  ntype calcenergyi(int i, int opt=0) const;

  //calculate total energy
  ntype totenergy() const;

  //implements periodic boundary conditions
  void pbc(int i);

  // saves all particles positions
  void save_mgl_snapshot(long int t);

public:

  //read parameters from file
  void read_params(const std::string& filename);

  //get n_iter
  int get_n_iter() const;

  // prepares initial configuration
  void prepare_initial_conf();

  //initializing random number generator
  void init_rng() const;

  static void run();
};



//MC simulation
class mcsim: public sim
{

  using bc=sim;
  using bc::calcenergyi, bc::parts, bc::pars,
        bc::totenergy, bc::save_mgl_snapshot, bc::pbc;


  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  // counters used to calculate acceptance rates
  long int tot_tra, tot_vol, tra_rej, vol_rej;

  //  trial move
  void alpha(int i);

  // acceptance move
  void acc(int i, ntype eno);

  // NTV move
  void move_NTV(int i);

  // adjust the deltra parameter looking at the acceptance rates
  void calc_acceptance_and_adjust();

  //start measurements
  void init_measures() const;

  //save measurements
  void save_measures(long int t);

  // restore all particle positions
  void restore_all_pars();

  // store all particle position
  void store_all_pars();

  // trial box move
  void alpha_box(ntype& DG, ntype& fact);

  // accept or reject box move
  void acc_box(ntype DG, ntype fact);

  void move_box();

 public:

  //run the simulation
  void run();
};


// MD simulation
class mdsim: public sim
{
  using ntype=double;
  using bc = sim;
  using bc::pars, bc::parts, bc::pbc;
  ntype Us, U;
  int i_print=0;

  //start measurements
  void init_measures() const;

  // calculate total kinetic energy
  ntype calcK() const;

  //save measurements
  void save_measures(long int t) const;

  // calculate force acting on each particle
  void calc_forces(ntype& Us, ntype& U);

  // Generate gaussian number by Box Muller algorithm
  static ntype gauss();

public:
  void prepare_initial_conf();

  void Measures(long int t);

  void Measures_LF(long int t, ntype dt);

  void VelocityVerlet(ntype dt, ntype &Us, ntype &U);

  void PositionVerlet(ntype dt, ntype &Us, ntype &U);

  void LeapFrog(ntype dt, ntype &Us, ntype &U);

  void run(int i_print_var);
};
#endif
