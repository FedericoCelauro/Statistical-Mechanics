#ifndef _SIMCLASS_
#define _SIMCLASS_

#include "./particle.hpp"

using ntype=double;


//generic simulation class
class sim
{
  using simp = simpars;

protected:
  simp pars;  // parameters of the simulation
  std::vector<NucleiLJ> parts; //parts contains all particles in the simulation

  // calculate energy of particle i
  // if opt=1 calculate energies only for i < j
  ntype calcenergyi(int i, int opt) const;

  //calculate total energy
  ntype totenergy() const;

  //implements periodic boundary conditions
  void pbc(int i) ;

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


// molecular dynamics sim
class mdsim: public sim
{
  using ntype=double;
  using bc = sim;
  using bc::pars, bc::parts, bc::pbc;
  ntype Us, U, U_id = 0;
  int i_print=0;


  //start measurements
  void init_measures() const;

  // calculate total kinetic energy
  ntype calcK() const;

  //save measurements
  void save_measures(long int t) const;

  // calculate force acting on each particle
  void calc_forces( int i_start,  int i_end, ntype& Us, ntype& U);

  // Generate gaussian number by Box Muller algorithm
  static ntype gauss();


public:
  void prepare_initial_conf();

  void Measures(long int t);

  void VelocityVerlet(int i_start, int i_end, ntype dt, ntype &Us, ntype &U);

  void MultipleTimeStep(ntype &Us, ntype &U);

  void NormalTimeStepIntegration(ntype &Us, ntype &U);

  void run(int i);
};

#endif
