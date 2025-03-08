#include "./sim.hpp"

int main(int argc, char **argv)
{

  mdsim md; //Lennard-Jones type particles

  //read the parameters
  const std::string paramsFile = "sim_parameters.txt";
  md.read_params(paramsFile);


  for(int i=0; i < md.get_n_iter(); i++)
  {
    md.init_rng(); // random number generator initialization
    md.prepare_initial_conf(); // creates the initial configuration
    md.run(i); // runs the MD simulation
  }
  return 0;
}
