#include "./sim_MTS.hpp"
#include <fstream>
#include <chrono>

int main(int argc, char **argv)
{

  mdsim md; //Lennard-Jones type particles

  //read the parameters
  const std::string paramsFile = "sim_parameters.txt";
  md.read_params(paramsFile);


  //open the output file
  std::fstream f;

  f.open("time_measures.txt", std::ios::out|std::ios::app);

    if (f.is_open())
    {
    f << "iteration, duration" << "\n";
    }else
    {
    std::cerr << "Error: can't open file 'time_measures.txt'" << std::endl;
    exit(1);
  }
  f.close();


  for(int i=0; i < md.get_n_iter(); i++)
  {
    md.init_rng(); // random number generator initialization
    md.prepare_initial_conf(); // creates the initial configuration

    auto start = std::chrono::high_resolution_clock::now();

    md.run(i); // runs the MC simulation
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "duration = "  << duration.count() << "ms" << std::endl;

    f.open("time_measures_1.txt", std::ios::out|std::ios::app);
    if (f.is_open())
    {
    // save duration time per iteration
    f << i << ", "  << duration.count()/1000.0 << "\n";
    }else
    {
    std::cerr << "Error: can't open file 'time_measures.txt'" << std::endl;
    exit(1);
    }
    f.close();

  }
  return 0;
}
