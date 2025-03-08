#include "./params.hpp"

#include <fstream>
#include "map"
#include <algorithm>
#include <regex>


//trim string custom function
inline std::string trimString(const std::string& s) {
  std::regex const pattern("^\\s+|\\s+$");
  return std::regex_replace(s, pattern, "");
}

//split custom function
std::vector<std::string> splitString(const std::string& s, const char delimiter) {
  std::vector<std::string> tokens;
  std::stringstream ss(s);
  std::string token;

  while (std::getline(ss, token, delimiter)) {
    token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
    tokens.push_back(token);
  }

  return tokens;
}


//simple constructor
simpars::simpars()
{
  simtype = 0; // 0 NTV, 1 NPT
  algotype = 0;

  nx = 8;  // number of particles along each direction
  ny = 8;
  nz = 8;

  sigma = 1.0;
  epsilon = 1.0;
  rho = 0.5;
  rc = 2.5;
  n = 0.3;

  seed = 0;
  n_iter = 1;
  mass1 = 1.0;
  mass2 = 100.0;
  Nratio = 0.01;

  adjstps = 200;
  maxadjstps = 2000;
  eqstps = 500;
  totsteps = 20000;

  save_mgl_snapshot = 1000;
  savemeasure = 20;
  outstps = 200;

  T = 2.0;
  P = 3.838; //se P*=\beta*P*v0, 1 < P* < 10 dove v0 è il volume di una particella
  deltra = 0.2;
  vmax = 10.0;
  dt = 0.002;
}

//constructor with input parameters
simpars::simpars(std::string const& filename)
{


  //map of all the parameters that can be accessed by keyword during the reading of the file
  std::map <std::string, int*> intparams = {{"nx", &nx}, {"ny", &ny}, {"nz", &nz}, {"Np", &Np}, {"simtype", &simtype}, {"seed", &seed}, {"algotype", &algotype},
                                              {"n_iter", &n_iter}};
  std::map <std::string, long int*> lintparams = { {"adjstps", &adjstps}, {"maxadjstps", &maxadjstps}, {"eqstps", &eqstps}, {"outstps", &outstps}, {"totsteps", &totsteps},
                                                  {"save_mgl_snapshot", &save_mgl_snapshot}, {"savemeasure", &savemeasure} };
  std::map <std::string, double*> doubleparams = {{"T", &T}, {"P", &P}, {"deltra", &deltra}, {"vmax", &vmax}, {"dt", &dt}, {"sigma", &sigma}, {"epsilon", &epsilon},
                                                 {"rho", &rho}, {"rc", &rc}, {"n", &n}, {"mass1", &mass1}, {"mass2", &mass2}, {"Lx", &L(0)}, {"ly", &L(1)}, {"Lz", &L(2)},
                                                    {"n", &n}, {"Nratio", &Nratio}};



  std::ifstream fin(filename);

  if(fin.is_open())
  {
    std:: string _;
    char delimiter;

    //reading what delimiter has been used in the file
    fin >> _ >> _ >> delimiter;

    //if whitespace
    if (delimiter == '0')
    {
      std::string pname, ptype, pvalue;

      //parameters reading
      while(fin >> pname >>  ptype >>  pvalue)
      {
        if (ptype == "int")
          *intparams[pname]  = stoi(pvalue);

        else if (ptype == "long")
          *lintparams[pname] = stol(pvalue);

        else if (ptype == "double")
          *doubleparams[pname] = stod(pvalue);

        else
          std::cerr << "Error: unknown type " << ptype << std::endl;
      }
      }
    else //delimiter is a character different than whitespace
    {

      std:: vector<std::string> paramvector; //0 = name, 1 = type, 2 = value
      //parameters reading
      while ( getline(fin , _))
      {
        if (_.empty()) // if _ is empty, we have an empty line
          continue;

        paramvector = splitString(_, delimiter);

        if (paramvector[1] == "int")
          *intparams[paramvector[0]]  = stoi(paramvector[2]);

        else if (paramvector[1] == "long")
          *lintparams[paramvector[0]] = stol(paramvector[2]);

        else if (paramvector[1] == "double")
          *doubleparams[paramvector[0]] = stod(paramvector[2]);

        else
          std::cerr << "Error: unknown type " << paramvector[1] << std::endl;
      }
    }
  }
  else
  {
    std::cerr << "Error: can't open file " << filename << std::endl;
    exit(1);
  }
  fin.close();
}


// << overload
std::ostream& operator<<(std::ostream& os, const simpars& params)
{
  os << "General parameters:\n\n";

  os << "nx: " << params.nx << "\n";
  os << "ny: " << params.ny << "\n";
  os << "nz: " << params.nz << "\n";
  os << "T: " << params.T << "\n";
  os << "P: " << params.P << "\n";
  os << "n " << params.n << "\n";
  os << "rho: " << params.rho << "\n";
  os << "rc: " << params.rc << "\n";
  os << "sigma: " << params.sigma << "\n";
  os << "epsilon: " << params.epsilon << "\n";
  os << "totsteps: " << params.totsteps << "\n";
  os << "save_mgl_snapshot: " << params.save_mgl_snapshot << "\n";
  os << "savemeasure: " << params.savemeasure << "\n";
  os << "outstps: " << params.outstps << "\n";
  os << "seed: " << params.seed << "\n";
  os << "n_iter: " << params.n_iter << "\n";

  os << "\nMC parameters:\n\n";

  os << "simptype: " << params.simtype << "\n";
  os << "deltra: " << params.deltra << "\n";
  os << "vmax: " << params.vmax << "\n";
  os << "adjstps: " << params.adjstps << "\n";
  os << "maxadjstps: " << params.maxadjstps << "\n";
  os << "eqstps: " << params.eqstps << "\n";

  os << "\nMD parameters:\n\n";

  os << "algotype: " << params.algotype << "\n";
  os << "mass1: " << params.mass1 << "\n";
  os << "mass2: " << params.mass2 << "\n";
  os << "Nratio: " << params.Nratio << "\n";
  os << "dt: " << params.dt << "\n";

  return os;
}

