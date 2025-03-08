#pragma once
#ifndef _PARAMS_
#define _PARAMS_

#include "./pvector.hpp"
#include "string"

//trim string custom function
inline std::string trimString(const std::string& s);

//split custom function
std::vector<std::string> splitString(const std::string& s, const char delimiter);


class simpars
{
  using ntype=double;

public:
  //general
  int nx, ny, nz; /* nx*ny*nz particelle */
  double T, P; // temperature and pressure
  int Np; // numero di particelle
  long int save_mgl_snapshot, savemeasure, outstps, totsteps; // Nsteps = simulations steps, outstps steps print something on current simulation status
  double rho, rc; // density
  int seed; // -1 means random
  double n;
  pvector<ntype, 3> L; // box
  double sigma, epsilon; // Lennard-Jones parameters
  int n_iter; //number of iteration of the whole code

  //Monte Carlo specific

  long int maxadjstps, eqstps, adjstps;
  int simtype; // simulation type (0 NTV, 1 NPT)
  double deltra, vmax; // parameter of MC moves

  //Molecular Dynamics specific

  double mass;
  int algotype; // integration algorithm (0 normal velocity Verlet, 1 multi time step algo)
  double Nratio; // (# particles mass1 / total) for multiple time step algo
  int M; // M = sqrt(m1/m2) is m1 < m2
  int N1; //N1 = Np*Nratio
  double dt, deltat; //integration steps

  //simple constructor
  simpars();

  //constructor with input parameters
  simpars(std::string const& filename);

  // << overload
  friend std::ostream& operator<<(std::ostream& os, const simpars& params);


};


#endif
