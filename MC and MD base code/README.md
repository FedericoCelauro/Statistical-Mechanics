# Code overview

This is the base code to simulate systems of interacting particles (for now only Lennard-Jones and soft spheres interactions) using either Monte Carlo or Molecular Dynamics approach. 

**To run the code** 
- modifying the sim.hpp file choose the type of interaction (standard is Lennard-Jones)
- modifying the main.cpp file choose to do either a MC or a MD simulation
- use the sim_parameters.txt file to modify the parameters.
  - if MC, simtype let you choose the ensemble (0 = NTV, 1 = NPT)
  - if MD, algotype let you choose the integration algorithm (0 = Velocity Verlet, 1 = Position Verlet, 2 = Leap Frog)
