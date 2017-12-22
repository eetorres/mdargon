# mdargon
A basic Molecular Dynamics (MD) simulation of argon using OpenGL

<img src="https://github.com/eetorres/mdargon/blob/master/screenshot.png?raw=true" alt="Main GUI"> 

Compile with:

 > make

Run the program as:

 > argon_gl

The program display an OpenGL window with the coordinates of the
argon particles and store the values of the energy:

 + Potential
 + Kinetics
 + Total

Every 10 time steps in the file: argon_energy.dat

To plot the results run the following command:

 > gnuplot -persist energy.plot


