# Solid-Hydrogen-Phonon-Dispersion
A simple program that calculates phonon dispersion of solid H2. Can be easily modified for other materials. 

Phonon_Dispersion.pdf succintly covers the theory basics and some implementation detail.

The code sets up a HCP lattice then calculates the phonon frequency on a chosen path in the Brilluoin zone (which is a hexagonal prism for HCP). The code is quite modular, so you could:

1. Change the lattice to CCP or any other lattices.
2. Choose a different path in the Brillouin zone. (Although you have to do a bit of geometry to work them out)
3. For the version 2 of the code, you could use your own interatomic potential. (Note: Unit of potential is in k_B T, position in angstrom, an example file Vpimd.dat is provided, where the second column is the centroid potential of mean force at 5.4 K and third column is the Silvera-Goldman potential)

Have fun. 
