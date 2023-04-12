# ece4370_final_project
 
Adapted from the matlab code provided by Professor Karan Mehta, which was originally written by Carl Poitras. 

Initialize the distribution using the setup_parameters.txt. This has the parameters:
nCladding: cladding index
nCore: core (effective) index
lambd: wavelength (um)
widthDomain: width of total thing simulating (um)
lengthDomain: length of total thing simulating (um)
w_B: width of the bus waveguide (um)
s: separation between the two waveguides (um)
w: width of the output waveguide (um)
taper_start: top of the taper for the output waveguide (um)
taper_end: bottom point of the taper for the output waveguide (um)
turn_start: startpoint for where the output waveguide turns
R: radius of the curve (um)
angle: angle the output waveguide curves to (degrees)
tip_width: the width of the tip of the taper of the output waveguide (um)
sig: width of the initial amplitude distribution (um)
del_x: discritization in x (um)
del_z: discritization in z (um)
alpha: mix between the two methods of simulation

An image showing how some of these parameters map onto the structure can be seen below. 0 x is in the middle, 0 z is the top edge. 
![Image showing the index disbribution.] (./structure_image.jpg)
