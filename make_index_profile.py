import numpy as np
import matplotlib.pyplot as plt 

um=1e-6
i = 1j

#Parameters of the simulation

#Free space wavelength
lambd = 1.5*um
k0 = 2*np.pi/lambd

#Indices of refraction
nCladding = 1
nCore = 1.1
nBar = (nCladding + nCore)/2

#Length of waveguide (z-direction)
inputWGLength = 100*um #400*um;

#Width of the domain (x-direction)
widthDomain = 15*um
#Length of simulation domain (z-direction)
lengthDomain = inputWGLength

#Width of bus waveguide (x-direction)
BusWidth = 0.22*um

#gap between waveguides: 
waveguide_separation = 0.5*um

#width of the output waveguide
output_width = 0.3*um

#end point of the output waveguide taper (z)
taper_end = inputWGLength/2

#start point of the taper (z)
taper_start = inputWGLength/4
#width of the tip of the taper
tip_width = .15*um


#where to start angling away (relative to the tip start)
#angle_away_start = 3*inputWGLength/4
turn_start = 3*inputWGLength/4

#angle to go at
#angle_away_angle = 10*np.pi/180
turn_angle = 10*np.pi/180

#radius of turn
radius = 100*um
#Width of initial Gaussian
sig = .5*1.2*um

#Discretization in z-direction
#deltaz = .25*um
#deltaz = .01*um
deltaz = .1*um

Nzpts = round(lengthDomain/deltaz)
#excluding z=0; 

#Discretization in x-direction
deltax = 0.01 * um
#Number of points in x-direction resulting from discretization
N = round(widthDomain / deltax)

#Create the x-axis and z-axis: 0 is at center of domain
x = np.linspace(-1*widthDomain/2, widthDomain/2, N)
z = np.linspace(0, deltaz*Nzpts, Nzpts+1)


def make_ind_profile(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end):
	# Index profile for WG 
	#n_Input_WG = nCladding*np.ones((1, N))
	n_Input_WG = nCladding*np.ones((Nzpts+1, N), dtype=complex) 


	#define the bus index stuff
	coreinds = np.where((x<=BusWidth/2) & (x>=-BusWidth/2))
	n_Input_WG[:,coreinds] = nCore



	#define the coupler index main part locations
	coreinds = np.where((x<=BusWidth/2 + waveguide_separation + output_width) & (x>=BusWidth/2 + waveguide_separation))
	post_taper = np.where(z <= taper_end)[0][-1]
	#n_Input_WG[post_taper,coreinds] = nCore
	n_Input_WG[post_taper:,coreinds] = nCore

	#now the taper
	taper_init_ind = np.where(z >= taper_start)[0][0]
	print(taper_init_ind)
	z_taper_inds = np.where((z >= taper_start) & (z <= taper_end))
	print(z_taper_inds[0])
	for z_ind in z_taper_inds[0]:
		taper_ind = np.where((x >= BusWidth/2 + waveguide_separation) & (x <= BusWidth/2 + waveguide_separation + (z[z_ind] -taper_start)*output_width/(taper_end - taper_start)))
		n_Input_WG[z_ind, taper_ind] = nCore 

	plt.imshow(np.real(n_Input_WG))
	plt.xlabel('x (\mum)')
	plt.ylabel('Refractive index')
	plt.show()
	return(n_Input_WG)


def make_ind_profile_angle_away(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, angle_away_start, angle_away_angle):
	#TODO: FIGURE OUT WHAT ANGLE IS SMALL ENOUGH FOR ADIABATIC
	# Index profile for WG 
	#n_Input_WG = nCladding*np.ones((1, N))
	n_Input_WG = nCladding*np.ones((Nzpts+1, N), dtype=complex) 
	border_line_bounds = []

	#define the bus index stuff
	coreinds = np.where((x<=BusWidth/2) & (x>=-BusWidth/2))
	n_Input_WG[:,coreinds] = nCore
	border_line_bounds.append([[BusWidth/2, BusWidth/2],[z[0], z[-1]]])
	border_line_bounds.append([[-1*BusWidth/2, -1*BusWidth/2],[z[0], z[-1]]])



	#define the coupler index main part locations
	coreinds = np.where((x<=BusWidth/2 + waveguide_separation + output_width) & (x>=BusWidth/2 + waveguide_separation))
	post_taper = np.where((z >= taper_end) & (z <= angle_away_start))[0]
	n_Input_WG[post_taper[0]:post_taper[-1]+1,coreinds] = nCore

	border_line_bounds.append([[BusWidth/2 + waveguide_separation + output_width, BusWidth/2 + waveguide_separation + output_width],[taper_end, angle_away_start]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation],[taper_start, angle_away_start]])

	#n_Input_WG[post_taper:,coreinds] = nCore

	#now the taper
	taper_init_ind = np.where(z >= taper_start)[0][0]
	print(taper_init_ind)
	z_taper_inds = np.where((z >= taper_start) & (z <= taper_end))
	print(z_taper_inds[0])
	for z_ind in z_taper_inds[0]:
		taper_ind = np.where((x >= BusWidth/2 + waveguide_separation) & (x <= BusWidth/2 + waveguide_separation + (z[z_ind] -taper_start)*output_width/(taper_end - taper_start)))
		n_Input_WG[z_ind, taper_ind] = nCore 
	border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation + (taper_end -taper_start)*output_width/(taper_end - taper_start)],[taper_start, taper_end]])


	#now angled stuff
	z_angled_inds = np.where(z>=angle_away_start)
	for z_ind in z_angled_inds[0]:
		ang_ind = np.where((x >= BusWidth/2 + waveguide_separation +np.tan(angle_away_angle)*(z[z_ind] - angle_away_start)) & (x <= BusWidth/2 + waveguide_separation +np.tan(angle_away_angle)*(z[z_ind]-angle_away_start) + output_width))
		#ang_ind = np.where(x >= BusWidth/2 + waveguide_separation +np.tan(angle_away_angle)*(z[z_ind] - angle_away_start))
		#ang_ind = np.where(x <= BusWidth/2 + waveguide_separation +np.tan(angle_away_angle)*(z[z_ind]-angle_away_start) + output_width)
		n_Input_WG[z_ind, ang_ind] = nCore 
		#print(ang_ind)

	border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation +np.tan(angle_away_angle)*(z[-1] - angle_away_start)],[angle_away_start, z[-1]]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation + output_width, BusWidth/2 + output_width + waveguide_separation +np.tan(angle_away_angle)*(z[-1] - angle_away_start)],[angle_away_start, z[-1]]])

	#plt.imshow(np.real(n_Input_WG))
	plt.imshow(np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)

	plt.xlabel('$x (\mu m)$')
	plt.ylabel('$z (\mu m)$')
	#for i in range(len(border_line_bounds)):
	#	plt.plot(border_line_bounds[i][0],border_line_bounds[i][1], c='k')
	plt.show()
	return n_Input_WG, border_line_bounds




def make_ind_profile_angle_away_with_radius(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius):
	#TODO: FIGURE OUT WHAT ANGLE IS SMALL ENOUGH FOR ADIABATIC
	# Index profile for WG 
	#radius will be the radius of the inner edge
	#n_Input_WG = nCladding*np.ones((1, N))
	n_Input_WG = nCladding*np.ones((Nzpts+1, N), dtype=complex) 
	border_line_bounds = []

	#define the bus index stuff
	coreinds = np.where((x<=BusWidth/2) & (x>=-BusWidth/2))
	n_Input_WG[:,coreinds] = nCore
	border_line_bounds.append([[BusWidth/2, BusWidth/2],[z[0], z[-1]]])
	border_line_bounds.append([[-1*BusWidth/2, -1*BusWidth/2],[z[0], z[-1]]])


	
	#define the coupler index main part locations
	coreinds = np.where((x<=BusWidth/2 + waveguide_separation + output_width) & (x>=BusWidth/2 + waveguide_separation))
	post_taper = np.where((z >= taper_end) & (z <= turn_start))[0]
	n_Input_WG[post_taper[0]:post_taper[-1]+1,coreinds] = nCore

	border_line_bounds.append([[BusWidth/2 + waveguide_separation + output_width, BusWidth/2 + waveguide_separation + output_width],[taper_end, turn_start]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation],[taper_start, turn_start]])
	
	#n_Input_WG[post_taper:,coreinds] = nCore


	#now the taper
	taper_init_ind = np.where(z >= taper_start)[0][0]
	print(taper_init_ind)
	z_taper_inds = np.where((z >= taper_start) & (z <= taper_end))
	print(z_taper_inds[0])
	for z_ind in z_taper_inds[0]:
		taper_ind = np.where((x >= BusWidth/2 + waveguide_separation) & (x <= BusWidth/2 + waveguide_separation + (z[z_ind] -taper_start)*output_width/(taper_end - taper_start)))
		n_Input_WG[z_ind, taper_ind] = nCore 
	border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation + (taper_end -taper_start)*output_width/(taper_end - taper_start)],[taper_start, taper_end]])


	#the curve
	z_curve_inds = np.where((z >= turn_start) & (z <= turn_start + (radius + output_width)*np.sin(turn_angle)))
	#print(z_curve_inds)
	for z_ind in z_curve_inds[0]:
		#x_part_ind = np.where((x >= BusWidth/2 + waveguide_separation + np.sqrt(radius**2 - z[z_ind]**2)) & (x <= BusWidth/2 + waveguide_separation + np.sqrt((radius + output_width)**2 - z[z_ind]**2)))
		x_part_ind = np.where((x >= BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_ind]-turn_start)**2)) & (x <= BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_ind]-turn_start)**2)))
		#print(x_part_ind)
		n_Input_WG[z_ind, x_part_ind] = nCore 
	
	border_line_bounds.append([BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_curve_inds]-turn_start)**2), z[z_curve_inds]])
	border_line_bounds.append([BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_curve_inds]-turn_start)**2), z[z_curve_inds]])
	plt.plot(BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_curve_inds]-turn_start)**2),z[z_curve_inds])
	#plt.plot(BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_curve_inds]-turn_start)**2),z[z_curve_inds])
	#plt.show()
	#will be some slight inaccuracy sicne the interface has a little hang over that shouldnt exist, but since large radius shouldnt really matter


	
	#now angled stuff after the curve
	angle_away_start = turn_start + (radius + output_width)*np.sin(turn_angle)
	z_angled_inds = np.where(z>=turn_start + (radius + output_width)*np.sin(turn_angle))
	x_start = BusWidth/2 + waveguide_separation + radius + output_width - np.cos(turn_angle)*(radius + output_width)
	for z_ind in z_angled_inds[0]:
		ang_ind = np.where((x >= x_start + np.tan(turn_angle)*(z[z_ind] - angle_away_start)) & (x <= x_start +np.tan(turn_angle)*(z[z_ind]-angle_away_start) + output_width))
		#ang_ind = np.where(x >= x_start + np.tan(turn_angle)*(z[z_ind] - angle_away_start))
		#ang_ind = np.where(x <= x_start +np.tan(turn_angle)*(z[z_ind]-angle_away_start) + output_width)

		n_Input_WG[z_ind, ang_ind] = nCore 
		#print(ang_ind)

	border_line_bounds.append([[x_start + np.tan(turn_angle)*(z[z_curve_inds[0][-1]] - angle_away_start), x_start + np.tan(turn_angle)*(z[-1] - angle_away_start)],[angle_away_start, z[-1]]])
	border_line_bounds.append([[x_start +np.tan(turn_angle)*(z[z_curve_inds[0][-1]]-angle_away_start) + output_width, x_start +np.tan(turn_angle)*(z[-1]-angle_away_start) + output_width],[angle_away_start, z[-1]]])
	


	#plt.imshow(np.real(n_Input_WG))
	plt.imshow(np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)
	#plt.imshow(np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)))
	

	plt.xlabel('$x (m)$')
	plt.ylabel('$z (m)$')
	for i in range(len(border_line_bounds)):
		plt.plot(border_line_bounds[i][0],border_line_bounds[i][1], c='k')
	plt.show()
	return n_Input_WG, border_line_bounds

def make_ind_profile_angle_away_with_radius_and_finite_tip(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius, tip_width):
	#TODO: FIGURE OUT WHAT ANGLE IS SMALL ENOUGH FOR ADIABATIC
	# Index profile for WG 
	#radius will be the radius of the inner edge
	#n_Input_WG = nCladding*np.ones((1, N))
	n_Input_WG = nCladding*np.ones((Nzpts+1, N), dtype=complex) 
	border_line_bounds = []

	#define the bus index stuff
	coreinds = np.where((x<=BusWidth/2) & (x>=-BusWidth/2))
	n_Input_WG[:,coreinds] = nCore
	border_line_bounds.append([[BusWidth/2, BusWidth/2],[z[0], z[-1]]])
	border_line_bounds.append([[-1*BusWidth/2, -1*BusWidth/2],[z[0], z[-1]]])


	
	#define the coupler index main part locations
	coreinds = np.where((x<=BusWidth/2 + waveguide_separation + output_width) & (x>=BusWidth/2 + waveguide_separation))
	post_taper = np.where((z >= taper_end) & (z <= turn_start))[0]
	n_Input_WG[post_taper[0]:post_taper[-1]+1,coreinds] = nCore

	border_line_bounds.append([[BusWidth/2 + waveguide_separation + output_width, BusWidth/2 + waveguide_separation + output_width],[taper_end, turn_start]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation],[taper_start, turn_start]])
	
	#n_Input_WG[post_taper:,coreinds] = nCore


	#now the taper
	#print(taper_init_ind)
	z_taper_inds = np.where((z >= taper_start) & (z <= taper_end))
	#print(z_taper_inds[0])
	for z_ind in z_taper_inds[0]:
		#taper_ind = np.where((x >= BusWidth/2 + waveguide_separation) & (x <= BusWidth/2 + waveguide_separation + (z[z_ind] -taper_start)*output_width/(taper_end - taper_start)))
		taper_ind = np.where((x >= BusWidth/2 + waveguide_separation + tip_width) & (x <= BusWidth/2 + waveguide_separation + tip_width + (z[z_ind] -taper_start)*(output_width - tip_width)/(taper_end - taper_start)))

		n_Input_WG[z_ind, taper_ind] = nCore 
	coreinds = np.where((x<=BusWidth/2 + waveguide_separation + tip_width) & (x>=BusWidth/2 + waveguide_separation))
	n_Input_WG[z_taper_inds[0][0]:z_taper_inds[0][-1],coreinds] = nCore

	#border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation + (taper_end -taper_start)*output_width/(taper_end - taper_start)],[taper_start, taper_end]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation + tip_width, BusWidth/2 + waveguide_separation + tip_width + (taper_end -taper_start)*(output_width - tip_width)/(taper_end - taper_start)],[taper_start, taper_end]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation + tip_width, BusWidth/2 + waveguide_separation], [taper_start, taper_start]])

	#the curve
	z_curve_inds = np.where((z >= turn_start) & (z <= turn_start + (radius + output_width)*np.sin(turn_angle)))
	#print(z_curve_inds)
	for z_ind in z_curve_inds[0]:
		#x_part_ind = np.where((x >= BusWidth/2 + waveguide_separation + np.sqrt(radius**2 - z[z_ind]**2)) & (x <= BusWidth/2 + waveguide_separation + np.sqrt((radius + output_width)**2 - z[z_ind]**2)))
		x_part_ind = np.where((x >= BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_ind]-turn_start)**2)) & (x <= BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_ind]-turn_start)**2)))
		#print(x_part_ind)
		n_Input_WG[z_ind, x_part_ind] = nCore 
	
	border_line_bounds.append([BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_curve_inds]-turn_start)**2), z[z_curve_inds]])
	border_line_bounds.append([BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_curve_inds]-turn_start)**2), z[z_curve_inds]])
	plt.plot(BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_curve_inds]-turn_start)**2),z[z_curve_inds])
	#plt.plot(BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_curve_inds]-turn_start)**2),z[z_curve_inds])
	#plt.show()
	#will be some slight inaccuracy sicne the interface has a little hang over that shouldnt exist, but since large radius shouldnt really matter


	
	#now angled stuff after the curve
	angle_away_start = turn_start + (radius + output_width)*np.sin(turn_angle)
	z_angled_inds = np.where(z>=turn_start + (radius + output_width)*np.sin(turn_angle))
	x_start = BusWidth/2 + waveguide_separation + radius + output_width - np.cos(turn_angle)*(radius + output_width)
	for z_ind in z_angled_inds[0]:
		ang_ind = np.where((x >= x_start + np.tan(turn_angle)*(z[z_ind] - angle_away_start)) & (x <= x_start +np.tan(turn_angle)*(z[z_ind]-angle_away_start) + output_width))
		#ang_ind = np.where(x >= x_start + np.tan(turn_angle)*(z[z_ind] - angle_away_start))
		#ang_ind = np.where(x <= x_start +np.tan(turn_angle)*(z[z_ind]-angle_away_start) + output_width)

		n_Input_WG[z_ind, ang_ind] = nCore 
		#print(ang_ind)

	border_line_bounds.append([[x_start + np.tan(turn_angle)*(z[z_curve_inds[0][-1]] - angle_away_start), x_start + np.tan(turn_angle)*(z[-1] - angle_away_start)],[angle_away_start, z[-1]]])
	border_line_bounds.append([[x_start +np.tan(turn_angle)*(z[z_curve_inds[0][-1]]-angle_away_start) + output_width, x_start +np.tan(turn_angle)*(z[-1]-angle_away_start) + output_width],[angle_away_start, z[-1]]])
	


	#plt.imshow(np.real(n_Input_WG))
	plt.imshow(np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)
	#plt.imshow(np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)))
	

	plt.xlabel('$x (m)$')
	plt.ylabel('$z (m)$')
	for i in range(len(border_line_bounds)):
		plt.plot(border_line_bounds[i][0],border_line_bounds[i][1], c='k')
	plt.show()
	return n_Input_WG, border_line_bounds

def make_ind_profile_with_finite_tip_width(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, tip_width):
	# Index profile for WG 
	#n_Input_WG = nCladding*np.ones((1, N))
	n_Input_WG = nCladding*np.ones((Nzpts+1, N), dtype=complex) 


	#define the bus index stuff
	coreinds = np.where((x<=BusWidth/2) & (x>=-BusWidth/2))
	n_Input_WG[:,coreinds] = nCore



	#define the coupler index main part locations
	coreinds = np.where((x<=BusWidth/2 + waveguide_separation + output_width) & (x>=BusWidth/2 + waveguide_separation))
	post_taper = np.where(z <= taper_end)[0][-1]
	#n_Input_WG[post_taper,coreinds] = nCore
	n_Input_WG[post_taper:,coreinds] = nCore

	#now the taper
	taper_init_ind = np.where(z >= taper_start)[0][0]
	#print(taper_init_ind)
	z_taper_inds = np.where((z >= taper_start) & (z <= taper_end))
	#print(z_taper_inds[0])
	for z_ind in z_taper_inds[0]:
		taper_ind = np.where((x >= BusWidth/2 + waveguide_separation + tip_width) & (x <= BusWidth/2 + waveguide_separation + tip_width + (z[z_ind] -taper_start)*(output_width-tip_width)/(taper_end - taper_start)))
		n_Input_WG[z_ind, taper_ind] = nCore 
	#the wide part of the tip
	coreinds = np.where((x<=BusWidth/2 + waveguide_separation + tip_width) & (x>=BusWidth/2 + waveguide_separation))
	n_Input_WG[taper_init_ind:,coreinds] = nCore

	plt.imshow(np.real(n_Input_WG), extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)
	plt.xlabel('$x (\mu m)$')
	#plt.ylabel('Refractive index')
	plt.ylabel('$z (\mu m)$')
	plt.show()
	return(n_Input_WG)




#make_ind_profile(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end)
#make_ind_profile_with_finite_tip_width(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, tip_width)
#make_ind_profile_angle_away(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, angle_away_start, angle_away_angle)
#make_ind_profile_angle_away_with_radius(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius)
make_ind_profile_angle_away_with_radius_and_finite_tip(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius, tip_width)






