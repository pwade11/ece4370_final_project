import numpy as np
import matplotlib.pyplot as plt 

#adding Z varying n
#Parameters of the simulation

#Free space wavelength
#lambd = 1.5*um


#Indices of refraction

#nCladding = 1
#nCore = 1.1
#NOT SURE WHAT TO USE SINCE THESE ARE NOT THE Neff WE WOULD NEED, BUT PUT THESE HERE FOR NOW, ASK PROF
#nCladding = 1.4446 #silicon dioxide at 1500 nm https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson
#nCore =  1.747	#aluminum oxide at 1500 nm https://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-o

#Length of waveguide (z-direction)
#inputWGLength = 100*um #400*um;

#Width of the domain (x-direction)
#widthDomain = 15*um
#Length of simulation domain (z-direction)
#lengthDomain = inputWGLength

#Width of bus waveguide (x-direction)
#BusWidth = 0.22*um
#BusWidth = 1.5*um
#BusWidth = 0.6*um

#gap between waveguides: 
#waveguide_separation = 0.2*um
#waveguide_separation = 0.125*um
#waveguide_separation = 0.5*um
#waveguide_separation = 1*um
#waveguide_separation = 1.5*um	#this so far the best for high transfer, but probably depends on the length set up
#waveguide_separation = 4.5*um	#this has like 0 couplign at all 


#width of the output waveguide
#output_width = 0.3*um
#output_width = 1.9*um
#output_width = 1*um

#end point of the output waveguide taper (z)
#taper_end = inputWGLength/2
#taper_end =0 
#start point of the taper (z)
#taper_start = inputWGLength/4
#taper_start =0 

#Width of initial Gaussian
#sig = .5*1.2*um
#sig = 0.4*1.2*um
#sig = 1.2*BusWidth/3 #keep it same ratio of bus width as the example


#Discretization in z-direction
#deltaz = .25*um
#need to decrease for the angled one to make sense
#deltaz = .1*um

#Discretization in x-direction
#deltax = 0.01 * um
#deltax = 0.001 * um


#where to start angling away (relative to the tip start)
#angle_away_start = 3*inputWGLength/4
#angle to go at
#angle_away_angle = 10*np.pi/180

#radius = 100*um
#turn_angle = 10*np.pi/180
#turn_start = 3*inputWGLength/4

#tip_width = .15*um


#load parameters from the file
with open("setup_parameters.txt", "r") as f:
	vals_ind = 0
	for line in f:
		tmp_val = np.float(line.split(' ')[-1])
		if(vals_ind == 0):
			nCladding = tmp_val
		elif(vals_ind == 1):
			nCore = tmp_val
		elif(vals_ind == 2):
			lambd = tmp_val*1e-6
		elif(vals_ind == 3):
			widthDomain = tmp_val*1e-6;
		elif(vals_ind == 4):
			lengthDomain = tmp_val*1e-6
		elif(vals_ind == 5):
			BusWidth = tmp_val*1e-6
		elif(vals_ind == 6):
			waveguide_separation = tmp_val*1e-6
		elif(vals_ind == 7):
			output_width = tmp_val*1e-6
		elif(vals_ind == 8):
			taper_start = tmp_val*1e-6
		elif(vals_ind == 9):
			taper_end = tmp_val*1e-6
		elif(vals_ind == 10):
			turn_start = tmp_val*1e-6
		elif(vals_ind == 11):
			radius = tmp_val*1e-6
		elif(vals_ind == 12):
			turn_angle = tmp_val*np.pi/180
		elif(vals_ind == 13):
			tip_width = tmp_val*1e-6
		elif(vals_ind == 14):
			sig = tmp_val*1e-6
		elif(vals_ind == 15):
			deltax = tmp_val*1e-6
		elif(vals_ind == 16):
			deltaz = tmp_val*1e-6
		elif(vals_ind == 17):
			alpha = tmp_val
		vals_ind += 1


k0 = 2*np.pi/lambd
nBar = (nCladding + nCore)/2

def look_at_effective_index(nCore, nCladding, k0, t = 500e-9, plot=False):
	#from vertical stack
	#cladding effective = nCladding since it is fully enclosed in same mat
	#core effective:
	beta = np.linspace(k0*nCladding, k0*nCore*.99, num=100000)
	h = np.sqrt((k0*nCore)**2 - beta**2)
	gamma = np.sqrt(beta**2 - (k0*nCladding)**2)
	#this is for TE
	lhs = np.tan(h*t)
	rhs_te = h*2*gamma/(h**2 -gamma**2)
	#for tm:
	rhs_tm = (h/(nCore**2))*(2*gamma/(nCladding**2))/((h/(nCore**2))**2  - (gamma/(nCladding**2))**2)
	if(plot):
		plt.scatter(beta/k0, lhs)
		plt.scatter(beta/k0, rhs_te)
		plt.scatter(beta/k0, rhs_tm)
		plt.ylim(-40,40)
		plt.show()

		plt.scatter(beta/k0, np.abs(rhs_te-lhs))
		plt.scatter(beta/k0, np.abs(rhs_tm-lhs))
		plt.ylim(-1, 10)
		print(beta[np.argmin(np.abs(rhs_te-lhs))]/k0)
		print(beta[np.argmin(np.abs(rhs_tm-lhs))]/k0)
		plt.show()
	
	n_eff_TE = beta[np.argmin(np.abs(rhs_te-lhs))]/k0
	n_eff_TM = beta[np.argmin(np.abs(rhs_tm-lhs))]/k0
	if(np.min(np.abs(rhs_te-lhs)) > 0.005):
		n_eff_TE = nCore #no confined waves
	if(np.min(np.abs(rhs_tm-lhs)) > 0.005):
		n_eff_TM = nCore #no confined waves
	return n_eff_TE, n_eff_TM

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



Nzpts = round(lengthDomain/deltaz)
#excluding z=0; 


#Number of points in x-direction resulting from discretization
N = round(widthDomain / deltax)

#Create the x-axis and z-axis: 0 is at center of domain
x = np.linspace(-1*widthDomain/2, widthDomain/2, N, dtype=complex)
z = np.linspace(0, deltaz*Nzpts, Nzpts+1, dtype=complex)

# Index profile for WG 
#n_Input_WG = nCladding*np.ones((1, N))
#n_Input_WG = nCladding*np.ones((Nzpts+1, N), dtype=complex) 

#print(n_Input_WG.shape)
#coreinds = np.nonzero(x) 
#print(x<=(3*inputWGWidth/2 + coupGap))
#print(x>=(inputWGWidth/2 + coupGap))
#coreinds = np.where((x<=inputWGWidth/2) & (x>=-inputWGWidth/2))
#print(np.min(coreinds))
#print(np.max(coreinds))
#find returns nonzero element indices; the boolean expressions 
#identifies values of x for which we are in the waveguide core. 

#n_Input_WG[:,coreinds] = nCore
#if(SecondWG):
#	core2inds = np.where((x<=(3*inputWGWidth/2 + coupGap)) & (x>=(inputWGWidth/2 + coupGap)))
#	n_Input_WG[:, core2inds] = nCore


n_eff_TE, n_eff_TM = look_at_effective_index(nCore, nCladding, k0)
#look_at_possible_confined_modes(k0, n_eff_TE, nCladding, BusWidth)
#look_at_possible_confined_modes(k0, n_eff_TM, nCladding, BusWidth)
print(n_eff_TE)
nCore = n_eff_TE

#nCore = n_eff_TE

#n_Input_WG = make_ind_profile(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end)
#n_Input_WG, border_line_bounds = make_ind_profile_angle_away(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, angle_away_start, angle_away_angle)
#n_Input_WG, border_line_bounds = make_ind_profile_angle_away_with_radius(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius)
n_Input_WG, border_line_bounds = make_ind_profile_angle_away_with_radius_and_finite_tip(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius, tip_width)



#compare with my c code output
n_data = []
with open("project_sim_n.csv",'r') as f:
	for line in f:
		n_data.append(np.array(line.split(',')).astype(np.float))
n_data = np.array(n_data)

plt.imshow(n_data,extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)
plt.show()
#plt.imshow(np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)

plt.imshow(n_data - np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)
plt.colorbar()
plt.show()
plt.imshow(np.abs(n_data - np.real(n_Input_WG))/np.real(n_Input_WG),extent=(np.min(x), np.max(x), np.max(z),np.min(z)),aspect=widthDomain/lengthDomain)
plt.colorbar()
plt.show()






