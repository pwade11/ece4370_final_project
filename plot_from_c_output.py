import numpy as np
import matplotlib.pyplot as plt 

#plot from the spec file and the c data
#file = "project_sim_output_on.csv"
#file = "project_sim_output_off.csv"
file = "project_sim_output.csv"




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




def make_border(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius, tip_width):
	border_line_bounds = []

	#define the bus index stuff
	coreinds = np.where((x<=BusWidth/2) & (x>=-BusWidth/2))
	border_line_bounds.append([[BusWidth/2, BusWidth/2],[z[0], z[-1]]])
	border_line_bounds.append([[-1*BusWidth/2, -1*BusWidth/2],[z[0], z[-1]]])


	
	#define the coupler index main part locations

	border_line_bounds.append([[BusWidth/2 + waveguide_separation + output_width, BusWidth/2 + waveguide_separation + output_width],[taper_end, turn_start]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation],[taper_start, turn_start]])
	


	#now the taper
	#border_line_bounds.append([[BusWidth/2 + waveguide_separation, BusWidth/2 + waveguide_separation + (taper_end -taper_start)*output_width/(taper_end - taper_start)],[taper_start, taper_end]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation + tip_width, BusWidth/2 + waveguide_separation + tip_width + (taper_end -taper_start)*(output_width - tip_width)/(taper_end - taper_start)],[taper_start, taper_end]])
	border_line_bounds.append([[BusWidth/2 + waveguide_separation + tip_width, BusWidth/2 + waveguide_separation], [taper_start, taper_start]])

	#the curve
	z_curve_inds = np.where((z >= turn_start) & (z <= turn_start + (radius + output_width)*np.sin(turn_angle)))

	
	border_line_bounds.append([BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_curve_inds]-turn_start)**2), z[z_curve_inds]])
	border_line_bounds.append([BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_curve_inds]-turn_start)**2), z[z_curve_inds]])
	#plt.plot(BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt(radius**2 - (z[z_curve_inds]-turn_start)**2),z[z_curve_inds])
	#plt.plot(BusWidth/2 + waveguide_separation + radius + output_width - np.sqrt((radius+output_width)**2 - (z[z_curve_inds]-turn_start)**2),z[z_curve_inds])
	#plt.show()
	#will be some slight inaccuracy sicne the interface has a little hang over that shouldnt exist, but since large radius shouldnt really matter


	
	#now angled stuff after the curve
	angle_away_start = turn_start + (radius + output_width)*np.sin(turn_angle)
	z_angled_inds = np.where(z>=turn_start + (radius + output_width)*np.sin(turn_angle))
	x_start = BusWidth/2 + waveguide_separation + radius + output_width - np.cos(turn_angle)*(radius + output_width)

	border_line_bounds.append([[x_start + np.tan(turn_angle)*(z[z_curve_inds[0][-1]] - angle_away_start), x_start + np.tan(turn_angle)*(z[-1] - angle_away_start)],[angle_away_start, z[-1]]])
	border_line_bounds.append([[x_start +np.tan(turn_angle)*(z[z_curve_inds[0][-1]]-angle_away_start) + output_width, x_start +np.tan(turn_angle)*(z[-1]-angle_away_start) + output_width],[angle_away_start, z[-1]]])

	#for i in range(len(border_line_bounds)):
	#	plt.plot(border_line_bounds[i][0],border_line_bounds[i][1], c='k')
	#plt.show()
	return border_line_bounds


Nzpts = round(lengthDomain/deltaz)
#excluding z=0; 


#Number of points in x-direction resulting from discretization
N = round(widthDomain / deltax)

#Create the x-axis and z-axis: 0 is at center of domain
x = np.linspace(-1*widthDomain/2, widthDomain/2, N, dtype=complex)
z = np.linspace(0, deltaz*Nzpts, Nzpts+1, dtype=complex)


border_line_bounds = make_border(N, nCladding, nCore, x, z, BusWidth, waveguide_separation, output_width, taper_start, taper_end, turn_start, turn_angle, radius, tip_width)


data = []
#with open("project_sim_output.csv",'r') as f:
with open(file,'r') as f:
	for line in f:
		data.append(np.array(line.split(',')).astype(np.float))
data = np.array(data)
plt.imshow(data, aspect=widthDomain/lengthDomain,extent=(np.min(x)/1e-6, np.max(x)/1e-6, np.max(z)/1e-6,np.min(z)/1e-6))
plt.colorbar()
for i in range(len(border_line_bounds)):
	plt.plot(np.array(border_line_bounds[i][0])/1e-6,np.array(border_line_bounds[i][1])/1e-6, c='k')
plt.title("Intensity distribution")
plt.ylabel("z ($\mu m$)")
plt.xlabel("x ($\mu m$)")

plt.show()


start_looking_ind = 50 #how far down from the start to begin analysis as our starting point

plt.imshow(data[start_looking_ind:], aspect=widthDomain/lengthDomain,extent=(np.min(x), np.max(x), np.max(z[start_looking_ind:]),np.min(z[start_looking_ind:])))
plt.colorbar()
for i in range(len(border_line_bounds)):
	plt.plot(border_line_bounds[i][0],border_line_bounds[i][1], c='k')
plt.show()

#look at the power at the input (down a bit so have time to settle in) vs output
bus_coreinds = np.where((x<=BusWidth/2) & (x>=-BusWidth/2))

angle_away_start = turn_start + (radius + output_width)*np.sin(turn_angle)
x_start = BusWidth/2 + waveguide_separation + radius + output_width - np.cos(turn_angle)*(radius + output_width)
ang_ind = np.where((x >= x_start + np.tan(turn_angle)*(z[-1] - angle_away_start)) & (x <= x_start +np.tan(turn_angle)*(z[-1]-angle_away_start) + output_width))
#this assumes the end is on the z axis, which may not always be the case, so need to look

max_bus_start = np.max(data[start_looking_ind,bus_coreinds])
print(max_bus_start)

max_output_end = np.max(data[-1, ang_ind])
print(max_output_end)

max_bus_end = np.max(data[-1,bus_coreinds])
print(max_bus_end)


print("power transfer (dB) = " + str(10*np.log10(max_output_end/max_bus_start)))

print("power out (dB) = " + str(10*np.log10(max_bus_end/max_bus_start)))





