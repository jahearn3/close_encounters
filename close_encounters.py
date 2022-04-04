# Close Encounters 
#
# Author: Joseph A'Hearn
# Created 10/24/2017
#
# This program calculates and plots
#    various quantities over the course of 
#    close encounters
#

import bodies as bd 
import data_loader as dl 
import constants_of_mercury6 as cm6 
import useful as uf 
import orbital_motion as om 
import plot_assistant as pa 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import xyz2geo as x2g

# STRUCTURE
# Functions to read data from a file
# Functions to sort out close encounters
# Functions to compute
# Functions to make plots
#    Close encounter plots
#    Energy plots
# Functions to write data to a file



# FUNCTIONS TO READ DATA FROM A FILE --------------------------------------------------------------------

def get_hill_and_body_radii(M_Sat, bodies, masses):
	print("Calculating Hill and body radii...")
	hill_radii = np.zeros(len(bodies))
	body_radii = np.zeros(len(bodies))
	for i in range(len(bodies)):
		a = dl.geo1mean(bodies[i])
		hill_radii[i] = ((a * (masses[i] / (3 * M_Sat))**(1/3)) * 1.0E-03) # converting to km
		print('   hill radius of ' + bodies[i] + ': ' + str(np.round(hill_radii[i], decimals=2)) + ' km')
		body_radii[i] = (((3 * masses[i] / (4 * np.pi * 500))**(1/3)) * 1.0E-03)
		print('   body radius of ' + bodies[i] + ': ' + str(np.round(body_radii[i], decimals=2)) + ' km')
	return hill_radii, body_radii

def get_xyz_data(bodies, data_points):
	# assemble a matrix with xy data from all bodies
	print("Importing position data...")
	x = np.zeros((len(bodies), data_points))
	y = np.zeros((len(bodies), data_points))
	z = np.zeros((len(bodies), data_points))
	for i in range(len(bodies)):
		x_i, y_i, z_i = dl.xyz_data(bodies[i])
		x[i] = x_i
		y[i] = y_i	
		z[i] = z_i
	return x, y, z

def get_uvw_data(bodies, data_points):
	# assemble a matrix with xy data from all bodies
	print("Importing velocity data...")
	u = np.zeros((len(bodies), data_points))
	v = np.zeros((len(bodies), data_points))
	w = np.zeros((len(bodies), data_points))
	for i in range(len(bodies)):
		u_i, v_i, w_i = dl.uvw_data(bodies[i])
		u[i] = u_i
		v[i] = v_i	
		w[i] = w_i
	return u, v, w

def get_geo_data(bodies, data_points):
	# assemble a matrix with xy data from all bodies
	print("Importing geometric elements...")
	a = np.zeros((len(bodies), data_points))
	e = np.zeros((len(bodies), data_points))
	m = np.zeros((len(bodies), data_points))
	p = np.zeros((len(bodies), data_points))
	for i in range(len(bodies)):
		a_i, e_i, m_i, p_i = dl.geo1245data(bodies[i])
		a[i] = a_i
		e[i] = e_i
		m[i] = m_i	
		p[i] = p_i
	return a, e, m, p

# FUNCTIONS TO SORT OUT CLOSE ENCOUNTERS ------------------------------------------------------------------------------------

def compare_space_between(bodies, x, y, t, hill_radii, body_radii, i, j, unite=True):
	AU = cm6.AU()
	space_between = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2) * AU * 1.0E-03 # units of km  
	#print('   closest approach: ' + str(np.round(np.amin(space_between), decimals=2)) + ' km')
	relevant = highlight_close_approaches(space_between, hill_radii[i], len(t))
	if(len(relevant) != 0):
		if(unite == True):
			begindex, endindex = distinguish_events(relevant, bodies, i, j)
			#print('indices after running distinguish_events:')
			#for k in range(len(begindex)):
			#	print(begindex[k], endindex[k])
			begidx, endidx = unite_near_events(begindex, endindex)
			#print('indices after running unite_near_events:')
			#for k in range(len(begidx)):
			#	print(begidx[k], endidx[k])
		else:
			begidx, endidx = distinguish_events(relevant, bodies, i, j)
		some_encounter = True
	else:
		#print('No close encounters between ' + str(bodies[i]) + ' and ' + str(bodies[j]))
		begidx = 0
		endidx = 0
		some_encounter = False
	return begidx, endidx, some_encounter

def highlight_close_approaches(space_between, hill_radius, data_points, relevant_distance=6):
	#print('Highlighting close approaches...')
	# make a list of the indices where the bodies are within some specific distance
	relevant = []
	#print('   ' + str(relevant_distance) + ' hill radii: ' + str(np.round(relevant_distance * hill_radius, decimals=2)) + ' km')
	k = 0
	while(k < data_points):
		#uf.print_progress(k, data_points - 1)
		if(space_between[k] < relevant_distance * hill_radius):
			# iterate upward from that datum
			while(space_between[k] < 2.5 * relevant_distance * hill_radius):
				k -= 1
			# jump back from the lava	
			k += 1
			# iterate back to the datum and forward
			while(space_between[k] < 2.5 * relevant_distance * hill_radius):
				relevant.append(k)
				k += 1
		k += 1
	return relevant 

def distinguish_events(relevant, bodies, i, j):
	# determining how many separate events there are and keeping track of the indices of their beginning and end
	begidx = [] 
	endidx = []
	begidx.append(relevant[0])
	for l in range(1, len(relevant)):
		# comparing indices for continuity
		# if the indices are not consecutive then we are dealing with a separate event (another close encounter) 
		if (relevant[l] != (relevant[l - 1] + 1)):		
			endidx.append(relevant[l-1])
			begidx.append(relevant[l])

	endidx.append(relevant[-1])
	#if(len(begidx) > 1):
	#	print(str(len(begidx)) + ' close encounters detected between ' + str(bodies[i]) + ' and ' + str(bodies[j]) + '.')
	#elif(len(begidx) == 1):
	#	print(str(len(begidx)) + ' close encounter detected between ' + str(bodies[i]) + ' and ' + str(bodies[j]) + '.')
	return begidx, endidx

def unite_near_events(begindex, endindex):
	begidx = [] 
	endidx = []
	l = 0
	begidx.append(begindex[0])
	for l in range(1, len(begindex)):
		if(begindex[l] - endindex[l-1] > 200):
			endidx.append(endindex[l-1])
			begidx.append(begindex[l])
	endidx.append(endindex[-1])
	#if(len(begidx) != len(begindex)):
	#	if(len(begidx) == 1):
	#		print('Close encounter tally adjusted to ' + str(len(begidx)) + ' distinct event.')
	#	else:
	#		print('Close encounter tally adjusted to ' + str(len(begidx)) + ' distinct events.')
	return begidx, endidx

def s_outcome(masses, s_matrix, i, j, begidx, endidx, pos_delta_s, neg_delta_s, zer_delta_s):
	delta_s_i = s_matrix[i][endidx] - s_matrix[i][begidx]
	delta_s_j = s_matrix[j][endidx] - s_matrix[j][begidx]
	delta_s_m = ((masses[i] * delta_s_i) + (masses[j] * delta_s_j)) / (masses[i] + masses[j])
	outcome = 0
	if(delta_s_i > 0):
		if(delta_s_j > 0):
			pos_delta_s.append(delta_s_m)
			outcome = 2
		else:
			zer_delta_s.append(delta_s_m)
	else:
		if(delta_s_j < 0):
			neg_delta_s.append(delta_s_m)
			outcome = 1
		else:
			zer_delta_s.append(delta_s_m)
	return pos_delta_s, neg_delta_s, zer_delta_s, delta_s_m, outcome

def initialize_array_set(N, nmax=300):
	return np.zeros((N, nmax)), np.zeros((N, nmax)), np.zeros((N, nmax)), np.zeros((N, nmax))

# FUNCTIONS TO COMPUTE ------------------------------------------------------------------------------------------------------

def compute_precession(bodies, a):
	# using equation 0.36 in Hedman, An Introduction to Planetary Ring Dynamics
	time_avg_derivative = np.zeros(len(bodies))
	for i in range(len(bodies)):
		time_avg_derivative[i] = 1.5 * np.mean(om.n_geo(bodies[i])) * 8.64E+04 * bd.J2() * ((bd.R_Sat() / np.mean(a[i]))**2)
	return time_avg_derivative # units of deg/day

def xy_path(i, j, begidx, endidx, x, y):
	# defining the fixed frame for this body
	theta = om.define_corotating_frame(x[i], y[i])			
	x_i_c = (x[i] * np.cos(theta)) + (y[i] * np.sin(theta))
	# convert to fixed frame
	x_j_c = (x[j] * np.cos(theta)) + (y[j] * np.sin(theta))
	y_f =  (-x[j] * np.sin(theta)) + (y[j] * np.cos(theta))
	x_f = x_j_c - x_i_c
	AU_km = cm6.AU() * 1.0E-03
	return x_f * AU_km, y_f * AU_km

def angular_momentum(i, j, begidx, endidx, t, x, y, z, u, v, w, masses):
	AU_m = cm6.AU()
	L_i = masses[i] * np.sqrt((((y[i] * w[i]) - (z[i] * v[i]))**2) + (((z[i] * u[i]) - (x[i] * w[i]))**2) + (((x[i] * v[i]) - (y[i] * u[i]))**2)) * (AU_m**2) / 8.64E+04
	L_j = masses[j] * np.sqrt((((y[j] * w[j]) - (z[j] * v[j]))**2) + (((z[j] * u[j]) - (x[j] * w[j]))**2) + (((x[j] * v[j]) - (y[j] * u[j]))**2)) * (AU_m**2) / 8.64E+04
	delta_L_i = L_i[begidx:endidx+1] - L_i[begidx]
	delta_L_j = L_j[begidx:endidx+1] - L_j[begidx] 
	return delta_L_i, delta_L_j, L_i, L_j

def energy(i, j, begidx, endidx, energy_matrix):
	E3 = energy_matrix[i]
	E4 = energy_matrix[j]
	delta_E3 = E3[begidx:endidx+1] - E3[begidx]
	delta_E4 = E4[begidx:endidx+1] - E4[begidx]
	return delta_E3, delta_E4, E3, E4 

def mean_longitude_difference(lmda1, lmda0):
	lmda_diff = lmda1 - lmda0
	for q in range(len(lmda_diff)):
		if(lmda_diff[q] > 180):
			lmda_diff[q] -= 360
		elif(lmda_diff[q] < -180):
			lmda_diff[q] += 360
	return lmda_diff

def angle_from_pericenter(lmda, curlypi):
	if((lmda - curlypi) < 0):
		lmda += 360
	if((lmda - curlypi) < 180):
		proximity = (lmda - curlypi)
	else:
		proximity = 360 - (lmda - curlypi)
	return proximity

def stats_at_closest_approach(i, j, begidx, endidx, lmda_Mim, lmda, curlypi, x, y, x_f, y_f):
	space_between = np.sqrt((x_f[begidx:endidx])**2 + (y_f[begidx:endidx])**2)
	prox_idx =  np.argmin(space_between)
	#print('   close approach: ' + str(np.round(np.amin(space_between), decimals=2)) + ' km')
	proximity_i = angle_from_pericenter(lmda[i][begidx+prox_idx], curlypi[i][begidx+prox_idx])
	proximity_j = angle_from_pericenter(lmda[j][begidx+prox_idx], curlypi[j][begidx+prox_idx])
	lmda_i_prox = lmda[i][begidx+prox_idx]
	lmda_Mim_prox = lmda_Mim[begidx+prox_idx]
	return proximity_i, proximity_j, x_f[begidx+prox_idx], y_f[begidx+prox_idx], lmda_Mim_prox, lmda_i_prox

def compute_s(a, a_cer, lmda, lambda_cer):
	return np.sqrt((((a - a_cer) * 1.0E-03) / 37)**2 + ((((lmda - lambda_cer + 180) % 360) - 180) / 30)**2)

def compute_phase_space_distance_for_all_bodies(bodies, a, a_cer, lmda, lambda_cer):
	print('Computing phase space distance from exact resonance...')
	s_matrix = np.zeros((len(bodies), len(a_cer)))
	for i in range(len(bodies)):
		s_matrix[i] = compute_s(a[i], a_cer, lmda[i], lambda_cer)
	#s_m = ((masses[i] * s_i) + (masses[j] * s_j)) / (masses[i] + masses[j])
	return s_matrix

def compute_energy(G, AU_m, bodies, masses, M_Sat, R_Sat, J2, J4, J6, M_Mim, x_Mim, y_Mim, z_Mim, x, y, z, u, v, w, i):
	v_i = np.sqrt((u[i]**2) + (v[i]**2) + (w[i]**2)) * AU_m / 8.64E+04
	r_i = np.sqrt((x[i]**2) + (y[i]**2) + (z[i]**2)) * AU_m
	j2term = J2 * ((R_Sat / r_i)**2) * uf.P2(z[i] * AU_m / r_i)
	j4term = J4 * ((R_Sat / r_i)**4) * uf.P4(z[i] * AU_m / r_i)
	j6term = J6 * ((R_Sat / r_i)**6) * uf.P6(z[i] * AU_m / r_i)
	U13 = -(G * M_Sat * masses[i] / r_i) * (1 - (j2term + j4term + j6term))
	U23 = -G * M_Mim * masses[i] / (np.sqrt(((x[i] - x_Mim)**2) + ((y[i] - y_Mim)**2) + ((z[i] - z_Mim)**2)) * AU_m)
	U = U13 + U23
	T = 0.5 * masses[i] * (v_i**2)
	# energy due to each body in the ring arc
	for k in range(len(bodies)):
		if(i != k):
			r_k = np.sqrt(((x[k] - x[i])**2) + ((y[k] - y[i])**2) + ((z[k] - z[i])**2)) * AU_m
			U += -(G * masses[i] * masses[k] / r_k)
	return T + U

def compute_energy_of_Mimas(G, AU_m, M_Sat, R_Sat, J2, J4, J6, M_Mim, x_Mim, y_Mim, z_Mim, u_Mim, v_Mim, w_Mim, x, y, z, u, v, w, bodies, masses):
	v_i = np.sqrt((u_Mim**2) + (v_Mim**2) + (w_Mim**2)) * AU_m / 8.64E+04
	r_i = np.sqrt((x_Mim**2) + (y_Mim**2) + (z_Mim**2)) * AU_m
	j2term = J2 * ((R_Sat / r_i)**2) * uf.P2(z_Mim * AU_m / r_i)
	j4term = J4 * ((R_Sat / r_i)**4) * uf.P4(z_Mim * AU_m / r_i)
	j6term = J6 * ((R_Sat / r_i)**6) * uf.P6(z_Mim * AU_m / r_i)
	U = -(G * M_Sat * M_Mim / r_i) * (1 - (j2term + j4term + j6term))
	T = 0.5 * M_Mim * (v_i**2)
	# energy due to each body in the ring arc
	for k in range(len(bodies)):
		r_k = np.sqrt(((x[k] - x_Mim)**2) + ((y[k] - y_Mim)**2) + ((z[k] - z_Mim)**2)) * AU_m
		U += -(G * masses[k] * M_Mim / r_k)
	return T + U 

def compute_energy_for_all_bodies(bodies, masses, M_Sat, M_Mim, x_Mim, y_Mim, z_Mim, u_Mim, v_Mim, w_Mim, x, y, z, u, v, w):
	G = cm6.G()
	AU_m = cm6.AU()
	R_Sat = bd.R_Sat()
	J2 = bd.J2()
	J4 = bd.J4()
	J6 = bd.J6()
	print('Computing energy...')
	energy_matrix = np.zeros((len(bodies) + 2, len(x_Mim)))
	energy_Mim = compute_energy_of_Mimas(G, AU_m, M_Sat, R_Sat, J2, J4, J6, M_Mim, x_Mim, y_Mim, z_Mim, u_Mim, v_Mim, w_Mim, x, y, z, u, v, w, bodies, masses)
	energy_matrix[len(bodies)] = energy_Mim
	energy_matrix[len(bodies) + 1] += energy_matrix[len(bodies)]
	for i in range(len(bodies)):
		energy_matrix[i] = compute_energy(G, AU_m, bodies, masses, M_Sat, R_Sat, J2, J4, J6, M_Mim, x_Mim, y_Mim, z_Mim, x, y, z, u, v, w, i)
		energy_matrix[len(bodies) + 1] += energy_matrix[i]
	return energy_matrix



# PLOTTING FUNCTIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# CLOSE ENCOUNTER PLOTS -------------------------------------------------------------------

def plot_xy_path(fig, bodies, i, j, t, begidx, endidx, hill_radius, body_radius, x_f, y_f, x, y, nrows, ncols, initrow, initcol, rs, cs):
	td = plt.gca().transData
	ax1 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	plt.gca().set_aspect('equal', adjustable='box')
	#fig, ax1 = pa.title_and_axes(fig, ax1, "Path of " + uf.decapitalize(bodies[j]) + "\nin " + uf.decapitalize(bodies[i]) + "'s frame", 'x [km]', 'y [km]', 
		#xmin=-hill_radius * 9, xmax=hill_radius * 9, ymin=-hill_radius * 9, ymax=hill_radius * 9)
	fig, ax1 = pa.title_and_axes(fig, ax1, '', 'x [km]', 'y [km]', 
		xmin=-hill_radius * 9, xmax=hill_radius * 9, ymin=-hill_radius * 9, ymax=hill_radius * 9)
	hill_circle = plt.Circle((0,0), hill_radius, color='blue', fill=False, linestyle='dashed')
	hyb_changeover_circle = plt.Circle((0,0), hill_radius * 3, color='green', fill=False, linestyle='dashed')
	body_circle = plt.Circle((0,0), body_radius, color='red', fill=True)
	ax1.add_artist(hill_circle)
	ax1.add_artist(hyb_changeover_circle)
	ax1.add_artist(body_circle)
	
	line = ax1.plot(x_f[begidx:endidx+1], y_f[begidx:endidx+1], color='c', label=str(bodies[j]) + " (days " 
		+ str("{:.2f}".format(t[begidx])) + " - " + str("{:.2f}".format(t[endidx])) + ")", linewidth=5)[0]
	uf.add_arrows(line, begidx, endidx)

	if(hill_radius < 5):
		ax1.set_xticks([-30, -15, 0, 15, 30])
		ax1.set_yticks([-30, -15, 0, 15, 30])
	else:
		ax1.set_xticks([-100, -50, 0, 50, 100])
		ax1.set_yticks([-100, -50, 0, 50, 100])
	
	# set color cycle 
	ax1.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(bodies))))
	AU = cm6.AU()
	for k in range(len(bodies)):
		if((k != i) and (k != j)):
			# compare space between during the close encounter but for the other bodies
			close_by = False
			space_between = np.sqrt((x[i] - x[k])**2 + (y[i] - y[k])**2) * AU * 1.0E-03 # units of km
			for l in range(endidx - begidx):
				if(space_between[begidx+l] < 100):
					close_by = True
			if(close_by == True):
				x_f, y_f = xy_path(i, k, begidx, endidx, x, y)
				ax1.plot(x_f[begidx:endidx+1], y_f[begidx:endidx+1], label=str(bodies[k]))
	ax1.text(-0.02, 1.02, 'Body 4', color='c', fontsize=24, transform=ax1.transAxes)
	ax1.text(0.205, 1.02, "'s path in", color='k', fontsize=24, transform=ax1.transAxes)
	ax1.text(0.530, 1.02, 'Body 3', color='r', fontsize=24, transform=ax1.transAxes)
	ax1.text(0.755, 1.02, "'s frame", color='k', fontsize=24, transform=ax1.transAxes)
	#ax1.legend(loc=2, prop={'size':11})
	#ax1.spines['top'].set_visible(False)
	#ax1.spines['right'].set_visible(False)
	return fig, ax1

def plot_phase_space_path(fig, bodies, i, j, begidx, endidx, a, e, lmda, t, a_cer, lambda_cer, s_matrix, nrows, ncols, initrow, initcol, rs, cs):
	ax2 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax2 = pa.title_and_axes(fig, ax2, 'Phase space', r'$a - a_{CER}$' + ' [km]', r'$\lambda - \lambda_{CER}$' + ' [deg]', xmin=-39, xmax=39, ymin=-31, ymax=31)
	# contour of the edge of the corotation site
	fringe = patches.Ellipse(xy=[0,0], width=74, height=60, angle=0, color='g', fill=False, linestyle='dashed')
	plt.draw()
	ax2.add_patch(fringe)
	# error bars to account for eccentricity
	#xerri = a[i][int((begidx + endidx) / 2)] * e[i][int((begidx + endidx) / 2)] * 1.0E-03
	#xerrj = a[j][int((begidx + endidx) / 2)] * e[j][int((begidx + endidx) / 2)] * 1.0E-03
	#print('xerri: ' + str(xerri))
	#print('xerrj: ' + str(xerrj))
	#
	#ax2.errorbar((a[j][int((begidx + endidx) / 2)] - a_cer[int((begidx + endidx) / 2)]) * 1.0E-03, ((lmda[j][int((begidx + endidx) / 2)] - lambda_cer[int((begidx + endidx) / 2)] + 180) % 360) - 180, xerr=xerrj, color='c', capsize=10, capthick=2)
	#ax2.errorbar((a[i][int((begidx + endidx) / 2)] - a_cer[int((begidx + endidx) / 2)]) * 1.0E-03, ((lmda[i][int((begidx + endidx) / 2)] - lambda_cer[int((begidx + endidx) / 2)] + 180) % 360) - 180, xerr=xerri, color='r', capsize=10, capthick=2)

	ax2.scatter((a[j][begidx:endidx] - a_cer[begidx:endidx]) * 1.0E-03, ((lmda[j][begidx:endidx] - lambda_cer[begidx:endidx] + 180) % 360) - 180, 
		label=uf.decapitalize(bodies[j]) + ' (' + r'$s_0 =$' + ' ' + str("{:.3f}".format(s_matrix[j][begidx])) + ', ' + r'$s_f =$' + ' ' + str("{:.3f}".format(s_matrix[j][endidx]))+ ')', color='c', s=180)
	ax2.scatter((a[i][begidx:endidx] - a_cer[begidx:endidx]) * 1.0E-03, ((lmda[i][begidx:endidx] - lambda_cer[begidx:endidx] + 180) % 360) - 180, 
		label=uf.decapitalize(bodies[i]) + ' (' + r'$s_0 =$' + ' ' + str("{:.3f}".format(s_matrix[i][begidx])) + ', ' + r'$s_f =$' + ' ' + str("{:.3f}".format(s_matrix[i][endidx]))+ ')', color='r', s=180)
	
	#ax2.legend(loc=0, prop={'size':14})
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	return fig, ax2

def plot_s(fig, bodies, i, j, begidx, endidx, t, s_matrix, nrows, ncols, initrow, initcol, rs, cs):
	ax3 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax3 = pa.title_and_axes(fig, ax3, 'Phase space distance \nfrom exact corotation', 't  [days]', 's', t[begidx], t[endidx], 0, 1.0)
	ax3.plot(t[begidx:endidx+1], s_matrix[j][begidx:endidx+1], color='c', label=uf.decapitalize(bodies[j]))
	ax3.plot(t[begidx:endidx+1], s_matrix[i][begidx:endidx+1], color='r', label=uf.decapitalize(bodies[i]))
	start, end = ax3.get_xlim()
	ax3.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	#ax3.legend(loc=0, prop={'size':14})
	return fig, ax3

def plot_s_j(fig, bodies, j, begidx, endidx, t, s_matrix, nrows, ncols, initrow, initcol, rs, cs):
	ax3a = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax3a = pa.title_and_axes(fig, ax3a, 'Phase space distance from CER', '', 's', t[begidx], t[endidx])
	ax3a.plot(t[begidx:endidx+1], s_matrix[j][begidx:endidx+1], color='c', label=uf.decapitalize(bodies[j]), linewidth=5)
	ymax = np.amax(s_matrix[j][begidx:endidx+1])
	start, end = ax3a.get_xlim()
	ax3a.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	#ax3a.text(0.5,0.9, 'Body 4', horizontalalignment='center', verticalalignment='center', fontsize=16, transform=ax3a.transAxes)
	#for xlabel_i in ax3a.axes.get_xticklabels():
		#xlabel_i.set_visible(False)
	#ax3a.axes.get_xaxis().set_ticks([])
	ax3a.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	#ax3a.legend(loc=0, prop={'size':14})
	plt.setp(ax3a.get_xticklabels(), visible=False) # make these tick labels invisible
	ax3a.spines['bottom'].set_visible(False)
	ax3a.spines['top'].set_visible(False)
	ax3a.spines['right'].set_visible(False)
	ax3a.tick_params(axis='x', colors='white')

	d = 0.015 # how big to make the diagonal lines in axes coordinates
	kwargs = dict(transform=ax3a.transAxes, color='k', clip_on=False)
	ax3a.plot((-d,+d), (+d,-d), **kwargs)
	ax3a.plot((-d,+d), (0,-2*d), **kwargs)

	return fig, ax3a

def plot_s_i(fig, bodies, i, begidx, endidx, t, s_matrix, nrows, ncols, initrow, initcol, rs, cs):
	ax3b = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax3b = pa.title_and_axes(fig, ax3b, '', 't  [days]', 's', t[begidx], t[endidx])
	ax3b.plot(t[begidx:endidx+1], s_matrix[i][begidx:endidx+1], color='r', label=uf.decapitalize(bodies[i]), linewidth=5)
	ymax = np.amax(s_matrix[i][begidx:endidx+1])
	start, end = ax3b.get_xlim()
	ax3b.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	#ax3b.text(0.5, 0.9, 'Body 3', horizontalalignment='center', verticalalignment='center', fontsize=16, transform=ax3b.transAxes)
	ax3b.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	#ax3.legend(loc=0, prop={'size':14})
	plt.setp(ax3b.get_xticklabels(), visible=False) # make these tick labels invisible
	ax3b.set_xlabel('')
	ax3b.spines['top'].set_visible(False)
	ax3b.spines['right'].set_visible(False)

	d = 0.015 # how big to make the diagonal lines in axes coordinates
	kwargs = dict(transform=ax3b.transAxes, color='k', clip_on=False)
	ax3b.plot((-d,+d), (1+d,1-d), **kwargs)
	ax3b.plot((-d,+d), (1+d+d,1), **kwargs)

	return fig, ax3b

def plot_e(fig, bodies, i, j, begidx, endidx, e, t, nrows, ncols, initrow, initcol, rs, cs):
	ax4 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	delta_e_i = (e[i] - e[i][begidx]) * 1.0E+06
	delta_e_j = (e[j] - e[j][begidx]) * 1.0E+06
	ymax = np.amax(np.fmax(np.absolute(delta_e_i[begidx:endidx+1]), np.absolute(delta_e_j[begidx:endidx+1]))) * 1.1
	fig, ax4 = pa.title_and_axes(fig, ax4, 'Change in Eccentricity', 't [days]', r'$e-e_0$' + ' [' + r'$10^{-6}$' + ']', t[begidx], t[endidx])
	ax4.plot(t[begidx:endidx+1], delta_e_j[begidx:endidx+1], color='c', label=uf.decapitalize(bodies[j]) + ' (' + r'$e_0=$' + str(np.round(e[j][begidx], decimals=7)) + ')', linewidth=5)
	ax4.plot(t[begidx:endidx+1], delta_e_i[begidx:endidx+1], color='r', label=uf.decapitalize(bodies[i]) + ' (' + r'$e_0=$' + str(np.round(e[i][begidx], decimals=7)) + ')', linewidth=5)
	start, end = ax4.get_xlim()
	ax4.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	#ax4.legend(loc=0, prop={'size':14})
	#plt.setp(ax4.get_xticklabels(), visible=False) # make these tick labels invisible
	#ax4.set_xlabel('')
	ax4.spines['top'].set_visible(False)
	ax4.spines['right'].set_visible(False)
	return fig, ax4

def plot_angular_momentum(fig, bodies, i, j, begidx, endidx, t, delta_L_i, delta_L_j, L_i, L_j, nrows, ncols, initrow, initcol, rs, cs):
	ymax = np.amax(np.fmax(np.absolute(delta_L_i), np.absolute(delta_L_j))) * 1.1
	ax5 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax5 = pa.title_and_axes(fig, ax5, 'Change in Angular Momentum', 't [days]', r'$L-L_0$' + ' [kg m' + r'$^2$' + '/s]', 
		t[begidx], t[endidx], -ymax, ymax)
	ax5.plot(t[begidx:begidx+len(delta_L_j)], delta_L_j, color='c', label=uf.decapitalize(bodies[j]) + ' (' + r'$L_0=$' + str("{:.6e}".format(L_j[begidx])) + ' kg m' + r'$^2$' + '/s' + ')', linewidth=5)
	ax5.plot(t[begidx:begidx+len(delta_L_i)], delta_L_i, color='r', label=uf.decapitalize(bodies[i]) + ' (' + r'$L_0=$' + str("{:.6e}".format(L_i[begidx])) + ' kg m' + r'$^2$' + '/s' + ')', linewidth=5)
	start, end = ax5.get_xlim()
	ax5.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	uf.modify_yaxis_label(ax5, ymax)
	ax5.legend(loc=0, prop={'size':14})
	#if(len(delta_L_i) != (1 + endidx - begidx)):
	#	print('Weird: len(delta_L_i) = ' + str(len(delta_L_i)) + ', 1 + endidx - begidx = ' + str(1 + endidx - begidx))
	plt.setp(ax5.get_xticklabels(), visible=False) # make these tick labels invisible
	ax5.set_xlabel('')
	ax5.spines['top'].set_visible(False)
	ax5.spines['right'].set_visible(False)
	return fig, ax5

def plot_delta_v(fig, bodies, i, j, begidx, endidx, t, x, y, z, u, v, w, masses):
	G = cm6.G()
	AU_m = cm6.AU()
	ax6 = plt.subplot2grid((18,22),(12, 6), rowspan=6, colspan=8)
	delta_v = np.sqrt(((u[i] - u[j])**2) + ((v[i] - v[j])**2) + ((w[i] - w[j])**2)) * AU_m / 8.64E+04
	v_esc = np.sqrt(2 * G * masses[i] / (np.sqrt((((x[i] - x[j]) * AU_m)**2) + (((y[i] - y[j]) * AU_m)**2) + (((z[i] - z[j]) * AU_m)**2))))
	fig, ax6 = pa.title_and_axes(fig, ax6, 'Relative velocity', 't [days]', r'$\Delta$' + r'$v$' + ' (m/s)', t[begidx], t[endidx], 0, 50)
	ax6.plot(t[begidx:endidx+1], delta_v[begidx:endidx+1], color='g', label=r'$\Delta$' + r'$v$')
	ax6.plot(t[begidx:endidx+1], v_esc[begidx:endidx+1],   color='r', label=r'$v_{esc}$')
	start, end = ax6.get_xlim()
	ax6.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax6.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	ax6.legend(loc=0, prop={'size':14})
	return fig, ax6

def plot_F_by_lambda_offset(fig, bodies, i, j, begidx, endidx, x, y, z, x_Mim, y_Mim, z_Mim, lmda, lmda_Mim, masses, M_Mim, nrows, ncols, initrow, initcol, rs, cs):
	G = cm6.G()
	AU_m = cm6.AU()
	ax6 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	fig, ax6 = pa.title_and_axes(fig, ax6, 'F(' + r'$\lambda$' + ')', r'$\lambda - \lambda_{encounter}$', r'$|F|$' + ' [N]', -180, 180)
	ax6.set_yscale('log')
	# set color cycle 
	ax6.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(bodies))))
	for k in range(len(bodies)):
		if(k != i):
			F_mag = G * masses[k] * masses[i] / ((np.sqrt(((x[i] - x[k])**2) + ((y[i] - y[k])**2) + ((z[i] - z[k])**2)) * AU_m)**2)
			lambda_offset = mean_longitude_difference(lmda[k], lmda[i])
			ax6.scatter(lambda_offset[begidx:endidx+1], F_mag[begidx:endidx+1], label=str(bodies[k]), s=15)
	# due to Mimas
	F_mag = G * M_Mim * masses[i] / ((np.sqrt(((x[i] - x_Mim)**2) + ((y[i] - y_Mim)**2) + ((z[i] - z_Mim)**2)) * AU_m)**2)
	lambda_offset = mean_longitude_difference(lmda_Mim, lmda[i])
	ax6.scatter(lambda_offset[begidx:endidx+1], F_mag[begidx:endidx+1], label='Mimas', s=15, color='r')
	ax6.set_xticks( [-180, -90, 0, 90, 180])
	ax6.legend(loc=0, prop={'size':11})
	return fig, ax6 

def plot_pericenter(fig, bodies, i, j, begidx, endidx, curlypi, t, precession, nrows, ncols, initrow, initcol, rs, cs):
	ax7 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	delta_curlypi_i = curlypi[i] - curlypi[i][begidx] 
	delta_curlypi_j = curlypi[j] - curlypi[j][begidx]
	fig, ax7 = pa.title_and_axes(fig, ax7, 'Change in Longitude of pericenter', 't [days]', r'$\varpi-\varpi_0$' + ' [deg]', t[begidx], t[endidx])
	ax7.plot(t[begidx:endidx+1], delta_curlypi_j[begidx:endidx+1], color='c', label=uf.decapitalize(bodies[j]) + ' (' + r'$\varpi_0$' + r'$\approx$' + str(int(curlypi[j][begidx])) + ')')
	ax7.plot(t[begidx:endidx+1], delta_curlypi_i[begidx:endidx+1], color='r', label=uf.decapitalize(bodies[i]) + ' (' + r'$\varpi_0$' + r'$\approx$' + str(int(curlypi[i][begidx])) + ')')
	ax7.plot(t[begidx:endidx+1], delta_curlypi_i[begidx] + (precession[i] * (t[begidx:endidx+1] - t[begidx])), color='g', linestyle='-', label=r'$\left< \frac{d\varpi}{dt} \right>$')
	ax7.legend(loc=0, prop={'size':14})
	start, end = ax7.get_xlim()
	ax7.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax7.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	#ax7.set_yticks([0, 90, 180, 270, 360])
	return fig, ax7

def plot_energy(fig, bodies, i, j, begidx, endidx, t, delta_E3, delta_E4, E3, E4, nrows, ncols, initrow, initcol, rs, cs):
	ax8 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	ymax = np.amax(np.fmax(np.absolute(delta_E3), np.absolute(delta_E4))) * 1.1
	fig, ax8 = pa.title_and_axes(fig, ax8, 'Change in Energy', 't [days]', r'$E-E_0$' + ' [J]', t[begidx], t[endidx], -ymax, ymax)
	ax8.plot(t[begidx:begidx+len(delta_E4)], delta_E4, color='c', label=uf.decapitalize(bodies[j]) + ' (' + r'$E_0=$' + str("{:.5e}".format(E4[begidx])) + ' J)', linewidth=5)
	ax8.plot(t[begidx:begidx+len(delta_E3)], delta_E3, color='r', label=uf.decapitalize(bodies[i]) + ' (' + r'$E_0=$' + str("{:.5e}".format(E3[begidx])) + ' J)', linewidth=5)
	start, end = ax8.get_xlim()
	ax8.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax8.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	uf.modify_yaxis_label(ax8, ymax)
	ax8.legend(loc=0, prop={'size':14})
	plt.setp(ax8.get_xticklabels(), visible=False) # make these tick labels invisible
	ax8.set_xlabel('')
	ax8.spines['top'].set_visible(False)
	ax8.spines['right'].set_visible(False)
	return fig, ax8

def plot_delta_r(fig, bodies, i, j, begidx, endidx, t, x, y, z, hill_radius, body_radius):
	AU_km = cm6.AU() * 1.0E-03
	ax9 = plt.subplot2grid((18,22),(12, 14), rowspan=6, colspan=8)
	delta_r = np.sqrt((((x[i] - x[j]) * AU_km)**2) + (((y[i] - y[j]) * AU_km)**2) + (((z[i] - z[j]) * AU_km)**2))
	fig, ax9 = pa.title_and_axes(fig, ax9, 'Relative distance', 't [days]', r'$\Delta$' + r'$r$' + ' (km)', t[begidx], t[endidx], 0, np.amax(delta_r[begidx:endidx]) * 1.1)
	ax9.plot(t[begidx:endidx+1], delta_r[begidx:endidx+1], color='g', label=r'$\Delta$' + r'$r$')
	ax9.axhline(y=hill_radius, color='orangered', linestyle='-', label='Hill radius')
	ax9.axhline(y=3 * hill_radius, color='goldenrod', linestyle='-', label='Hybrid changeover radius')
	ax9.axhline(y=body_radius, color='chartreuse', linestyle='-', label='Body radius')
	start, end = ax9.get_xlim()
	ax9.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax9.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	ax9.legend(loc=0, prop={'size':14})
	return fig, ax9

def plot_positions_of_other_bodies(fig, bodies, i, j, begidx, endidx, t, x, y, nrows, ncols, initrow, initcol, rs, cs):
	ax9 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	plt.gca().set_aspect('equal', adjustable='box')
	boxsize = 1.0E+05 
	fig, ax9 = pa.title_and_axes(fig, ax9, "Whereabouts of other bodies", 'x [km]', 'y [km]', 
		xmin=-boxsize, xmax=boxsize, ymin=-boxsize, ymax=boxsize)

	# scatter point at origin (location of close encounter)
	ax9.scatter(0, 0, color='r', s=10)

	# set color cycle 
	ax9.set_prop_cycle('color', plt.cm.rainbow(np.linspace(0, 1, len(bodies))))

	for k in range(len(bodies)):
		if((k != i) and (k != j)):
			x_f, y_f = xy_path(i, k, begidx, endidx, x, y)
			ax9.scatter(x_f[(begidx):(endidx+1)], y_f[(begidx):(endidx+1)], label=str(bodies[k]), s=10)
	
	ax9.legend(loc=2, prop={'size':11})
	return fig, ax9

def plot_a(fig, bodies, i, j, begidx, endidx, a, a_cer, t, nrows, ncols, initrow, initcol, rs, cs):
	ax10 = plt.subplot2grid((nrows,ncols),(initrow, initcol), rowspan=rs, colspan=cs)
	#ymax = np.amax(np.fmax(e[i][begidx:endidx], e[j][begidx:endidx])) * 1.5
	delta_a_i = (a[i] - a[i][begidx]) * 1.0E-03
	delta_a_j = (a[j] - a[j][begidx]) * 1.0E-03
	fig, ax10 = pa.title_and_axes(fig, ax10, 'Change in Semi-major axis', 't [days]', r'$a - a_0$' + ' [km]', t[begidx], t[endidx])
	ax10.plot(t[begidx:endidx+1], delta_a_j[begidx:endidx+1], color='c', label=uf.decapitalize(bodies[j]) + ' (' + r'$a_0=$' + str("{:.9e}".format(a[j][begidx] * 1.0E-03)) + ' km)', linewidth=5)
	ax10.plot(t[begidx:endidx+1], delta_a_i[begidx:endidx+1], color='r', label=uf.decapitalize(bodies[i]) + ' (' + r'$a_0=$' + str("{:.9e}".format(a[i][begidx] * 1.0E-03)) + ' km)', linewidth=5)
	start, end = ax10.get_xlim()
	ax10.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax10.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	ax10.legend(loc=0, prop={'size':14})
	#plt.setp(ax10.get_xticklabels(), visible=False) # make these tick labels invisible
	#ax10.set_xlabel('')
	ax10.spines['top'].set_visible(False)
	ax10.spines['right'].set_visible(False)
	return fig, ax10


def close_encounter_plots(bodies, i, j, a, e, lmda, curlypi, t, a_cer, lambda_cer, begidx, endidx, hill_radii, body_radii, x, y, z, u, v, w, masses, M_Sat, M_Mim, x_Mim, y_Mim, z_Mim, x_f, y_f, delta_L_i, delta_L_j, L_i, L_j, delta_E3, delta_E4, E3, E4, s_matrix, lmda_Mim, black=False):
	# plots for each close encounter
	print('Making plots...')
	nrows = 18
	ncols = 22
	for m in range(len(begidx)):
		if (black == True):
			fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
		else:
			fig = plt.figure(figsize=(ncols,nrows))
		# aerial view of close encounter
		fig, ax1 = plot_xy_path(fig, bodies, i, j, t, begidx[m], endidx[m], hill_radii[i], body_radii[i], x_f, y_f, x, y, nrows, ncols, 0, 0, 6, 6)

		# phase space plot 
		fig, ax2 = plot_phase_space_path(fig, bodies, i, j, begidx[m], endidx[m], a, e, lmda, t, a_cer, lambda_cer, s_matrix, nrows, ncols, 6, 0, 6, 6)

		# plot of s 
		fig, ax3 = plot_s(fig, bodies, i, j, begidx[m], endidx[m], t, s_matrix, nrows, ncols, 12, 0, 6, 6)

		# plot of eccentricity
		fig, ax4 = plot_e(fig, bodies, i, j, begidx[m], endidx[m], e, t, nrows, ncols, 0, 6, 6, 8)
		
		# plot of angular momentum
		fig, ax5 = plot_angular_momentum(fig, bodies, i, j, begidx[m], endidx[m], t, delta_L_i[m], delta_L_j[m], L_i, L_j, nrows, ncols, 6, 6, 6, 8)

		# plot of relative velocity
		#fig, ax6 = plot_delta_v(fig, bodies, i, j, begidx[m], endidx[m], t, x, y, z, u, v, w, masses)
		fig, ax6 = plot_F_by_lambda_offset(fig, bodies, i, j, begidx[m], endidx[m], x, y, z, x_Mim, y_Mim, z_Mim, lmda, lmda_Mim, masses, M_Mim, nrows, ncols, 12, 6, 6, 8)

		# plot of pericenter 
		precession = compute_precession(bodies, a)
		fig, ax7 = plot_pericenter(fig, bodies, i, j, begidx[m], endidx[m], curlypi, t, precession, nrows, ncols, 0, 14, 6, 8)

		# plot of energy
		fig, ax8 = plot_energy(fig, bodies, i, j, begidx[m], endidx[m], t, delta_E3[m], delta_E4[m], E3, E4, nrows, ncols, 6, 14, 6, 8)

		# plot of relative distance
		#fig, ax9 = plot_delta_r(fig, bodies, i, j, begidx[m], endidx[m], t, x, y, z, hill_radii[i], body_radii[i])

		# plot of wherabouts of other bodies
		fig, ax9 = plot_positions_of_other_bodies(fig, bodies, i, j, begidx[m], endidx[m], t, x, y, nrows, ncols, 12, 14, 6, 8)

		plt.tight_layout()
		#pa.save_and_clear_plot(fig, ax1, filename='close_encounter_' + str(m) + '_between_' + str(bodies[i]) + '_and_' + str(bodies[j]), black=black, dpi=400)	
		ax1.figure.savefig('close_encounter_' + str(m) + '_between_' + str(bodies[i]) + '_and_' + str(bodies[j]) + '.png', dpi=100)	
		plt.clf()
		#uf.print_progress(m + 1, len(begidx))		

def close_encounter_plots_for_paper(bodies, i, j, a, e, lmda, curlypi, t, a_cer, lambda_cer, begidx, endidx, hill_radii, body_radii, x, y, z, u, v, w, masses, M_Sat, M_Mim, x_Mim, y_Mim, z_Mim, x_f, y_f, delta_L_i, delta_L_j, L_i, L_j, delta_E3, delta_E4, E3, E4, s_matrix, lmda_Mim, black=False):
	#print('Making plots...')
	nrows = 16
	ncols = 23
	for m in range(len(begidx)):
		#uf.print_progress(m + 1, len(begidx))
		plot_for_paper = False
		if((m == 1) and (i == 0) and (j == 10)): # this one is just so that the plots for the other three are good 
			plot_for_paper = True
		if((m == 2) and (i == 0) and (j == 10)):
			plot_for_paper = True
		if((m == 1) and (i == 5) and (j == 8)):
			plot_for_paper = True
		if((m == 0) and (i == 6) and (j == 9)):
			plot_for_paper = True
		if(plot_for_paper == True):
			if (black == True):
				fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
			else:
				fig = plt.figure(figsize=(ncols,nrows))
			# aerial view of close encounter
			fig, ax1 = plot_xy_path(fig, bodies, i, j, t, begidx[m], endidx[m], hill_radii[i], body_radii[i], x_f, y_f, x, y, nrows, ncols, 0, 0, 7, 7)
	
			# phase space plot 
			fig, ax2 = plot_phase_space_path(fig, bodies, i, j, begidx[m], endidx[m], a, e, lmda, t, a_cer, lambda_cer, s_matrix, nrows, ncols, 7, 0, 8, 7)
	
			# plot of s 
			fig, ax3a = plot_s_j(fig, bodies, j, begidx[m], endidx[m], t, s_matrix, nrows, ncols, 0, 7, 5, 8)
			fig, ax3b = plot_s_i(fig, bodies, i, begidx[m], endidx[m], t, s_matrix, nrows, ncols, 5, 7, 5, 8)
	
			# plot of semi-major axis
			fig, ax10 = plot_a(fig, bodies, i, j, begidx[m], endidx[m], a, a_cer, t, nrows, ncols, 10, 15, 5, 8)
	
			# plot of eccentricity
			fig, ax4 = plot_e(fig, bodies, i, j, begidx[m], endidx[m], e, t, nrows, ncols, 10, 7, 5, 8)
			
			# plot of angular momentum
			fig, ax5 = plot_angular_momentum(fig, bodies, i, j, begidx[m], endidx[m], t, delta_L_i[m], delta_L_j[m], L_i, L_j, nrows, ncols, 0, 15, 5, 8)
	
			# plot of energy
			fig, ax8 = plot_energy(fig, bodies, i, j, begidx[m], endidx[m], t, delta_E3[m], delta_E4[m], E3, E4, nrows, ncols, 5, 15, 5, 8)
	
			plt.tight_layout()
			#pa.save_and_clear_plot(fig, ax1, filename='close_encounter_' + str(m) + '_between_' + str(bodies[i]) + '_and_' + str(bodies[j]), black=black, dpi=400)	
			ax1.figure.savefig('close_encounter_' + str(m) + '_between_' + str(bodies[i]) + '_and_' + str(bodies[j]) + '.png', dpi=300)	
			plt.clf()
			print('New plot!')
			

# ENERGY PLOTS -----------------------------------------------------------------------------------

def plot_E_ij(fig, begidx, endidx, t, energy_matrix, bodies, i, j):
	ax1 = plt.subplot2grid((18,22),(0, 0), rowspan=9, colspan=11)
	ymax = np.amax(np.fmax(np.absolute(energy_matrix[i][begidx:endidx+1] - energy_matrix[i][begidx]), np.absolute(energy_matrix[j][begidx:endidx+1] - energy_matrix[j][begidx]))) * 1.1
	fig, ax1 = pa.title_and_axes(fig, ax1, r'$\Delta E$' + ' for interacting bodies', 't  [days]', r'$E-E_0$' + ' [J]', t[begidx], t[endidx], -ymax, ymax)
	ax1.plot(t[begidx:endidx+1], energy_matrix[i][begidx:endidx+1] - energy_matrix[i][begidx], color='c', label=bodies[i] + ' (' + r'$E_0=$' + str("{:.3e}".format(energy_matrix[i][begidx])) + ' J)')
	ax1.plot(t[begidx:endidx+1], energy_matrix[j][begidx:endidx+1] - energy_matrix[j][begidx], color='r', label=bodies[j] + ' (' + r'$E_0=$' + str("{:.3e}".format(energy_matrix[j][begidx])) + ' J)')
	start, end = ax1.get_xlim()
	ax1.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	uf.modify_yaxis_label(ax1, ymax)
	ax1.legend(loc=0, prop={'size':14})
	return fig, ax1

def plot_E_al(fig, begidx, endidx, t, energy_matrix, bodies, i, j):
	ax2 = plt.subplot2grid((18,22),(9, 0), rowspan=9, colspan=11)
	ymax = 0
	for k in range(len(bodies)):
		if((k != i) and (k != j)):
			ymaxk = np.amax(np.absolute(energy_matrix[k][begidx:endidx+1] - energy_matrix[k][begidx]))
			if(ymaxk > ymax):
				ymax = ymaxk
	ymax *= 1.1
	fig, ax2 = pa.title_and_axes(fig, ax2, r'$\Delta E$' + ' for other bodies', 't  [days]', r'$E-E_0$' + ' [J]', t[begidx], t[endidx], -ymax, ymax)
	ax2.set_prop_cycle('color', plt.cm.nipy_spectral(np.linspace(0.1, 0.9, len(bodies) - 2)))
	for k in range(len(bodies)):
		if((k != i) and (k != j)):
			ax2.plot(t[begidx:endidx+1], energy_matrix[k][begidx:endidx+1] - energy_matrix[k][begidx], label=bodies[k] + ' (' + r'$E_0=$' + str("{:.3e}".format(energy_matrix[k][begidx])) + ' J)')
	start, end = ax2.get_xlim()
	ax2.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	uf.modify_yaxis_label(ax2, ymax)
	ax2.legend(loc=0, prop={'size':14})
	return fig, ax2

def plot_E_Mim(fig, begidx, endidx, t, E_Mim):
	ax3 = plt.subplot2grid((18,22),(0, 11), rowspan=9, colspan=11)
	ymax = np.amax(np.absolute(E_Mim[begidx:endidx+1] - E_Mim[begidx])) * 1.1
	fig, ax3 = pa.title_and_axes(fig, ax3, r'$\Delta E$' + ' for Mimas', 't  [days]', r'$E-E_0$' + ' [J]', t[begidx], t[endidx], -ymax, ymax)
	ax3.plot(t[begidx:endidx+1], E_Mim[begidx:endidx+1] - E_Mim[begidx], color='g', label='Mimas (' + r'$E_0=$' + str("{:.3e}".format(E_Mim[begidx])) + ' J)')
	start, end = ax3.get_xlim()
	ax3.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	uf.modify_yaxis_label(ax3, ymax)
	ax3.legend(loc=0, prop={'size':14})
	return fig, ax3

def plot_E_Tot(fig, begidx, endidx, t, E_Tot):
	ax4 = plt.subplot2grid((18,22),(9, 11), rowspan=9, colspan=11)
	ymax = np.amax(np.absolute(E_Tot[begidx:endidx+1] - E_Tot[begidx])) * 1.1
	fig, ax4 = pa.title_and_axes(fig, ax4, r'$\Delta E$' + ' for all bodies', 't  [days]', r'$E-E_0$' + ' [J]', t[begidx], t[endidx], -ymax, ymax)
	ax4.plot(t[begidx:endidx+1], E_Tot[begidx:endidx+1] - E_Tot[begidx], color='c', label='sum of all bodies (' + r'$E_0=$' + str("{:.3e}".format(E_Tot[begidx])) + ' J)')
	start, end = ax4.get_xlim()
	ax4.set_xticks( [start, start + ((end - start) / 3), start + ((end - start) * 2 / 3), end])
	ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	uf.modify_yaxis_label(ax4, ymax)
	ax4.legend(loc=0, prop={'size':14})
	return fig, ax4

def energy_plots(bodies, i, j, begidx, endidx, energy_matrix, t, black=False):
	print('Making energy plots...')
	for m in range(len(begidx)):
		if (black == True):
			fig = plt.figure(facecolor='black', figsize=(22,18))
		else:
			fig = plt.figure(figsize=(22,18))
		# aerial view of close encounter
		fig, ax1 = plot_E_ij(fig, begidx[m], endidx[m], t, energy_matrix, bodies, i, j)
		fig, ax2 = plot_E_al(fig, begidx[m], endidx[m], t, energy_matrix, bodies, i, j) 
		fig, ax3 = plot_E_Mim(fig, begidx[m], endidx[m], t, energy_matrix[len(bodies)])
		fig, ax4 = plot_E_Tot(fig, begidx[m], endidx[m], t, energy_matrix[len(bodies) + 1])
		plt.tight_layout()
		ax1.figure.savefig('energy_during_close_encounter_' + str(m) + '_between_' + str(bodies[i]) + '_and_' + str(bodies[j]) + '.png', dpi=100)	
		plt.clf()
		#uf.print_progress(m + 1, len(begidx))

# WRITE TO FILE ----------------------------------------------------------------------------------------------------------------------------------------------

def stats_for_color_maps(i, t0, tf, x_i, x_j, y_i, y_j, delta_L_i, delta_L_j, prox_i, prox_j, xj, yj, hr, delta_E3, delta_E4, E3, E4, lmda_Mim_prox, lmda_i_prox, delta_s_m, outcome):
	with open('stats' + str(i) + '.txt', 'a') as myfile:
		myfile.write("%.2f" % t0)
		myfile.write("\t\t%.2f" % tf)
		myfile.write("\t\t%.6f" % x_i)
		myfile.write("\t\t%.6f" % x_j)
		myfile.write("\t\t%.6f" % y_i)
		myfile.write("\t\t%.6f" % y_j)
		myfile.write("\t\t%.5e" % delta_L_i)
		myfile.write("\t\t%.5e" % delta_L_j)
		myfile.write("\t\t%.6f" % prox_i)
		myfile.write("\t\t%.6f" % prox_j)
		myfile.write("\t\t%.8f" % xj)
		myfile.write("\t\t%.8f" % yj)
		myfile.write("\t\t%.8f" % delta_E3)
		myfile.write("\t\t%.8f" % delta_E4)
		myfile.write("\t\t%.8f" % E3)
		myfile.write("\t\t%.8f" % E4)
		myfile.write("\t\t%.8f" % lmda_Mim_prox)
		myfile.write("\t\t%.8f" % lmda_i_prox)
		myfile.write("\t\t%.8f" % delta_s_m)
		myfile.write("\t\t%.0f" % outcome)
		myfile.write("\t\t%.8f" % hr + "\n")
	with open('stats_tot.txt', 'a') as myfile:
		myfile.write("%.2f" % t0)
		myfile.write("\t\t%.2f" % tf)
		myfile.write("\t\t%.6f" % x_i)
		myfile.write("\t\t%.6f" % x_j)
		myfile.write("\t\t%.6f" % y_i)
		myfile.write("\t\t%.6f" % y_j)
		myfile.write("\t\t%.5e" % delta_L_i)
		myfile.write("\t\t%.5e" % delta_L_j)
		myfile.write("\t\t%.6f" % prox_i)
		myfile.write("\t\t%.6f" % prox_j)
		myfile.write("\t\t%.8f" % xj)
		myfile.write("\t\t%.8f" % yj)
		myfile.write("\t\t%.8f" % delta_E3)
		myfile.write("\t\t%.8f" % delta_E4)
		myfile.write("\t\t%.8f" % E3)
		myfile.write("\t\t%.8f" % E4)
		myfile.write("\t\t%.8f" % lmda_Mim_prox)
		myfile.write("\t\t%.8f" % lmda_i_prox)
		myfile.write("\t\t%.8f" % delta_s_m)
		myfile.write("\t\t%.0f" % outcome)
		myfile.write("\t\t%.8f" % hr + "\n")

def simulation_results(pos_delta_s, neg_delta_s, zer_delta_s):
	with open('simulation_results.txt', 'a') as myfile:
		myfile.write("%.0f" % (len(pos_delta_s) + len(neg_delta_s) + len(zer_delta_s)))
		myfile.write("\t\t\t%.0f" % len(pos_delta_s))
		myfile.write("\t\t\t%.0f" % len(neg_delta_s))
		myfile.write("\t\t\t%.0f" % len(zer_delta_s) + "\n")


# MAIN ----------------------------------------------------------------------------------------------------------------------------------------

def main(plots=False, energy_change_plots=False):
	# importing data from simulation files
	M_Sat = bd.M_Sat()
	M_Mim = bd.M_Mim()
	bodies = bd.pocket_body_list()
	masses = bd.pocket_mass_list()
	hill_radii, body_radii = get_hill_and_body_radii(M_Sat, bodies, masses)
	t, a_cer, lambda_cer = dl.exact_cer014data('CERFXGEO.txt')
	data_points = len(t)
	x_Mim, y_Mim, z_Mim = dl.xyz_data('MIMAS')
	u_Mim, v_Mim, w_Mim = dl.uvw_data('MIMAS')
	tm, lmda_Mim = dl.geo04data('MIMAS')
	x, y, z = get_xyz_data(bodies, data_points)
	u, v, w = get_uvw_data(bodies, data_points)
	a, e, lmda, curlypi = get_geo_data(bodies, data_points)

	# computing data from imported data
	energy_matrix = compute_energy_for_all_bodies(bodies, masses, M_Sat, M_Mim, x_Mim, y_Mim, z_Mim, u_Mim, v_Mim, w_Mim, x, y, z, u, v, w) 
	s_matrix = compute_phase_space_distance_for_all_bodies(bodies, a, a_cer, lmda, lambda_cer)

	# initializing arrays to count outcomes 
	pos_delta_s = []
	neg_delta_s = []
	zer_delta_s = []

	for i in range(len(bodies)):
		for j in range(len(bodies)):
			if(i != j):
				begidx, endidx, some_encounter = compare_space_between(bodies, x, y, t, hill_radii, body_radii, i, j)
				if(some_encounter == True):
					delta_L_i, delta_L_j, delta_E3, delta_E4 = initialize_array_set(len(begidx))
					for m in range(len(begidx)):
						x_f, y_f = xy_path(i, j, begidx[m], endidx[m], x, y)
						delta_L_i_m, delta_L_j_m, L_i, L_j = angular_momentum(i, j, begidx[m], endidx[m], t, x, y, z, u, v, w, masses)
						encounter_length = len(delta_L_i_m)
						for k in range(encounter_length):
							delta_L_i[m][k] = delta_L_i_m[k]
							delta_L_j[m][k] = delta_L_j_m[k]
						delta_E3_m, delta_E4_m, E3, E4 = energy(i, j, begidx[m], endidx[m], energy_matrix)
						for k in range(encounter_length):
							delta_E3[m][k] = delta_E3_m[k]
							delta_E4[m][k] = delta_E4_m[k]
						proximity_i, proximity_j, xj, yj, lmda_Mim_prox, lmda_i_prox = stats_at_closest_approach(i, j, begidx[m], endidx[m], lmda_Mim, lmda, curlypi, x, y, x_f, y_f)
						pos_delta_s, neg_delta_s, zer_delta_s, delta_s_m, outcome = s_outcome(masses, s_matrix, i, j, begidx[m], endidx[m], pos_delta_s, neg_delta_s, zer_delta_s)
						stats_for_color_maps(i, t[begidx[m]], t[endidx[m]], (a[i][begidx[m]] - a_cer[begidx[m]]) * 1.0E-03, (a[j][begidx[m]] - a_cer[begidx[m]]) * 1.0E-03, 
							((lmda[i][begidx[m]] - lambda_cer[begidx[m]] + 180) % 360) - 180, ((lmda[j][begidx[m]] - lambda_cer[begidx[m]] + 180) % 360) - 180, delta_L_i[m][encounter_length-1], delta_L_j[m][encounter_length-1], 
							proximity_i, proximity_j, xj, yj, hill_radii[i], delta_E3[m][encounter_length-1], delta_E4[m][encounter_length-1], E3[endidx[m]], E4[endidx[m]], lmda_Mim_prox, lmda_i_prox, delta_s_m, outcome)
					if(plots == True):
						close_encounter_plots_for_paper(bodies, i, j, a, e, lmda, curlypi, t, a_cer, lambda_cer, begidx, endidx, hill_radii, body_radii, x, y, z, u, v, w, masses, M_Sat, M_Mim, x_Mim, y_Mim, z_Mim, x_f, y_f, delta_L_i, delta_L_j, L_i, L_j, delta_E3, delta_E4, E3, E4, s_matrix, lmda_Mim)
						#close_encounter_plots(bodies, i, j, a, e, lmda, curlypi, t, a_cer, lambda_cer, begidx, endidx, hill_radii, body_radii, x, y, z, u, v, w, masses, M_Sat, M_Mim, x_Mim, y_Mim, z_Mim, x_f, y_f, delta_L_i, delta_L_j, L_i, L_j, delta_E3, delta_E4, E3, E4, s_matrix, lmda_Mim)
					if(energy_change_plots == True):
						energy_plots(bodies, i, j, begidx, endidx, energy_matrix, t)
	print('Total number of close encounters: ' + str(len(pos_delta_s) + len(neg_delta_s) + len(zer_delta_s)))					
	print('Close encounters where both bodies move out: ' + str(len(pos_delta_s)))
	print('Close encounters where both bodies move in: '  + str(len(neg_delta_s)))
	print('Close encounters where one body moves in and the other moves out: ' + str(len(zer_delta_s)))
	#simulation_results(pos_delta_s, neg_delta_s, zer_delta_s)
# RUN --------------------------------------------------------------------------------

main(plots=True, energy_change_plots=False)	
