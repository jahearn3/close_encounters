# Tricolor Phase Space Plot
#
# Author: Joseph A'Hearn
# Created 12/13/2018
#
# This program gathers data from the suite of simulations
#   and makes phase space plots colored by outcome

import numpy as np 
import sys 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import plot_assistant as pa 


def folder_name(i, third_body):
	if(third_body == True):
		if(i == 0):
			folder = 'alphecca'
		elif(i < 10):
			folder = 'suite_with_massive_body/suite000' + str(i)
		elif(i < 100):
			folder = 'suite_with_massive_body/suite00' + str(i)
		elif(i < 1000):
			folder = 'suite_with_massive_body/suite0' + str(i)
		else:
			folder = 'suite_with_massive_body/suite' + str(i)
	else:
		if(i == 0):
			folder = 'suite_without_massive_body/suitebalphecca'
		elif(i < 10):
			folder = 'suite_without_massive_body/suiteb000' + str(i)
		elif(i < 100):
			folder = 'suite_without_massive_body/suiteb00' + str(i)
		elif(i < 1000):
			folder = 'suite_without_massive_body/suiteb0' + str(i)
		else:
			folder = 'suite_without_massive_body/suiteb' + str(i)
	return folder

def load_data(i, folder):
	data = np.loadtxt(folder + '/stats_tot.txt', skiprows=1)
	x = data[:,2] 
	y = data[:,4] 
	outcome = data[:,19]
	return x, y, outcome

def import_data(N_simulations):
	x0 = []
	y0 = []
	x1 = []
	y1 = []
	x2 = []
	y2 = []
	for i in range(N_simulations):
		folder = folder_name(i, True)
		xa, ya, oa = load_data(i, folder)
		for j in range(len(xa)):
			if(oa[j] == 0):
				x0.append(xa[j])
				y0.append(ya[j])
			elif(oa[j] == 1):
				x1.append(xa[j])
				y1.append(ya[j])
			elif(oa[j] == 2):
				x2.append(xa[j])
				y2.append(ya[j])

		folder = folder_name(i, False)
		xb, yb, ob = load_data(i, folder)
		for j in range(len(xb)):
			if(ob[j] == 0):
				x0.append(xb[j])
				y0.append(yb[j])
			elif(ob[j] == 1):
				x1.append(xb[j])
				y1.append(yb[j])
			elif(ob[j] == 2):
				x2.append(xb[j])
				y2.append(yb[j])
	return x0, y0, x1, y1, x2, y2

def plot_data(x0, y0, x1, y1, x2, y2, black=False):
	print('Making plots...')
	nrows = 14
	ncols = 12
	rs = 7
	cs = 6

	xspan = 40
	yspan = 32

	fig, ax = pa.initialize_plot()
	fig, ax = pa.title_and_axes(fig, ax, 'title', r'$a - a_{CER}$' + ' [km]', r'$\lambda - \lambda_{CER}$' + ' [deg]', xmin=-xspan, xmax=xspan, ymin=-yspan, ymax=yspan)
	pa.save_and_clear_plot(fig, ax, filename='figure')

	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(ncols,nrows))
	else:
		fig = plt.figure(figsize=(ncols,nrows))

	ax1 = plt.subplot2grid((nrows,ncols),(0, 0), rowspan=rs, colspan=cs)
	fig, ax1 = pa.title_and_axes(fig, ax1, 'title', r'$a - a_{CER}$' + ' [km]', r'$\lambda - \lambda_{CER}$' + ' [deg]', xmin=-xspan, xmax=xspan, ymin=-yspan, ymax=yspan)
	ax1.set_title('Phase space positions in outcomes\nwhere one moves in and one moves out', fontsize=20)
	# contour of the edge of the corotation site
	fringe = patches.Ellipse(xy=[0,0], width=74, height=60, angle=0, color='g', fill=False, linestyle='dashed')
	plt.draw()
	ax1.add_patch(fringe)
	ax1.scatter(x0, y0, c='r', s=10)

	ax2 = plt.subplot2grid((nrows,ncols),(0, 6), rowspan=rs, colspan=cs)
	fig, ax2 = pa.title_and_axes(fig, ax2, 'title', r'$a - a_{CER}$' + ' [km]', r'$\lambda - \lambda_{CER}$' + ' [deg]', xmin=-xspan, xmax=xspan, ymin=-yspan, ymax=yspan)
	ax2.set_title('Phase space positions in outcomes\nwhere both move in', fontsize=20)
	# contour of the edge of the corotation site
	fringe = patches.Ellipse(xy=[0,0], width=74, height=60, angle=0, color='g', fill=False, linestyle='dashed')
	plt.draw()
	ax2.add_patch(fringe)
	ax2.scatter(x1, y1, c='k', s=10)

	ax3 = plt.subplot2grid((nrows,ncols),(7, 0), rowspan=rs, colspan=cs)
	fig, ax3 = pa.title_and_axes(fig, ax3, 'title', r'$a - a_{CER}$' + ' [km]', r'$\lambda - \lambda_{CER}$' + ' [deg]', xmin=-xspan, xmax=xspan, ymin=-yspan, ymax=yspan)
	ax3.set_title('Phase space positions in outcomes\nwhere both move out', fontsize=20)
	# contour of the edge of the corotation site
	fringe = patches.Ellipse(xy=[0,0], width=74, height=60, angle=0, color='g', fill=False, linestyle='dashed')
	plt.draw()
	ax3.add_patch(fringe)
	ax3.scatter(x2, y2, c='c', s=10)

	ax4 = plt.subplot2grid((nrows,ncols),(7, 6), rowspan=rs, colspan=cs)
	fig, ax4 = pa.title_and_axes(fig, ax4, 'title', r'$a - a_{CER}$' + ' [km]', r'$\lambda - \lambda_{CER}$' + ' [deg]', xmin=-xspan, xmax=xspan, ymin=-yspan, ymax=yspan)
	ax4.set_title('Phase space positions\nfor all outcomes', fontsize=20)
	# contour of the edge of the corotation site
	fringe = patches.Ellipse(xy=[0,0], width=74, height=60, angle=0, color='g', fill=False, linestyle='dashed')
	plt.draw()
	ax4.add_patch(fringe)
	ax4.scatter(x0, y0, c='r', s=5)
	ax4.scatter(x1, y1, c='k', s=5)
	ax4.scatter(x2, y2, c='c', s=5)

	plt.tight_layout()
	ax1.figure.savefig('encounters_in_phase_space.png', dpi=300)	
	plt.clf()

def main():
	if(len(sys.argv) < 2):
		N_simulations = 1
	else:
		N_simulations = int(sys.argv[1])

	x0, y0, x1, y1, x2, y2 = import_data(N_simulations)
	plot_data(x0, y0, x1, y1, x2, y2)

main()