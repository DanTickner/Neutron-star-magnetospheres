'''
Plot_Pulsar_Period_as_Function_of_Time.py
Pulsar period as a function of time
'''

import numpy as np
import os
import csv
import matplotlib.pyplot as plt

#----- Define variables -----

P_0 = 33e-3 # Initial period in s
B   = 7.6e8 # Magnetic field in tesla.
alpha_angle = np.pi/2.0
M = 1.4 * 1.989e30 # Mass in kilograms.
R = 10 * 1e3 # Radius in metres.
I = 0.4 * M * R*R # Moment of inertia of a uniform solid sphere.
c = 299792458.0
mu_0 = 1.25663706127e-6 # Permeability of free space in N A^-2.

P    = [] # Period as a function of time, to be calculated.
Pdot = [] # Period derivative as a function of time, to be calculated.

tau_min = 0 # Minimum pulsar age to plot, in seconds.
tau_max = 1e6 * 365 * 24 * 3600 # Maximum pulsar age to plot, in seconds.
n_points = 100 # Number of time values to evaluate the expression for.

csv_folder = ""
csv_filename = "20240917_Crab_pulsar_P_Pdot_evolution.csv"
csv_filepath = os.path.join( csv_folder, csv_filename )

fig_folder = "../Images"
fig_filename = "20240917-Crab-pulsar-P-Pdot-evolution"
fig_filepath = os.path.join( fig_folder, fig_filename )


tau = np.linspace( tau_min, tau_max, n_points )


#----- Calculate P and Pdot as a function of time -----

factor = ( 16.0 * np.pi**3 / ( 3.0 * mu_0 * c**3 ) ) * ( B**2 * R**6 * np.sin(alpha_angle)**2 ) / ( I * P_0**2 )

for t in tau:
	P.append( P_0 * np.sqrt( 1.0 + factor * t ) )
	Pdot.append( P_0 * 0.5 * factor / np.sqrt( 1.0 + factor * t ) )


#----- Export to CSV in case further analysis needed -----
with open( csv_filename, "w", newline="" ) as csv_file:
	
	writer = csv.writer( csv_file )
	
	writer.writerow( [ "age (sec)", "P (sec)", "Pdot" ] )
	
	for i in range( len( tau ) ):
		
		writer.writerow( [ tau[i], P[i], Pdot[i] ] )


print( f"CSV file saved:\t{ csv_filepath }" )


#----- Plot on a graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.set_xlabel( "Period (s)" )
ax.set_ylabel( r"Period derivative (s s$^{-1}$)" )

ax.set_xscale( "log" )
ax.set_yscale( "log" )


ax.scatter(
	P[1:-1],
	Pdot[1:-1],
	s=5,
	color="grey",
	marker="*"
	)

ax.scatter(
	P[0],
	Pdot[0],
	s=20,
	color="r",
	marker="*"
	)

ax.scatter(
	P[-1],
	Pdot[-1],
	s=20,
	color="b",
	marker="*"
	)

fig.tight_layout()
plt.savefig( fig_filepath )
print( f"Saved figure: { fig_filepath }" )