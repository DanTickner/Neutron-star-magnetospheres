'''
Geneate_magnetic_field_values_CSV_dipole.py

Evaluate the magnetic field of a dipole at given values of r and theta.
Calculate the correspoding Cartesian coordinates x and z, and the
corresponding magnetic field components in Cartesian components B_x and B_z.
Output these values to CSV file.
'''

import numpy as np
import csv


#----- Define functions -----
def B_x_function( r, theta ):
	# x-component of magnetic field of a dipole in terms of spherical polar coordinates (r,theta).
	return 3.0 * np.cos(theta) * np.sin(theta) * r**-3

def B_z_function( r, theta ):
	# z-component of magnetic field of a dipole in terms of spherical polar coordinates (r,theta).
	return ( 3.0 * np.cos(theta)**2 - 1.0 ) * r**-3

#----- Define variables -----
r_min          = 1
r_max          = 2
theta_min      = 0
theta_max      = np.pi
n_points_r     = 1000
n_points_theta = 1000
csv_filename   = "CSV/Magnetic_field_values_dipole_1e6.csv"


#----- Generate lists of gridpoints -----
r_list     = np.linspace( r_min    , r_max    , n_points_r     )
theta_list = np.linspace( theta_min, theta_max, n_points_theta )

x_grid   = np.zeros( shape=(n_points_r,n_points_theta) )
z_grid   = np.zeros( shape=(n_points_r,n_points_theta) )
B_x_grid = np.zeros( shape=(n_points_r,n_points_theta) )
B_z_grid = np.zeros( shape=(n_points_r,n_points_theta) )


#----- Create CSV file and output field values -----
with open( csv_filename, "w", newline="" ) as csv_file:
	writer = csv.writer( csv_file )
	
	writer.writerow( [ "r", "theta", "x", "z", "B_x", "B_z" ] )
	
	for r in r_list:
		for theta in theta_list:
			x = r * np.sin(theta)
			z = r * np.cos(theta)
			B_x = B_x_function( r, theta )
			B_z = B_z_function( r, theta )
			
			writer.writerow( [ r, theta, x, z, B_x, B_z ] )


#----- Finished -----
print( f"CSV file created:\t{ csv_filename }" )