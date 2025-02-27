'''
Plot_ATNF_P0_and_P1_v02
Plot the period in seconds (P0) and period derivative in seconds per second (P1)
of pulsars from the ATNF pulsar catalogue.
Use the same variable names as the catalogue to avoid added confusion.
https://www.atnf.csiro.au/research/pulsar/psrcat/

V02: Include "Type" parameter in order to make separate data labels in the P-Pdot diagram.
Note that "Type" is called "PSR".
'''

import csv
import math # use math.nan for missing values
import matplotlib.pyplot as plt

#----- Define variables -----
csv_filename = "20240912_P0_and_P1_and_Type_and_Assoc_and_Binary.csv"
fig_filename = "20240912-ATNF-P-Pdot-02"

got_headers = False
got_units   = False

headers   = []
units     = []
P0        = []
P1        = []
is_binary = []
is_AXP    = []


i_P0  = 0 # Column numbers for data, to be determined from the CSV file.
i_P1  = 0
i_PSR = 0

n = 3.0 # Braking index. n=3 for a dipole.
alpha_angle = math.pi/2.0 # Angle of magnetic to rotation axis.
R = 10 * 1e3 # Radius in metres
M = 1.4 * 1.989e30 # Mass in kilograms.
I = 0.4 * M * R*R # Moment of inertia of a uniform solid sphere.
c = 299792458.0
mu_0 = 1.25663706127e-6 # Permeability of free space in N A^-2.

tau_yr   = [ 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10 ] # Characteristic age in years to plot.
tau_text = [ "1 kyr", "10 kyr", "100 kyr", "1 Myr", "10 Myr", "100 Myr", "1 Gyr", "10 Gyr" ] # Text for graph labels.
tau_s    = [ tau * 3600 * 24 * 365 for tau in tau_yr ] # For calculations.

B_T    = [ 1e6, 1e7, 1e8, 1e9, 1e10 ] # Polar magnetic field in tesla to plot.
B_text = [ r"$10^6$ T", r"$10^7$ T", r"$10^8$ T", r"$10^9$ T", r"$10^{10}$ T" ]

text_rotation_tau = 12 # Angle at which the tau labels should be written. Trial-and-error.
text_vertical_multiplier_tau = 2
text_horizontal_offset_tau = 0

text_rotation_B = -12
text_vertical_multiplier_B = 2
text_horizontal_offset_B = 0

P1_for_line_of_constant_tau = [ [ 0, 0 ] for tau in tau_s ] # For the lines of constant tau. To be calculated.
P1_for_line_of_constant_B   = [ [ 0, 0 ] for B   in B_T   ] # For the lines of constant B. To be calculated.

xmin = 5e-4
xmax = 5e1


#----- Read CSV file and extract data -----
with open( csv_filename ) as csv_file:
	
	data = csv.reader( csv_file )
	
	for row in data:
		
		if row == []:
			continue
		if not ";" in row[0]:
			continue
		
		# "ASSOC" uses commas and this produces multi-element lists when reading with Python.
		# Join all elements into single string first, then split into list with semicolons as separators.
		row = "".join(row).split( ";" )
		
		if not got_headers:
			headers = row
			got_headers = True
			i_P0     = headers.index( "P0"     )
			i_P1     = headers.index( "P1"     )
			i_BINARY = headers.index( "BINARY" )
			i_PSR    = headers.index( "PSR"    )
		
		elif not got_units:
			units = row
			got_units = True
		
		else:
			
			try:
				P0.append( float( row[ i_P0 ] ) )
			except ValueError:
				P0.append( math.nan )
			
			try:
				P1.append( float( row[ i_P1 ] ) )
			except ValueError:
				P1.append( math.nan )
			
			is_binary.append( row[ i_BINARY ] != "*" )
			
			is_AXP.append( "AXP" in row[ i_PSR ] )
			

#----- Print brief data summary to screen -----
print( f"Number of pulsars:\t{ len( P0 ) }" )
print( f"Number of binaries:\t{ is_binary.count( True ) }" )
print( f"Number of AXPs:\t{ is_AXP.count( True ) }" )

num_binary_and_AXP = 0
for i in range( len( P0 ) ):
	if is_binary[i] and is_AXP[i]:
		num_binary_and_AXP += 1

print( f"Number of binaries which are also AXPs:\t{ num_binary_and_AXP }" )
# If the above is zero, as it should, we can label them separately on the graph without issue.


#----- Reorganise into separate lists for normal, AXP and binary -----
# We could have done this as we read the CSV file, but it was better to check for crossovers first.

normal_P0 = []
normal_P1 = []
binary_P0 = []
binary_P1 = []
AXP_P0    = []
AXP_P1    = []

for i in range( len( P0 ) ):
	if is_AXP[i]:
		AXP_P0.append( P0[i] )
		AXP_P1.append( P1[i] )
	elif is_binary[i]:
		binary_P0.append( P0[i] )
		binary_P1.append( P1[i] )
	else:
		normal_P0.append( P0[i] )
		normal_P1.append( P1[i] )


#----- Plot data -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.set_xlabel( "Period (s)" )
ax.set_ylabel( r"Period derivative (s s$^{-1}$)" )

ax.set_xscale( "log" )
ax.set_yscale( "log" )

ax.scatter(
	AXP_P0,
	AXP_P1,
	s=2,
	color="b",
	marker=".",
	label="AXP"
	)

ax.scatter(
	binary_P0,
	binary_P1,
	s=2,
	color="r",
	marker=".",
	label="Binary"
	)

ax.scatter(
	normal_P0,
	normal_P1,
	s=2,
	color="k",
	marker=".",
	label="Normal"
	)



#----- Add lines of constant dipole age -----
P_min = min(P0)
P_max = max(P0)

for i, tau in enumerate( tau_s ):
	
	P1_for_line_of_constant_tau[i][0] = 1.0/(n-1.0) * P_min / tau
	P1_for_line_of_constant_tau[i][1] = 1.0/(n-1.0) * P_max / tau
	
	ax.plot(
		[ P_min, P_max ],
		P1_for_line_of_constant_tau[i],
		color="grey",
		lw=0.5,
		ls="--"
		)

	plt.text(
		P_min + text_horizontal_offset_tau,
		P1_for_line_of_constant_tau[i][0] * text_vertical_multiplier_tau,
		tau_text[i],
		ha="center",
		va="center",
		rotation = text_rotation_tau,
		size=14,
		color="grey"
		)

#----- Add lines of constant polar magnetic field -----
for i, B in enumerate( B_T ):
	
	P1_for_line_of_constant_B[i][0] = 16.0*math.pi**3/(6.0*mu_0*c**3) * R**6 * math.sin(alpha_angle)**2 * B**2 / ( I * P_min )
	P1_for_line_of_constant_B[i][1] = 16.0*math.pi**3/(6.0*mu_0*c**3) * R**6 * math.sin(alpha_angle)**2 * B**2 / ( I * P_max )
	
	ax.plot(
		[ P_min, P_max ],
		P1_for_line_of_constant_B[i],
		color="darkgrey",
		lw=0.5,
		ls=":"
		)
	
	plt.text(
		P_max + text_horizontal_offset_B,
		P1_for_line_of_constant_B[i][1] * text_vertical_multiplier_B,
		B_text[i],
		ha="center",
		va="center",
		rotation = text_rotation_B,
		size=14,
		color="darkgrey"
		)



#----- Final formatting and save figure -----
ax.set_xlim( xmin, xmax )

# Specify location of legend so it doesn't coincide with a magnetic field line label.
# This uses the graph coordinates of the lower left corner,
# so that [0,0] is the origin and [1,1] is the upper right regardless of the values plotted.
ax.legend( loc=[ 0.53, 0.03 ] )

fig.tight_layout()
plt.savefig( "../Figures/" + fig_filename )
print( f"Saved figure: ../Figures/{ fig_filename }" )