'''
Test_time_and_rotation.py

Plot angular velocity, rotational period and total angle rotated as a function of time.
Include both code-units and SI units as axis labels.
'''


import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
SI_or_codeunits = "SI" # "SI" or "Code units". Toggle to give the graph you want.

csv_filename      = "../CSV/20240403_Test_time_and_rotation.csv"
fig_filename_base = "../Figures/20240403/20240403_Test_time_and_rotation"

if   SI_or_codeunits == "SI":
	fig_filename_base += "_SI"
	print( "Plotting in SI units." )
elif SI_or_codeunits == "Code units":
	fig_filename_base += "_codeunits"
	print( "Plotting in code units." )
else:
	raise ValueError( f"SI_or_codeunits = { SI_or_codeunits } not understood. Choose either 'SI' or 'Code units'." )

fig_filename_Omega = fig_filename_base + "_Omega"
fig_filename_P     = fig_filename_base + "_P"
fig_filename_angle = fig_filename_base + "_angle"

T_index  = []	# Number of values not known a priori.
T        = []
T_SI     = []
Omega    = []
Omega_SI = []
P        = []
P_SI     = []
angle    = []




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })




#----- Read CSV -----
with open( csv_filename, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	# Get correct columns
	i_T_index  = headers.index( "T_index"             )
	i_T        = headers.index( "T"                   )
	i_T_SI     = headers.index( "T_SI"                )
	i_Omega    = headers.index( "Omega"               )
	i_Omega_SI = headers.index( "Omega_SI"            )
	i_P        = headers.index( "P"                   )
	i_P_SI     = headers.index( "P_SI"                )
	i_angle    = headers.index( "star_rotation_angle" )
	
	for row in data:
		T_index .append( int  ( row[ i_T_index  ] ) )
		T       .append( float( row[ i_T        ] ) )
		T_SI    .append( float( row[ i_T_SI     ] ) )
		Omega   .append( float( row[ i_Omega    ] ) )
		Omega_SI.append( float( row[ i_Omega_SI ] ) )
		P       .append( float( row[ i_P        ] ) )
		P_SI    .append( float( row[ i_P_SI     ] ) )
		angle   .append( float( row[ i_angle    ] ) )




#----- Plot Omega as a function of T -----
fig_Omega = plt.figure( figsize=(8,5) )
ax_Omega  = fig_Omega.add_subplot( 111 )

if   SI_or_codeunits == "SI":
	ax_Omega.plot( T_SI, Omega_SI )
	ax_Omega.set_ylabel( r"Angular velocity $\Omega$ (rad s$^{-1}$)" )
	ax_Omega.set_xlabel( r"Time $T$ (s)" )
	ax_Omega.set_xlim( min(T_SI), max(T_SI) )
	ax_Omega.set_ylim( 0, max( Omega_SI ) * 1.05 )

elif SI_or_codeunits == "Code units":
	ax_Omega.plot( T, Omega )
	ax_Omega.set_ylabel( r"Angular velocity $\Omega$ (code units)" )
	ax_Omega.set_xlabel( r"Time $\tilde{T}$ (code units)" )
	ax_Omega.set_xlim( min(T), max(T) )
	ax_Omega.set_ylim( 0, max( Omega ) * 1.05 )

fig_Omega.tight_layout()
plt.savefig( fig_filename_Omega )
print( f"Saved figure:\t{ fig_filename_Omega }" )




#----- Plot P as a function of T -----
fig_P = plt.figure( figsize=(8,5) )
ax_P  = fig_P.add_subplot( 111 )

if   SI_or_codeunits == "SI":
	ax_P.plot( T_SI, P_SI )
	ax_P.set_ylabel( r"Rotation period $P$ (s)" )
	ax_P.set_xlabel( r"Time $T$ (s)" )
	ax_P.set_xlim( min(T_SI), max(T_SI) )
	ax_P.set_ylim( 0, max( P_SI ) * 1.05 )

elif SI_or_codeunits == "Code units":
	ax_P.plot( T, P )
	ax_P.set_ylabel( r"Rotation period $P$ (code units)" )
	ax_P.set_xlabel( r"Time $\tilde{T}$ (code units)" )
	ax_P.set_xlim( min(T), max(T) )
	ax_P.set_ylim( 0, max( P ) * 1.05 )

fig_P.tight_layout()
plt.savefig( fig_filename_P )
print( f"Saved figure:\t{ fig_filename_P }" )




#----- Plot angle as a function of T -----
fig_angle = plt.figure( figsize=(8,5) )
ax_angle  = fig_angle.add_subplot( 111 )

if   SI_or_codeunits == "SI":
	ax_angle.plot( T_SI, angle )
	ax_angle.set_ylabel( r"Total angle rotated (rad)" )
	ax_angle.set_xlabel( r"Time $T$ (s)" )
	ax_angle.set_xlim( min(T_SI), max(T_SI) )
	ax_angle.set_ylim( 0, max( angle ) * 1.05 )

elif SI_or_codeunits == "Code units":
	ax_angle.plot( T, angle )
	ax_angle.set_ylabel( r"Total angle rotated (rad)" )
	ax_angle.set_xlabel( r"Time $\tilde{T}$ (code units)" )
	ax_angle.set_xlim( min(T), max(T) )
	ax_angle.set_ylim( 0, max( angle ) * 1.05 )

fig_angle.tight_layout()
plt.savefig( fig_filename_angle )
print( f"Saved figure:\t{ fig_filename_angle }" )