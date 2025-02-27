# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 11:57:11 2024

@author: Dan

Plot_NS_Neutron_to_Proton_Ratio.py

Read the CSV file
NS_neutron_to_proton_ratio.csv
Which was created manually by running the code
Decimal_search_NS_neutron_number_density.py
for a few representative mass densities, and plot mass density vs neutron number density and neutron-proton ratio.
"""

import csv
import matplotlib.pyplot as plt

#----- Define variables -----
csv_filename = "NS_neutron_to_proton_ratio.csv"
fig_filename = "../Images/NS-neutron-to-proton-ratio.png"

rho        = [] # Mass density
n_n        = [] # Neutron number density
n_n_to_n_p = [] # Neutron-proton ratio


#----- Read CSV -----
with open( csv_filename, "r" ) as csv_file:
	data = csv.reader( csv_file )
	
	headers = next( data )
	
	for row in data:
		rho       .append( float( row[0] ) )
		n_n       .append( float( row[1] ) )
		n_n_to_n_p.append( float( row[2] ) )


#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax1 = fig.add_subplot( 111 )
ax2 = ax1.twinx()

ax1.plot( rho, n_n       , color="blue", ls="-"  )
ax2.plot( rho, n_n_to_n_p, color="red" , ls="--" )

ax1.set_ylabel( "Neutron number density $n_n$ (m$^{-3}$)" )
ax2.set_ylabel( "Neutron-proton ratio $n_n/n_p$" )
ax1.set_xlabel( "Mass density $\\rho$ (kg m$^{-3})$" )

fig.tight_layout()
plt.savefig( fig_filename )