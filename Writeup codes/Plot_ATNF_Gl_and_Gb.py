'''
Plot_ATNF_Gl_and_Gb.py
Plot the Galactic longitude in degrees (Gl) and Galactic latitude in degrees (Gb) of pulsars from the ATNF pulsar catalogue.
Use the same variable names as the catalogue to avoid added confusion.
https://www.atnf.csiro.au/research/pulsar/psrcat/
'''

import csv
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u

#----- Define variables -----
csv_filename = "20240911_Gl_and_Gb.csv"
fig_filename = "ATNF_galactic_pulsar_distribution"

got_headers = False
got_units   = False

headers = []
units   = []
Gl      = []
Gb      = []

i_Gl = 0 # Column numbers for data, to be determined from the CSV file.
i_Gb = 0


#----- Read CSV file and extract data -----
with open( csv_filename ) as csv_file:
	
	data = csv.reader( csv_file )
	
	for row in data:
		
		if row == []:
			continue
		if not ";" in row[0]:
			continue
		
		row = row[0].split( ";" )
		
		if not got_headers:
			headers = row
			got_headers = True
			i_Gl = headers.index( "Gl" )
			i_Gb = headers.index( "Gb" )
		elif not got_units:
			units = row
			got_units = True
		else:
			Gl.append( float( row[ i_Gl ] ) )
			Gb.append( float( row[ i_Gb ] ) )


#----- Plot data -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111, projection="aitoff" ) # 'projection="aitoff"' makes a global, not rectangular, graph.

coordinates = coord.SkyCoord( Gl, Gb, frame="galactic", unit=u.deg )

ax.scatter(
		coordinates.l.wrap_at( "180d" ).radian,
		coordinates.b.radian,
		s=0.5,
		color="k",
		marker="."
		)

ax.grid( True )

fig.tight_layout()
plt.savefig( "../Figures/" + "20240911-ATNF-galactic-pulsar-distribution" )