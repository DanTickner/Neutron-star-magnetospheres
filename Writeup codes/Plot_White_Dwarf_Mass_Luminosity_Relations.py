'''
Plot_White_Dwarf_Mass_Luminosity_Relations.py
Nauenberg 1972 https://ui.adsabs.harvard.edu/abs/1972ApJ...175..417N/abstract
Ryan-Norton 2010 Eq. (6.6) https://www.cambridge.org/gb/universitypress/subjects/physics/astrophysics/stellar-evolution-and-nucleosynthesis?format=HB&isbn=9780521196093
'''

import numpy as np
import matplotlib.pyplot as plt

#----- Define variables -----
M_min = 0.5 # Minimum mass in Solar units.
M_max = 1.5 # Maximum mass in Solar units.

n_points = 1000

R_odot = 6.9634e8 # Solar radius in metres.

fig_filepath = "../Images/20241003-White-Dwarf-Mass-Luminosity-Relations"

#----- Build arrays -----
M = np.linspace( M_min, M_max, n_points )

R_Nauenberg   = []
R_Ryan_Norton = []

for m in M:
	m_Ch = m / 1.44 # Mass in Chandrasekhar units.
	R_Nauenberg.append( 7.83e6 * np.sqrt( m_Ch**(-2.0/3.0) - m_Ch**(2.0/3.0) ) )
	R_Ryan_Norton.append( ( R_odot / 74.0 ) * m**(-1.0/3.0) )


#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.set_xlabel( r"Mass ($M_\odot$)" )
ax.set_ylabel( r"Radius (m)" )

ax.plot( M, R_Nauenberg  , label="Nauenberg"  , color="k", ls="-"  )
ax.plot( M, R_Ryan_Norton, label="Ryan-Norton", color="b", ls="--" )

ax.set_ylim( 0, max( max(R_Nauenberg), max(R_Ryan_Norton) )*1.05 )
ax.set_xlim( M_min, M_max )

ax.legend( loc="best" )

fig.tight_layout()

plt.savefig( fig_filepath )

print( f"Image saved:\t{ fig_filepath }" )