'''
20250225_Aliasing_visualisation.py

Aliasing model
Parfrey PhD thesis, section 3.6
'''
import numpy as np
import matplotlib.pyplot as plt

#----- Define variables -----
N = 10 # Number of gridpoints
omega = 0.4

#----- Build lists -----
x = np.linspace( -np.pi, np.pi, N )

ftest = [ 1 for x_i in x ] # fake data to delete later

q_1 = np.pi
q_2 = 2.0 * np.pi



mode_1 = np.cos( ( omega + q_1 * N ) * x )
mode_2 = np.cos( ( omega + q_2 * N ) * x )

#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure( figsize=(8,5) )
ax  = fig.add_subplot( 111 )

ax.plot( x, ftest )

ax.plot( x, mode_1 )
ax.plot( x, mode_2 )

for n in range( N ):
	ax.scatter( [ x[n] ], [ 0 ], color="k" )