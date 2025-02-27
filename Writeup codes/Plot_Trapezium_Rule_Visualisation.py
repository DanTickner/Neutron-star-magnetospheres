'''
Plot_Trapezium_Rule_Visualisation.py

Integrate an arbitrary function f(x) between x=a and x=b by the trapezium rule,
using N integration steps.
Plot the function and the trapezia in the same figure, to demonstrate how the
trapezium rule approximates the shape of the region enclosed by the function.
'''
import matplotlib.pyplot as plt
import numpy             as np


#----- Define function to plot -----
def f(x):
	return np.sin(x)


#----- Define variables -----
x_min = 0.25*np.pi
x_max = 0.65*np.pi
a     = 0.3 *np.pi		# Value to integrate from
b     = 0.6 *np.pi		# Value to integrate to
n     = 100				# Number of points to plot
N     = 4				# Number of integration steps

x = np.linspace( x_min, x_max, n )
y = f( x )

h = ( b - a ) / N
x_trapezium = [ a + i * h for i in range(N+1) ]
y_trapezium = [ f(x) for x in x_trapezium ]


#----- Setup graphs -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":20 })


#----- Plot graph with one step -----
fig1 = plt.figure( figsize=(7,5) )
ax1  = fig1.add_subplot( 111 )

ax1.plot( x, y )

ax1.plot( [ a, b ], [ f(a), f(b) ], "-o", color="black" )
ax1.plot( [ a, a ], [ 0   , f(a) ], "--", color="black"  )
ax1.plot( [ b, b ], [ 0   , f(b) ], "--", color="black"  )

ax1.set_xlim( x_min, x_max )
ax1.set_ylim( 0.5, max(y)*1.1 )

xticks_values = [ a, b ]
xticks_labels = [ r"$a$", r"$b$" ]
yticks_values = [ f(a), f(b) ]
yticks_labels = [ f"$f(a)$", r"$f(b)$" ]
plt.xticks( xticks_values, xticks_labels )
plt.yticks( yticks_values, yticks_labels )

section = np.linspace( a, b, 2 )
plt.fill_between( section, f(section), color="lightgrey" )

fig1.tight_layout()
plt.savefig( "Trapezium_Rule_Visualisation_1" )


#----- Plot graph with N steps -----
fig2 = plt.figure( figsize=(7,5) )
ax2  = fig2.add_subplot( 111 )

ax2.plot( x, y )

ax2.plot( x_trapezium, y_trapezium, "-o", color="black" )

for i in range( N+1 ):
	ax2.plot( [ x_trapezium[i], x_trapezium[i] ], [ 0, y_trapezium[i] ], "--", color="black" )

for i in range( N ):
	section = np.linspace( x_trapezium[i], x_trapezium[i+1], 2 )
	plt.fill_between( section, f(section), color="lightgrey" )

ax2.set_xlim( x_min, x_max )
ax2.set_ylim( 0.5, max(y)*1.1 )

xticks_values = [ a, b ]
xticks_labels = [ r"$a$", r"$b$" ]
yticks_values = [ f(a), f(b) ]
yticks_labels = [ f"$f(a)$", r"$f(b)$" ]
plt.xticks( xticks_values, xticks_labels )
plt.yticks( yticks_values, yticks_labels )

fig2.tight_layout()
plt.savefig( "Trapezium_Rule_Visualisation_N" )