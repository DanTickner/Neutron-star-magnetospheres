# Neutron-star-magnetospheres
C++ code to produce time-dependent simulations of the magnetic field surrounding a twisted, rotating neutron star.


All information regarding theoretical background, mathematical methods, numerical intepretations, and code purposes and running instructions, is to be found within the attached thesis.

The three files which represent the numerical method are:
Header_Initial_Conditions_xx.h
Header_Time_Evolution_xx.h
Time_Evolution_xx.cpp

Results will be saved to the subfolders "CSV" and "Logs". If these folders do not exist, the code will run as normal but no saved files will appear.

The subfolder "Python" contains some Python codes which may be used to analyse the results of a given simulation. Of these, the most versatile are
Plot_History_Multiple.py
Plot_Profiles.py
whose results are explained in the Appendix of the thesis.

The subfolder "Figures" is not necessary but proves a useful location for saving the output figures from the Python codes. It is not necessary to add internal subfolders by date as I have done, but I found that it helps with organisation.

The remaining codes are examples demonstrating the numerical methods developed during the project, verifying calculations that have been carried out, testing the accuracy of the numerical simulation, and so on. They are somewhat disorganised and many are probably not useful, but are included for completeness and for reference. The most useful in terms of the thesis are in the subfolders
Test codes
Writeup codes

The codes within the folder
Encoded mathematical functions
should be ready-to-use methods that can be used for your own applications.

I do not claim that any of these extra example codes are up-to-date, or that the file structure is correct so that their dependent header files are in the right place - you will probably need to modify the
#include
lines within these codes to suit.
