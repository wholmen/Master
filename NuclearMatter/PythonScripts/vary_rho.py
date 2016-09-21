from pylab import *

infile = open("../Results/Solver_Nh14_varyRho_Ns20.txt",'r')
infile.readline(); infile.readline(); infile.readline(); # Skipping three first lines filled with text

n = []; Nh = []; E = []; Ns = []; 
for line in infile:
	splitline = line.split();
	n.append(float(splitline[0]))
	E.append(float(splitline[-1]))
	Ns.append(float(splitline[1]))
	Nh.append(float(splitline[2]))

plot(n,E,'-o',color="orange",label=r'$N_p = 14$',linewidth=1.5)


infile = open("../Results/Solver_Nh36_varyRho_Ns20.txt",'r')
infile.readline(); infile.readline(); infile.readline(); # Skipping three first lines filled with text

n = []; Nh = []; E = []; Ns = []; 
for line in infile:
	splitline = line.split();
	n.append(float(splitline[0]))
	E.append(float(splitline[-1]))
	Ns.append(float(splitline[1]))
	Nh.append(float(splitline[2]))

plot(n,E,'-^',color="red",label=r'$N_p = 36$',linewidth=1.5)

"""
infile = open("../Results/Solver_Nh54_varyRho_Ns20.txt",'r')
infile.readline(); infile.readline(); infile.readline(); # Skipping three first lines filled with text

n = []; Nh = []; E = []; Ns = []; 
for line in infile:
	splitline = line.split();
	n.append(float(splitline[0]))
	E.append(float(splitline[-1]))
	Ns.append(float(splitline[1]))
	Nh.append(float(splitline[2]))

plot(n,E,'-o',color="blue",label=r'$N_p = 54, $',linewidth=1.5)
"""

legend(loc="best")
title("A plot showing energy per particle as a function of density for Ns = %g"%Ns[0])
xlabel(r'$\rho$ [fm$-3$]')
ylabel(r'$E/A$ [MeV] ',rotation=90)
savefig("../Results/Figures/Vary_rho_Nh14.png")

show()