from pylab import *

infile = open("../Results/ThermodynamicLimit_Nh14_n0.2.txt",'r')

infile.readline(); infile.readline(); infile.readline(); # Skipping three first lines filled with text

E = []; Ns = []; 
for line in infile:
	splitline = line.split();
	Ns.append( float(splitline[1]))
	E.append(float(splitline[-2]))


plot(Ns,E,color="orange",label=r'$\Delta E$',linewidth=2.5)

legend(loc="best")
title("Approaching the thermodynamic limit")
xlabel(r'$N_s$')
ylabel(r'$\Delta E$ [MeV] ',rotation=90)
savefig("../Results/Figures/Thermodynamic_limit.png")

show()