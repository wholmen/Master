from pylab import *

infile = open("../Results/Solver_Nshells3to35Nfilled2.txt",'r')

infile.readline(); infile.readline(); infile.readline(); # Skipping three first lines filled with text

E = []; Ns = []; 
for line in infile:
	splitline = line.split();
	Ns.append( float(splitline[1]))
	E.append(float(splitline[-1]))

infile = open("../Results/Solver_Nshells36to40Nfilled2.txt",'r')

infile.readline(); infile.readline(); infile.readline(); # Skipping three first lines filled with text

for line in infile:
	splitline = line.split();
	Ns.append( float(splitline[1]))
	E.append(float(splitline[-1]))


plot(Ns,E,color="orange",label=r'$\Delta E$',linewidth=1.5)

legend(loc="best")
title("Approaching the thermodynamic limit")
xlabel(r'$N_s$')
ylabel(r'$\Delta E$',rotation=0)
savefig("../Results/Figures/Thermodynamic_limit.png")

show()