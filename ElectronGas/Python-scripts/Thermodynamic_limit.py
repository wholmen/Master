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


plot(Ns,E,label="dE")

legend(loc="best")
title("Approaching the thermodynamic limit")
xlabel("Ns")
ylabel("E")
savefig("../Results/Figures/Thermodynamic_limit.png")

show()