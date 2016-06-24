from pylab import *

infile = open("../Results/Comparison_FCI_MBPT_CCD.txt",'r')

for i in range(5):
	infile.readline()

g = []; ReferenceEnergy = []; FCI = []; CI = []; MBPT2 = []; MBPT2exact = []; MBPT3 = []; MBPT4 = []; CCD = [];

for line in infile:
	splitline = line.split()

	g.append( float(splitline[0]) )
	ReferenceEnergy.append( float(splitline[1]) )
	FCI.append( float(splitline[2]) )
	CI.append( float(splitline[3]) )
	MBPT2.append( float(splitline[4]) )
	MBPT2exact.append( float(splitline[5]) )
	MBPT3.append( float(splitline[6]) )
	MBPT4.append( float(splitline[7]) )
	CCD.append( float(splitline[8]) )


# Generate a comparison of FCI, CI, MBPT2, MBPT3, MBPT4, CCD total energy
plot(g,array(ReferenceEnergy)+array(FCI),   color="blue"  ,label="FCI")
plot(g,array(ReferenceEnergy)+array(CI),    color="green" ,label="CI")
plot(g,array(ReferenceEnergy)+array(MBPT2), color="red"   ,label="MBPT2")
plot(g,array(ReferenceEnergy)+array(MBPT3), color="orange",label="MBPT3")
plot(g,array(ReferenceEnergy)+array(MBPT4), color="yellow",label="MBPT4")
plot(g,array(ReferenceEnergy)+array(CCD),   color="purple",label="CCD")

legend(loc="best")
title("A comparison of Ground state energy for 4p4h Pairing model \n calculated by various methods. Relaxing factor on CCD is 0.5")
xlabel("Interaction strength, g")
ylabel("Ground state Energy, dE")
savefig("../Results/Figures/Pairing4p4h_CompareE_AllMethods.png")

show()


# Generate a comparison of FCI, CI, MBPT2, MBPT3, MBPT4, CCD corrolation energy
plot(g,FCI,   color="blue"  ,label="FCI")
plot(g,CI,    color="green" ,label="CI")
plot(g,MBPT2, color="red"   ,label="MBPT2")
plot(g,MBPT3, color="orange",label="MBPT3")
plot(g,MBPT4, color="yellow",label="MBPT4")
plot(g,CCD,   color="purple",label="CCD")

legend(loc="best")
title("A comparison of corrolation energy for 4p4h Pairing model \n calculated by various methods. Relaxing factor on CCD is 0.5")
xlabel("Interaction strength, g")
ylabel("Corrolation Energy, dE")
savefig("../Results/Figures/Pairing4p4h_CompareDE_AllMethods.png")

show()