# Master
I repository for my Master thesis on Coupled Cluster theory

The folder contains following subfolders: 

ElectronGasBlockStructure:
	This folder contains a program for computing ground state energy for electron gas using Coupled Cluster Doubles set up as matrix-matrix multiplications.
	It will use the physical symmetries of the equations to reduce the computation cost of matrices. 

	basis_set.ccp contains basis set for electron gas. Computes one-body and two-body interaction
	ccdintermediates.cpp is a full solver computing ccd equations when provided with basis_set
	solver.cpp is the matrix-matrix block solver
	block.cpp will assist solver.cpp. Each block is an instance of block.cpp

Electron_Gas:
	This folder contains a program for computing ground state energy for electron gas using Coupled Cluster Doubles.
	It uses the naive approach with loops.

	This program will not be used to produce results. It has been created for educational purposes

Pairing_Model:
	his folder contains a program for computing ground state energy for the pairing model using Coupled Cluster Doubles set up as matrix-matrix multiplications.
	It will use the physical symmetries of the equations to reduce the computation cost of matrices. 

	basis_set.ccp contains basis set for pairing model. Computes one-body and two-body interaction
	ccdintermediates.cpp is a full solver computing ccd equations when provided with basis_set
	ccdnaive.cpp is a full solver computing ccd equations when provided with basis_set
	mbptnaive.cpp is a full solver computing mbpt equations when provided with basis_set
	fci.cpp is a naive solver for 4p4h full configuration interaction

	Results: 
		contains results for Pairing model

Report:
	I write the thesis here. 