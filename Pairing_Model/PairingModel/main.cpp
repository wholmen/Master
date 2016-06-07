#include <iostream>
#include <armadillo>
#include <pairingbasis.h>

//#include <ccdnaive.h>
//#include <ccdintermediates.h>
//#include <mbptnaive.h>

#include <fci.h>
#include <../../Solvers/Solvers/ccdnaive.h>
#include <../../Solvers/Solvers/ccdintermediates.h>
#include <../../Solvers/Solvers/mbptnaive.h>

using namespace std;
using namespace arma;


void TestCCDNaive();
void TestMBPTNaive();
void TestCCDIntermediates();
void TestFCI();
void TestCCDBlocks();

void CompareMethods();
void CompareResultsIntermediates();
void CompareIterationsIntermediates();
void CompareMBPT2_CCD1order();
void CompareIterationsIntermediatesWithWeight();

int main()
{

    //TestCCDNaive();
    //TestMBPTNaive();
    //TestCCDIntermediates();
    //TestFCI();

    // Produce results
    CompareMBPT2_CCD1order();
    CompareIterationsIntermediates();
    CompareResultsIntermediates();
    CompareMethods();
    CompareIterationsIntermediatesWithWeight();
}


void CompareMBPT2_CCD1order(){
    // As the first order in t0 should compute the same energy as mbpt to second order, I can use this to benchmark the naive implementation of CCD equations

    int ng = 201; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    double g = 1.0; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_Results_MBPT2_CCD1order.txt");
    myfile << "A comparison of results computed with 1. order CCD equations and MBPT 2. order. These should be equal by construct. Benchmark for CCD Naive." << endl;
    myfile << "The system consist of " << basis.Nshells << " shells, and " << basis.Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << " values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies in following order" << endl;
    myfile << "g, Reference energy, Naive CCD 1. order, MBPT 2. order, MBPT 2. order exact" << endl;

    for (int i=0; i<ng; i++){

        basis.g = G(i);
        double mbpt2exact = -G(i)*G(i) / 4.0 *(1/(4+G(i)) + 1/(6+G(i)) + 1/(2+G(i)) + 1/(4+G(i)));
        CCDNaive ccdn = CCDNaive(basis);
        MBPTNaive mbpt = MBPTNaive(basis);

        myfile << basis.g << " " <<  basis.ReferenceEnergy() << " " << ccdn.CCD(0) << " " << mbpt.MBPT2() << " " << mbpt2exact << endl;
    }
    myfile.close();

}

void CompareIterationsIntermediates(){
    int ng = 21; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    double g = 1.0; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_Iterations_Intermediate_Naive_CCD.txt");
    myfile << "A comparison of results for every iteration computed with CCD equations, using naive approach and intermediates. This works as a benchmark for the inclusion of Intermediates." << endl;
    myfile << "The system consist of " << basis.Nshells << " shells, and " << basis.Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << " values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies for all iterations in following order" << endl;
    myfile << "Line 1: g, Reference energy, Number of iterations, energies for all iterations by Naive implementation of CCD" << endl;
    myfile << "Line 2: g, Reference energy, Number of iterations, energies for all iterations by Intermediate implementation of CCD";

    for (int i=0; i<ng; i++){

        basis.g = G(i);

        CCDNaive ccdn = CCDNaive(basis);
        CCDIntermediates ccdi = CCDIntermediates(basis);

        vec energies = ccdn.CCD_ReturnAllIterations();
        myfile << endl << basis.g << " " << basis.ReferenceEnergy() << " " << ccdn.NIterations;
        for (int n=0; n<ccdn.NIterations; n++){
            myfile << " " << energies(n);
        }
        energies = ccdi.CCD_ReturnAllIterations();
        myfile << endl << basis.g << " " << basis.ReferenceEnergy() << " " << ccdi.NIterations;
        for (int n=0; n<ccdn.NIterations; n++){
            myfile << " " << energies(n);
        }

    }
    myfile.close();
}

void CompareIterationsIntermediatesWithWeight(){
    int ng = 21; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    double g = 1.0; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_Iterations_with_weight_Intermediate_Naive_CCD.txt");
    myfile << "A comparison of results for every iteration computed with CCD equations, using naive approach and intermediates. I have added a weight of 0.5 to deal with divergence. This works as a benchmark for the inclusion of Intermediates." << endl;
    myfile << "The system consist of " << basis.Nshells << " shells, and " << basis.Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << " values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies for all iterations in following order" << endl;
    myfile << "Line 1: g, Reference energy, Number of iterations, energies for all iterations by Naive implementation of CCD" << endl;
    myfile << "Line 2: g, Reference energy, Number of iterations, energies for all iterations by Intermediate implementation of CCD";

    for (int i=0; i<ng; i++){

        basis.g = G(i);

        CCDNaive ccdn = CCDNaive(basis); ccdn.weight = 0.5;
        CCDIntermediates ccdi = CCDIntermediates(basis);  ccdi.weight = 0.5;

        vec energies = ccdn.CCD_ReturnAllIterations();
        myfile << endl << basis.g << " " << basis.ReferenceEnergy() << " " << ccdn.NIterations;
        for (int n=0; n<ccdn.NIterations; n++){
            myfile << " " << energies(n);
        }
        energies = ccdi.CCD_ReturnAllIterations();
        myfile << endl << basis.g << " " << basis.ReferenceEnergy() << " " << ccdi.NIterations;
        for (int n=0; n<ccdn.NIterations; n++){
            myfile << " " << energies(n);
        }

    }
    myfile.close();
}

void CompareResultsIntermediates(){
    int ng = 201; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    double g = 1.0; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_Results_Intermediate_Naive_CCD.txt");
    myfile << "A comparison of results computed with CCD equations, using naive approach and intermediates. This works as a benchmark for the inclusion of Intermediates." << endl;
    myfile << "The system consist of " << basis.Nshells << " shells, and " << basis.Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << " values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies in following order" << endl;
    myfile << "g, Reference energy, Naive CCD, Intermediate CCD" << endl;

    for (int i=0; i<ng; i++){

        basis.g = G(i);

        CCDNaive ccdn = CCDNaive(basis);
        CCDIntermediates ccdi = CCDIntermediates(basis);

        myfile << basis.g << " " <<  basis.ReferenceEnergy() << " " << ccdn.CCD(20) << " " << ccdi.CCD(20) << endl;
    }
    myfile.close();
}

void CompareMethods(){
    int ng = 201; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    double g = 1.0; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_FCI_MBPT_CCD.txt");
    myfile << "A Comparison of FCI, CI, MBPT to second, third and fourth order and CCD equations. Performed with intermediates." << endl;
    myfile << "The system consist of " << basis.Nshells << " shells, and " << basis.Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << " values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies in following order" << endl;
    myfile << "g, Reference energy; FCI, CI, MBPT2, MBPT2 exact, MBPT3, MBPT4, CCD Intermediates." << endl;

    for (int i=0; i<ng; i++){

        basis.g = G(i);

        double mbpt2exact = -G(i)*G(i) / 4.0 *(1/(4+G(i)) + 1/(6+G(i)) + 1/(2+G(i)) + 1/(4+G(i)));
        FCI fci = FCI(basis);
        CCDIntermediates ccdi = CCDIntermediates(basis);
        MBPTNaive mbpt = MBPTNaive(basis);
        myfile << basis.g << " " << basis.ReferenceEnergy() << " " << fci.CalculateFCI() << " " << fci.CalculateCI() << " " << mbpt.MBPT2() << " " << mbpt2exact << " " << mbpt.MBPT3() << " " << mbpt.MBPT4() << " " << ccdi.CCD(20) << endl;
    }
    myfile.close();
}

void TestMBPTNaive(){

    double g = -0.5; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    double dEexact = -g*g / 4.0 *(1/(4+g) + 1/(6+g) + 1/(2+g) + 1/(4+g));

    MBPTNaive solve = MBPTNaive(basis);
    cout << "MBPT2: " << solve.MBPT2() << " MHJ MBPT2: " << solve.MBPT2_MHJ() << "Exact: " << dEexact << endl;
}

void TestCCDNaive(){

    double g = -0.5; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    CCDNaive solve = CCDNaive(basis);
    cout << "CCD Naive: " << solve.CCD(20) << endl;

}

void TestCCDIntermediates(){
    double g = -0.5; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    CCDIntermediates solve = CCDIntermediates(basis);
    cout << "CCD Intermediates" << solve.CCD(20) << endl;
}

void TestFCI(){
    double g = -0.5; int Nshells = 4; int Nshellsfilled = 2; double delta = 1.0;

    PairingBasis basis = PairingBasis(Nshells, Nshellsfilled, g, delta);

    FCI solve = FCI(basis);
    cout << "FCI: " << solve.CalculateFCI() << endl;

}
