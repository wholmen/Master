#include <iostream>
#include <armadillo>
#include <basis_set.h>

#include <ccdnaive.h>
#include <ccdintermediates.h>
#include <ccdblocks.h>
#include <mbptnaive.h>
#include <fci.h>

using namespace std;
using namespace arma;


void TestCCDNaive();
void TestMBPTNaive();
void TestCCDIntermediates();
void TestFCI();
void CompareMethods();
void CompareResultsIntermediates();
void CompareIterationsIntermediates();
void CompareMBPT2_CCD1order();

int main()
{
    //CompareMBPT2_CCD1order();
    //CompareIterationsIntermediates();
    //CompareResultsIntermediates();
    //CompareMethods();
    TestCCDNaive();
    TestMBPTNaive();
    //TestCCDIntermediates();
    //TestFCI();
}


void CompareMBPT2_CCD1order(){
    // As the first order in t0 should compute the same energy as mbpt to second order, I can use this to benchmark the naive implementation of CCD equations

    int ng = 201; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    int Nholes = 4; int Nshells = 4; double delta = 1.0;
    basis_set basis = basis_set(Nholes,Nshells, 1, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_Results_MBPT2_CCD1order.txt");
    myfile << "A comparison of results computed with 1. order CCD equations and MBPT 2. order. These should be equal by construct. Benchmark for CCD Naive." << endl;
    myfile << "The system consist of " << Nshells << " shells, and " << Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << "values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies in following order" << endl;
    myfile << "g, Reference energy, Naive CCD 1. order, MBPT 2. order" << endl;

    for (int i=0; i<ng; i++){

        basis.g = G(i);

        CCDNaive ccdn = CCDNaive(basis);
        MBPTNaive mbpt = MBPTNaive(basis);

        myfile << basis.g << " " <<  basis.ReferenceEnergy() << " " << ccdn.CCD(0) << " " << mbpt.MBPT2() << endl;
    }
    myfile.close();

}

void CompareIterationsIntermediates(){
    int ng = 11; double gval = 0.5;
    vec G = linspace<vec>(-gval,gval,ng);

    int Nholes = 4; int Nshells = 4; double delta = 1.0;
    basis_set basis = basis_set(Nholes,Nshells, 1, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_Iterations_Intermediate_Naive_CCD.txt");
    myfile << "A comparison of results for every iteration computed with CCD equations, using naive approach and intermediates. This works as a benchmark for the inclusion of Intermediates." << endl;
    myfile << "The system consist of " << Nshells << " shells, and " << Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << "values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
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

void CompareResultsIntermediates(){
    int ng = 201; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    int Nholes = 4; int Nshells = 4; double delta = 1.0;
    basis_set basis = basis_set(Nholes,Nshells, 1, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_Results_Intermediate_Naive_CCD.txt");
    myfile << "A comparison of results computed with CCD equations, using naive approach and intermediates. This works as a benchmark for the inclusion of Intermediates." << endl;
    myfile << "The system consist of " << Nshells << " shells, and " << Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << "values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies in following order" << endl;
    myfile << "g, Reference energy, Naive CCD, Intermediate CCD" << endl;

    for (int i=0; i<ng; i++){

        basis.g = G(i);

        CCDNaive ccdn = CCDNaive(basis);
        CCDIntermediates ccdi = CCDIntermediates(basis);

        myfile << basis.g << " " <<  basis.ReferenceEnergy() << " " << ccdn.CCD(10) << " " << ccdi.CCD(10) << endl;
    }
    myfile.close();
}

void CompareMethods(){
    int ng = 201; double gval = 1.0;
    vec G = linspace<vec>(-gval,gval,ng);

    int Nholes = 4; int Nshells = 4; double delta = 1.0;
    basis_set basis = basis_set(Nholes,Nshells, 1, delta);

    ofstream myfile;
    myfile.open("../Results/Comparison_FCI_MBPT_CCD.txt");
    myfile << "A Comparison of FCI, CI, MBPT to second, third and fourth order and CCD equations. Performed with intermediates." << endl;
    myfile << "The system consist of " << Nshells << " shells, and " << Nholes << " hole states. " << endl;
    myfile << "I have used " << ng << "values of g, equally spaced, ranging from " << -gval << " to " << gval << endl;
    myfile << "Below are the reference energy and calculated correlation energies in following order" << endl;
    myfile << "g, Reference energy; FCI, CI, MBPT2, MBPT3, MBPT4, CCD Intermediates." << endl;

    for (int i=0; i<ng; i++){

        basis.g = G(i);

        FCI fci = FCI(basis);
        CCDIntermediates ccdi = CCDIntermediates(basis);
        MBPTNaive mbpt = MBPTNaive(basis);
        myfile << basis.g << " " << basis.ReferenceEnergy() << " " << fci.CalculateFCI() << " " << fci.CalculateCI() << " " << mbpt.MBPT2() << " " << mbpt.MBPT3() << " " << mbpt.MBPT4() << " " << ccdi.CCD(10) << endl;
    }
    myfile.close();
}

void TestMBPTNaive(){

    double g = -0.5;
    basis_set basis = basis_set(4,4,g,1);

    double dEexact = -g*g / 4.0 *(1/(4+g) + 1/(6+g) + 1/(2+g) + 1/(4+g));

    MBPTNaive solve = MBPTNaive(basis);
    cout << "MBPT2: " << solve.MBPT2() << " MHJ MBPT2: " << solve.MBPT2_MHJ() << "Exact: " << dEexact << endl;
}

void TestCCDNaive(){

    double g = -0.5;
    basis_set basis = basis_set(4,4,g,1);

    CCDNaive solve = CCDNaive(basis);
    cout << "CCD Naive: " << solve.CCD(10) << endl;

}

void TestCCDIntermediates(){
    double g = 0.712;
    basis_set basis = basis_set(4,4,g,1);

    CCDIntermediates solve = CCDIntermediates(basis);
    cout << "CCD Intermediates" << solve.CCD(10) << endl;
}

void TestFCI(){
    double g = 0.712;
    basis_set basis = basis_set(4,4,g,1);

    FCI solve = FCI(basis);
    cout << "FCI: " << solve.CalculateFCI() << endl;

}
