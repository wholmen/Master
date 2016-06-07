#include <iostream>
#include <armadillo>
#include <math.h>
#include <solver.h>
#include <electronbasis.h>
#include <block.h>

#include <../../Solvers/Solvers/ccdintermediates.h>
#include <../../Solvers/Solvers/ccdnaive.h>



using namespace arma;
using namespace std;

void TestSolver1();
void ResultsNh14Ns54_Intermediates();
void AllIterationsNh14Ns54_Intermediates();

int main()
{
    TestSolver1();
    //AllIterationsNh14Ns54_Intermediates();
    //ResultsNh14Ns54_Intermediates();
}

void TestSolver1(){
    int Nshells = 4;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    Solver solve = Solver(basis);

    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using block structure to update amplitudes." << endl;
    cout << "Nstates: " << basis.Nstates << " Nholes: " << basis.Nholes << endl;
    cout << endl << "CCD Energy:" << solve.CCD(20)  / basis.Nholes << endl;

    //CCDIntermediates solve3 = CCDIntermediates(basis);
    //cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using intermediates to update amplitudes." << endl;
    //cout << endl << "CCD Energy:" << solve3.CCD(20) << endl;

    /*
    Solver solve4 = Solver(basis);
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using naive approach to update amplitudes." << endl;
    cout << endl << "CCD Energy:" << solve4.CCD(true) << endl;
    */
    /*
    Solver solve2 = Solver(basis);

    cout << endl << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using naive implementation to update amplitudes." << endl;
    cout << endl << "CCD Energy:" << solve2.CCD(true) << endl;
    */
}



void AllIterationsNh14Ns54_Intermediates(){
    int Nshells = 4;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/CCDIntermediates_AllIterations_Nh14_Ns54.txt");
    myfile << "A coupled cluster study of electron gas using intermediates. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "rs, Reference Energy, Niterations, all iterations CCDIntermediate results" << endl;

    CCDIntermediates solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " ";
    vec energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << "   ";
    for (int i=0; i<solve.NIterations; i++) myfile << energies(i) << " ";
    myfile << endl;

    basis.rs = 0.5;
    solve = CCDIntermediates(basis);
    myfile << basis.rs << "   " << basis.ReferenceEnergy() << "   ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << "   ";
    for (int i=0; i<solve.NIterations; i++) myfile << energies(i) << " ";
    myfile << endl;

    basis.rs = 2.0;
    solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << "   ";
    for (int i=0; i<solve.NIterations; i++) myfile << energies(i) << " ";
    myfile << endl;

    myfile.close();
}

void ResultsNh14Ns54_Intermediates(){

    int Nshells = 4;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/CCDIntermediates_Nh14_Ns54.txt");
    myfile << "A coupled cluster study of electron gas using intermediates. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "rs, Reference Energy, CCDIntermediate result" << endl;

    CCDIntermediates solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " " << solve.CCD(20) << endl;

    basis.rs = 0.5;
    solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " " << solve.CCD(20) << endl;

    basis.rs = 2.0;
    solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " " << solve.CCD(20) << endl;

    myfile.close();
}
