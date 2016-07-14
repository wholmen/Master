#include <iostream>
#include <armadillo>
#include <math.h>
#include <solver.h>
#include <electronbasis.h>
#include <block.h>
#include <iomanip>

#include <../../Solvers/Solvers/ccdintermediates.h>
#include <../../Solvers/Solvers/ccdnaive.h>


using namespace arma;
using namespace std;

void PrintMagicNumbers();
void TestSolver1();
void ResultsNh14Ns54_Intermediates();
void ResultsNshells4Nfilled2_Blocks();
void ResultsNshells5Nfilled2_Blocks();
void ResultsNshells6Nfilled2_Blocks();
void ResultsNshells7Nfilled2_Blocks();
void ResultsNshells8Nfilled2_Blocks();

void AllIterationsNh14Ns54_Intermediates();
void AllIterationsNh14Ns54_Blocks();


int main()
{
    TestSolver1();
    //PrintMagicNumbers();

    /*
    ResultsNshells4Nfilled2_Blocks();
    ResultsNshells5Nfilled2_Blocks();
    ResultsNshells6Nfilled2_Blocks();
    ResultsNshells7Nfilled2_Blocks();
    ResultsNshells8Nfilled2_Blocks();

    AllIterationsNh14Ns54_Blocks();

    AllIterationsNh14Ns54_Intermediates();
    ResultsNh14Ns54_Intermediates();
    */
}

void PrintMagicNumbers(){
    int NfilledShells = 1;
    double rs = 1.0;

    for (int i=2; i<12; i++){

        ElectronBasis basis = ElectronBasis(i, NfilledShells, rs);
        cout << "Shell level: " << i-1 << "  Nstates: " << basis.Nstates << endl;
    }
}

void TestSolver1(){
    int Nshells = 7;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    Solver solve = Solver(basis);
    solve.weight = 0.3;
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using block structure to update amplitudes." << endl;
    cout << "Nstates: " << basis.Nstates << " Nholes: " << basis.Nholes << endl;
    cout << endl << "CCD Energy:" << setprecision(12) << solve.CCD(200) << endl;

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




void ResultsNshells4Nfilled2_Blocks(){

    // Computing results using block matrix-matrix multiplications with rs as 0.5, 1.0 and 2.0.
    // Changing weight from 1.0, 0.7 and 0.3
    // Using 54 basis states where 14 are hole states.

    cout << "Function 'ResultsNshells4Nfilled2_Blocks()' started" << endl;

    int Nshells = 4;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells4Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;

    Solver solve = Solver(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;


    myfile.close();

    cout << "Function 'ResultsNshells4Nfilled2_Blocks()' ended" << endl;
}

void ResultsNshells5Nfilled2_Blocks(){

    // Computing results using block matrix-matrix multiplications with rs as 0.5, 1.0 and 2.0.
    // Changing weight from 1.0, 0.7 and 0.3
    // Using 54 basis states where 14 are hole states.

    cout << "Function 'ResultsNshells5Nfilled2_Blocks()' started" << endl;

    int Nshells = 5;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells5Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;

    Solver solve = Solver(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;
    myfile.close();

    cout << "Function 'ResultsNshells5Nfilled2_Blocks()' ended" << endl;
}


void ResultsNshells6Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells6Nfilled2_Blocks()' started" << endl;

    int Nshells = 6;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells6Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;


    Solver solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    myfile.close();

    cout << "Function 'ResultsNshells6Nfilled2_Blocks()' ended" << endl;
}

void ResultsNshells7Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells7Nfilled2_Blocks()' started" << endl;

    int Nshells = 7;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells7Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;


    Solver solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    myfile.close();

    cout << "Function 'ResultsNshells7Nfilled2_Blocks()' ended" << endl;
}

void ResultsNshells8Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells8Nfilled2_Blocks()' started" << endl;

    int Nshells = 8;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells8Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;

    Solver solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(200);
    myfile << " " << solve.NIterations << endl;

    myfile.close();

    cout << "Function 'ResultsNshells8Nfilled2_Blocks()' ended" << endl;
}

void AllIterationsNh14Ns54_Blocks(){
    cout << "Function 'AllIterationsNh14Ns54_Blocks()' started" << endl;

    int Nshells = 4;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_AllIterations_Nh14_Ns54.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "rs, weight, Reference Energy, Niterations, all iterations CCDIntermediate results" << endl;

    Solver solve = Solver(basis);
    solve.weight = 1.0;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    vec energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;


    // New rs
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 1.0;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;


    // New rs
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 1.0;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = Solver(basis);
    solve.weight = 0.7;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;


    myfile.close();

    cout << "Function 'AllIterationsNh14Ns54_Blocks()' ended" << endl;
}

void AllIterationsNh14Ns54_Intermediates(){

    cout << "Function 'AllIterationsNh14Ns54_Intermediates()' started" << endl;

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
    solve.weight = 1.0;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    vec energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.7;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.3;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;


    // New rs
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = CCDIntermediates(basis);
    solve.weight = 1.0;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.7;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.3;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;


    // New rs
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = CCDIntermediates(basis);
    solve.weight = 1.0;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.7;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.3;
    myfile << basis.rs << " " << solve.weight << " " << basis.ReferenceEnergy() << " ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << " ";
    for (int i=0; i<solve.NIterations; i++) myfile << setprecision(12) << energies(i) << " ";
    myfile << endl;

    myfile.close();

    cout << "Function 'AllIterationsNh14Ns54_Intermediates()' ended" << endl;
}

void ResultsNh14Ns54_Intermediates(){

    cout << "Function 'ResultsNh14Ns54_Intermediates()' started" << endl;

    int Nshells = 4;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/CCDIntermediates_Nh14_Ns54.txt");
    myfile << "A coupled cluster study of electron gas using intermediates. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "rs, weight, Reference Energy, CCDIntermediate result, Iterations" << endl;

    CCDIntermediates solve = CCDIntermediates(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = CCDIntermediates(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;


    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = CCDIntermediates(basis);
    solve.weight = 1.0;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.7;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;

    solve = CCDIntermediates(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(50);
    myfile << " " << solve.NIterations << endl;


    myfile.close();
    cout << "Function 'ResultsNh14Ns54_Intermediates()' ended" << endl;
}
