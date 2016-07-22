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
void ResultsNshells9Nfilled2_Blocks();
void ResultsNshells10Nfilled2_Blocks();
void ResultsNshells11Nfilled2_Blocks();
void ResultsNshells12Nfilled2_Blocks();
void ResultsNshells13Nfilled2_Blocks();

void ResultsNshells3to40Nfilled2_Blocks();

void AllIterationsNh14Ns54_Intermediates();
void AllIterationsNh14Ns54_Blocks();


int main()
{
    // Function to test solver
    //TestSolver1();

    ResultsNshells3to40Nfilled2_Blocks();

    // Producing results for rs in [0.5,1.0,2.0] for different amount of states.
    /*
    ResultsNshells4Nfilled2_Blocks();
    ResultsNshells5Nfilled2_Blocks();
    ResultsNshells6Nfilled2_Blocks();
    ResultsNshells7Nfilled2_Blocks();
    ResultsNshells8Nfilled2_Blocks();
    ResultsNshells9Nfilled2_Blocks();
    ResultsNshells10Nfilled2_Blocks();
    ResultsNshells11Nfilled2_Blocks();
    ResultsNshells12Nfilled2_Blocks();
    ResultsNshells13Nfilled2_Blocks();
    */
    // Functions to produce a comparison between blocks and intermediates
    /*
    AllIterationsNh14Ns54_Blocks();
    AllIterationsNh14Ns54_Intermediates();
    ResultsNh14Ns54_Intermediates();
    */

    // Function to compute all magic numbers
    //PrintMagicNumbers();
}

void PrintMagicNumbers(){

    int NfilledShells = 1;
    double rs = 1.0;

    for (int i=2; i<14; i++){

        ElectronBasis basis = ElectronBasis(i, NfilledShells, rs);
        cout << "Shell level: " << i-1 << "  Nstates: " << basis.Nstates << endl;
    }
}



void TestSolver1(){
    clock_t start, finish;
    int Nshells = 40;
    int NfilledShells = 2;
    double rs = 1.0;

    start = clock();
    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);
    finish = clock(); double ptime = (double(finish-start)/CLOCKS_PER_SEC);
    cout << "Electronbasis needed " << ptime << " seconds." << endl;

    start = clock();
    Solver solve = Solver(basis);
    finish = clock(); ptime = (double(finish-start)/CLOCKS_PER_SEC);
    cout << "Setting up solver needed " << ptime << " seconds." << endl;

    solve.weight = 0.3;
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using block structure to update amplitudes." << endl;
    cout << "Nstates: " << basis.Nstates << " Nholes: " << basis.Nholes << endl;

    start = clock();
    cout << endl << "CCD Energy:" << setprecision(12) << solve.CCD(200) << endl;
    finish = clock(); ptime = (double(finish-start)/CLOCKS_PER_SEC);
    cout << "Iterations needed " << ptime << " seconds." << endl;

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

void ResultsNshells9Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells9Nfilled2_Blocks()' started" << endl;

    int Nshells = 9;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells9Nfilled2.txt");
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

    cout << "Function 'ResultsNshells9Nfilled2_Blocks()' ended" << endl;
}

void ResultsNshells10Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells10Nfilled2_Blocks()' started" << endl;

    int Nshells = 10;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells10Nfilled2.txt");
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

    cout << "Function 'ResultsNshells10Nfilled2_Blocks()' ended" << endl;
}

void ResultsNshells11Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells11Nfilled2_Blocks()' started" << endl;

    int Nshells = 11;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells11Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;

    Solver solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    myfile.close();

    cout << "Function 'ResultsNshells11Nfilled2_Blocks()' ended" << endl;
}

void ResultsNshells12Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells12Nfilled2_Blocks()' started" << endl;

    int Nshells = 12;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells12Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;

    Solver solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    myfile.close();

    cout << "Function 'ResultsNshells12Nfilled2_Blocks()' ended" << endl;
}

void ResultsNshells13Nfilled2_Blocks(){

    cout << "Function 'ResultsNshells13Nfilled2_Blocks()' started" << endl;

    int Nshells = 13;
    int NfilledShells = 2;
    double rs = 1.0;

    ElectronBasis basis = ElectronBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells13Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using blocks. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "weight, rs, Reference Energy, CCDIntermediate result, Iterations" << endl;

    Solver solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 0.5; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    // Changing rs for new results.
    rs = 2.0; basis = ElectronBasis(Nshells, NfilledShells, rs);

    solve = Solver(basis);
    solve.weight = 0.3;
    myfile << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << setprecision(12) << solve.CCD(300);
    myfile << " " << solve.NIterations << endl;

    myfile.close();

    cout << "Function 'ResultsNshells13Nfilled2_Blocks()' ended" << endl;
}



void ResultsNshells3to40Nfilled2_Blocks(){

    ofstream myfile;
    myfile.open("../Results/Solver_Nshells3to40Nfilled2.txt");
    myfile << "A coupled cluster study of electron gas using a parallel code with blocks implementation for different amount of states" << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "Shells, States, Occupied States, weight, rs, Reference Energy, Iterations, time spent initiate, time spent iterating, CCD Energy" << endl;


    int Nfilled = 2; double rs = 1.0;
    for (int Nshells=3; Nshells<41; Nshells++){

        ElectronBasis basis = ElectronBasis(Nshells,Nfilled,rs);

        double time_init_0 = omp_get_wtime();
        Solver solve = Solver(basis);
        double time_init_1 = omp_get_wtime();

        solve.weight = 0.3;

        double time_iter_0 = omp_get_wtime();
        double E = solve.CCD(200);
        double time_iter_1 = omp_get_wtime();

        myfile << setprecision(12) << Nshells << " " << basis.Nstates << " " << basis.Nholes << " " << solve.weight << " " << basis.rs << " " << basis.ReferenceEnergy() << " " << solve.NIterations << " " << time_init_1 - time_init_0 << " " << time_iter_1 - time_iter_0 << " " << E << endl;
    }

}

// Functions to produce a comparison between blocks and intermediates

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
