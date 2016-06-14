#include <iostream>
#include <nuclearbasis.h>
#include <solver.h>
#include <../../Solvers/Solvers/ccdintermediates.h>

using namespace std;

void AllIterationsNh14Ns54_Intermediates();
void ResultsNh14Ns54_Intermediates();
void TestBlocks();
void TestIntermediates();

int main()
{

    //TestIntermediates();
    TestBlocks();
    //AllIterationsNh14Ns54_Intermediates();
    //ResultsNh14Ns54_Intermediates();

    /*NuclearBasis basis = NuclearBasis(3,2,0.5);

    cout << "Reference Energy: " << basis.ReferenceEnergy() << endl;
    CCDIntermediates solver = CCDIntermediates(basis);
    cout << solver.CCD(5); */
}

void TestIntermediates(){
    int Nshells = 3; int NfilledShells = 2; double rs = 1;

    NuclearBasis basis = NuclearBasis(Nshells,NfilledShells,rs);
    cout << "Reference Energy: " << basis.ReferenceEnergy() << endl;
    CCDIntermediates solver = CCDIntermediates(basis);
    cout << solver.CCD(10) << endl;
}

void TestBlocks(){
    int Nshells = 2; int NfilledShells = 1; double rs = 0.4;

    NuclearBasis basis = NuclearBasis(Nshells,NfilledShells,rs);
    cout << "Reference Energy: " << basis.ReferenceEnergy() << endl;
    Solver solver = Solver(basis);
    cout << solver.CCD(10) << endl;
}

void AllIterationsNh14Ns54_Intermediates(){
    int Nshells = 2;
    int NfilledShells = 1;
    double rs = 1.0;

    NuclearBasis basis = NuclearBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/CCDIntermediates_AllIterations_Nshells4_Nfilled2.txt");
    myfile << "A coupled cluster study of nuclear matter using intermediates. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "rs, Reference Energy, Niterations, all iterations CCDIntermediate results" << endl;

    CCDIntermediates solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " ";
    vec energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << "   ";
    for (int i=0; i<solve.NIterations; i++) myfile << energies(i) << " ";
    myfile << endl;

    rs = 0.5; basis = NuclearBasis(Nshells, NfilledShells, rs);
    solve = CCDIntermediates(basis);
    myfile << basis.rs << "   " << basis.ReferenceEnergy() << "   ";
    energies = solve.CCD_ReturnAllIterations();
    myfile << solve.NIterations << "   ";
    for (int i=0; i<solve.NIterations; i++) myfile << energies(i) << " ";
    myfile << endl;

    rs = 2.0; basis = NuclearBasis(Nshells, NfilledShells, rs);
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

    NuclearBasis basis = NuclearBasis(Nshells, NfilledShells, rs);

    ofstream myfile;
    myfile.open("../Results/CCDIntermediates_Nshells4_Nfilled2.txt");
    myfile << "A coupled cluster study of nuclear matter using intermediates. " << endl;
    myfile << "The system consist of " << basis.Nstates << " states and " << basis.Nholes << " holes. These values are chosen to compare results with Audun Skau Hansen's master thesis." << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "rs, Reference Energy, CCDIntermediate result" << endl;

    CCDIntermediates solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " " << solve.CCD(20) << endl;

    rs = 0.5; basis = NuclearBasis(Nshells, NfilledShells, rs);
    solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " " << solve.CCD(20) << endl;

    rs = 2.0; basis = NuclearBasis(Nshells, NfilledShells, rs);
    solve = CCDIntermediates(basis);
    myfile << basis.rs << " " << basis.ReferenceEnergy() << " " << solve.CCD(20) << endl;

    myfile.close();
}
