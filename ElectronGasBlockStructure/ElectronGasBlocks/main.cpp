#include <iostream>
#include <armadillo>
#include <math.h>
#include <solver.h>
#include <basis_set.h>
#include <block.h>
#include <ccdintermediates.h>

using namespace arma;
using namespace std;

void TestSolver1();

int main()
{
    TestSolver1();
}

void TestSolver1(){
    int Nshells = 6;
    int NfilledShells = 2;
    double rs = 1.0;

    basis_set basis = basis_set(Nshells, NfilledShells, rs);

    Solver solve = Solver(basis);

    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using block structure to update amplitudes." << endl;
    cout << "Nstates: " << basis.Nstates << " Nholes: " << basis.Nholes << endl;
    cout << endl << "CCD Energy:" << solve.CCD(false)  / basis.Nholes << endl;
    /*
    CCDIntermediates solve3 = CCDIntermediates(basis);
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using intermediates to update amplitudes." << endl;
    cout << endl << "CCD Energy:" << solve3.CCD(10) << endl;

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
