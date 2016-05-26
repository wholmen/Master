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
    int ns = 3;
    int np = 14;
    double rs = 1.0;

    basis_set basis = basis_set(np, ns, rs);

    Solver solve = Solver(basis);

    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using block structure to update amplitudes." << endl;
    cout << endl << "CCD Energy:" << solve.CCD(false) << endl;

    CCDIntermediates solve3 = CCDIntermediates(basis);
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using intermediates to update amplitudes." << endl;
    cout << endl << "CCD Energy:" << solve3.CCD(10) << endl;

    /*
    Solver solve2 = Solver(basis);

    cout << endl << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << "Using naive implementation to update amplitudes." << endl;
    cout << endl << "CCD Energy:" << solve2.CCD(true) << endl;
    */
}
