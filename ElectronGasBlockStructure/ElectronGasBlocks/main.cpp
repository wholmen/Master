#include <iostream>
#include <armadillo>
#include <math.h>
#include <solver.h>
#include <basis_set.h>
#include <block.h>

using namespace arma;
using namespace std;

void TestSolver1();

int main()
{
    TestSolver1();
}

void TestSolver1(){
    int ns = 2;
    int np = 2;
    double rs = 1.0;

    basis_set basis = basis_set(np, ns, rs);

    Solver solve = Solver(basis);


    cout << solve.CCD();


    Block block = Block(basis, np, 12);

    rowvec I0, I1, I2, I3, I4;

    I0 = zeros<rowvec>(3); I0(0) = 1; I0(1) = 0;
    I1 = zeros<rowvec>(3); I1(0) = 2; I1(1) = 0;
    I2 = zeros<rowvec>(3); I2(0) = 0; I2(1) = 1;
    I3 = zeros<rowvec>(3); I3(0) = 1; I3(1) = 0;
    I4 = zeros<rowvec>(3); I4(0) = 4; I4(1) = 1;

    block.AddStates(I0,I0);
    block.AddStates(I1,I1);
    block.AddStates(I2,I2);
    block.AddStates(I3,I3);
    block.AddStates(I4,I4);
}
