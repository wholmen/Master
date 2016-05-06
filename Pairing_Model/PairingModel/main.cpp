#include <iostream>
#include <armadillo>
#include <basis_set.h>
#include <naivesolvers.h>
#include <ccdblocks2.h>

using namespace std;
using namespace arma;


void TestNaiveSolver();
void TestBlocks2();


int main()
{
    TestBlocks2();
}


void TestBlocks2(){
    double g = 0.712;
    basis_set basis = basis_set(4,4,g,1);

    CCDBlocks2 solver = CCDBlocks2(basis);


}

void TestNaiveSolver(){
    double g = 0.712;
    basis_set basis = basis_set(4,4,g,1);

    cout << "g: " << g << endl;
    cout << "Reference energy: " << basis.ReferenceEnergy() << endl;
    //cout << "2 - g: " << 2-g << endl;

    NaiveSolvers mbpt2 = NaiveSolvers(basis);
    double e1 = mbpt2.MBPT2_MHJ();
    double e2 = mbpt2.MBPT2();
    double e3 = mbpt2.MBPT3();
    double e4 = mbpt2.MBPT4();

    double eccd = mbpt2.CCD();

    cout << "MHJ MBPT2: " << e1 << "  ; MBPT2 energy: " << e2 << endl << "  ; MBPT3 energy: " << e3 << "  ; MBPT4 energy: " << e4 << endl;
    cout << "CC energy; " << eccd << endl;
}
