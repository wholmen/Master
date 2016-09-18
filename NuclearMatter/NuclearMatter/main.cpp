#include <iostream>
#include <nuclearbasis.h>
#include <solver.h>
#include <../../Solvers/Solvers/ccdintermediates.h>
#include <iomanip>
#include <omp.h>

using namespace std;

void AllIterationsNh14Ns54_Intermediates();
void ResultsNh14Ns54_Intermediates();
void TestBlocks();
void TestIntermediates();

int main()
{

    TestBlocks();

}


void TestBlocks(){
    int Nshells = 10; int NfilledShells = 2; double n = 0.05;

    NuclearBasis basis = NuclearBasis(Nshells,NfilledShells,n,true);
    cout << "Reference Energy: " << basis.ReferenceEnergy() << endl;
    Solver solver = Solver(basis);
    solver.weight = 1.0;
    cout << "Nparticles: " << solver.Nholes << " Nstates: " << solver.Nstates << " Energy:" << solver.CCD(200) << endl;
}
