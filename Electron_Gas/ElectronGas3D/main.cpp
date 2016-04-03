#include <iostream>
#include <armadillo>
#include <basis_set.h>
#include <naivesolvers.h>

using namespace std;

void TestFromMHJ();

int main()
{
    int ns = 3;
    int np = 2;
    double rs = 1.0;

    basis_set basis = basis_set(np, ns, rs);
    NaiveSolvers solver = NaiveSolvers(basis);

    cout << solver.CCD() << endl;
    cout << solver.MBPT4() << endl;

    //TestFromMHJ();
}

void TestFromMHJ(){
    int number_of_shells = 5;
    int number_of_particles = 14;
    double rs = 1.0;

    basis_set basis = basis_set(number_of_particles, number_of_shells, rs);

    cout << "Number of shells: " << number_of_shells << endl;
    cout << "Number of particles: " << number_of_particles << endl;
    cout << "Value of rs: " << rs << endl;
    cout << "Number of states: " << basis.nstates << endl;
    cout << endl << "Reference energy: " << basis.ReferenceEnergy() << " hartrees." << endl;
    cout << "Ref.E. per particle: " << basis.ReferenceEnergy() / number_of_particles << " hartrees" << endl;
}
