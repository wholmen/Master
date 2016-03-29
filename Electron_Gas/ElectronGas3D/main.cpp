#include <iostream>
#include <armadillo>
#include <basis_set.h>

using namespace std;

void TestFromMHJ();

int main()
{
    //basis_set Helium = basis_set(14,2,0.5);

    //Helium.states.print();
    /*
    cout << Helium.Absolute_Difference2(3,37) << endl; // Lovende resultat
    cout << Helium.KDelta_integer(1,1) << endl; // Lovende resultat
    cout << Helium.KDelta_spin(0,0) << endl; // Lovende resultat
    cout << Helium.KDelta_k(2,1) << endl;
    cout << "sum: " << Helium.KDelta_sum(1,1,3,0) << endl;

    rowvec a = zeros<rowvec>(5);rowvec b = zeros<rowvec>(5); b(2) = 1;
    cout << Helium.KDelta_array(a,b) << endl; // Lovende resultat


    //cout << Helium.OneBodyOperator(3,4) << endl;
    cout << Helium.ReferenceEnergy() << endl;
    cout << Helium.ReferenceEnergy2() << endl;

    */
    //Helium.states.print();

    TestFromMHJ();
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
    cout << endl << "Reference energy: " << basis.ReferenceEnergy2() << " hartrees." << endl;
    cout << "Ref.E. per particle: " << basis.ReferenceEnergy2() / number_of_particles << " hartrees" << endl;
}
