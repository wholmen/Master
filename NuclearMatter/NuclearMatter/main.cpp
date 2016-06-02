#include <iostream>
#include <basis_set.h>
#include <ccdblocks.h>
#include <ccdintermediates.h>

using namespace std;

int main()
{
    basis_set basis = basis_set(3,1,0.5);

    cout << "Reference Energy: " << basis.ReferenceEnergy() << endl;
    CCDIntermediates solver = CCDIntermediates(basis);
    cout << solver.CCD(5);
}

