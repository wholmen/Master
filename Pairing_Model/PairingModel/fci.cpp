#include "fci.h"

FCI::FCI()
{

}

FCI::FCI(PairingBasis BASIS){
    basis = BASIS;
    g = basis.g;
}

double FCI::CalculateFCI(){

    mat H = zeros<mat>(6,6);

    H(0,0)= 2-g ; H(0,1)=-g/2.; H(0,2)=-g/2.; H(0,3)=-g/2.; H(0,4)=-g/2.; H(0,5)= 0   ;
    H(1,0)=-g/2.; H(1,1)= 4-g ; H(1,2)=-g/2.; H(1,3)=-g/2.; H(1,4)= 0;    H(1,5)=-g/2.;
    H(2,0)=-g/2.; H(2,1)=-g/2.; H(2,2)= 6-g ; H(2,3)= 0   ; H(2,4)=-g/2.; H(2,5)=-g/2.;
    H(3,0)=-g/2.; H(3,1)=-g/2.; H(3,2)= 0   ; H(3,3)= 6-g ; H(3,4)=-g/2.; H(3,5)=-g/2.;
    H(4,0)=-g/2.; H(4,1)= 0   ; H(4,2)=-g/2.; H(4,3)=-g/2.; H(4,4)= 8-g ; H(4,5)=-g/2.;
    H(5,0)= 0   ; H(5,1)=-g/2.; H(5,2)=-g/2.; H(5,3)=-g/2.; H(5,4)=-g/2.; H(5,5)= 10-g;

    vec eigval = eig_sym(H);

    double E = eigval.min();

    double ReferenceEnergy = basis.ReferenceEnergy();

    double CorrelationEnergy = E - ReferenceEnergy;

    return CorrelationEnergy;;
}

double FCI::CalculateCI(){

    mat H = zeros<mat>(5,5);

    H(0,0)= 2-g ; H(0,1)=-g/2.; H(0,2)=-g/2.; H(0,3)=-g/2.; H(0,4)=-g/2.;
    H(1,0)=-g/2.; H(1,1)= 4-g ; H(1,2)=-g/2.; H(1,3)=-g/2.; H(1,4)= 0;
    H(2,0)=-g/2.; H(2,1)=-g/2.; H(2,2)= 6-g ; H(2,3)= 0   ; H(2,4)=-g/2.;
    H(3,0)=-g/2.; H(3,1)=-g/2.; H(3,2)= 0   ; H(3,3)= 6-g ; H(3,4)=-g/2.;
    H(4,0)=-g/2.; H(4,1)= 0   ; H(4,2)=-g/2.; H(4,3)=-g/2.; H(4,4)= 8-g ;

    vec eigval = eig_sym(H);

    double E = eigval.min();

    double ReferenceEnergy = basis.ReferenceEnergy();

    double CorrelationEnergy = E - ReferenceEnergy;

    return CorrelationEnergy;
}
