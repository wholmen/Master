#include "electronbasis.h"

ElectronBasis::ElectronBasis() {}

ElectronBasis::ElectronBasis(int Nshells_input, int NshellsFilled_input, double rs_input){


    Nshells = Nshells_input; NshellsFilled = NshellsFilled_input; rs = rs_input;

    States = zeros<mat>(0,5);

    Nstates = 0; Nholes = 0;

    for (int shell=0; shell<Nshells; shell++){
        // shell is the value x^2 + y^2 + z^2 must match

        for (int nx=-shell; nx <= shell; nx++){
            for (int ny=-shell; ny <= shell; ny++){
                for (int nz=-shell; nz <= shell; nz++){
                    // Looping through all possible combinations of nx, ny, nz below n_shell

                    // Energy integer
                    int e = nx*nx + ny*ny + nz*nz;

                    // Removing duplicates from previous values for shell
                    if (e == shell){

                        // Adding spin
                        for (int ms = -1; ms <= 1; ms+=2){

                            States.insert_rows(Nstates,1);

                            States(Nstates,0) = e;
                            States(Nstates,1) = nx;
                            States(Nstates,2) = ny;
                            States(Nstates,3) = nz;
                            States(Nstates,4) = ms;

                            Nstates++;
                            if (shell < NshellsFilled) Nholes ++;
                        }
                    }
                }
            }
        }
    }

    Nparticles = Nstates - Nholes;
    // Various calculations of variables needed
    L3 = (4*pi*rs*rs*rs*Nholes) / 3.0;
    L2 = pow(L3, 2.0/3 );
    L1 = pow(L3, 1.0/3 );

    kstep = 2*pi / L1;

    // Calculating the actual one-body energy for given e=nx^2+ny^2+nz^2
    for (int i=0; i<Nstates; i++){
        States(i,0) = States(i,0) * 2*pi*pi / L2;
    }
}

double ElectronBasis::ReferenceEnergy(){

    double Energy = 0.0;

    if (Nholes <= Nstates){

        for (int p=0; p<Nholes; p++){
            Energy += OneBodyOperator(p,p);

            for (int q=0; q<Nholes; q++){
                if (p != q) {
                    Energy += 0.5*TwoBodyOperator(p,q,p,q);
                }
            }
        }
    }
    else {
        cout << "Number of particles cannot exceed number of available states.";
    }
    return Energy;
}

double ElectronBasis::ei(int q){
    double interaction = 0;
    for (int i=0; i<Nholes; i++){
        interaction += TwoBodyOperator(q,i,q,i);
    }
    return OneBodyOperator(q,q) + interaction;
}

double ElectronBasis::epsilon(int i, int j, int a, int b){
    // Function to compute the sum of h(i) + h(j) - h(a) - h(b)
    return ei(i) + ei(j) - ei(a) - ei(b);
}


double ElectronBasis::OneBodyOperator(int p, int q){
    // Need to return the one-body energy value for p if p == q
    return ( p == q ) ? States(p,0) : 0;
}

double ElectronBasis::TwoBodyOperator(int p, int q, int r, int s){
    // Two body interaction for electron gas

    double asym, asym1, asym2;
    asym = 0.0; asym1 = 0.0; asym2 = 0.0;

    if (KDelta_sum(p,q,r,s) == 1){
        asym = 1.0/L3;

        if ( KDelta_spin(p,r)*KDelta_spin(q,s) == 1 ){

            if (KDelta_k(p,r) != 1){
                asym1 = L2 / (pi*Absolute_Difference2(r,p));
            }
        }
        if ( KDelta_spin(p,s)*KDelta_spin(q,r) == 1 ){

            if (KDelta_k(p,s) != 1){
                asym2 = L2 / (pi*Absolute_Difference2(s,p));
            }
        }
    }
    return asym * (asym1 - asym2);
}

int ElectronBasis::KDelta_integer(int a, int b){
    // Returns 1 if a and b are equal
    return ( a == b ) ? 1:0;
}

int ElectronBasis::KDelta_array(rowvec a, rowvec b){
    // Returns 1 if all elements are equal, else return 0
    int value = 1;
    for (int i=1; i<=3; i++){
        value = ( a(i) == b(i) ) ? value*1 : value*0;
    }
    return value;
}

int ElectronBasis::KDelta_k(int a, int b){
    // Kroenecker delta for wave numbers ka and kb
    return KDelta_array( States.row(a), States.row(b));
}

int ElectronBasis::KDelta_spin(int a, int b){
    // Kroenecker delta for spin integer ms_a and ms_b
    return KDelta_integer( States(a,4), States(b,4));
}

int ElectronBasis::KDelta_sum(int a, int b, int c, int d){
    // Kroenecker delta comparing a+b and c+d
    return KDelta_array( States.row(a)+States.row(b) , States.row(c)+States.row(d));
}

double ElectronBasis::Absolute_Difference2(int a, int b){
    // Returning difference squared between a and b
    double diff = 0.0;
    for (int i=1; i<=3; i++){
        diff += (States(a,i) - States(b,i)) * (States(a,i) - States(b,i));
    }
    return diff;
}





