#include "basis_set.h"

basis_set::basis_set() {}

basis_set::basis_set(int nparticles, int nshells, double RS){
    Nparticles = nparticles; Nshells = nshells; rs = RS;

    states = zeros<mat>(0,5);

    nstates = 0;

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

                            states.insert_rows(nstates,1);

                            states(nstates,0) = e;
                            states(nstates,1) = nx;
                            states(nstates,2) = ny;
                            states(nstates,3) = nz;
                            states(nstates,4) = ms;

                            nstates++;
                        }
                    }
                }
            }
        }
    }

    // Various calculations of variables needed
    L3 = (4*pi*rs*rs*rs*Nparticles) / 3.0;
    L2 = pow(L3, 2.0/3 );
    L1 = pow(L3, 1.0/3 );

    kstep = 2*pi / L1;

    // Calculating the actual one-body energy for given e=nx^2+ny^2+nz^2
    for (int i=0; i<nstates; i++){
        states(i,0) = states(i,0) * 2*pi*pi / L2;
    }
}

double basis_set::ReferenceEnergy(){

    double Energy = 0.0;

    if (Nparticles <= nstates){

        for (int p=0; p<Nparticles; p++){
            Energy += OneBodyOperator(p,p);

            for (int q=0; q<Nparticles; q++){
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

double basis_set::epsilon(int i, int j, int a, int b){
    // Function to compute the sum of h(i) + h(j) - h(a) - h(b)
    return OneBodyOperator(i,i) + OneBodyOperator(j,j) - OneBodyOperator(a,a) - OneBodyOperator(b,b);
}

double basis_set::epsilon4(int i, int j, int k, int l, int a, int b, int c, int d){
    // Function to compute the sum of h(i) + h(j) - h(a) - h(b)
    return OneBodyOperator(i,i) + OneBodyOperator(j,j) + OneBodyOperator(k,k) + OneBodyOperator(l,l) - OneBodyOperator(a,a) - OneBodyOperator(b,b) - OneBodyOperator(c,c) - OneBodyOperator(d,d);
}


double basis_set::OneBodyOperator(int p, int q){
    // Need to return the one-body energy value for p if p == q
    return ( p == q ) ? states(p,0) : 0;
}

double basis_set::TwoBodyOperator(int p, int q, int r, int s){
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

int basis_set::KDelta_integer(int a, int b){
    // Returns 1 if a and b are equal
    return ( a == b ) ? 1:0;
}

int basis_set::KDelta_array(rowvec a, rowvec b){
    // Returns 1 if all elements are equal, else return 0
    int value = 1;
    for (int i=1; i<=3; i++){
        value = ( a(i) == b(i) ) ? value*1 : value*0;
    }
    return value;
}

int basis_set::KDelta_k(int a, int b){
    // Kroenecker delta for wave numbers ka and kb
    return KDelta_array( states.row(a), states.row(b));
}

int basis_set::KDelta_spin(int a, int b){
    // Kroenecker delta for spin integer ms_a and ms_b
    return KDelta_integer( states(a,4), states(b,4));
}

int basis_set::KDelta_sum(int a, int b, int c, int d){
    // Kroenecker delta comparing a+b and c+d
    return KDelta_array( states.row(a)+states.row(b) , states.row(c)+states.row(d));
}

double basis_set::Absolute_Difference2(int a, int b){
    // Returning difference squared between a and b
    double diff = 0.0;
    for (int i=1; i<=3; i++){
        diff += (states(a,i) - states(b,i)) * (states(a,i) - states(b,i));
    }
    return diff;
}










