#include "basis_set.h"

basis_set::basis_set() { }

basis_set::basis_set(int nparticles, int nshells, double G, double Delta){
    Nparticles = nparticles;
    Nshells = nshells;
    g = G;
    delta = Delta;

    // Setting up all states in basis
    states = zeros<mat>(0,2);
    nstates = 0;

    for (int p=1; p <= Nshells; p++){

        for (int ms=-1; ms <= 1; ms += 2){

            states.insert_rows(nstates,1);

            states(nstates,0) = p;
            states(nstates,1) = ms;

            nstates ++;
        }
    }
}

double basis_set::ReferenceEnergy(){
    double E0 = 0;
    if (Nparticles < nstates){
        for (int p=0; p<Nparticles; p++){
            E0 += OneBodyOperator(p,p);

            for (int q=0; q<Nparticles; q++){
                E0 += 0.5 * TwoBodyOperator(p,q,p,q);
            }
        }
    }
    else{
        cout << "Nparticles cannot exceed Nstates." << endl;
    }
    return E0;
}

double basis_set::OneBodyOperator(int p, int q){
    double obo = 0.0;
    if (KD_state(p,q) == 1){
        obo = delta*( states(p,0) - 1 );
    }
    return obo;
}

double basis_set::TwoBodyOperator(int p, int q, int r, int s){
    double tbo = 0.0;

    if ( KD_state(p,q)*KD_state(r,s) == 1 ){

        if ( KD_spin(p,q) == 0 ){
            if (KD_spin(r,s) == 0){

                if ( KD_spin(p,r) == 1){
                    tbo = -g/2.0;
                }
                if ( KD_spin(p,r) == 0){
                    tbo = g/2.0;
                }
            }
        }
    }
    return tbo;
}

double basis_set::epsilon(int i, int j, int a, int b){
    // Function to compute the sum of h(i) + h(j) - h(a) - h(b)
    return OneBodyOperator(i,i) + OneBodyOperator(j,j) - OneBodyOperator(a,a) - OneBodyOperator(b,b);
}

double basis_set::epsilon4(int i, int j, int k, int l, int a, int b, int c, int d){
    // Function to compute the sum of h(i) + h(j) - h(a) - h(b)
    return OneBodyOperator(i,i) + OneBodyOperator(j,j) + OneBodyOperator(k,k) + OneBodyOperator(l,l) - OneBodyOperator(a,a) - OneBodyOperator(b,b) - OneBodyOperator(c,c) - OneBodyOperator(d,d);
}

int basis_set::KD_integer(int a, int b){
    // Kroenecker delta for integers a and b;
    return (a == b) ? 1 : 0;
}

int basis_set::KD_state(int a, int b){
    // Kroenecker delta for states a and b;
    return KD_integer( states(a,0), states(b,0) );
}

int basis_set::KD_spin(int a, int b){
    // Kroenecker delta for spin value of state a and b
    return KD_integer( states(a,1), states(b,1) );
}
