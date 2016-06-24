#include "pairingbasis.h"

PairingBasis::PairingBasis() { }

PairingBasis::PairingBasis(int Nshells_input, int NshellsFilled_input, double g_input, double delta_input){

    Nshells = Nshells_input;
    NshellsFilled = NshellsFilled_input;
    g = g_input;
    delta = delta_input;

    // Setting up all states in basis
    States = zeros<mat>(0,2);
    Nstates = 0; Nholes = 0;

    for (int p=1; p <= Nshells; p++){

        for (int ms=-1; ms <= 1; ms += 2){

            States.insert_rows(Nstates,1);

            States(Nstates,0) = p;
            States(Nstates,1) = ms;

            Nstates ++;
            if (p <= NshellsFilled) Nholes++;
        }
    }
    Nparticles = Nstates - Nholes;
    SetUpEpsilon();
}

double PairingBasis::ReferenceEnergy(){
    double E0 = 0;
    if (Nholes < Nstates){
        for (int p=0; p<Nholes; p++){
            E0 += OneBodyOperator(p,p);

            for (int q=0; q<Nholes; q++){
                E0 += 0.5 * TwoBodyOperator(p,q,p,q);
            }
        }
    }
    else{
        cout << "Nholes cannot exceed Nstates." << endl;
    }
    return E0;
}

double PairingBasis::OneBodyOperator(int p, int q){
    double obo = 0.0;
    if (KD_state(p,q) == 1){
        obo = delta*( States(p,0) - 1 );
    }
    return obo;
}

double PairingBasis::TwoBodyOperator(int p, int q, int r, int s){
    double tbo = 0.0;

    if ( KD_state(p,q)*KD_state(r,s) == 1 ){

        if ( KD_spin(p,q) == 0 ){
            if (KD_spin(r,s) == 0){

                if ( KD_spin(p,r) == 1 && KD_spin(q,s) == 1 ){
                    tbo = -g/2.0;
                }
                if ( KD_spin(p,s) == 1 && KD_spin(q,r) == 1 ){
                    tbo = g/2.0;
                }
            }
        }
    }
    //return tbo;
    tbo = 0;
    if ( States(p,0) == States(q,0) && States(r,0) == States(s,0) ){

        if ( States(p,1) != States(q,1) && States(r,1) != States(s,1) ){

            if (States(p,1) == States(r,1) && States(q,1) == States(s,1)) tbo =-g/2.0;
            if (States(p,1) == States(s,1) && States(q,1) == States(r,1)) tbo = g/2.0;
        }
    }
    return tbo;
}


double PairingBasis::ei(int q){
    double interaction = 0;
    for (int i=0; i<Nholes; i++){
        interaction += TwoBodyOperator(q,i,q,i);
    }
    return OneBodyOperator(q,q) + interaction;
}

double PairingBasis::epsilon(int i, int j, int a, int b){
    // Function to compute the sum of h(i) + h(j) - h(a) - h(b)
    return ei(i) + ei(j) - ei(a) - ei(b);
}

void PairingBasis::SetUpEpsilon(){
    EpsilonMatrix = zeros<vec>(Nparticles*Nparticles*Nholes*Nholes);

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    int a = aa + Nholes; int b = bb + Nholes;
                    EpsilonMatrix(aa + bb*Nparticles + i*Nparticles*Nparticles + j*Nparticles*Nparticles*Nholes) = epsilon(i,j,a,b);
                }
            }
        }
    }
}

double PairingBasis::epsilon4(int i, int j, int k, int l, int a, int b, int c, int d){
    return ei(i) + ei(j) + ei(k) + ei(l) - ei(a) - ei(b) - ei(c) - ei(d);
}

int PairingBasis::KD_integer(int a, int b){
    // Kroenecker delta for integers a and b;
    return (a == b) ? 1 : 0;
}

int PairingBasis::KD_state(int a, int b){
    // Kroenecker delta for states a and b;
    return KD_integer( States(a,0), States(b,0) );
}

int PairingBasis::KD_spin(int a, int b){
    // Kroenecker delta for spin value of state a and b
    return KD_integer( States(a,1), States(b,1) );
}
