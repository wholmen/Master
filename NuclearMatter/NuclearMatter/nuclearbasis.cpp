#include "nuclearbasis.h"

NuclearBasis::NuclearBasis() {}

NuclearBasis::NuclearBasis(int Nshells_input, int NFilledShells_input, double rs_input, bool OnlyNeutrons){

    Nshells = Nshells_input; NfilledShells = NFilledShells_input; rs = rs_input;

    v0R = 200; v0T = 178; v0S = 91.85; kappaR = 1.487; kappaT = 0.649; kappaS = 0.465;

    States = zeros<mat>(0,6);
    Nstates = 0; Nholes = 0;

    int IsospinMax = (OnlyNeutrons) ? 0:1; // If OnlyNeutrons = true, then there will only be one

    for (int shell=0; shell<Nshells; shell++){
        // Shell is defined by x^2 + y^2 + z^2 as these quantum numbers give energy.

        for (int nx=-shell; nx <= shell; nx++){
            for (int ny=-shell; ny <= shell; ny++){
                for (int nz=-shell; nz <= shell; nz++){
                    // Looping over all possible combinations of nx, ny and nz all below shell n.

                    // Energy integer
                    int e = nx*nx + ny*ny + nz*nz;

                    // Only look at values on this exact shell, and none below.
                    if (e == shell){

                        // Looping over spin
                        for (int ms=-1; ms<=1; ms+=2){

                            // Looping over isospin
                            for (int mt=-1; mt <= IsospinMax; mt+=2){

                                States.insert_rows(Nstates,1);
                                States(Nstates,0) = e;
                                States(Nstates,1) = nx;
                                States(Nstates,2) = ny;
                                States(Nstates,3) = nz;
                                States(Nstates,4) = ms;
                                States(Nstates,5) = mt;

                                Nstates ++;
                                if (shell < NfilledShells) Nholes++;
                            }
                        }
                    }
                }
            }
        }
    }
    // Various calculations of variables needed
    L3 = (4*pi*rs*rs*rs*Nholes) / 3.0;
    L2 = pow(L3, 2.0/3 );
    L1 = pow(L3, 1.0/3 );

    kstep = 2*pi / L1;

    // Calculating the actual one-body energy for given e=nx^2+ny^2+nz^2
    for (int i=0; i<Nstates; i++){
        States(i,0) = States(i,0) * 2*pi*pi / L2;
    }

    Nstates2 = Nstates*Nstates; Nstates3 = Nstates2*Nstates;
    Nparticles = Nstates - Nholes;
    Nparticles2 = Nparticles*Nparticles; Nholes2 = Nholes*Nholes;
}

double NuclearBasis::ReferenceEnergy(){

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
        cout << "Number of particles cannot exceed number of available States.";
    }
    return Energy;
}

double NuclearBasis::OneBodyOperator(int p, int q){
    // Need to return the one-body energy value for p if p == q
    return ( p == q ) ? States(p,0) : 0;
}

double NuclearBasis::Minnesota(int p, int q, int r, int s){
    // Two body interaction for electron gas

    double V = 0.0;

    if ( KD_sum(p,q,r,s) ) { // Testing wether kp + kq = kr + ks. Kroenecker delta that persist through whole calculation

        // Setting up vectors containing wave numbers k for all states
        vec kp = zeros<vec>(3);
        vec kq = zeros<vec>(3);
        vec kr = zeros<vec>(3);
        vec ks = zeros<vec>(3);

        for (int i=0; i<3; i++){
            kp(i) = States(p,i+1);
            kq(i) = States(q,i+1);
            kr(i) = States(r,i+1);
            ks(i) = States(s,i+1);
        }
        kp = kp*2*pi/L1;
        kq = kq*2*pi/L1;
        kr = kr*2*pi/L1;
        ks = ks*2*pi/L1;

        vec p1 = 0.5*(kp - kq);
        vec p2 = 0.5*(kr - ks);

        double qsquared = dot(p1-p2,p1-p2);

        double VR =  v0R * pow(pi/kappaR, 1.5) / L3 * exp( -qsquared/(4*kappaR) );
        double VT = -v0T * pow(pi/kappaT, 1.5) / L3 * exp( -qsquared/(4*kappaT) );
        double VS = -v0S * pow(pi/kappaS, 1.5) / L3 * exp( -qsquared/(4*kappaS) );

        double ExchangeSpin = KD_spin(p,s)*KD_spin(q,r); // Spin exchanged
        double ConserveSpin = KD_spin(p,r)*KD_spin(q,s); // Spin conserved

        double ExchangeIsoSpin = KD_isospin(p,s)*KD_isospin(q,r); // Isospin exchanged
        double ConserveIsoSpin = KD_isospin(p,r)*KD_isospin(q,s); // Isospin conserved

        V =  0.50*(VR + 0.5*VT + 0.5*VS) * ConserveSpin * ConserveIsoSpin
           + 0.25*(VT - VS) * ExchangeSpin * ConserveIsoSpin
           - 0.50*(VR + 0.5*VT + 0.5*VS) * ExchangeSpin * ExchangeIsoSpin
           - 0.25*(VT - VS) * ConserveSpin * ExchangeIsoSpin;
    }
    return V;
}

double NuclearBasis::TwoBodyOperator(int p, int q, int r, int s){

    return Minnesota(p,q,r,s);
}

double NuclearBasis::ei(int q, vec &v){
    double interaction = 0;
    for (int i=0; i<Nholes; i++){
        interaction += v(q + i*Nstates + q*Nstates2 + i*Nstates3);
    }
    return OneBodyOperator(q,q) + interaction;
}

double NuclearBasis::epsilon(int i, int j, int a, int b, vec &v){
    // Function to compute the sum of h(i) + h(j) - h(a) - h(b)
    return ei(i,v) + ei(j,v) - ei(a,v) - ei(b,v);
}





// Various kroenecker deltas
int NuclearBasis::KD_integer(int a, int b){
    // Returns 1 if a and b are equal
    return ( a == b ) ? 1:0;
}

int NuclearBasis::KD_array(rowvec a, rowvec b){
    // Returns 1 if all elements are equal, else return 0
    int value = 1;
    for (int i=1; i<=3; i++){
        value = ( a(i) == b(i) ) ? value*1 : value*0;
    }
    return value;
}

int NuclearBasis::KD_spin(int a, int b){
    // Kroenecker delta for spin integer ms_a and ms_b
    return KD_integer( States(a,4), States(b,4));
}

int NuclearBasis::KD_isospin(int a, int b){
    // Kroenecker delta for isospin integer ms_a and ms_b
    return KD_integer( States(a,5), States(b,5));
}

int NuclearBasis::KD_sum(int a, int b, int c, int d){
    // Kroenecker delta comparing a+b and c+d
    return KD_array( States.row(a)+States.row(b) , States.row(c)+States.row(d));
}

int NuclearBasis::KD_k(int a, int b){
    // Kroenecker delta for wave numbers ka and kb
    return KD_array( States.row(a), States.row(b));
}
