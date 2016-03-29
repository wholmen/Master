#include "basis_set2d.h"

basis_set2D::basis_set2D()
{

}

basis_set2D::basis_set2D(int nparticles, int nshells, double RS){
    Nparticles = nparticles; Nshells = nshells; rs = RS;

    states = zeros<mat>(0,4);

    nstates = 0;

    for (int shell=0; shell<Nshells; shell++){
        // shell is the value x^2 + y^2 must match

        for (int nx=-shell; nx <= shell; nx++){
            for (int ny=-shell; ny <= shell; ny++){
                // Looping through all possible combinations of nx, ny, nz below n_shell

                // Energy integer
                int e = nx*nx + ny*ny;

                // Removing duplicates from previous values for shell
                if (e == shell){

                    // Adding spin
                    for (int ms = -1; ms <= 1; ms+=2){

                        states.insert_rows(nstates,1);

                        states(nstates,0) = e;
                        states(nstates,1) = nx;
                        states(nstates,2) = ny;
                        states(nstates,3) = ms;

                        nstates++;
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

    // Calculating the actual one-body energy for given e=nx^2+ny^2
    for (int i=0; i<nstates; i++){
        states(i,0) = states(i,0) * 2*pi*pi / L2;
    }
}


double basis_set2D::ReferenceEnergy(){

    double Energy = 0.0;

    for (int p=0; p<Nparticles; p++){

        for (int q=0; q<Nparticles; q++){
            Energy += OneBodyOperator(p,q);
        }
        for (int q=0; q<Nparticles; q++){
            for (int r=0; r<Nparticles; r++){
                for (int s=0; s<Nparticles; s++){
                    if (p != q && r != s){
                        Energy += 0.5*TwoBodyOperator(p,q,r,s);
                    }
                }
            }
        }
    }
    return Energy;
}

double basis_set2D::ReferenceEnergy2(){

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


double basis_set2D::OneBodyOperator(int p, int q){
    // Need to return the one-body energy value for p if p == q
    return ( p == q ) ? states(p,0) : 0;
}

double basis_set2D::TwoBodyOperator(int p, int q, int r, int s){
    // Two body interaction for electron gas

    double asym, asym1, asym2;
    asym = 0.0; asym1 = 0.0; asym2 = 0.0;

    if (KDelta_sum(p,q,r,s) == 1){
        asym = 1.0/L3; //4*pi/L3;

        if ( KDelta_spin(p,r)*KDelta_spin(q,s) == 1 ){

            if (KDelta_k(p,r) != 1){
                asym1 = 1.0 / (Absolute_Difference2(r,p));
                //asym1 = L2 / (pi*Absolute_Difference2(r,p));
            }
        }
        if ( KDelta_spin(p,s)*KDelta_spin(q,r) == 1 ){

            if (KDelta_k(p,s) != 1){
                asym2 = 1.0 / (Absolute_Difference2(s,p));
                //asym2 = L2 / (pi*Absolute_Difference2(s,p));
            }

        }
    }
    return asym * (asym1 - asym2);
}

int basis_set2D::KDelta_integer(int a, int b){
    // Returns 1 if a and b are equal
    return ( a == b ) ? 1:0;
}

int basis_set2D::KDelta_array(rowvec a, rowvec b){
    // Returns 1 if all elements are equal, else return 0
    int value = 1;
    for (int i=1; i<=2; i++){
        value = ( a(i) == b(i) ) ? value*1 : value*0;
    }
    return value;
}

int basis_set2D::KDelta_k(int a, int b){
    // Kroenecker delta for wave numbers ka and kb
    return KDelta_array( states.row(a), states.row(b));
}

int basis_set2D::KDelta_spin(int a, int b){
    // Kroenecker delta for spin integer ms_a and ms_b
    return KDelta_integer( states(a,3), states(b,3));
}

int basis_set2D::KDelta_sum(int a, int b, int c, int d){
    // Kroenecker delta comparing a+b and c+d
    return KDelta_array( states.row(a)+states.row(b) , states.row(c)+states.row(d));
}

double basis_set2D::Absolute_Difference2(int a, int b){
    // Returning difference squared between a and b
    double diff = 0.0;
    for (int i=1; i<=2; i++){
        diff += (states(a,i) - states(b,i)) * (states(a,i) - states(b,i));
    }
    return diff;
}










