#include "ccdblocks2.h"

CCDBlocks2::CCDBlocks2()
{

}

CCDBlocks2::CCDBlocks2(basis_set BASIS){

    basis = BASIS;

    Nholes = basis.Nparticles;
    Nstates = basis.nstates;
    Nparticles = Nstates - Nholes;

    //cout << Nholes << " " << Nstates << "  " << Nparticles << endl;

    T = zeros<mat>(0,3); int i=0;
    X = zeros<mat>(0,3);
    t = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);
    t0 = zeros<mat>(Nparticles*Nparticles,Nholes*Nholes);

    Nmax = basis.Nshells;
    m = 2 * floor(sqrt(Nmax));
    M = 2*m + 1;

    //basis.states.print();


    for (int p=0; p<Nstates; p++){
        for (int q=0; q<Nstates; q++){

            if (p != q) { // Pauli exclusion satisfied

                if ( basis.states(p,0) == basis.states(q,0)){ // No broken pairs

                    T.insert_rows(i,1);
                    X.insert_rows(i,1);

                    int P = basis.states(p,0) + basis.states(q,0);
                    int Sz = basis.states(p,1) + basis.states(q,1);

                    T(i,0) = p; T(i,1) = q; T(i,2) = Index(P, Sz);

                    P = basis.states(p,0) - basis.states(q,0);
                    Sz = basis.states(p,1) - basis.states(q,1);

                    X(i,0) = p; X(i,1) = q; X(i,2) = Index(P,Sz);

                    //cout << "p: " << T(i,0) << "   q: " << T(i,1) << "   index: " << T(i,2) << "   X-matrix index: " << X(i,2) << endl;
                    i++;
                }
            }

        }
    }
    // Now the channel, T, has been created. The cross-channel, X, has also been created. I do not know what I need X for yet.
    // The channel T consist of all the two-body configurations that can exist. No other configurations will give a non-zero contribution.
    // Therefore I will now loop over all elements in channel T when I calculate two-body elements.

    Nchannels = T.n_rows;
    //cout << "Nchannels: " << Nchannels << " i: " << i << endl;

    for (int P=0; P<Nchannels; P++){
        for (int R=0; R<Nchannels; R++){

            //cout << "here I take twobody product of " << T(P,0) << " " << T(P,1) << " and " << T(R,0) << " " << T(R,1) << endl;
        }
    }
}


int CCDBlocks2::Index(int P, int Sz){
    return 2*(P+m)*M + 2*(Sz+1);
}

double CCDBlocks2::Abs(double a, double b){
    return sqrt( pow(a-b, 2) );
}

double CCDBlocks2::CCD(){
    int n = 0;
    double error = 0.00001;

    update_t();

    double Eold = 0.0;
    double Enew = CorrolationEnergy();

    while ( n < 10 && Abs(Enew,Eold) > error){

        update_t();
        Eold = Enew;
        Enew = CorrolationEnergy();

        n ++;
    }
    cout << "Coupled Cluster method used " << n << " iterations." << endl;
    return Enew;
}

double CCDBlocks2::CorrolationEnergy(){
    double E = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    E += basis.TwoBodyOperator(i,j,a,b) * t(aa+bb*Nparticles, i+j*Nholes);
                }
            }
        }
    }
    return E / 4.0;
}

void CCDBlocks2::update_t(){
    t0 = t;

    for (int P=0; P<Nchannels; P++){
        for (int R=0; R<Nchannels; R++){
            double tau = 0;

            int p = T(P,0); int q = T(P,1);
            int r = T(R,0); int s = T(R,1);

            t() += basis.TwoBodyOperator(a,b,i,j) +
        }
    }

}









