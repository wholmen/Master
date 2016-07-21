#include <iostream>
#include <omp.h>
#include <iomanip>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    /*
    for (int NTHREADS=1; NTHREADS <= 4; NTHREADS++){

        int Nsteps = 10000000000;
        double step = 1.0 / double(Nsteps);
        double pi = 0;

        double time0 = omp_get_wtime();



        omp_set_num_threads(NTHREADS);
        #pragma omp parallel
        {

            double x;
            int id = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            //cout << nthreads << endl;
            double sum = 0;
            for (int i=id; i<Nsteps; i+=nthreads){
                x = (i + 0.5)*step;
                sum += 4 / (1.0+x*x);
            }

            #pragma omp atomic
            pi += sum*step;

        }
        double time1 = omp_get_wtime();

        cout << pi << " Nthreads: " << NTHREADS << " time: " << time1-time0 << endl;
    }*/

    /*
    // First test
    for (int NTHREADS=1; NTHREADS <= 4; NTHREADS++){
        int Nholes = 300;
        mat Holes = zeros<mat>(0,2);

        double time0 = omp_get_wtime();
        omp_set_num_threads(NTHREADS);

#pragma omp parallel
        {
            int id = omp_get_thread_num();
            int nthreads = omp_get_num_threads();

            mat smallHole = zeros<mat>(0,2);

            int n=0;
            for (int i=id; i<Nholes; i += nthreads){
                for (int j=0; j<Nholes; j++){

                    if (i != j){ // Pauli principle demands that the particles must be unequal

                        // Adding a new two-hole-state configuration to matrix. (i, j, Identifier)
                        smallHole.insert_rows(n,1);
                        smallHole(n, 0) = i;
                        smallHole(n, 1) = j;
                        n++;
                    }
                }
            }
#pragma omp critical
            Holes.insert_rows(0,smallHole);

        }
        double time1 = omp_get_wtime();
        //Holes.print();
        cout << " Nthreads: " << NTHREADS << " time: " << time1-time0 << endl;
    }
    */



    // Second test
    for (int NTHREADS=1; NTHREADS <= 4; NTHREADS++){
        int Nholes = 10000;
        mat Holes = zeros<mat>(0,2);

        double time0 = omp_get_wtime();
        omp_set_num_threads(NTHREADS);

#pragma omp parallel
        {
            int id = omp_get_thread_num();
            int nthreads = omp_get_num_threads();

            int size = floor( Nholes/nthreads) * (Nholes-1);
            if ( id < Nholes%nthreads) size += Nholes - 1;

            mat smallHole = zeros<mat>( size, 2 );

            int n=0;
            for (int i=id; i<Nholes; i += nthreads){
                for (int j=0; j<Nholes; j++){

                    if (i != j){ // Pauli principle demands that the particles must be unequal

                        // Adding a new two-hole-state configuration to matrix. (i, j, Identifier)
                        smallHole(n, 0) = i;
                        smallHole(n, 1) = j;
                        n ++;
                    }
                }
            }
#pragma omp critical
            Holes.insert_rows(0,smallHole);

        }
        double time1 = omp_get_wtime();
        //Holes.print();
        cout << " Nthreads: " << NTHREADS << " time: " << time1-time0 << endl;
    }



    /*
    int Nholes = 4000;
    mat Holes = zeros<mat>( Nholes*Nholes ,3);

    double time0 = omp_get_wtime();

    omp_set_num_threads(1);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int threads = omp_get_num_threads();

        for (int i=id; i<Nholes; i += threads){
            for (int j=0; j<Nholes; j++){

                if (i != j){ // Pauli principle demands that the particles must be unequal

                    // Adding a new two-hole-state configuration to matrix. (i, j, Identifier)
                    Holes(i + j*Nholes, 0) = i;
                    Holes(i + j*Nholes, 1) = j;

                }
            }
        }
    }
    double time1 = omp_get_wtime();
    cout << time1 - time0 << endl;

    time0 = omp_get_wtime();
    for (double i=0; i<Holes.n_rows; i+=100) if ( Holes(i,0) == Holes(i,1) ) Holes.shed_row(i);
    time1 = omp_get_wtime();

    cout << time1 - time0 << endl;

    //Holes.print();
    */
}

