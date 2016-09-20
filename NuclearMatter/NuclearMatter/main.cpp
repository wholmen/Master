#include <iostream>
#include <nuclearbasis.h>
#include <solver.h>
#include <../../Solvers/Solvers/ccdintermediates.h>
#include <iomanip>
#include <omp.h>

using namespace std;


void TestBlocks();

void Results_Nh14_ThermodynamicLimit();
void Results_Nh14_Ns20_vary_n();
void Results_Nh36_Ns20_vary_n();
void Results_Nh54_Ns20_vary_n();

int main()
{

    //TestBlocks();

    //IResults_Nh14_Ns20_vary_n();
    //Results_Nh36_Ns20_vary_n();
    Results_Nh54_Ns20_vary_n();
    Results_Nh14_ThermodynamicLimit();
}


void TestBlocks(){
    int Nshells = 8; int NfilledShells = 2; double n = 0.2;

    NuclearBasis basis = NuclearBasis(Nshells,NfilledShells,n,true);
    cout << "Reference Energy: " << basis.ReferenceEnergy() << endl;
    Solver solver = Solver(basis);
    solver.weight = 0.3;
    cout << "Nparticles: " << solver.Nholes << " Nstates: " << solver.Nstates << " Energy:" << solver.CCD(200) << endl;
}

void Results_Nh14_ThermodynamicLimit(){
    ofstream myfile;
    myfile.open("../Results/ThermodynamicLimit_Nh14_n0.2.txt");
    myfile << "A coupled cluster study of electron gas using a parallel code with blocks implementation for different amount of states" << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "density [n], States, Occupied States, weight, Reference Energy, Iterations, time spent initiate, time spent iterating, CCD correlation Energy, total energy per particle" << endl;


    int Nfilled = 2; double n = 0.20;
    for (int Nshells=3; Nshells<35; Nshells++){

        if (Nshells != 8 && Nshells != 16 && Nshells != 24 && Nshells != 29 && Nshells != 32){

            NuclearBasis basis = NuclearBasis(Nshells,Nfilled,n,true);

            double time_init_0 = omp_get_wtime();
            Solver solve = Solver(basis);
            double time_init_1 = omp_get_wtime();

            solve.weight = 0.3; solve.tolerance = 1e-6;

            double time_iter_0 = omp_get_wtime();
            double E = solve.CCD(200);
            double time_iter_1 = omp_get_wtime();

            double EN = (basis.ReferenceEnergy()+E) / basis.Nholes;
            myfile << setprecision(6) << basis.n << " " << basis.Nstates << " " << basis.Nholes << " " << solve.weight << " " << basis.ReferenceEnergy() << " " << solve.NIterations << " " << time_init_1 - time_init_0 << " " << time_iter_1 - time_iter_0 << " " << E << " " << EN << endl;
        }
    }
}


void Results_Nh14_Ns20_vary_n(){

    int Nshells = 20; int NfilledShells = 2;

    vec n = zeros<vec>(14);

    for (int i=0; i<14; i++){
        n(i) = 0.025 + i*0.025;
    }

    ofstream myfile;
    myfile.open("../Results/Solver_Nh14_varyRho_Ns20.txt");
    myfile << "A coupled cluster study of neutrom matter using a parallel code with blocks implementation for different amount of states" << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "density [n], States, Occupied States, weight, Reference Energy, Iterations, time spent initiate, time spent iterating, CCD correlation Energy, total energy per particle" << endl;


    for (int i=0; i<14; i++){

        double rho = n(i);
        NuclearBasis basis = NuclearBasis(Nshells,NfilledShells, rho, true);

        double time_init_0 = omp_get_wtime();
        Solver solve = Solver(basis);
        double time_init_1 = omp_get_wtime();

        solve.weight = 0.3; solve.tolerance = 1e-6;

        double time_iter_0 = omp_get_wtime();
        double E = solve.CCD(200);
        double time_iter_1 = omp_get_wtime();

        double EN = (basis.ReferenceEnergy()+E) / basis.Nholes;
        myfile << setprecision(6) << basis.n << " " << basis.Nstates << " " << basis.Nholes << " " << solve.weight << " " << basis.ReferenceEnergy() << " " << solve.NIterations << " " << time_init_1 - time_init_0 << " " << time_iter_1 - time_iter_0 << " " << E << " " << EN << endl;
    }

}


void Results_Nh36_Ns20_vary_n(){

    int Nshells = 20; int NfilledShells = 3;

    vec n = zeros<vec>(14);

    for (int i=0; i<14; i++){
        n(i) = 0.025 + i*0.025;
    }

    ofstream myfile;
    myfile.open("../Results/Solver_Nh36_varyRho_Ns20.txt");
    myfile << "A coupled cluster study of neutrom matter using a parallel code with blocks implementation for different amount of states" << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "density [n], States, Occupied States, weight, Reference Energy, Iterations, time spent initiate, time spent iterating, CCD correlation Energy, total energy per particle" << endl;


    for (int i=0; i<14; i++){

        double rho = n(i);
        NuclearBasis basis = NuclearBasis(Nshells,NfilledShells, rho, true);

        double time_init_0 = omp_get_wtime();
        Solver solve = Solver(basis);
        double time_init_1 = omp_get_wtime();

        solve.weight = 0.3; solve.tolerance = 1e-6;

        double time_iter_0 = omp_get_wtime();
        double E = solve.CCD(200);
        double time_iter_1 = omp_get_wtime();

        double EN = (basis.ReferenceEnergy()+E) / basis.Nholes;
        myfile << setprecision(6) << basis.n << " " << basis.Nstates << " " << basis.Nholes << " " << solve.weight << " " << basis.ReferenceEnergy() << " " << solve.NIterations << " " << time_init_1 - time_init_0 << " " << time_iter_1 - time_iter_0 << " " << E << " " << EN << endl;
    }
}

void Results_Nh54_Ns20_vary_n(){

    int Nshells = 15; int NfilledShells = 4;

    vec n = zeros<vec>(14);

    for (int i=0; i<14; i++){
        n(i) = 0.025 + i*0.025;
    }

    ofstream myfile;
    myfile.open("../Results/Solver_Nh54_varyRho_Ns20.txt");
    myfile << "A coupled cluster study of neutrom matter using a parallel code with blocks implementation for different amount of states" << endl;
    myfile << "The results are presented in the following order: " << endl;
    myfile << "density [n], States, Occupied States, weight, Reference Energy, Iterations, time spent initiate, time spent iterating, CCD correlation Energy, total energy per particle" << endl;


    for (int i=0; i<14; i++){

        double rho = n(i);
        NuclearBasis basis = NuclearBasis(Nshells,NfilledShells, rho, true);

        double time_init_0 = omp_get_wtime();
        Solver solve = Solver(basis);
        double time_init_1 = omp_get_wtime();

        solve.weight = 0.3; solve.tolerance = 1e-6;

        double time_iter_0 = omp_get_wtime();
        double E = solve.CCD(200);
        double time_iter_1 = omp_get_wtime();

        double EN = (basis.ReferenceEnergy()+E) / basis.Nholes;
        myfile << setprecision(6) << basis.n << " " << basis.Nstates << " " << basis.Nholes << " " << solve.weight << " " << basis.ReferenceEnergy() << " " << solve.NIterations << " " << time_init_1 - time_init_0 << " " << time_iter_1 - time_iter_0 << " " << E << " " << EN << endl;
    }
}


