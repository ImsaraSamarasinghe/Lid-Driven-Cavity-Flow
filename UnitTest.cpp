#include "SolverCG.h"
#include <cmath>
#include <mpi.h>
using namespace std;

#define BOOST_TEST_MODULE LidDrivenCavityTest
#include <boost/test/included/unit_test.hpp>
#define IDX(I,J) (J*(Nx)+I)

struct MPIFixture {
    public:
        explicit MPIFixture(){
            int initialized;
            MPI_Initialized(&initialized);
            if (!initialized) {
                argc = boost::unit_test::framework::master_test_suite().argc;
                argv = boost::unit_test::framework::master_test_suite().argv;
                cout << "Initialising MPI" << endl;
                MPI_Init(&argc,&argv);
            }
        }
        
        ~MPIFixture() {
            int finalized;
            MPI_Finalized(&finalized);
            if (!finalized) {
                cout << "Finalising MPI" << endl;
                MPI_Finalize();
            }
        }
        int argc;
        char **argv;
};
BOOST_GLOBAL_FIXTURE(MPIFixture);
			
BOOST_AUTO_TEST_CASE(InitialiseTest)
{	
	int Nx = 9;
	int Ny = 9;
	double Lx = 1.0;
	double Ly = 1.0;
	double dx = Lx / (Nx-1);
   double dy = Ly / (Ny-1);
   double tol = 1000.0; // tolerance
	double* stream = new double [Nx*Ny]; // array to store the analytical solution
	double* stream_cal = new double [Nx*Ny];
	double* vorticity = new double [Nx*Ny];
	
	SolverCG Solver(Nx,Ny,dx,dy); // create object of SolverCG class
	MPIFixture obj;
	// fill vorticity and stream arrays
	const int k=3;
	const int l=3;
	for (int i = 0; i < Nx; ++i){
		for (int j = 0; j < Ny; ++j) {
			vorticity[IDX(i,j)] = -M_PI*M_PI*(k*k+l*l)*sin(M_PI*k*i*dx)*sin(M_PI*l*j*dy);
			stream[IDX(i,j)] = sin(M_PI*k*i*dx)*sin(M_PI*l*j*dy);
		}
	}
	
	//cg->Solve(vorticity,stream_cal,argc,argv);
	Solver.Solve(vorticity,stream_cal,obj.argc,obj.argv); // run Solve function
	// check arrays
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j){
			BOOST_CHECK_CLOSE(stream[IDX(i,j)],stream_cal[IDX(i,j)],tol); // run the check
		}
	}
	delete[] vorticity;
	delete[] stream;
	delete[] stream_cal;
}
