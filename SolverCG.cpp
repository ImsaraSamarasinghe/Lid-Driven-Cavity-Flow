#include <iostream>
#include <algorithm>
#include <cstring>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

#include <cblas.h>

#include "SolverCG.h"

#define IDX(I,J) ((J)*Nx + (I))

SolverCG::SolverCG(int pNx, int pNy, double pdx, double pdy)
{
    dx = pdx;
    dy = pdy;
    Nx = pNx;
    Ny = pNy;
    int n = Nx*Ny;
    r = new double[n];
    p = new double[n];
    z = new double[n];
    t = new double[n]; //temp
}


SolverCG::~SolverCG()
{
    delete[] r;
    delete[] p;
    delete[] z;
    delete[] t;
}


void SolverCG::Solve(double* b, double* x,int argc,char **argv) {
    unsigned int n = Nx*Ny; // note: nis the total number of elements
    int k;

    double tol = 0.001;
    int size,rank; // initialising rank and size for MPI
    
    // MPI intialise
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    
    // initialise variables
    /* Of type int, global and local */
    double alpha_local_num;
    double alpha_local_den;
    double alpha_global_num;
    double alpha_global_den;
    double alpha_global;
    
    double beta_local_num;
    double beta_local_den;
    double beta_global_num;
    double beta_global_den;
    double beta_global;
    
    double eps_local;
    double eps_global;
    //

	// Include a test to check the appropriate number of processors

    // Calculate the number of columns each process will handle
    int BigRows = Nx / size; // Number of rows each process gets
    int SmallRows = Nx % size; // Number of processes that get an extra row
    
    // calculate the start and end of each sub-domain based on the rank
    int start_row, end_row;
    if (rank < SmallRows) {
        // Processes with rank less than SmallRows handle an extra column
        start_row = rank * (BigRows + 1);
        end_row = start_row + BigRows;
    } else {
        // Processes with rank greater than or equal to SmallRows handle BigRows columns
        start_row = rank * BigRows + SmallRows;
        end_row = start_row + BigRows - 1;
    }
    int local_rows = end_row-start_row+1; // number of rows in local arrays
    
    // --- Initialise arrays ---
    double* b_local = new double [Nx*local_rows]; // Nx x local_rows
    double* x_local = new double [Nx*local_rows]; // Nx x local_rows
    double* r_local = new double [Nx*local_rows]; // Nx x local_rows
    double* p_local = new double [Nx*local_rows]; // Nx x local_rows
    double* z_local = new double [Nx*local_rows]; // Nx x local_rows
    double* t_local = new double [Nx*local_rows]; // Nx x local_rows
    // --------------------------
    
    // ---- Initialise arrays for sending data ----
    /*These arrays are used for sending data between the row-partioned
    subdomains. Since each sub-domain only requires the next row they
    are 1D arrays*/
    double* xsendT = new double [Nx]; // Nx
    double* xsendB = new double [Nx]; // Nx
    double* xrecvT = new double [Nx]; // Nx
    double* xrecvB = new double [Nx]; // Nx
    
    double* psendT = new double [Nx]; // Nx
    double* psendB = new double [Nx]; // Nx
    double* precvT = new double [Nx]; // Nx
    double* precvB = new double [Nx]; // Nx
    // ----------------------------------------------
    
    // synchronize
    MPI_Barrier(MPI_COMM_WORLD);
    
    // broad cast initial inputs
    MPI_Bcast(b,n,MPI_DOUBLE,0,MPI_COMM_WORLD); // send b from root
    MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD); // send b from root

    /*Fill the local arrays with values from global arrays*/
    for (int i=0;i<Nx;++i){
        for (int j=start_row;j<=end_row;++j){
            b_local[(j-start_row)*Nx+i] = b[j*Nx+i];
            x_local[(j-start_row)*Nx+i] = x[j*Nx+i];
            r_local[(j-start_row)*Nx+i] = r[j*Nx+i];
            p_local[(j-start_row)*Nx+i] = p[j*Nx+i];
            z_local[(j-start_row)*Nx+i] = z[j*Nx+i];
            t_local[(j-start_row)*Nx+i] = t[j*Nx+i];
        }
    }
    
    //###### Finding eps global #######
    eps_local = cblas_dnrm2(Nx*local_rows, b_local, 1); // find eps_local from b_local
    eps_local = eps_local*eps_local;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&eps_local,&eps_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); // sum to eps_global from all eps_locals
    eps_global = sqrt(eps_global);
    
    // remains unchanged - keep performing to x global in each rank (should end simultaneously then)
    if (eps_global < tol*tol) {
        std::fill(x, x+n, 0.0);
        cout << "Norm is " << eps_global << endl;
        MPI_Finalize();
        return;
    }
    
    MPI_Barrier(MPI_COMM_WORLD); // ######## synchronize here
    Send_and_Receive(x_local,xsendT,xsendB,xrecvT,xrecvB,rank,size,local_rows);// find the neigbouring rows
    ApplyOperator(x_local, t_local, xrecvT, xrecvB, rank, size,local_rows); // send the neighbouring rows 

    cblas_dcopy(Nx*local_rows, b_local, 1, r_local, 1);        // r_0 = b (i.e. b)
    ImposeBC(r_local,rank,size,local_rows);
    
    MPI_Barrier(MPI_COMM_WORLD);
    cblas_daxpy(Nx*local_rows, -1.0, t_local, 1, r_local, 1);
    Precondition(r_local, z_local, rank, size, local_rows); // change this function
    cblas_dcopy(Nx*local_rows, z_local, 1, p_local, 1);        // p_0 = r_0
    MPI_Barrier(MPI_COMM_WORLD); // ####### synchronise here
    //## set last row of p to zero
    if (rank==size-1){
    	for (int i=0;i<Nx;++i){
    		p_local[IDX(i,local_rows-1)] = 0;
 		}
    }
    // ### printing p_local and t_local ###
    //printSubArrayToFile(p_local,rank,0,Nx,local_rows,"p_local_outside_loop");
    //printSubArrayToFile(t_local,rank,0,Nx,local_rows,"t_local_outside_loop");
    k = 0;
    
    do {
        k++;

        Send_and_Receive(p_local,psendT,psendB,precvT,precvB,rank,size,local_rows); // find adjacent rows
        ApplyOperator(p_local, t_local, precvT, precvB, rank, size, local_rows); // apply operator
        // ## printing P_local and t_local
        //printSubArrayToFile(p_local,rank,k,Nx,local_rows,"p_local_inLoop");
        //printSubArrayToFile(t_local,rank,k,Nx,local_rows,"t_local_inLoop");
        /*// ##### FINDING ALPHA ######
        alpha_local_den = cblas_ddot(Nx*local_rows, t_local, 1, p_local, 1);  // alpha = p_k^T A p_k
        MPI_Barrier(MPI_COMM_WORLD); // ##### synchronize
        MPI_Allreduce(&alpha_local,&alpha_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); // find global alpha by summing all local alphas
        // ###### END ALPHA ########*/
        
        /*Finding Alpha*/
        alpha_local_den = cblas_ddot(Nx*local_rows,t_local,1,p_local,1);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&alpha_local_den,&alpha_global_den,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        alpha_local_num = cblas_ddot(Nx*local_rows,r_local,1,z_local,1);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&alpha_local_num,&alpha_global_num,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        alpha_global = alpha_global_num/alpha_global_den;
        // end find alpha
        
        // #### FINDING BETA - INITIAL ######
        beta_local_den  = cblas_ddot(Nx*local_rows, r_local, 1, z_local, 1);  // z_k^T r_k
        MPI_Barrier(MPI_COMM_WORLD); // ##### synchronize
        MPI_Allreduce(&beta_local_den,&beta_global_den,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			
        cblas_daxpy(Nx*local_rows,  alpha_global, p_local, 1, x_local, 1);  // x_{k+1} = x_k + alpha_k p_k
        cblas_daxpy(Nx*local_rows, -alpha_global, t_local, 1, r_local, 1); // r_{k+1} = r_k - alpha_k A p_k
        
        // ### finding eps ###
        MPI_Barrier(MPI_COMM_WORLD); // ##### synchronize
        eps_local = cblas_dnrm2(Nx*local_rows, r_local, 1);
        eps_local = eps_local*eps_local;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&eps_local,&eps_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); // find the global error in each process
        eps_global = sqrt(eps_global);
        // #######
        
        MPI_Barrier(MPI_COMM_WORLD);
        if (eps_global < tol*tol) { // each processor will break if eps_global < tol*tol
            break;
        }
        
        Precondition(r_local, z_local, rank, size, local_rows);
        
        // ##### FINDING BETA  - FINAL ######
        beta_local_num = cblas_ddot(Nx*local_rows, r_local, 1, z_local, 1);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&beta_local_num,&beta_global_num,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        beta_global = beta_global_num/beta_global_den;
        // ############
        
        cblas_dcopy(Nx*local_rows, z_local, 1, t_local, 1);
        cblas_daxpy(Nx*local_rows, beta_global, p_local, 1, t_local, 1);
        cblas_dcopy(Nx*local_rows, t_local, 1, p_local, 1);
        MPI_Barrier(MPI_COMM_WORLD);
        cout<<"-------------->"<<"Rank: "<<rank<<" Error: "<<eps_global<<" k:"<<k<<" alpha_global: "<<alpha_global<<"<----------------------"<<endl;
    } while (k < 5000); // Set a maximum number of iterations k =5000
    
    // recombine b_local and x_local to be printed ### only works for exactly divisible case
    MPI_Allgather(b_local,Nx*local_rows,MPI_DOUBLE,b,Nx*local_rows,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgather(x_local,Nx*local_rows,MPI_DOUBLE,x,Nx*local_rows,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (k == 5000) {
        cout << "FAILED TO CONVERGE" << endl;
        exit(-1);
    }
    
    delete[] b_local;
    delete[] x_local;
    delete[] r_local;
    delete[] p_local;
    delete[] z_local;
    delete[] t_local;
    
    delete[] xsendB;
    delete[] xsendT;
    delete[] xrecvB;
    delete[] xrecvT;
    
    delete[] psendB;
    delete[] psendT;
    delete[] precvB;
    delete[] precvT;
    
    if (rank==0) {
        cout << "Converged in " << k << " iterations. eps = " << eps_global << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


void SolverCG::ApplyOperator(double* in, double* out,double* top,double* bot,int& rank, int& size, int& local_rows) {
    // Assume ordered with y-direction fastest (column-by-column)
    // NOTE: i now means down a column and j means along a row (0<=i<=Ny-1 & 0<=j<=local_cols-1)
    int i,j;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
/*    int jm1 = 0, jp1 = 2; // original code
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx - 1; ++i) {
            out[IDX(i,j)] = ( -     in[IDX(i-1, j)]
                              + 2.0*in[IDX(i,   j)]
                              -     in[IDX(i+1, j)])*dx2i
                          + ( -     in[IDX(i, jm1)]
                              + 2.0*in[IDX(i,   j)]
                              -     in[IDX(i, jp1)])*dy2i;
        }
        jm1++;
        jp1++;
    
    }*/
    // check right most sub_domain
    if (size>1) {// only perform with more than 1 process
        if (rank==0){
            for (i = 1;i<Nx-1;++i){
                for (j = 1;j<local_rows-1;++j) { // all other values in rank 0
                    out[IDX(i,j)] = (-in[IDX(i-1, j)]+ 2.0*in[IDX(i,j)]-in[IDX(i+1, j)])*dx2i+ (-in[IDX(i, j-1)]+ 2.0*in[IDX(i,   j)]-in[IDX(i, j+1)])*dy2i;
                }
                // bottom boundary for rank 0 
                j = local_rows-1;
                out[IDX(i,j)] = (-in[IDX(i-1, j)]+ 2.0*in[IDX(i,j)]-in[IDX(i+1, j)])*dx2i+ (-in[IDX(i, j-1)]+ 2.0*in[IDX(i,   j)]-bot[IDX(i,0)])*dy2i;
            }
        }
        //check topmost boundary
        else if (rank==size-1){
            for (i = 1;i<Nx-1;++i){
                for (j = 1;j<local_rows-1;++j) { // all other values in rank 0
                    out[IDX(i,j)] = (-in[IDX(i-1,j)]+ 2.0*in[IDX(i,j)]-in[IDX(i+1,j)])*dx2i+(-in[IDX(i,j-1)]+ 2.0*in[IDX(i,j)]-in[IDX(i,j+1)])*dy2i;
                }
                // top boundary for last rank
                j = 0;
                out[IDX(i,j)] = (-in[IDX(i-1, j)]+ 2.0*in[IDX(i,j)]-in[IDX(i+1, j)])*dx2i+ (-top[IDX(i,0)]+ 2.0*in[IDX(i,   j)]-in[IDX(i, j+1)])*dy2i;
            }
        }
        // all other subdomains
        else {
            for (i = 1;i<Nx-1;++i){
                for (j = 1;j<local_rows-1;++j) { // all other values in rank 0
                    out[IDX(i,j)] = (-in[IDX(i-1,j)]+ 2.0*in[IDX(i,j)]-in[IDX(i+1,j)])*dx2i+(-in[IDX(i,j-1)]+ 2.0*in[IDX(i,j)]-in[IDX(i,j+1)])*dy2i;
                }
                // top boundary for middle ranks - j=0 
                j = 0;
                out[IDX(i,j)] = (-in[IDX(i-1, j)]+ 2.0*in[IDX(i,j)]-in[IDX(i+1, j)])*dx2i+ (-top[IDX(i,0)]+ 2.0*in[IDX(i,   j)]-in[IDX(i, j+1)])*dy2i;
                //bottom boundary - j=local_rows-1
                j=local_rows-1;
                out[IDX(i,j)] = (-in[IDX(i-1, j)]+ 2.0*in[IDX(i,j)]-in[IDX(i+1, j)])*dx2i+ (-in[IDX(i,j-1)]+ 2.0*in[IDX(i,   j)]-bot[IDX(i,0)])*dy2i;
            }
        }
    }
    else { // if rank=0 and size=0 then perfrom the original algorithm
        //int jm1 = 0, jp1 = 2; // original code
        for (int j = 1; j < local_rows - 1; ++j) {
            for (int i = 1; i < Nx - 1; ++i) {
                out[IDX(i,j)] = ( -     in[IDX(i-1, j)]
                                  + 2.0*in[IDX(i,   j)]
                                  -     in[IDX(i+1, j)])*dx2i
                              + ( -     in[IDX(i, j-1)]
                                  + 2.0*in[IDX(i,   j)]
                                  -     in[IDX(i, j+1)])*dy2i;
            }
            //jm1++;
            //jp1++;
        }
    }
}

void SolverCG::Send_and_Receive(double* localMatrix,double* sendMatrixT,double* sendMatrixB, double* recvMatrixT, double* recvMatrixB,int& rank,int& size,int& local_rows){
    // send and receive boundaries between sub-domains - rows
    if (size!=1) { // sending data if only one process
		 if (size>1){
		     if (rank>0){ //only sending to top rank
		         // create send matrix top
		         for (int i=0;i<Nx;++i){
		             sendMatrixT[i] = localMatrix[(0)*Nx+i];
		         }
		         MPI_Send(sendMatrixT,Nx,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);
		         MPI_Recv(recvMatrixT,Nx,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		     }
		     if (rank<size-1){ //only sending to the bottom
		         // create the send matrix bottom
		         for (int i=0;i<Nx;++i){
		             sendMatrixB[i] = localMatrix[(local_rows-1)*Nx+i];
		         }
		         MPI_Send(sendMatrixB,Nx,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD);
		         MPI_Recv(recvMatrixB,Nx,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		     }
		 }
	 }
	 else {
	 	cout <<"No Transmission ####ONLY ONE PROCESS####"<<endl;
 	}
}

void SolverCG::Precondition(double* in, double* out,int& rank,int& size,int& local_rows) {
    int i, j;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    double factor = 2.0*(dx2i + dy2i);
    for (i = 0; i < Nx ; ++i) { // do for all elements
        for (j = 0; j < local_rows ; ++j) {
            out[IDX(i,j)] = in[IDX(i,j)]/factor;
        }
    }
    // Apply top and bottom and top boundaries
    if (rank==0) {
        for (int i=0;i<Nx;++i){
            out[IDX(i, 0)] = in[IDX(i,0)]; // top
        }
        if (size==1){ // perform top boundary if only 1 process
            for (int i=0;i<Nx;++i) {
                out[IDX(i, local_rows-1)] = in[IDX(i, local_rows-1)]; // bottom
            }
        }
    }
    else if(rank==size-1){ // top boundary of last rank
        for (int i=0;i<Nx;++i) {
            out[IDX(i, local_rows-1)] = in[IDX(i, local_rows-1)]; // bottom
        }
    }
    
    // left and right most boundaries for all ranks
    for (j = 0; j < local_rows; ++j) {
        out[IDX(0, j)] = in[IDX(0, j)];
        out[IDX(Nx - 1, j)] = in[IDX(Nx - 1, j)];
    }
}

void SolverCG::ImposeBC(double* inout,int& rank,int& size,int& local_rows) {
        // Boundaries
    if (rank==0){// bottom boundary
        for (int i = 0; i < Nx; ++i) {
            inout[IDX(i, 0)] = 0.0;
        }
        if (size==1){ // top boundary too if only one process
            for (int i = 0; i < Nx; ++i) {
                inout[IDX(i, local_rows-1)] = 0.0;
            }
        }
    }
    else if (rank==size-1){ // only the top boundary of last rank
        for (int i = 0; i < Nx; ++i) {
            inout[IDX(i, local_rows-1)] = 0.0;
        }
    }
    // left and right boundaries for all processes.
    for (int j = 0; j < local_rows; ++j) {
        inout[IDX(0, j)] = 0.0;
        inout[IDX(Nx - 1, j)] = 0.0;
    }

}

void SolverCG::printSubArrayToFile(double* subArray, int rank,int k, int n,int local_rows, string name) {
    ostringstream filename;
    filename <<"/home/is420/lid-driven-cavity/output/"<< name << setfill('0') << setw(2) << rank << setw(2)<<k<< ".txt";
    ofstream file(filename.str());
    file << "The rank is: " << rank <<" k: "<<k<<endl;
    if (file.is_open()) {
        for (int j = 0; j < local_rows; ++j) {
            for (int i = 0; i < n; ++i) {
                    file << subArray[j*n+i] << "\t";
            }
            file << endl;
        }
        file.close();
        cout << "Process " << rank << " wrote subarray to " << filename.str() << endl;
    } else {
        cerr << "Process " << rank << " failed to open file " << filename.str() << endl;
    }
}



