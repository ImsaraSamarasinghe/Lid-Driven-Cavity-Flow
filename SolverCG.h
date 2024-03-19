#pragma once

class SolverCG
{
public:
    SolverCG(int pNx, int pNy, double pdx, double pdy);
    ~SolverCG();

    void Solve(double* b, double* x, int argc, char** argv);

private:
    double dx;
    double dy;
    int Nx;
    int Ny;
    double* r;
    double* p;
    double* z;
    double* t;

    void ApplyOperator(double* in, double* out,double* top,double* bot,int& rank, int& size, int& local_rows);
    void Precondition(double* in, double* out,int& rank,int& size,int& local_rows);
    void ImposeBC(double* inout,int& rank,int& size,int& local_rows);
    void Send_and_Receive(double* localMatrix,double* sendMatrixT,double* sendMatrixB, double* recvMatrixT, double* recvMatrixB,int& rank,int& size,int& local_rows);
    void printSubArrayToFile(double* subArray, int rank,int k,int n,int local_rows, string name);
};

