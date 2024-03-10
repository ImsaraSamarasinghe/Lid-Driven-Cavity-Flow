CXX = g++
CFLAGS = -std=c++11 -Wall -O2
default:solver

SolverCG.o: SolverCG.cpp SolverCG.h
	$(CXX) $(CFLAGS) -o SolverCG.o -c SolverCG.cpp -lblas

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp LidDrivenCavity.h
	$(CXX) $(CFLAGS) -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp -lboost_program_options

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h SolverCG.h
	$(CXX) $(CFLAGS) -o LidDrivenCavity.o -c LidDrivenCavity.cpp -lblas

solver: LidDrivenCavity.o LidDrivenCavitySolver.o SolverCG.o
	$(CXX) -o solver LidDrivenCavity.o LidDrivenCavitySolver.o SolverCG.o -lblas -lboost_program_options




