CXX = mpicxx
CFLAGS = -std=c++11 -Wall -O2
BOOST_LIBS = -lboost_unit_test_framework

SolverCG.o: SolverCG.cpp SolverCG.h
	$(CXX) $(CFLAGS) -o SolverCG.o -c SolverCG.cpp -lblas

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp LidDrivenCavity.h
	$(CXX) $(CFLAGS) -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp -lboost_program_options

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h SolverCG.h
	$(CXX) $(CFLAGS) -o LidDrivenCavity.o -c LidDrivenCavity.cpp -lblas

solver: LidDrivenCavity.o LidDrivenCavitySolver.o SolverCG.o
	$(CXX) -o solver LidDrivenCavity.o LidDrivenCavitySolver.o SolverCG.o -lblas -lboost_program_options

UnitTest.o: UnitTest.cpp
	$(CXX) $(CFLAGS) -o UnitTest.o -c UnitTest.cpp $(BOOST_LIBS) -lblas -lboost_program_options

test: UnitTest.o SolverCG.o
	$(CXX) -o test UnitTest.o SolverCG.o $(BOOST_LIBS) -lblas -lboost_program_options


clean:
	rm -f solver
	rm -f *.o
	rm -f test
	
doc:
	doxygen Doxyfile
	




