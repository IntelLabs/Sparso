CXX=icpc

test/CSR_test: test/CSR_test.c CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -fopenmp -g -I. -o $@ $^
