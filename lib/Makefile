CXX=g++

test/CSR0_test: test/CSR0_test.c CSR0.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -fopenmp -g -I. -o $@ $^

test/CSR1_test: test/CSR1_test.c CSR1.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -fopenmp -g -I. -o $@ $^
	
# 0 based CSR library
lib0: CSR0.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr0.so  -o libcsr0.so -fopenmp -g -I. -fPIC $^
	
# 1 based CSR library
lib1: CSR1.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr1.so  -o libcsr1.so -fopenmp -g -I. -fPIC $^
	
