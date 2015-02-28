CXX=g++

test/CSR_test: test/CSR_test.c CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -fopenmp -g -I. -o $@ $^
	
lib: CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr.so  -o libcsr.so -fopenmp -g -I. -fPIC $^
	
	
