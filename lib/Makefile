CXX=g++

CFLAGS = -I. -fopenmp
ifeq (yes, $(DBG))
  CFLAGS += -g -O0 -DNDEBUG
else
  CFLAGS += -O3 -DNDEBUG
endif

test/CSR0_test: test/CSR0_test.c CSR0.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) $(CFLAGS) -o $@ $^
	
test/CSR1_test: test/CSR1_test.c CSR1.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) $(CFLAGS) -o $@ $^
	
# 0 based CSR library
lib0: CSR0.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr0.so  -o libcsr0.so -fopenmp -g -I. -fPIC $^
	
# 1 based CSR library
lib1: CSR1.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr1.so  -o libcsr1.so -fopenmp -g -I. -fPIC $^
	
