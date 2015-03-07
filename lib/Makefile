CXX=g++

CFLAGS = -I. -fopenmp
ifeq (yes, $(DBG))
  CFLAGS += -g -O0 -DNDEBUG
else
  CFLAGS += -O3 -DNDEBUG
endif

test/CSR_test: test/CSR_test.c CSR.cpp CSR_Interface.cpp Vector.cpp mm_io.cpp
	$(CXX) $(CFLAGS) -o $@ $^

test/CSR_small_test: test/CSR_small_test.c CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) $(CFLAGS) -o $@ $^
	
test/CSR1_test: test/CSR_1based_indexing_test.c CSR1.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) $(CFLAGS) -o $@ $^
	
# 0 based CSR library
lib0: CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr0.so  -o libcsr0.so -fopenmp -g -I. -fPIC $^
	
# 1 based CSR library
lib1: CSR_1based_indexing.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr1.so  -o libcsr1.so -fopenmp -g -I. -fPIC $^
	
