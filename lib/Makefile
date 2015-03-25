CXX=icpc

CFLAGS = -I. -fopenmp --std=c++11
ifeq (yes, $(DBG))
  CFLAGS += -g -O0
else
  CFLAGS += -O3 -DNDEBUG
endif

test/CSR_test: test/CSR_test.c CSR.cpp CSR_Interface.cpp Vector.cpp mm_io.cpp reordering/RCM.cpp
	$(CXX) $(CFLAGS) -o $@ $^

test/CSR_small_test: test/CSR_small_test.c CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) $(CFLAGS) -o $@ $^
	
test/CSR1_test: test/CSR_1based_indexing_test.c CSR_1based_indexing.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) $(CFLAGS) -o $@ $^
	
# 0 based CSR library
lib: CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr.so  -o libcsr.so -fopenmp -g -I. -fPIC $^
	
# 1 based CSR library
lib1: CSR_1based_indexing.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr1.so  -o libcsr1.so -fopenmp -g -I. -fPIC $^
	
