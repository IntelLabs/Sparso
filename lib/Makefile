CXX=g++

CFLAGS = -I. -fopenmp
ifeq (yes, $(DBG))
  CFLAGS += -g -O0 -DNDEBUG
else
  CFLAGS += -O3 -DNDEBUG
endif

test/CSR_test: test/CSR_test.c CSR.cpp CSR_Interface.cpp Vector.cpp mm_io.cpp
	$(CXX) $(CFLAGS) -o $@ $^
	
lib: CSR.cpp CSR_Interface.cpp Vector.cpp
	$(CXX) -shared -Wl,-soname,libcsr.so  -o libcsr.so -fopenmp -g -I. -fPIC $^
