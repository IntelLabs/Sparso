CXX=icpc

CFLAGS = -I. -fopenmp --std=c++11
CFLAGS += -I/opt/apps/intel13/boost/1.51.0/include # for stampede

ifeq (yes, $(DBG))
  CFLAGS += -g -O0
else
  CFLAGS += -O3 -DNDEBUG -xHost
endif
ifeq (yes, $(BOOST))
  CFLAGS += -DUSE_BOOST
endif

test/CSR_test: test/CSR_test.c CSR.cpp CSR_Interface.cpp Vector.cpp mm_io.cpp reordering/RCM.cpp
	$(CXX) $(CFLAGS) -o $@ $^

# 0 based CSR library
lib: CSR.cpp CSR_Interface.cpp Vector.cpp reordering/RCM.cpp mm_io.cpp 
	$(CXX) $(CFLAGS) -shared -Wl,-soname,libcsr.so  -o libcsr.so -fopenmp -g -I. -fPIC $^
		
all:	test/CSR_test lib
libs:	lib

clean:
	rm -f test/CSR_test libcsr.so

.PHONY: clean
