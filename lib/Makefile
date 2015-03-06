CXX=icpc

CFLAGS = -I. -openmp
ifeq (yes, $(DBG))
  CFLAGS += -g -O0 -DNDEBUG
else
  CFLAGS += -O3 -DNDEBUG
endif

test/CSR_test: test/CSR_test.c CSR.cpp CSR_Interface.cpp Vector.cpp mm_io.cpp
	$(CXX) $(CFLAGS) -o $@ $^
