CXX=icpc

SPMP_DIR=$(HOME)/SpMP

CFLAGS = -I. -I$(SPMP_DIR)/.. -fopenmp --std=c++11

ifeq (yes, $(DBG))
  CFLAGS += -g -O0
else
  CFLAGS += -O3 -DNDEBUG -xHost
endif
ifeq (yes, $(SEP))
  SEP_BASE_DIR ?= /opt/intel/sep # typically set by sep_vars.sh
  CFLAGS += -DSEP -I$(SEP_BASE_DIR)/../include $(SEP_BASE_DIR)/libprog_api.so
endif

test/CSR_test: CSR_Interface.cpp Vector.cpp mm_io.cpp SpMV.cpp PageRank.cpp test/CSR_test.c $(SPMP_DIR)/libspmp.a
	$(CXX) -o $@ $^ $(CFLAGS)

test/sparse_orig: CSR_Interface.cpp Vector.cpp mm_io.cpp SpMV.cpp PageRank.cpp test/sparse_orig.cc $(SPMP_DIR)/libspmp.a
	$(CXX) -o $@ $^ $(CFLAGS)

# 0 based CSR library
lib: CSR_Interface.cpp Vector.cpp mm_io.cpp SpMV.cpp PageRank.cpp $(SPMP_DIR)/libspmp.a
	$(CXX) $(CFLAGS) -shared -Wl,-soname,libcsr.so  -o libcsr.so -fopenmp -g -I. -fPIC $^

$(SPMP_DIR)/libspmp.a:
	cd $(SPMP_DIR); $(MAKE)
		
all:	test/CSR_test lib
libs:	lib

clean:
	rm -f test/CSR_test libcsr.so

.PHONY: clean $(SPMP_DIR)/libspmp.a
