CXX=icpc

SPMP_DIR=$(HOME)/SpMP

CFLAGS = -I. -I$(SPMP_DIR)/.. -fopenmp --std=c++11 -Wall -I$(MKLROOT)/include

ifeq (yes, $(DBG))
  CFLAGS += -g -O0
else
  CFLAGS += -O3 -DNDEBUG -xHost
endif
ifeq (yes, $(SEP))
  SEP_BASE_DIR ?= /opt/intel/sep # typically set by sep_vars.sh
  CFLAGS += -DSEP -I$(SEP_BASE_DIR)/../include $(SEP_BASE_DIR)/libprog_api.so
endif

CPP_SRCS = $(wildcard *.cpp)
C_SRCS = $(wildcard *.c)
OBJS = $(CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)
LIBS = $(SPMP_DIR)/libspmp.a
LINKFLAGS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_rt.so -Wl,--end-group

lib: $(OBJS) $(LIBS)
	$(CXX) $(CFLAGS) -shared -Wl,-soname,libcsr.so -o libcsr.so -fopenmp -g -I. -fPIC $^ $(LINKFLAGS)

test/CSR_test: CSR_Interface.cpp Vector.cpp mm_io.cpp SpMV.cpp PageRank.cpp test/CSR_test.c $(SPMP_DIR)/libspmp.a
	$(CXX) -o $@ $^ $(CFLAGS)

test/sparse_orig: CSR_Interface.cpp Vector.cpp mm_io.cpp SpMV.cpp PageRank.cpp test/sparse_orig.cc $(SPMP_DIR)/libspmp.a
	$(CXX) -o $@ $^ $(CFLAGS)

$(SPMP_DIR)/libspmp.a:
	cd $(SPMP_DIR); $(MAKE)
		
all:	test/CSR_test lib
libs:	lib

# By some reason, -O3 doesn't work for SpGEMM (FIXME)
SpGEMM.o: SpGEMM.cpp
	$(CXX) $(CFLAGS) -O0 -fPIC -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -fPIC -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -fPIC -c $< -o $@

clean:
	rm -f test/CSR_test libcsr.so *.o

.PHONY: clean $(SPMP_DIR)/libspmp.a
