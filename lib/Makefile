# Copyright (c) 2015, Intel Corporation
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Intel Corporation nor the names of its contributors
#       may be used to endorse or promote products derived from this software
#       without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
CXX=icpc

SPMP_DIR=../deps/SpMP

CFLAGS = -I. -I$(SPMP_DIR)/.. -qopenmp --std=c++11 -Wall -I$(MKLROOT)/include

ifeq (yes, $(DBG))
  CFLAGS += -g -O0
else
  CFLAGS += -O3 -DNDEBUG -xHost #-DMKL_DIRECT_CALL
endif
ifeq (yes, $(SEP))
  SEP_BASE_DIR ?= /opt/intel/sep # typically set by sep_vars.sh
  CFLAGS += -DSEP -I$(SEP_BASE_DIR)/../include $(SEP_BASE_DIR)/libprog_api.so
endif

CPP_SRCS = $(wildcard *.cpp)
C_SRCS = $(wildcard *.c)
OBJS = $(CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)
LIBS = $(SPMP_DIR)/libspmp.a
LINKFLAGS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_rt.so -Wl,--end-group # single dynamic library
#LINKFLAGS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group # static
#LINKFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread # dynamic

lib: $(OBJS) $(LIBS)
	$(CXX) $(CFLAGS) -shared -Wl,-soname,libcsr.so -o libcsr.so -fopenmp -g -I. -fPIC $^ $(LINKFLAGS)

test/CSR_test: CSR_Interface.cpp Vector.cpp mm_io.cpp SpMV.cpp PageRank.cpp test/CSR_test.c $(SPMP_DIR)/libspmp.a
	$(CXX) -o $@ $^ $(CFLAGS)

test/sparse_orig: CSR_Interface.cpp Vector.cpp mm_io.cpp SpMV.cpp PageRank.cpp test/sparse_orig.cc $(SPMP_DIR)/libspmp.a
	$(CXX) -o $@ $^ $(CFLAGS)

test/BFSBipartite_test: $(OBJS) test/BFSBipartite_test.cpp $(SPMP_DIR)/libspmp.a
	$(CXX) -o $@ $^ $(CFLAGS) -mkl

$(SPMP_DIR)/libspmp.a:
	cd $(SPMP_DIR); $(MAKE)
		
all:	test/CSR_test lib
libs:	lib

%.o: %.cpp
	$(CXX) $(CFLAGS) -fPIC -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -fPIC -c $< -o $@

clean:
	rm -f test/CSR_test libcsr.so *.o
	cd $(SPMP_DIR); make clean

.PHONY: clean $(SPMP_DIR)/libspmp.a
