SPMP_DIR=$(HOME)/SpMP
CSR_LIB_DIR=$(HOME)/SparseAccelerator/lib

CXXFLAGS = -I$(SPMP_DIR)/.. -I$(CSR_LIB_DIR) $(SPMP_DIR)/libspmp.a $(CSR_LIB_DIR)/libcsr.so -mkl
CXX = icpc

ifeq (yes, $(DBG))
  CXXFLAGS += -O0 -g
else
	CXXFLAGS += -O3 -DNDEBUG -xHost
endif

dss_test: dss_test.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -f dss_test
