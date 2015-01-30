// TODO: we need to ask Jongsoo to expose the header of his library sources 
#include "OptimizeProblem.hpp"
#include "SparseAccelerator.h"
#include <vector.h>

typedef enum {
	COO,
	CSR,
	CSC,
	DIA,
	SKY,
	BSR
} SPARSE_MAT_STORAGE;

typedef enum {
	BACKSLASH,
	SPMV
} SPARSE_OPERATOR;

// A sparse matrix needs a few data structures to represent. 
// For example, a CSC-formatted matrix needs m(#rows), n(#columns),
// colptr, rowval, and nzval. See sparsematrix.jl
#define MAX_REPRESENTATION_STRUCTS 5

class MKNOB {
public:
	unsigned           version;
	boolean            constant_valued;
	SPARSE_MAT_STORAGE sparse_format;
	void*              matrix_representation[MAX_REPRESENTATION_STRUCTS];
}

class FKNOB {
public:
	std::vector<MKNOB*>   inputs;

protected: 
	std::vector<unsigned> versions; // versions of the inputs
}

class Solver_FKNOB : public FKNOB{
public:
	Solver_FKNOB() { opt_data = NULL; }

Private: 
	OptimizationData*     opt_data; // including dep graph, schedule, etc.
}

// In calling this function, end the parameter list with a NULL
void* create_matrix_knob(boolean constant_valued, SPARSE_MAT_STORAGE sparse_format, ...)
{
	MKNOB * mknob = new MKNOB();
	mknob->version = 0;
	mknob->constant_valued = constant_valued;
	mknob->sparse_format = sparse_format;
	
	va_list args;
	va_start(args, sparse_format);
	void * arg;
	int i = 0;
	for (arg = va_arg(args, void*); arg != NULL; arg = va_arg(args, void*)) {
		mknob->matrix_representation[i++] = arg; 
	}
	va_end(args);

	return mknob;
}

void associate_func_matrix_knobs(FKNOB* fknob, va_list mknobs)
{
	assert(fknob != NULL);
			
	MKNOB *mknob;
	for (mknob = va_arg(mknobs, MKNOB*); mknob != NULL; mknob = va_arg(mknobs, MKNOB*)) {
		fknob->inputs.push_back(mknob);
		fknob->versions.push_back(mknob->version);	
	}
}

// Pass in 1 or more matrix knob, end with NULL
void* create_function_knob(SPARSE_OPERATOR operator, ...)
{
	FKNOB *fknob = NULL; 
	va_list args;
	va_start(args, operator);
	
	switch (operator) {
	case BACKSLASH: {
		fknob = new Solver_FKNOB();
		associate_func_matrix_knobs(fknob, args);
		break;
	}
	case SPMV:
	default:
		assert(false); // not implemented yet.
	}
	va_end(args);
	
	return fknob;
}
