#ifndef SPARSE_ACCELERATOR_H
#define SPARSE_ACCELERATOR_H
 
extern "C" void* create_matrix_knob(boolean constant_valued, SPARSE_MAT_STORAGE sparse_format, ...);
extern "C" void* create_function_knob(SPARSE_OPERATOR operator, ...);

#endif
