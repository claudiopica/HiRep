#include "suN_types.h"
#include "linear_algebra.h"
#include "dirac.h"
#include "helper_functions.h"
#include "memory.h"
#include "hr_complex.h"
#include "error.h"
#include "logger.h"
#include "geometry.h"
#include "setup.h"
#include "inverters.h"

bool test_hermiticity(spinor_operator, spinor_operator, char*);
bool is_hermitian_on_GPU(spinor_field*, spinor_field*, 
		         spinor_field*, spinor_field*, spinor_operator);
bool is_hermitian_on_CPU(spinor_field*, spinor_field*, 
		         spinor_field*, spinor_field*, spinor_operator);
bool result_spinor_fields_not_identically_zero_gpu(spinor_field*, spinor_field*);
bool result_spinor_fields_not_identically_zero_cpu(spinor_field*, spinor_field*);
bool gpu_and_cpu_copies_identical(spinor_field*, spinor_field*);
void free_spinors(spinor_field**, spinor_field**, spinor_field**, spinor_field**);
void Q_operator(spinor_field*, spinor_field*);
void Q_operator_cpu(spinor_field*, spinor_field*);
void I_operator(spinor_field*, spinor_field*);
void I_operator_cpu(spinor_field*, spinor_field*);
