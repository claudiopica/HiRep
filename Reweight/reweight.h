#ifndef REWEIGHT_H
#define REWERIGHT_H

typedef struct {
	int hits;
	int steps;
	double precision;
	double old_mass;
	double new_mass;
	double old_theta_t;
	double old_theta_x;
	double old_theta_y;
	double old_theta_z;
	double new_theta_t;
	double new_theta_x;
	double new_theta_y;
	double new_theta_z;
	input_record_t read[14];
} input_reweight;

#define init_input_reweight(varname) \
{ \
	.read = {\
		{"number of gaussian vectors", "hits = %d", INT_T, &(varname).hits},\
		{"number of determinant substeps", "steps = %d", INT_T, &(varname).steps},\
		{"inverter precision", "precision = %lf", DOUBLE_T, &(varname).precision},\
		{"old bare mass", "old_mass = %lf", DOUBLE_T, &(varname).old_mass},\
		{"new bare mass", "new_mass = %lf", DOUBLE_T, &(varname).new_mass},\
		{"old twisting angle t-direction", "old_theta_t = %lf", DOUBLE_T, &(varname).old_theta_t},\
		{"old twisting angle x-direction", "old_theta_x = %lf", DOUBLE_T, &(varname).old_theta_x},\
		{"old twisting angle y-direction", "old_theta_y = %lf", DOUBLE_T, &(varname).old_theta_y},\
		{"old twisting angle z-direction", "old_theta_z = %lf", DOUBLE_T, &(varname).old_theta_z},\
		{"new twisting angle t-direction", "new_theta_t = %lf", DOUBLE_T, &(varname).new_theta_t},\
		{"new twisting angle x-direction", "new_theta_x = %lf", DOUBLE_T, &(varname).new_theta_x},\
		{"new twisting angle y-direction", "new_theta_y = %lf", DOUBLE_T, &(varname).new_theta_y},\
		{"new twisting angle z-direction", "new_theta_z = %lf", DOUBLE_T, &(varname).new_theta_z},\
		{NULL, NULL, INT_T, NULL}\
	}\
}

input_reweight rw_var = init_input_reweight(rw_var);

#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif

#endif