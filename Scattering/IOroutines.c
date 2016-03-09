
// Followed by the list of IO functions from meson_scattering.c
/* Mesons parameters */
typedef struct _input_scatt {
	char mstring[256];
	double precision;
	int nhits;
	int tsrc;
	char outdir[256], p1[16], p2[16];

	/* for the reading function */
	input_record_t read[8];

} input_scatt;

#define init_input_scatt(varname) \
{ \
	.read={\
		{"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring},\
		{"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
		{"number of inversions per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},\
		{"Source time:", "mes:tsrc = %d", INT_T, &(varname).tsrc},\
		{"Output directory:", "mes:outdir = %s", STRING_T, &(varname).outdir},\
		{"Momentum 1:", "mes:p1 = %s", STRING_T, &(varname).p1},\
		{"Momentum 2:", "mes:p2 = %s", STRING_T, &(varname).p2},\
		{NULL, NULL, INT_T, NULL}\
	}\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char prop_filename[256]="";
char source_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "meson_scattering.out";
int Nsource;
double M;

enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_scatt mes_var = init_input_scatt(mes_var);

typedef struct {
	char string[256];
	int t, x, y, z;
	int nc, nf;
	double b, m;
	int n;
	int type;
} filename_t;


int parse_cnfg_filename(char* filename, filename_t* fn) {
	int hm;
	char *tmp = NULL;
	char *basename;

	basename = filename;
	while ((tmp = strchr(basename, '/')) != NULL) {
		basename = tmp+1;
	}            

#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif
	hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%dr" repr_name "%*[Nn]f%db%lfm%lfn%d",
			&(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&(fn->m),&(fn->n));
	if(hm==9) {
		fn->m=-fn->m; /* invert sign of mass */
		fn->type=DYNAMICAL_CNFG;
		return DYNAMICAL_CNFG;
	}
#undef repr_name

	double kappa;
	hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%d%*[Nn]f%db%lfk%lfn%d",
			&(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&kappa,&(fn->n));
	if(hm==9) {
		fn->m = .5/kappa-4.;
		fn->type=DYNAMICAL_CNFG;
		return DYNAMICAL_CNFG;
	}

	hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%db%lfn%d",
			&(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
	if(hm==7) {
		fn->type=QUENCHED_CNFG;
		return QUENCHED_CNFG;
	}

	hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d",
			&(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
	if(hm==7) {
		fn->type=QUENCHED_CNFG;
		return QUENCHED_CNFG;
	}

	fn->type=UNKNOWN_CNFG;
	return UNKNOWN_CNFG;
}


void read_cmdline(int argc, char* argv[]) {
	int i, ai=0, ao=0, ac=0, al=0, am=0,ap=0,as=0;
	FILE *list=NULL;

	for (i=1;i<argc;i++) {
		if (strcmp(argv[i],"-i")==0) ai=i+1;
		else if (strcmp(argv[i],"-o")==0) ao=i+1;
		else if (strcmp(argv[i],"-c")==0) ac=i+1;
		else if (strcmp(argv[i],"-p")==0) ap=i+1;
		else if (strcmp(argv[i],"-l")==0) al=i+1;
		else if (strcmp(argv[i],"-s")==0) as=i+1;
		else if (strcmp(argv[i],"-m")==0) am=i;
	}


	if (am != 0) {
		print_compiling_info();
		exit(0);
	}


	if (ao!=0) strcpy(output_filename,argv[ao]);
	if (ai!=0) strcpy(input_filename,argv[ai]);

	error ((ap != 0 && as ==0),1,"parse_cmdline [discc]",
			"Syntax: disc { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m] -p <propagator name> -s <source_name> ");

	error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [meson_scattering.c]",
			"Syntax: disc { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m] -p <propagator name> -s <source_name> ");

	if(ap != 0) {
		strcpy(prop_filename,argv[ap]);
	} 
	if(as != 0) {
		strcpy(source_filename,argv[as]);
	} 

	if(ac != 0) {
		strcpy(cnfg_filename,argv[ac]);
		strcpy(list_filename,"");
	} else if(al != 0) {
		strcpy(list_filename,argv[al]);
		error((list=fopen(list_filename,"r"))==NULL,1,"parse_cmdline [meson_scattering.c]" ,
				"Failed to open list file\n");
		error(fscanf(list,"%s",cnfg_filename)==0,1,"parse_cmdline [meson_scattering.c]" ,
				"Empty list file\n");
		fclose(list);
	}

}
