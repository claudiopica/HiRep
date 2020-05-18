#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "type.h"
#include "utils.h"

void usage(char *code) // revamp
{
  fprintf(stderr, "Usage: %s -i <def correlator filename>  -f <1pt file name> [-J <jackknife bin size>] [-D <directories file>] [-r <op string>] [-G] [-x <inversion point>] [-d <diagonalization point>] [-e <number of output states>] [-T <lt>]\n", code);
  fprintf(stderr, " -i <def correlator filename> : name of the file containg the 2pt defining string\n   -f <1pt file name> : name of the 1pt file\n   -D <directories file> :  File containg the list of the directories containg the data files\n   -J <jacknife bin size> : Size of the Jacknife bin\n   -r <op string> : Defines the operators to be removed from the basis\n   \tsyntax n1,n2,n3-n4,with ni <= nj if i<j \n   -G : If enabled uses the automatic procedure to create the maximal\n   \toperatorial basis starting from the <op string> if specified\n   -x <inversion point> : defines the distance at which is evaluated the inverse correlator\n\t default is the minimal distance\n   -d <diagonalization point> : defines the distance at which is evaluated the diagonalized correlator\n   \tdefault is min dist + 1   \n   -e <number of states> : defines the number of eigenvalues that will be analysed in the code\n\t default is 6\n   -T <Lt> : Temporal extension of the lattice, if specified will fit with PBC otherwise will assume OBC\n");
};

static void select_opplist(par *apar, char *input);

void read_arg(par *apar, int argc, char **argv)
{
  
  int c;

  while ((c = getopt(argc, argv, "i:f:D:J:r:Gx:d:e:vT:")) != -1)
    switch (c)
    {
    case 'i':
      apar->iflag = 1;
      sprintf(apar->cor_def_filename, "%s", optarg);
      break;
    case 'f':
      apar->fflag = 1;
      sprintf(apar->vev_name, "%s", optarg);
      break;
    case 'v':
      apar->vflag = 1;
      break;
    case 'D':
      apar->Dflag = 1;
      sprintf(apar->list_dir, "%s", optarg);
      break;
    case 'J':
      apar->Jflag = 1;
      apar->binwidth = atoi(optarg);
      break;
    case 'r':
      apar->rflag = 1;
      sprintf(apar->opstring, "%s", optarg);
      break;
    case 'G':
      apar->Gflag = 1;
      break;
    case 'x':
      apar->xflag = 1;
      apar->p_inv = atoi(optarg);
      break;
    case 'd':
      apar->dflag = 1;
      apar->p_diag = atoi(optarg);
      break;
    case 'e':
      apar->eflag = 1;
      apar->n_states = atoi(optarg);
      break;
    case 'T':
      apar->Tflag = 1;
      apar->lt = atoi(optarg);
      break;
    case '?':
      fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      usage(argv[0]);
      exit(1);
    }

  if (apar->iflag == 0)
  {
    fprintf(stderr, "Missing -i argument.\n");
    usage(argv[0]);
    exit(1);
  }

  if (apar->fflag == 0)
  {
    fprintf(stderr, "Missing -f argument.\n");
    usage(argv[0]);
    exit(1);
  }

  // check that the directories in the input file contain necessary files
  // and set the analysis parameters in apar : like number of bins, nt, Number of operators etc...
  check_list_dir_and_data_file(apar);

  // check that the file containing the 2pt defs is fine
  // and set the analysis 2pt correlator map in apar (associate to each dt a list of pairs)
  read_2pt_def(apar);

  // set the map in apar  between the T position and memory index of the 1pt
  read_1pt_map(apar);

  if (apar->rflag == 1)
  {
    std::cout << "[INFO][read_arg] Dropping the operators identified by the string: " << apar->opstring << std::endl;
    select_opplist(apar, apar->opstring);
  }

  if (apar->n_states > apar->numop)
  {
    apar->n_states = apar->numop;
    std::cout << "[INFO][read_arg] N of excited states changed to the number of operators:" << apar->numop << std::endl;
  }

  if (apar->n_states < 1)
  {
    std::cerr << "Number of states to be analized (-e) must be greater than zero." << std::endl;
    usage(argv[0]);
    exit(1);
  }

  if (apar->p_inv == apar->p_diag || apar->p_diag < 0 || apar->p_diag >= apar->ndt || apar->p_inv < 0 || apar->p_inv >= apar->ndt)
  {
    std::cerr << "Inversion(t=" << apar->p_inv << ") or Diagonlization(t=" << apar->p_diag << ") point not correct." << std::endl;
    std::cerr << "It should be true that 0 <= " << apar->p_inv << " < " << apar->p_diag << " < " << apar->ndt << std::endl;
    usage(argv[0]);
    exit(1);
  }
}

static void select_opplist(par *apar, char *input)
{
  std::string str2;
  size_t found = 0, n1_start = 0, n1_end = 0;
  int op_start, op_end;
  bool control = true;

  std::string str(input);

  while (control)
  {
    found = str.find(",", n1_start);
    n1_end = found;
    if (found == std::string::npos)
    {
      control = false;
      n1_end = str.length();
    }
    str2 = str.substr(n1_start, n1_end - n1_start);

    if (std::sscanf(str2.c_str(), "%d-%d", &op_start, &op_end) == 2)
    {
      if (op_end < op_start)
      {
        std::cerr << "[ERROR][select_opplist]Wrong specification of the operator in the list of removed operators" << std::endl;
        exit(1);
      }
      for (int i = op_start; i <= op_end; i++)
        deactivate_op(apar, i);
    }
    else if (std::sscanf(str2.c_str(), "%d", &op_start) == 1)
      deactivate_op(apar, op_start);
    else
    {
      std::cerr << "[ERROR][select_opplist]Wrong specification of the operator in the list of removed operators" << std::endl;
      exit(1);
    }
    n1_start = n1_end + 1;
  }
}
