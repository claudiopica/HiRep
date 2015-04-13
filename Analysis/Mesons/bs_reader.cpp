#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <cstdlib>
#include <cstring>
#include "bs_type.h"
#include "bs_ctrl.h"


using namespace std;

extern int Lt;
extern int Ls;
extern double svol;
extern int blsize;
extern int nsamples1;
extern int nsamples2;
extern string channel;
extern string inputfilename;
extern string cutfilename;
extern int effmass_method;

void read_cut(const char *cutFile)
{
    ifstream CutFile(cutFile);

    while(1){
      stringstream line;
      stringbuf sb;
      CutFile.get(sb);
      line << sb.str();
      if(line.good()){
	string Channel;
	int left, right;
	line>>Channel;
        if(Channel.compare("")!=0)
	  if(Channel[0]!='#'){        
	    pair<map<string,int>::iterator,bool> ret; 
	    if(!(line >> left)) {
	      cerr<<"[read_cut]:Missing entry for "<<Channel<<" for left cut"<<endl;
	      exit(1);
	    }
	    ret=eval_ctrl::left_cut.insert(pair<string,int>(Channel,left));
	    if(ret.second==false) {
	      cerr<<"[read_cut]:Duplicate entry for "<<Channel<<" for left cut"<<endl;
	      exit(1);
	    }
	    if(!(line >> right)) {
	      cerr<<"[read_cut]:Missing entry for "<<Channel<<" for right cut"<<endl;
	      exit(1);
	    }
	    ret=eval_ctrl::right_cut.insert(pair<string,int>(Channel,right));
	    if(ret.second==false) {
	      cerr<<"[read_cut]:Duplicate entry for "<<Channel<<" for right cut"<<endl;
	      exit(1);
	    }
	    if(left<0) {
	      cerr<<"[read_cut]:Negative value for "<<Channel<<" for left cut"<<endl;
	      exit(1);
	    }
	    if(right<left) {
	      cerr<<"[read_cut]:Right cut value less than left cut value for "<<Channel<<""<<endl;
	      exit(1);
	    }
	  }
	}
	if(CutFile.eof()) break;
	CutFile.clear();
	CutFile.seekg(1,ios_base::cur);
      }
    CutFile.close();
}



static void sym(int n, double *data) {
  for(int i=1; i<(n+1)/2; ++i) {
    data[i]+=data[n-i];
    data[i]*=0.5;
    data[n-i]=data[i];
  }
}

static void antisym(int n, double *data) {
  for(int i=1; i<(n+1)/2; ++i) {
    data[i]-=data[n-i];
    data[i]*=0.5;
    data[n-i]=-data[i];
  }
  if(!(n&1)) data[n/2]=0.; //put to zero middle point
}


void read_input(const char * inputFile){  
  ifstream mesinput(inputFile);
  double corr[Lt];
  
  while(1) {
    stringstream line;
    stringbuf sb;
    mesinput.get(sb);
    if(!mesinput.good()) break;
    mesinput.seekg(1,ios_base::cur);
    line << sb.str();
    if(line.good()){
      string id;
      string Word;
      line>>Word; /* trajectory index */
      line>>id; /* channel id */
      int i=0;
      while(1) { 
    	  double d;
        line>>d;
        if(!line) break;
      	if(i<Lt) {
      	  corr[i]=d;
      	  i++;
      	} else {
      	  cerr << "[read_input]: Correlator is longer than Lt"<<endl;
      	  exit(1);
      	}
      }
      if(i<Lt) {
      	cerr << "[read_input]: Correlator is shorter than Lt"<<endl;
      	exit(1);
      }

      if(eval_ctrl::channel_map.find(id)!=eval_ctrl::channel_map.end()) {
      	//symmetrize or antisymmetrize data
      	if(id=="g5_g0g5_re") { antisym(Lt,corr);}
      	else { sym(Lt, corr); }
      	for(int i=0;i<eval_ctrl::channel_map[id]->length;++i) {
      	  eval_ctrl::channel_map[id]->d[i].push_back(corr[i]);
      	}
      }
    }
  }
  mesinput.close();

}


int read_cmdline(int argc, char*argv[]){
  if(argc!=8 && argc!=10 && argc!=3){
    cerr<<"[read_cmdline]: Missing Parameter\n\tUsage "<<argv[0]<<" <channel> <method> <inputfile> <cutfile> <Lt> <Ls> <blocksize> [<bootstrap1> <bootstrap2>]"<<endl;
    cerr<<"[read_cmdline]: Usage "<<argv[0]<<" info <channel>"<<endl;
    exit(1);
  }
  
  if(argc==3){
    if( strcmp(argv[1],"info") == 0 ){
      channel=argv[2];
      Lt=2;
      return 1;
    }
    else {
      cerr<<"[read_cmdline]: Missing Parameter\n\tUsage "<<argv[0]<<" <channel> <method> <inputfile> <cutfile> <Lt> <Ls> <blocksize> [<bootstrap1> <bootstrap2>]"<<endl;
      cerr<<"[read_cmdline]: Usage "<<argv[0]<<" info <channel>"<<endl;
      exit(1);
    }
  }
  
  channel=argv[1];
  
  effmass_method=atoi(argv[2]);
  if(effmass_method<0 || effmass_method>2) {
    cerr<<"[read_cmdline]: method can be only 0,1,2" <<endl;
    exit(1);
  }

  inputfilename=argv[3];
  ifstream input(inputfilename.c_str());
  if(!input){
    cerr<<"[read_cmdline]: Missing File "<<inputfilename<<endl;
    exit(1);
  }
  input.close();

  cutfilename=argv[4];
  ifstream cut(cutfilename.c_str());
  if(!cut){
    cerr<<"[read_cmdline]: Missing File "<<cutfilename<<endl;
    exit(1);
  }
  cut.close();

  Lt=atoi(argv[5]);
  if(Lt<=0) {
    cerr<<"[read_cmdline]: Wrong value for Lt"<<endl;
    exit(1);
  }
  Ls=atoi(argv[6]);
  if(Ls<=0) {
    cerr<<"[read_cmdline]: Wrong value for Ls"<<endl;
    exit(1);
  }
  
  svol=sqrt(double(Ls*Ls*Ls));

  blsize=atoi(argv[7]);
  if(blsize<=0) {
    cerr<<"[read_cmdline]: Wrong value for blocksize"<<endl;
    exit(1);
  }
  
  if(argc==9){
    nsamples1=atoi(argv[8]);
    if(nsamples1<=0) {
      cerr<<"[read_cmdline]: Wrong value for bootstrap1"<<endl;
      exit(1);
    }
    nsamples2=atoi(argv[9]);
    if(nsamples2<=0) {
      cerr<<"[read_cmdline]: Wrong value for bootstrap2"<<endl;
      exit(1);
    }
  }
  
  return 0;
}
