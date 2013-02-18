#include "global.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "communications.h"
#include <stdio.h>

void polyakov() {
  int x[4], i4d, i3d, size3d;
  int loc[4]={T,X,Y,Z};
  suNg *p;
  suNg *lp;
  suNg *bp;
  suNg tmp;
  complex poly;
  double adjpoly;
  double dtmp;
#ifdef WITH_MPI  
  int np[4]={NP_T,NP_X,NP_Y,NP_Z};
  int mpiret;
#endif /* WITH_MPI */
 
  for(int mu=0; mu<4; mu++){					// Loop over directions mu

    size3d=T*X*Y*Z/loc[mu];
    p=malloc(sizeof(suNg)*size3d);
    lp=malloc(sizeof(suNg)*size3d);
    bp=malloc(sizeof(suNg)*size3d);
    
    /* local wilson lines */
    
    i3d=0;
    for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)				// Loop over volume at slice mu
    for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
    for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
      _suNg_unit(lp[i3d]);
      for(x[mu]=0; x[mu]<loc[mu]; x[mu]++) {										// Loop over mu-slices
        i4d=ipt(x[0],x[1],x[2],x[3]);
        _suNg_times_suNg(tmp,lp[i3d],*pu_gauge(i4d,mu));									// Do product in lp
        lp[i3d] = tmp;
      }																		//
      i3d++;
    }																		//
  
  
    
    /* global wilson lines */
    
  #ifdef WITH_MPI  
    if(COORD[mu]!=0) {
      MPI_Status status;
      mpiret=MPI_Recv(bp, /* buffer */
          size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
          MPI_DOUBLE, /* basic datatype */
          proc_dn(CID,mu), /* cid of origin */
          COORD[mu]-1, /* tag of communication */
          cart_comm, /* use the cartesian communicator */
          &status);
  #ifndef NDEBUG
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        if (status.MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status.MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              0, 
              status.MPI_SOURCE, 
              status.MPI_TAG, 
              mesg);
        }
        error(1,1,"WL_Hamiltonian_gauge [wilsonloops.c]","Cannot receive buf_gtf[1]");
      }
  #endif /* NDEBUG */
    }
  #endif /* WITH_MPI */
  
  
    if(COORD[mu]==0) {												// If mu=0 copy result (lp) to p
      i3d=0;
      for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)
      for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
      for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
        p[i3d]=lp[i3d];
        i3d++;
      }
    } else {
      i3d=0;														// If mu=0 copy result (lp) to p
      for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)
      for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
      for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
        _suNg_times_suNg(p[i3d], bp[i3d], lp[i3d]);
        i3d++;
      }
    }
    
    
  #ifdef WITH_MPI  
    if(COORD[mu]!=np[mu]-1) {
      mpiret=MPI_Send(p, /* buffer */
          size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
          MPI_DOUBLE, /* basic datatype */
          proc_up(CID,mu), /* cid of destination */
          COORD[mu], /* tag of communication */
          cart_comm /* use the cartesian communicator */
          );
  #ifndef NDEBUG
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        error(1,1,"WL_Hamiltonian_gauge [wilsonloops.c]","Cannot send buf_gtf[0]");
      }
  #endif /* NDEBUG */
    }
  #endif /* WITH_MPI */
  
  
  
    /* broadcast wilson lines */
  
  #ifdef WITH_MPI
    int mpiret;
    if(COORD[mu]==np[mu]-1) {
      MPI_Request comm_req[np[mu]-1];
      int destCID=CID;
      for(int t=0; t<np[mu]-1; t++) {
        destCID=proc_up(destCID,mu);
        mpiret=MPI_Isend(p, /* buffer */
            size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
            MPI_DOUBLE, /* basic datatype */
            destCID, /* cid of destination */
            t, /* tag of communication */
            cart_comm, /* use the cartesian communicator */
            &(comm_req[t]) /* handle to communication request */
            );
  #ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot start send polyakov");
        }
  #endif /* NDEBUG */
      }
      
      MPI_Status status[np[mu]-1];
      mpiret=MPI_Waitall(np[mu]-1, comm_req, status);
  #ifndef NDEBUG
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen, k;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        for (k=0; k<np[mu]-1; ++k) {
          if (status[k].MPI_ERROR != MPI_SUCCESS) {
            MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
            lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                k, 
                status[k].MPI_SOURCE, 
                status[k].MPI_TAG, 
                mesg);
          }
        }
        error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot complete communications");
      }
  #endif /* NDEBUG */
  
    } else {
  
      int sCOORD[4], sCID;
      sCOORD[0]=COORD[0];sCOORD[1]=COORD[1];sCOORD[2]=COORD[2];sCOORD[3]=COORD[3];
      sCOORD[mu]=np[mu]-1;
      mpiret=MPI_Cart_rank(cart_comm, sCOORD, &sCID);
  #ifndef NDEBUG
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot retrieve source CID");
      }
  #endif /* NDEBUG */
      MPI_Status status;
      mpiret=MPI_Recv(p, /* buffer */
          size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
          MPI_DOUBLE, /* basic datatype */
          sCID, /* cid of destination */
          COORD[mu], /* tag of communication */
          cart_comm, /* use the cartesian communicator */
          &status);
  #ifndef NDEBUG
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        if (status.MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status.MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              0, 
              status.MPI_SOURCE, 
              status.MPI_TAG, 
              mesg);
        }
        error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot receive polyakov");
      }
  #endif /* NDEBUG */
    }
  #endif /* WITH_MPI */
  
  
    
    /* trace and average */
    poly.re=poly.im=0.0;
    adjpoly=0.;
  
    i3d=0;
    for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)	// Loop over all sites in 3d volume orthogonal to mu
    for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
    for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
      _suNg_trace_re(dtmp,p[i3d]);
/*      lprintf("LOC_POLYAKOV",0,"%d %d %d %d %d %1.8e\n",mu,x[(mu+1)%4],x[(mu+2)%4],x[(mu+3)%4],i3d,dtmp); */
      poly.re += dtmp;
      adjpoly +=dtmp*dtmp;
      //_suNg_trace_im(dtmp,p[i3d]);
        dtmp=0.;
        poly.im += dtmp;
      adjpoly +=dtmp*dtmp - 1;
      i3d++;
    }
  
    global_sum((double*)&poly,2);
    global_sum((double*)&adjpoly,1);
  
    
    poly.re /= NG*GLB_X*GLB_Y*GLB_Z*GLB_T/loc[mu];
    poly.im /= NG*GLB_X*GLB_Y*GLB_Z*GLB_T/loc[mu];
    adjpoly /= NG*GLB_X*GLB_Y*GLB_Z*GLB_T/loc[mu];
    
    lprintf("FUND_POLYAKOV",0,"Polyakov direction %d = %1.8e %1.8e\n",mu,poly.re,poly.im);
    lprintf("ADJ_POLYAKOV",0,"Polyakov direction %d = %1.8e\n",mu,adjpoly);
    
    free(p);
    free(lp);
    free(bp);
    
  }
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void polyakov_in_time() {
	int x[4], i4d, i3d, size3d;
	int loc[4]={T,X,Y,Z};
	suNg *p;
	suNg *lp;
	suNg *bp;
	suNg tmp;
	complex poly;
	double adjpoly;
	double dtmp;
	double *poly_t;
	poly_t=malloc(sizeof(double)*T);
#ifdef WITH_MPI  
	int np[4]={NP_T,NP_X,NP_Y,NP_Z};
	int mpiret;
#endif /* WITH_MPI */
	
	for(int mu=1; mu<4; mu++){					// Loop over directions mu
		
		size3d=T*X*Y*Z/loc[mu];
		p=malloc(sizeof(suNg)*size3d);
		lp=malloc(sizeof(suNg)*size3d);
		bp=malloc(sizeof(suNg)*size3d);
		
		
		/* local wilson lines */
		
		i3d=0;
		for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)				// Loop over volume at slice mu
			for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
				for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
					_suNg_unit(lp[i3d]);
					for(x[mu]=0; x[mu]<loc[mu]; x[mu]++) {										// Loop over mu-slices
						i4d=ipt(x[0],x[1],x[2],x[3]);
						_suNg_times_suNg(tmp,lp[i3d],*pu_gauge(i4d,mu));									// Do product in lp
						lp[i3d] = tmp;
					}																		//
					i3d++;
				}																		//
		
		
		
		/* global wilson lines */
		
#ifdef WITH_MPI  
		if(COORD[mu]!=0) {
			MPI_Status status;
			mpiret=MPI_Recv(bp, /* buffer */
							size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
							MPI_DOUBLE, /* basic datatype */
							proc_dn(CID,mu), /* cid of origin */
							COORD[mu]-1, /* tag of communication */
							cart_comm, /* use the cartesian communicator */
							&status);
#ifndef NDEBUG
			if (mpiret != MPI_SUCCESS) {
				char mesg[MPI_MAX_ERROR_STRING];
				int mesglen;
				MPI_Error_string(mpiret,mesg,&mesglen);
				lprintf("MPI",0,"ERROR: %s\n",mesg);
				if (status.MPI_ERROR != MPI_SUCCESS) {
					MPI_Error_string(status.MPI_ERROR,mesg,&mesglen);
					lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
							0, 
							status.MPI_SOURCE, 
							status.MPI_TAG, 
							mesg);
				}
				error(1,1,"WL_Hamiltonian_gauge [wilsonloops.c]","Cannot receive buf_gtf[1]");
			}
#endif /* NDEBUG */
		}
#endif /* WITH_MPI */
		
		
		if(COORD[mu]==0) {												// If mu=0 copy result (lp) to p
			i3d=0;
			for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)
				for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
					for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
						p[i3d]=lp[i3d];
						i3d++;
					}
		} else {
			i3d=0;														// If mu=0 copy result (lp) to p
			for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)
				for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
					for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
						_suNg_times_suNg(p[i3d], bp[i3d], lp[i3d]);
						i3d++;
					}
		}
		
		
#ifdef WITH_MPI  
		if(COORD[mu]!=np[mu]-1) {
			mpiret=MPI_Send(p, /* buffer */
							size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
							MPI_DOUBLE, /* basic datatype */
							proc_up(CID,mu), /* cid of destination */
							COORD[mu], /* tag of communication */
							cart_comm /* use the cartesian communicator */
							);
#ifndef NDEBUG
			if (mpiret != MPI_SUCCESS) {
				char mesg[MPI_MAX_ERROR_STRING];
				int mesglen;
				MPI_Error_string(mpiret,mesg,&mesglen);
				lprintf("MPI",0,"ERROR: %s\n",mesg);
				error(1,1,"WL_Hamiltonian_gauge [wilsonloops.c]","Cannot send buf_gtf[0]");
			}
#endif /* NDEBUG */
		}
#endif /* WITH_MPI */
		
		
		
		/* broadcast wilson lines */
		
#ifdef WITH_MPI
		int mpiret;
		if(COORD[mu]==np[mu]-1) {
			MPI_Request comm_req[np[mu]-1];
			int destCID=CID;
			for(int t=0; t<np[mu]-1; t++) {
				destCID=proc_up(destCID,mu);
				mpiret=MPI_Isend(p, /* buffer */
								 size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
								 MPI_DOUBLE, /* basic datatype */
								 destCID, /* cid of destination */
								 t, /* tag of communication */
								 cart_comm, /* use the cartesian communicator */
								 &(comm_req[t]) /* handle to communication request */
								 );
#ifndef NDEBUG
				if (mpiret != MPI_SUCCESS) {
					char mesg[MPI_MAX_ERROR_STRING];
					int mesglen;
					MPI_Error_string(mpiret,mesg,&mesglen);
					lprintf("MPI",0,"ERROR: %s\n",mesg);
					error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot start send polyakov");
				}
#endif /* NDEBUG */
			}
			
			MPI_Status status[np[mu]-1];
			mpiret=MPI_Waitall(np[mu]-1, comm_req, status);
#ifndef NDEBUG
			if (mpiret != MPI_SUCCESS) {
				char mesg[MPI_MAX_ERROR_STRING];
				int mesglen, k;
				MPI_Error_string(mpiret,mesg,&mesglen);
				lprintf("MPI",0,"ERROR: %s\n",mesg);
				for (k=0; k<np[mu]-1; ++k) {
					if (status[k].MPI_ERROR != MPI_SUCCESS) {
						MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
						lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
								k, 
								status[k].MPI_SOURCE, 
								status[k].MPI_TAG, 
								mesg);
					}
				}
				error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot complete communications");
			}
#endif /* NDEBUG */
			
		} else {
			
			int sCOORD[4], sCID;
			sCOORD[0]=COORD[0];sCOORD[1]=COORD[1];sCOORD[2]=COORD[2];sCOORD[3]=COORD[3];
			sCOORD[mu]=np[mu]-1;
			mpiret=MPI_Cart_rank(cart_comm, sCOORD, &sCID);
#ifndef NDEBUG
			if (mpiret != MPI_SUCCESS) {
				char mesg[MPI_MAX_ERROR_STRING];
				int mesglen;
				MPI_Error_string(mpiret,mesg,&mesglen);
				lprintf("MPI",0,"ERROR: %s\n",mesg);
				error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot retrieve source CID");
			}
#endif /* NDEBUG */
			MPI_Status status;
			mpiret=MPI_Recv(p, /* buffer */
							size3d*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
							MPI_DOUBLE, /* basic datatype */
							sCID, /* cid of destination */
							COORD[mu], /* tag of communication */
							cart_comm, /* use the cartesian communicator */
							&status);
#ifndef NDEBUG
			if (mpiret != MPI_SUCCESS) {
				char mesg[MPI_MAX_ERROR_STRING];
				int mesglen;
				MPI_Error_string(mpiret,mesg,&mesglen);
				lprintf("MPI",0,"ERROR: %s\n",mesg);
				if (status.MPI_ERROR != MPI_SUCCESS) {
					MPI_Error_string(status.MPI_ERROR,mesg,&mesglen);
					lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
							0, 
							status.MPI_SOURCE, 
							status.MPI_TAG, 
							mesg);
				}
				error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot receive polyakov");
			}
#endif /* NDEBUG */
		}
#endif /* WITH_MPI */
		
		
		
		/* trace and average */
		poly.re=poly.im=0.0;
		adjpoly=0.;
		
		for(x[0]=0; x[0]<T; x[0]++){								// initialize to zero
			poly_t[x[0]]=0;
		}
		
		
		i3d=0;
		for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)	// Loop over all sites in 3d volume orthogonal to mu
			for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
				for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
					_suNg_trace_re(dtmp,p[i3d]);
//					lprintf("LOC_POLYAKOV",0,"%d %d %d %d %d %1.8e\n",mu,x[(mu+1)%4],x[(mu+2)%4],x[(mu+3)%4],i3d,dtmp); 
					poly.re += dtmp;
					adjpoly +=dtmp*dtmp;
					poly_t[x[0]] += dtmp;
					//_suNg_trace_im(dtmp,p[i3d]);
					dtmp=0.;
					poly.im += dtmp;
					adjpoly +=dtmp*dtmp - 1;
					i3d++;
				}
		
		global_sum((double*)&poly,2);   //MPI communication
		global_sum((double*)&adjpoly,1);
		
		for(x[0]=0; x[0]<T; x[0]++){
			//global_sum((double*)&poly_t[x[0]],1); //MPI communication
			poly_t[x[0]] /= NG*GLB_X*GLB_Y*GLB_Z/loc[mu];
		}
		
		poly.re /= NG*GLB_X*GLB_Y*GLB_Z*GLB_T/loc[mu];
		poly.im /= NG*GLB_X*GLB_Y*GLB_Z*GLB_T/loc[mu];
		adjpoly /= NG*GLB_X*GLB_Y*GLB_Z*GLB_T/loc[mu];
		
		lprintf("FUND_POLYAKOV",0,"Polyakov direction %d = %1.8e %1.8e\n",mu,poly.re,poly.im);
		lprintf("ADJ_POLYAKOV",0,"Polyakov direction %d = %1.8e\n",mu,adjpoly);
		
		for(x[0]=0; x[0]<T; x[0]++){
			lprintf("FUND_POLYAKOV",0,"Polyakov ( d = %d t = %d ) : %1.8e\n",mu,x[0],poly_t[x[0]]);
		}
		
		
		free(p);
		free(lp);
		free(bp);
		
	}
	free(poly_t);
}











////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void polyakov_vol() {
	int x[4], i4d, i3d, size3d;
	int loc[4]={T,X,Y,Z};
	suNg *p;
	suNg *lp;
	suNg *bp;
	suNg tmp;
	double dtmp;
//	FILE * pFile;
//	pFile = fopen ("POLY_LOCAL.txt","w");
//	 fprintf (pFile, "Testing");
	
	for(int mu=1; mu<4; mu++){					// Loop over directions mu
		
		size3d=T*X*Y*Z/loc[mu];
		p=malloc(sizeof(suNg)*size3d);
		lp=malloc(sizeof(suNg)*size3d);
		bp=malloc(sizeof(suNg)*size3d);
		
		
		/* local wilson lines */
		
		i3d=0;
		for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)				// Loop over volume at slice mu
			for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
				for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
					_suNg_unit(lp[i3d]);
					for(x[mu]=0; x[mu]<loc[mu]; x[mu]++) {										// Loop over mu-slices
						i4d=ipt(x[0],x[1],x[2],x[3]);
						_suNg_times_suNg(tmp,lp[i3d],*pu_gauge(i4d,mu));									// Do product in lp
						lp[i3d] = tmp;
					}																		//
					i3d++;
				}																		//
		
		
		
		/* global wilson lines */
		

		
		
		if(COORD[mu]==0) {												// If mu=0 copy result (lp) to p
			i3d=0;
			for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)
				for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
					for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
						p[i3d]=lp[i3d];
						i3d++;
					}
		} else {
			i3d=0;														// If mu=0 copy result (lp) to p
			for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)
				for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
					for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
						_suNg_times_suNg(p[i3d], bp[i3d], lp[i3d]);
						i3d++;
					}
		}
		

		
		i3d=0;
		       //		 fprintf (pFile, "\nDirection: %d , Coordinate order: %d %d %d (Loop nesting)\n",mu,(mu+1)%4,(mu+2)%4,(mu+3)%4);
		      lprintf("LOC_POLYAKOV",0,"\nDirection: %d , Coordinate order: %d %d %d (Loop nesting)\n",mu,(mu+1)%4,(mu+2)%4,(mu+3)%4); 
		for(x[(mu+1)%4]=0; x[(mu+1)%4]<loc[(mu+1)%4]; x[(mu+1)%4]++)	// Loop over all sites in 3d volume orthogonal to mu
			for(x[(mu+2)%4]=0; x[(mu+2)%4]<loc[(mu+2)%4]; x[(mu+2)%4]++)
				for(x[(mu+3)%4]=0; x[(mu+3)%4]<loc[(mu+3)%4]; x[(mu+3)%4]++) {
					_suNg_trace_re(dtmp,p[i3d]);
//					 fprintf (pFile, " %1.8e",dtmp);
					 lprintf ("LOC_POLYAKOV_NUMBERS",0, " %1.8e",dtmp);
			//		      lprintf("LOC_POLYAKOV_NUMBERS",0,"( d = %d ) %d %d %d %1.8e ",mu,x[(mu+1)%4],x[(mu+2)%4],x[(mu+3)%4],dtmp); 
					i3d++;
				}
		
	//	fclose (pFile);
		free(p);
		free(lp);
		free(bp);
		
	}
}










////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plaq_matrix(suNg* rec, int ix,int mu,int nu)
{
	int iy,iz;
	suNg *v1,*v2,*v3,*v4,w1,w2;//w3;
	
	iy=iup(ix,mu);
	iz=iup(ix,nu);
	
	v1=pu_gauge(ix,mu);
	v2=pu_gauge(iy,nu);
	v3=pu_gauge(iz,mu);
	v4=pu_gauge(ix,nu);
	
	_suNg_times_suNg(w1,(*v1),(*v2));
	_suNg_times_suNg(w2,(*v4),(*v3));
	_suNg_times_suNg_dagger((*rec),w1,w2);      
	
	//rec=w3;
}


/*
void Tmonopoles() {
	int ix,t,x,y,z;
	
	lprintf("T-Monopoles",0,"Cubes orthogonal to t-direction: Nesting order t,x,y,z\n"); 
	
	suNg *p1,*p2,*p3;  //   You cannot assign while multiplying (_suNg_times_suNg(w1,(*v1),(*v2));) ... needs two temporaries.
	p1=malloc(sizeof(suNg));
	p2=malloc(sizeof(suNg));
	p3=malloc(sizeof(suNg));
	double p;
	int Mono=0;
	
	for (t=0; t<T; t++) {
		for (x=0; x<X; x++) {
			for (y=0; y<Y; y++) {
				for (z=0; z<Z; z++) {
					
					ix=ipt(t,x,y,z);
					
					
					int  xp1=iup(ix,1);
					int  xp2=iup(ix,2);
					int  xp3=iup(ix,3);
					
					// cube orthogonal to t
					plaq_matrix(p1,ix,2,1);
					plaq_matrix(p2,ix,3,1);
					_suNg_times_suNg((*p3),(*p1),(*p2));
					plaq_matrix(p1,ix,3,2);
					_suNg_times_suNg((*p2),(*p3),(*p1)); // now result is in p2
					
					plaq_matrix(p3,xp1,2,3);
					_suNg_times_suNg((*p1),(*p3),(*p2)); // now result is in p1
					
					plaq_matrix(p3,xp2,1,3);
					_suNg_times_suNg((*p2),(*p3),(*p1)); // now result is in p2
					
					plaq_matrix(p3,xp3,1,2);
					_suNg_times_suNg((*p1),(*p2),(*p3)); // now result is in p1
					
					_suNg_trace_re(p,(*p1));
					
					if (p>0) {
						Mono=1;
					}else {
						Mono=0;
					}
					
					lprintf("TMonopolesOut",0,"%d",Mono);
					
					
				}	
			}	
		}
	}
	
}
*/

void monopoles(int dir) {
	int ix,t,x,y,z;
	suNg *p1,*p2,*p3;  //   You cannot assign while multiplying (_suNg_times_suNg(w1,(*v1),(*v2));) ... needs two temporaries.
	double p;
	int Mono=0;
	int d1=1;
	int d2=2;
	int d3=3;
	
	p1=malloc(sizeof(suNg));
	p2=malloc(sizeof(suNg));
	p3=malloc(sizeof(suNg));
	
	lprintf("Monopoles",0,"Cubes orthogonal to %d-direction: Nesting order t,x,y,z\n",dir);
	lprintf("MonopolesOut",0,"(dir=%d)",dir); 
	
	if (dir==1) {d1=0;}
	if (dir==2) {d2=0;}
	if (dir==3) {d3=0;}
	
	for (t=0; t<T; t++) {
		for (x=0; x<X; x++) {
			for (y=0; y<Y; y++) {
				for (z=0; z<Z; z++) {
					
					ix=ipt(t,x,y,z);
					
					int  xp1=iup(ix,d1);
					int  xp2=iup(ix,d2);
					int  xp3=iup(ix,d3);
					
					// cube orthogonal to t
					plaq_matrix(p1,ix,d2,d1);
					plaq_matrix(p2,ix,d3,d1);
					_suNg_times_suNg((*p3),(*p1),(*p2));
					plaq_matrix(p1,ix,d3,d2);
					_suNg_times_suNg((*p2),(*p3),(*p1)); // now result is in p2
					
					plaq_matrix(p3,xp1,d2,d3);
					_suNg_times_suNg((*p1),(*p3),(*p2)); // now result is in p1
					
					plaq_matrix(p3,xp2,d1,d3);
					_suNg_times_suNg((*p2),(*p3),(*p1)); // now result is in p2
					
					plaq_matrix(p3,xp3,d1,d2);
					_suNg_times_suNg((*p1),(*p2),(*p3)); // now result is in p1
					
					_suNg_trace_re(p,(*p1));
					
					if (p>0) {
						Mono=1;
					}else {
						Mono=0;
					}
					lprintf("MonopolesOut",0,"%d",Mono);
					
					
				}	
			}	
		}
	}
	
	free(p1);
	free(p2);
	free(p3);
}



double plaq_sign(int ix,int mu,int nu)
{
  int iy,iz;
  suNg *v1,*v2,*v3,*v4;
  suNg w1,w2,w3;
  double tr;
  
  iy=iup(ix,mu);
  iz=iup(ix,nu);
  
  v1=pu_gauge(ix,mu);
  v2=pu_gauge(iy,nu);
  v3=pu_gauge(iz,mu);
  v4=pu_gauge(ix,nu);
  
  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,(*v4),(*v3));
  _suNg_times_suNg_dagger(w3,w1,w2);      
  _suNg_trace_re(tr,w3);
  return (tr<0) ? -1 : 1;
}



void twists(){
  double prod,sum;
  int ix;
  int d1,d2,d3,d4;
  int ids[4];
  int dm[4] = {T,X,Y,Z};
  for (d1=0;d1<3;d1++){
    for (d2=d1+1;d2<4;d2++){
      //This order d1,d2,d3,d4, s.t. \epsilon_{d1,d2,d3,d4} is positive
      d3 = (d2+d1+1) % 4;
      while ( (d3==d1) || (d3==d2) ) d3=(d3+1)%4;
      d4 = (d3+1) % 4;
      while ( (d4==d1) || (d4==d2) ) d4=(d4+1)%4;
      sum = 0;
      for (ids[d3] = 0;ids[d3]<dm[d3];ids[d3]++){
	for (ids[d4] = 0;ids[d4]<dm[d4];ids[d4]++){
	  prod = 1;
	  for (ids[d1] = 0; ids[d1] < dm[d1]; ids[d1]++){
	    for (ids[d2] = 0; ids[d2] < dm[d2]; ids[d2]++){
	      ix=ipt(ids[0],ids[1],ids[2],ids[3]);
	      prod *= plaq_sign(ix,d1,d2);
	    }
	  }
	  sum += prod;
	}
      }
      lprintf("TWIST",0,"mu: %d nu: %d twist: %1.5g\n",d1,d2,sum/dm[d3]/dm[d4]);
    }
  }
}

#undef SIGN
      
	






