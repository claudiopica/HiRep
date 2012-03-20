#ifndef GLOBAL_SUM_GPU_C
#define GLOBAL_SUM_GPU_C

#include "global.h"
#include "gpu.h"

#define GSUM_BLOCK_SIZE 256     // No more than 1024 on Tesla
#define BLOCK_SIZE_REM 32

template<typename REAL>
__global__ void copy_with_zero_padding(REAL* out, REAL* in, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  if (i<N){
    out[i]=in[i];
  }
  else{
    out[i]=0.0;
  }
}

//reduction of the last elements
template<typename REAL>
__device__ void warpReduce(volatile REAL* sdata, int tid){
  sdata[tid]+=sdata[tid+32];
  sdata[tid]+=sdata[tid+16];
  sdata[tid]+=sdata[tid+8];
  sdata[tid]+=sdata[tid+4];
  sdata[tid]+=sdata[tid+2];
  sdata[tid]+=sdata[tid+1];
}

////////////////////// Optimized Global Sum Kernel //////////////////////////
// Requires multiple kernel calls in order to produce global synchronization 
// (You might think that this gives some overhead, but you are wrong.)
template<unsigned int blockSize, typename REAL, typename REALOUT>
__global__ void global_sum_gpu_opt(REAL * g_idata, REALOUT * g_odata, unsigned int n) {   

  __shared__ REALOUT sdata[blockSize];
  unsigned int tid = threadIdx.x; 
  unsigned int i = blockIdx.x*(blockSize*2) + tid; 
  unsigned int gridSize = blockSize*2*gridDim.x; 

  sdata[tid] = 0;
	
  while (i < n) 
    { 
      sdata[tid] += g_idata[i] + g_idata[i+blockSize]; 
      i += gridSize;						// gridSize loop maintains coalescing
    } 
  __syncthreads();
	
  // For Tesla cards we need not consider blockSize higher than 1024
  if (blockSize >= 1024) {				// Evaluated at compile time
    if (tid < 512) { 
      sdata[tid] += sdata[tid + 512]; 
    } 
    __syncthreads(); 
  } 
  if (blockSize >= 512) {				// Evaluated at compile time
    if (tid < 256) { 
      sdata[tid] += sdata[tid + 256]; 
    } 
    __syncthreads(); 
  } 
  if (blockSize >= 256) {				// Evaluated at compile time 
    if (tid < 128) { 
      sdata[tid] += sdata[tid + 128]; 
    } 
    __syncthreads(); 
  } 
  if (blockSize >= 128) { 				// Evaluated at compile time
    if (tid <	64) { 
      sdata[tid] += sdata[tid +	64]; 
    }
    __syncthreads(); 
  }
  if (tid < 32) warpReduce(sdata, tid); 	// Total unrolling and no sync in last WARP
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template<unsigned int blockSize, typename REAL>
__global__ void global_sum_gpu_remainder(REAL * in, unsigned int n) {
  __shared__ volatile REAL sdata[blockSize];
  unsigned int tid = threadIdx.x;
  
	sdata[tid]=in[tid];
	sdata[tid+blockDim.x]=in[tid+blockDim.x];
	
  if (blockSize >=64) {sdata[tid]+=sdata[tid+32];}
  if (blockSize >=32) {sdata[tid]+=sdata[tid+16];}
  if (blockSize >=16) {sdata[tid]+=sdata[tid+8];}
  if (blockSize >=8) {sdata[tid]+=sdata[tid+4];}
  if (blockSize >=4) {sdata[tid]+=sdata[tid+2];}
  if (blockSize >=2) {sdata[tid]+=sdata[tid+1];}
  
  
  __syncthreads();
  if (tid == 0)  in[0] = sdata[0];
}


//Reduction of the last elements
template<typename COMPLEX>
__device__ void warpReduce_complex(volatile COMPLEX* sdata, int tid){
  _complex_add(sdata[tid],sdata[tid],sdata[tid+32]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+16]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+8]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+4]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+2]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+1]);
}

////////////////////////////////////////////////////// Optimized Global Sum Kernel
// Requires multiple kernel calls in order to produce global synchronization 
// (You might think that this gives some overhead, but you are wrong.)
template<unsigned int blockSize, typename COMPLEX>
__global__ void global_complex_sum_gpu_opt(COMPLEX * g_idata, complex * g_odata, unsigned int n) {   

  __shared__ complex sdata[blockSize];
  unsigned int tid = threadIdx.x; 
  unsigned int i = blockIdx.x*(blockSize*2) + tid; 
  unsigned int gridSize = blockSize*2*gridDim.x; sdata[tid].re = 0,sdata[tid].im = 0;
	
  while (i < n) 
    { 
      sdata[tid].re += g_idata[i].re + g_idata[i+blockSize].re; 
      sdata[tid].im += g_idata[i].im + g_idata[i+blockSize].im; 
      i += gridSize;	
    } 
  __syncthreads();
	
  // For Tesla cards we need not consider blockSize higher than 1024
  if (blockSize >= 1024) {				// Evaluated at compile time
    if (tid < 512) { 
      sdata[tid].re += sdata[tid + 512].re; 
      sdata[tid].im += sdata[tid + 512].im; 
    } 
    __syncthreads(); 
  } 
  if (blockSize >= 512) {				// Evaluated at compile time
    if (tid < 256) { 
      sdata[tid].re += sdata[tid + 256].re; 
      sdata[tid].im += sdata[tid + 256].im; 
    } 
    __syncthreads(); 
  } 
  if (blockSize >= 256) {				// Evaluated at compile time 
    if (tid < 128) { 
      sdata[tid].re += sdata[tid + 128].re; 
      sdata[tid].im += sdata[tid + 128].im; 
    } 
    __syncthreads(); 
  } 
  if (blockSize >= 128) { 				// Evaluated at compile time
    if (tid <	64) { 
      sdata[tid].re += sdata[tid + 64].re; 
      sdata[tid].im += sdata[tid + 64].im; 
    }
    __syncthreads(); 
  }
  if (tid < 32) warpReduce_complex(sdata, tid); 	// Total unrolling and no sync in last WARP 
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
	
}
template<unsigned int blockSize, typename COMPLEX>
__global__ void global_complex_sum_gpu_remainder(COMPLEX * in, unsigned int n) {   
  __shared__ volatile COMPLEX sdata[blockSize];
  unsigned int tid = threadIdx.x;
	
	sdata[tid].re=in[tid].re;
	sdata[tid].im=in[tid].im;
	sdata[tid+blockDim.x].re=in[tid+blockDim.x].re;
	sdata[tid+blockDim.x].im=in[tid+blockDim.x].im;
	
	
  if (blockSize >=64) {_complex_add(sdata[tid],sdata[tid],sdata[tid+32]);}  
  if (blockSize >=32) {_complex_add(sdata[tid],sdata[tid],sdata[tid+16]);}
  if (blockSize >=16) {_complex_add(sdata[tid],sdata[tid],sdata[tid+8]);}
  if (blockSize >=8) {_complex_add(sdata[tid],sdata[tid],sdata[tid+4]);}
  if (blockSize >=4) {_complex_add(sdata[tid],sdata[tid],sdata[tid+2]);}
  if (blockSize >=2) {_complex_add(sdata[tid],sdata[tid],sdata[tid+1]);}
  
  
  __syncthreads();
  if (tid == 0) {in[0].re = sdata[0].re;in[0].im = sdata[0].im;};
}



unsigned int next_pow2( unsigned int n )
{
  register unsigned int val = n;
  val--;
  val = (val >> 1) | val;
  val = (val >> 2) | val;
  val = (val >> 4) | val;
  val = (val >> 8) | val;
  val = (val >> 16) | val;
  return ++val; // Val is now the next highest power of 2.
} 

void global_reduction_sum(double* resField, unsigned int Npow2){
  unsigned int new_gridDim 		= 	Npow2/GSUM_BLOCK_SIZE ;//(Npow2/(2*GSUM_BLOCK_SIZE) > 32767) ? 32768 : Npow2/GSUM_BLOCK_SIZE ;
  unsigned int new_N 			= 	Npow2;
  
  while(new_gridDim>1){
    //    printf(" \t new_gridDim/2 = %d \n", new_gridDim/2);
    
    global_sum_gpu_opt<GSUM_BLOCK_SIZE><<<new_gridDim/2,GSUM_BLOCK_SIZE>>>(resField,resField,new_N);
    new_N=new_gridDim/2;
    new_gridDim=new_N/(GSUM_BLOCK_SIZE);
    CudaCheckError(); 
  }
  if (new_N > 64){
	  //	printf(" \t Intermediate step:\n \t new_N/(BLOCK_SIZE_REM*2) = %d \n", new_N/(BLOCK_SIZE_REM*2));	
 //   	printf(" \t new_N = %d \n", new_N);
    global_sum_gpu_opt<BLOCK_SIZE_REM><<<new_N/(BLOCK_SIZE_REM*2),BLOCK_SIZE_REM>>>(resField,resField,new_N);	
    new_N=new_N/(BLOCK_SIZE_REM*2);
    CudaCheckError(); 
  }

	  switch (new_N) {
		  case 64:
			  global_sum_gpu_remainder<64><<<1,32>>>(resField,new_N);
			  CudaCheckError();
			  break;
		  case 32:
			  global_sum_gpu_remainder<32><<<1,16>>>(resField,new_N);
			  CudaCheckError();
			  break;
		  case 16:
			  global_sum_gpu_remainder<16><<<1,8>>>(resField,new_N);
			  CudaCheckError();
			  break;
		  case 8:
			  global_sum_gpu_remainder<8><<<1,4>>>(resField,new_N);
			  CudaCheckError();
			  break;
		  case 4:
			  global_sum_gpu_remainder<4><<<1,2>>>(resField,new_N);
			  CudaCheckError();
			  break;
		  case 2:
			  global_sum_gpu_remainder<2><<<1,1>>>(resField,new_N);
			  CudaCheckError();
			  break;
	   default:
		   break;

	  }
}

void global_reduction_complex_sum(complex* resField, unsigned int Npow2){
  unsigned int new_gridDim 		= 	Npow2/GSUM_BLOCK_SIZE;
  unsigned int new_N 			= 	Npow2;
  
  while(new_gridDim>1){
    global_complex_sum_gpu_opt<GSUM_BLOCK_SIZE><<<new_gridDim/2,GSUM_BLOCK_SIZE>>>(resField,resField,new_N);
    new_N=new_gridDim/2;
    new_gridDim=new_N/(GSUM_BLOCK_SIZE);
    CudaCheckError(); 
  }
  if (new_N > 64){
    global_complex_sum_gpu_opt<BLOCK_SIZE_REM><<<new_N/(BLOCK_SIZE_REM*2),BLOCK_SIZE_REM>>>(resField,resField,new_N);	
    new_N=new_N/(BLOCK_SIZE_REM*2);
    CudaCheckError(); 
  }
	
	switch (new_N) {
		case 64:
			global_complex_sum_gpu_remainder<64><<<1,32>>>(resField,new_N);
			CudaCheckError();
			break;
		case 32:
			global_complex_sum_gpu_remainder<32><<<1,16>>>(resField,new_N);
			CudaCheckError();
			break;
		case 16:
			global_complex_sum_gpu_remainder<16><<<1,8>>>(resField,new_N);
			CudaCheckError();
			break;
		case 8:
			global_complex_sum_gpu_remainder<8><<<1,4>>>(resField,new_N);
			CudaCheckError();
			break;
		case 4:
			global_complex_sum_gpu_remainder<4><<<1,2>>>(resField,new_N);
			CudaCheckError();
			break;
		case 2:
			global_complex_sum_gpu_remainder<2><<<1,1>>>(resField,new_N);
			CudaCheckError();
			break;
		default:
			break;
			
	}
	
}


//OLD VERSIONS
//Global sum calculation
template<typename REAL>
__global__ void global_sum_gpu(REAL* in, double* out, int N){
  __shared__ double sdata[GSUM_BLOCK_SIZE];
  int tid = threadIdx.x;
  int i;
  //  REAL res;
  sdata[tid]=in[tid];
  for (i=tid+GSUM_BLOCK_SIZE;i<N;i+= GSUM_BLOCK_SIZE){ // Sum over all blocks 
    sdata[tid]+=(double)in[i];
  }
  __syncthreads();
  
  for (i = GSUM_BLOCK_SIZE/2;i>32;i/=2){
    if (tid < i) sdata[tid]+=sdata[tid+i];
    __syncthreads();
  }
  if (tid<32){ //Unroll all the threads in a WARP
    warpReduce(sdata,tid);
  }
  __syncthreads();
  if (tid == 0)  out[0] = sdata[0];
}


//Global sum calculation  - OLD VERSION
template<typename COMPLEX>
__global__ void global_sum_complex_gpu(COMPLEX* in, complex* out, int N){
  __shared__ complex sdata[GSUM_BLOCK_SIZE];
  int tid = threadIdx.x;
  int i;
  //  REAL res;
  sdata[tid].re=in[tid].re;
  sdata[tid].im=in[tid].im;
  for (i=tid+GSUM_BLOCK_SIZE;i<N;i+= GSUM_BLOCK_SIZE){ // Sum over all blocks
    _complex_add(sdata[tid],sdata[tid],in[i]);    
  }
  __syncthreads();
  
  for (i = GSUM_BLOCK_SIZE/2;i>32;i/=2){
    if (tid < i){
      _complex_add(sdata[tid],sdata[tid],sdata[tid+i]);
    }
    __syncthreads();
  }
  if (tid<32){ //Unroll all the threads in a WARP
    warpReduce_complex(sdata,tid);
  }
  __syncthreads();
  if (tid == 0)  out[0] = sdata[0];
}


double global_sum_gpu(double* vector, int n){
	unsigned int new_n = next_pow2(n);
	int grid_size = new_n / BLOCK_SIZE;
	double* padded_vector;
	double res,res2;
	cudaMalloc((void **) &padded_vector,new_n*sizeof(double));
	copy_with_zero_padding<<<grid_size,BLOCK_SIZE>>>(padded_vector, vector, n);
	global_reduction_sum(padded_vector,new_n);
	cudaMemcpy(&res,padded_vector,sizeof(res),cudaMemcpyDeviceToHost);
	cudaFree(padded_vector);
	return res;
}

#endif
