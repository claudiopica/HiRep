// Functions for easy timing of CUDA execution
//
// How to use the functions:
//
//		cudaEvent_t gpu_event = gpuTimerStart();
//		
//		... Call a kernel ...
//
//		float intermediateTime=gpuTimerMeasure(gpu_event);
//
//		... Call another kernel ...
//
//		float elapsedTime=gpuTimerStop(gpu_event);
#include "gpu.h"


cudaEvent_t gpuTimerStart(){
	cudaEvent_t e;	
	cudaEventCreate(&e);	
	cudaEventRecord(e,0);
	return e;
}

float gpuTimerMeasure(cudaEvent_t e){
	float t;
	cudaEvent_t stop;	
	cudaEventCreate(&stop);  
	cudaEventRecord(stop,0);						
	cudaEventSynchronize(stop);						
	cudaEventElapsedTime(&t, e, stop);
	cudaEventDestroy(stop);	
	return t;
}

float gpuTimerStop(cudaEvent_t e){
	float t;
	cudaEvent_t stop;	
	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);						
	cudaEventSynchronize(stop);						
	cudaEventElapsedTime(&t, e, stop);
	cudaEventDestroy(stop);	
	cudaEventDestroy(e);	
	return t;
}

