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
#ifdef WITH_GPU
#include "gpu.h"
#include "utils.h"


gpu_timer gpuTimerStart(){
	cudaEvent_t e;	
	cudaEventCreate(&e);	
	cudaEventRecord(e,0);
	return e;
}

float gpuTimerMeasure(gpu_timer e){
	float t;
	cudaEvent_t stop;	
	cudaEventCreate(&stop);  
	cudaEventRecord(stop,0);						
	cudaEventSynchronize(stop);						
	cudaEventElapsedTime(&t, e, stop);
	cudaEventDestroy(stop);	
	return t;
}

float gpuTimerStop(gpu_timer e){
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

#endif //WITH_GPU
