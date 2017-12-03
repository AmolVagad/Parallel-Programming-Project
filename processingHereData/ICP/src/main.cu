// HOST Code to compute ICP for localization 


#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_
#include <stdio.h>
#include "ICP_Compute.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include<vector>
#include<time.h>
#include<sys/time.h>
#include<ctime>
#include "dlib/optimization/optimization.h"
#include "dlib/optimization/find_optimal_parameters_abstract.h"
#include "dlib/optimization/optimization_bobyqa.h"
#include "dlib/optimization/find_optimal_parameters.h" 
#include "octree_code/octree.h"
#include "globals.h"
#include "icp_kernel.cu"
extern "C"

///// For initial testing purposes carrying out rotation and translation operation on cuda//////////////////




dlib::matrix<double> PerformRotation(dlib::matrix<double> R_h,dlib::matrix<double> t_h, dlib::matrix<double> Point_h,  dlib::matrix<double> Rotated_Point_h)
{
	int n = 3;
	int size = n*sizeof(double);

	// Declare the device variables 
	
	dlib::matrix<double> R_d;
	dlib::matrix<double> t_d;
	dlib::matrix<double> Point_d;
	dlib::matrix<double> Rotated_Point_d;
	
	// Declare the matrix to store the retured result 
	dlib::matrix<double> New_Point_h(3,1);
	
	// Allocate memory on the device
	
	cudaMalloc((void **)&R_d,size));
	cudaMalloc((void **)&t_d,size));
	cudaMalloc((void **)&Point_d,size));
	
	// Copy from host to device 
	
	cudaMemcpy(R_d,R_h,size, cudaMemcpyHostToDevice);
	cudaMemcpy(t_d,t_h,size, cudaMemcpyHostToDevice);
	cudaMemcpy(Point_d,Point_h,size, cudaMemcpyHostToDevice);
	
	
	// Allocate device memory for result 
	cudaMalloc((void **)&Rotated_Point_d,size));
	
	
	// Kernel Call 
	 
	 
	 
	 // Setup the execution configuration
    	int blocks_w = N.width/TILE_WIDTH ;
    	int blocks_h = M.height /TILE_WIDTH;
    
   
	if(M.width % TILE_WIDTH)

	{

		blocks_w ++;

	}
	if(M.height % TILE_WIDTH)

	{

		blocks_h ++;

	}
         
	dim3 dimGrid(blocks_w, blocks_h, 1);

	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH,1);

    // Launch the device computation threads!

     PerformRotationKernel<<<dimGrid,dimBlock>>>(P_d, t_d, Point_d, New_Point_d);
	
	
	
	// Transfer Rotated Point from device to host
     cudaMemcpy(Rotated_Point_h, Rotated_Point_d, size, cudaMemcpyDeviceToHost);
       // Free device memory for all
     cudaFree(P_d); cudaFree(t_d); cudaFree (Point_d);cudaFree (Rotated_Point_d);
     
     return Rotated_Point_h;
	
	
}
	
