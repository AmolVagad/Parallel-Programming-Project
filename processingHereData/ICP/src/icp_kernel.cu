// Device code for ICP computation 
// Currently working only on performing rotation and translation using cuda 


#ifndef _ICP_KERNEL_H_
#define _ICP_KERNEL_H_
#include <stdio.h>

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


#define TILE_WIDTH 4 

__global__  void PerformRotationKernel(Matrix R, Matrix t, Matrix Point,Matrix New_Point)
{
	


	
	

	// Create Matrices in the shared memory 
	__shared__ float R_s[3][3];
	__shared__ float t_s[3][1];
	__shared__ float Point_s[TILE_WIDTH][3];
	
	// Thread allocation 
	int bx = blockIdx.x; int by = blockIdx.y;
	int tx = threadIdx.x ; int ty = threadIdx.y;

	// Generate rows and cols of Result point 
	int row = by*TILE_WIDTH + ty;
	int col = bx*TILE_WIDTH + tx;
	float result = 0.00;

	//Loading R, t and point into the shared memory 
	 
	for ( int m = 0; m <= (R.width /TILE_WIDTH)  ; ++m)
	{
		if(row < R.height && (m*TILE_WIDTH + tx) < R.width)   //Checking the boundary conditions for matrix M               
			R_s[ty][tx] = R.elements[row*R.width + m*TILE_WIDTH + tx];
		else
			R_s[ty][tx] = 0;
	}
	
	for (int n = 0; n <= (Point.width/TILE_WIDTH); ++n)
	{
		if((n*TILE_WIDTH + ty) < Point.height && col < Point.width)   //Checking the boundary conditions for matrix Point
			Point_s[ty][tx] = Point.elements[(n*TILE_WIDTH + ty)*Point.width + col];
		else 
			Point_s[ty][tx] = 0;
		if((n*TILE_WIDTH + ty) < t.height && col < t.width)   //Checking the boundary conditions for matrix t
			t_s[ty][tx] = t.elements[(n*TILE_WIDTH + ty)*t.width + col];
		else 
			t_s[ty][tx] = 0;
	
	}
	
		__syncthreads();         // To ensure all elements of tile are loaded and consumed 

		for(int k = 0; k < TILE_WIDTH; ++k)
		{	
			result += R_s[ty][k]*Point_s[k][tx] + t_s[k][tx];
	        
	
		}
		__syncthreads();        // To ensure all elements of tile are loaded and consumed 


        
	if(row < R.height && col < Point.width)                    //Checking the boundary conditions for matrix new_point
		 New_Point.elements[row*Point.width + col] = result;

	

}

#endif // #ifndef _ICP_KERNEL_H_
