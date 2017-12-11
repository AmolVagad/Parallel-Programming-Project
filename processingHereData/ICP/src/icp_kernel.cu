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

__global__  void PerformRotationKernel( Matrix t, Matrix Point,Matrix New_Point)
{

	// Create Matrices in the shared memory 
	
	__shared__ float t_s[3][1];
	__shared__ float Point_s[TILE_WIDTH][3];
	
	// Thread allocation 
	int bx = blockIdx.x; int by = blockIdx.y;
	int tx = threadIdx.x ; int ty = threadIdx.y;

	// Generate rows and cols of Result point 
	int row = by*TILE_WIDTH + ty;
	int col = bx*TILE_WIDTH + tx;
	float result = 0.00;

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

	// carrying out the actual multiplication 
	 if(ty < TILE_WIDTH && tx < TILE_WIDTH)
	 {
	 	for(int i = 0; i < 3; i++)
		{
    			for(int j = 0; j < 3; j++)
    			{
				result += R_constant[i*3 + j]*Point_s[i+ty][j]+ t_s[i][j];
    			}
		}
	 }


         __syncthreads();        // To ensure all elements of tile are loaded and consumed 

	if(row < 3 && col < t.width)                    //Checking the boundary conditions for matrix new_point
	 New_Point.elements[row*t.width + col] = result;

	

}

/// Kernel function to find the error : 

/*

__global__  void FindTotalErrorInCloud_Kernel(Vector index, Vector bin_index_x,Vector bin_index_y,Vector bin_index_z, float icp_error , Vector x_coord, Vector y_coord, Vector z_coord)
{

	int i = blockdim.x*blockIdx.x + threadIdx.x;
	float error = 0.0f;
	int x_Idx = bin_index_x.elements[i];
	int y_Idx = bin_index_y.elements[i];
	int z_Idx = bin_index_z.elements[i];
	int j = index.elements[i];

	icp_error +=sqrt(pow((transformed_data.x_coord[i] - octree_icp(x_Idx, y_Idx, z_Idx)[3*j]),2) + pow((transformed_data.y_coord[i] - octree_icp(x_Idx, y_Idx, z_Idx)[3*j + 1]),2) + pow((transformed_data.z_coord[i] - octree_icp(x_Idx, y_Idx, z_Idx)[3*j + 2]),2)); 



}

*/















#endif // #ifndef _ICP_KERNEL_H_
