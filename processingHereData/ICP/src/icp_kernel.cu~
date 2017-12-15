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


#define TILE_WIDTH 256 

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
//Kernel Function to find the bin of the point
__global__ void find_bin_x_kernel(double * x_coord_d, int numPts, int * bin_x_d)
{
	int t = blockIdx.x*blockDim.x + threadIdx.x;
	if(t < numPts)
	{
		int bin_x_temp = floor(((x_coord_d[t]  - x_min_d)/range_x_d)*bin_size_d);
		bin_x_d[t] = max(min(bin_x_temp, bin_size_d - 1), 0);
	}
}

__global__ void find_bin_y_kernel(double * y_coord_d, int numPts, int * bin_y_d)
{
	int t = blockIdx.x*blockDim.x + threadIdx.x;
	if(t < numPts)
	{
		int bin_y_temp = floor(((y_coord_d[t]  - y_min_d)/range_y_d)*bin_size_d);
		bin_y_d[t] = max(min(bin_y_temp, bin_size_d - 1), 0);
	}
}

__global__ void find_bin_z_kernel(double * z_coord_d, int numPts, int * bin_z_d)
{
	int t = blockIdx.x*blockDim.x + threadIdx.x;
	if(t < numPts)
	{
		int bin_z_temp = floor(((z_coord_d[t]  - z_min_d)/range_z_d)*bin_size_d);
		bin_z_d[t] = max(min(bin_z_temp, bin_size_d - 1), 0);
	}
}

__global__ void find_closest_point_i(double point_x, double point_y, double point_z, double * x_coord_model_d, double * y_coord_model_d, double * z_coord_model_d, int * bin_index_d, double * distance_d, int size_data)
{
	__shared__ double distance_s[2*TILE_WIDTH];
	__shared__ unsigned int bin_smallest_index[TILE_WIDTH];
	unsigned int t = threadIdx.x;
	unsigned int start = 2*blockDim.x*blockIdx.x;

	if(start + t < size_data)
	{
		distance_s[t] = sqrt(pow(x_coord_model_d[start + t] - point_x,2) + pow(y_coord_model_d[start + t] - point_y,2) + pow(z_coord_model_d[start + t] - point_z,2));	
		bin_smallest_index[t] = 0; 	
	}
	else
		distance_s[t] = 65535;
	if(start + blockDim.x + t < size_data)
	{
		distance_s[blockDim.x + t] = sqrt(pow(x_coord_model_d[start + blockDim.x + t] - point_x,2) + pow(y_coord_model_d[start + blockDim.x + t] - point_y,2) + pow(z_coord_model_d[start + blockDim.x + t] - point_z,2));
	}
	else
		distance_s[t] = 65535;

	for(unsigned int stride = blockDim.x; stride >= 1; stride >>= 1)
	{
		__syncthreads();
		if(t < stride)
			if(distance_s[t] > distance_s[stride + t])
			{
				bin_smallest_index[t] = start + stride + t;
				distance_s[t] = distance_s[stride + t];
			}
			else
			{
				bin_smallest_index[t] = start + t;	
			}
	}

	if(t == 0)
	{
		distance_d[blockIdx.x] = distance_s[t];
		bin_index_d[blockIdx.x] = bin_smallest_index[t];	
	}

	
}

__global__ void find_closest_point_2(double * distance_d, int * bin_index_d, int size_data)
{
	__shared__ double distance_s[2*TILE_WIDTH];
	__shared__ unsigned int bin_smallest_index[TILE_WIDTH];
	unsigned int t = threadIdx.x;
	unsigned int start = 2*blockDim.x*blockIdx.x;
	
	if(start + t < size_data)
	{
		distance_s[t] = distance_d[start + t];
		bin_smallest_index[t] = bin_index_d[start + t];
	}
	else
	{
		distance_s[t] = 65535;
		bin_smallest_index[t] = 0;
	}
	if(start + blockDim.x + t < size_data)
	{
		distance_s[blockDim.x + t] = distance_d[start + blockDim.x + t];
		bin_smallest_index[blockDim.x + t] = bin_index_d[start + blockDim.x + t];
	}
	else
	{
		distance_s[blockDim.x + t] = 65535;
		bin_smallest_index[blockDim.x + t] = 0;
	}

	for(unsigned int stride = blockDim.x; stride >= 1; stride >>= 1)
	{
		__syncthreads();
		if(t < stride)
			if(distance_s[t] > distance_s[stride + t])
			{
				bin_smallest_index[t] = start + stride + t;
				distance_s[t] = distance_s[stride + t];
			}
			else
			{
				bin_smallest_index[t] = start + t;	
			}
	}

	if(t == 0)
	{
		distance_d[blockIdx.x] = distance_s[t];
		bin_index_d[blockIdx.x] = bin_smallest_index[t];	
	}


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
