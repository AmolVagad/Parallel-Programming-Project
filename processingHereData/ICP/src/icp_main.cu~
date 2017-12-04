// HOST Code to compute ICP for localization 



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
#include "icp_kernel.cu"
extern "C"
using namespace std;
///// For initial testing purposes carrying out rotation and translation operation on cuda//////////////////






// Point cloud data structure 


//////////////////////////////////////////////////////////
//Matrix column_vector;
//column_vector.height = 3;
//column_vector.width = 1;

struct point_cloud_data{

	std::vector <double> x_coord;
	std::vector <double> y_coord;
	std::vector <double> z_coord;
	std::vector <int> index;
	std::vector <int> bin_index_x;
	std::vector <int> bin_index_y;
	std::vector <int> bin_index_z;
	int size ;	
};

// Creating variables to store the measurement and model data

point_cloud_data measurement_data;
point_cloud_data model_data;

int bin_size = 4;







// Function to carry out Rotation of given point on the device 

double* PerformRotationOnDevice(const Matrix R_h,const Matrix t_h, const Matrix Point_h, Matrix Rotated_Point_h)
{
	int n = 3;
	int size = n*sizeof(double);

	// Declare the device variables 
	
	Matrix R_d, t_d,Point_d, Rotated_Point_d;
	
	// Allocate memory on the device
	
	cudaMalloc((void **)&R_d,size);
	cudaMalloc((void **)&t_d,size);
	cudaMalloc((void **)&Point_d,size);
	
	// Copy from host to device 
	
	cudaMemcpy(R_d,R_h,size, cudaMemcpyHostToDevice);
	cudaMemcpy(t_d,t_h,size, cudaMemcpyHostToDevice);
	cudaMemcpy(Point_d,Point_h,size, cudaMemcpyHostToDevice);
	
	
	// Allocate device memory for result 
	cudaMalloc((void **)&Rotated_Point_d,size);
	
	
	// Kernel Call 
	 
 // Setup the execution configuration
    	int blocks_w = R_h.width/TILE_WIDTH ;
    	int blocks_h = Point_h.height /TILE_WIDTH;
    
   
	if(t_h.width % TILE_WIDTH)
		blocks_w ++;

	if(R_h.height % TILE_WIDTH)
		blocks_h ++;
		         
	dim3 dimGrid(blocks_w, blocks_h, 1);

	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH,1);

    // Launch the device computation threads!

     PerformRotationKernel<<<dimGrid,dimBlock>>>(R_d, t_d, Point_d, Rotated_Point_d);
		
	// Transfer Rotated Point from device to host
     cudaMemcpy(Rotated_Point_h, Rotated_Point_d, size, cudaMemcpyDeviceToHost);
       // Free device memory for all
     cudaFree(R_d); cudaFree(t_d); cudaFree (Point_d);cudaFree (Rotated_Point_d);
     
     return Rotated_Point_h;	
	
}






// Function that calls transformation function and stores the transformed values 

void PerformTransformationToAllPoints(Matrix R,Matrix t, point_cloud_data * data, point_cloud_data * transformed_data, int skips)
{
	for(int i  = 0; i < data->size; i++)
	{
		Matrix point, rotated_point;
		point.height = rotated_point.height =3;
		point.width = rotated_point.width = 1;
		point = {data->x_coord.at(i), data->y_coord.at(i), data->z_coord.at(i)};
		rotated_point.elements =PerformRotationOnDevice(R, t, point, rotated_point);
		transformed_data->x_coord.push_back(rotated_point.elements[0]);
		transformed_data->y_coord.push_back(rotated_point.elements[1]);
		transformed_data->z_coord.push_back(rotated_point.elements[2]);		 
	}
	transformed_data->size = transformed_data->x_coord.size();
	
}





// The main function 

int main()
{
	

	
	
	ifstream infile1;
  	infile1.open ("icp_model.csv");
	//ifstream infile2;
  	//infile2.open ("icp_sensor_scan.csv");
	char* pEnd;
	string x,y,z;


	// Reading data from the model map data csv file 

	
	 while(!infile1.eof()){
		getline(infile1,x, ',');
		getline(infile1,y, ',');
		getline(infile1,z);
		//getline(infile,index);
		model_data.x_coord.push_back(strtod(x.c_str(),&pEnd));
		model_data.y_coord.push_back(strtod(y.c_str(),&pEnd));
		model_data.z_coord.push_back(strtod(z.c_str(),&pEnd));
		measurement_data.index.push_back(-1);
		measurement_data.bin_index_x.push_back(-1);
		measurement_data.bin_index_y.push_back(-1);
		measurement_data.bin_index_z.push_back(-1);
	
	}
	
	//Remove the last elements
	model_data.x_coord.pop_back();
	model_data.y_coord.pop_back();
	model_data.z_coord.pop_back();
	model_data.size = model_data.size - 1;

	// Calculating the min and max values of x,y,z
	double max_x =  *max_element(model_data.x_coord.begin(),model_data.x_coord.end()) ;
	double max_y =  *max_element(model_data.y_coord.begin(),model_data.y_coord.end()) ;
	double max_z =  *max_element(model_data.z_coord.begin(),model_data.z_coord.end()) ;
	min_x =  *min_element(model_data.x_coord.begin(),model_data.x_coord.end()) ;
	min_y =  *min_element(model_data.y_coord.begin(),model_data.y_coord.end()) ;
	min_z =  *min_element(model_data.z_coord.begin(),model_data.z_coord.end()) ;
	

	//cout<<"Min x value "<<min_x<<endl;
	// Calculating the range
	range_x = max_x - min_x; 
	range_y = max_y - min_y; 
	range_z = max_z - min_z;
		
	//cout<<"Range x value "<<range_x<<endl;

	model_data.size = model_data.x_coord.size();
	//cout<<"model data value "<<model_data.size<<endl;
	
	
	// Storing the data into Octrees 
	for(int i= 0; i < model_data.size; i++)
	{
		int index_x = 0;
		int index_y = 0;
		int index_z = 0;
		
		index_x = floor(((model_data.x_coord.at(i)  - min_x)/range_x)*bin_size);
		index_y= floor(((model_data.y_coord.at(i)  - min_y)/range_y)*bin_size);
		index_z = floor(((model_data.z_coord.at(i)  - min_z)/range_z)*bin_size);
		
		// Boundary conditon 
		index_x = min(index_x, bin_size - 1);
		index_y = min(index_y, bin_size - 1);
		index_z = min(index_z, bin_size - 1);
		
		
		octree_icp(index_x, index_y, index_z).push_back(model_data.x_coord.at(i));
		octree_icp(index_x, index_y, index_z).push_back(model_data.y_coord.at(i));
		octree_icp(index_x, index_y, index_z).push_back(model_data.z_coord.at(i));
		
	}
	
	//cout<<"Octree at 0 "<<octree_icp(0,0,0)[0]<<endl;		
		
	

	//Rotational function test
	double theta = 0.03;
	double point_x = 0.003;
	double point_y = 0.005;
	double point_z = 0.0;
	Matrix R;
	Matrix t;
	R.width = 3;
	R.height = 3;
	t.width = 1;
	t.height = 3;

	R.elements = {cos(theta), -sin(theta), 0,
	    sin(theta), cos(theta), 0,
	    0, 0, 1};

	t.elements[0] = point_x;
	t.elements[1]= point_y;
	t.elements[2] = point_z;
	

	// Generate mesasurement datra by rorating the model data
	PerformTransformationToAllPoints(R, t, &model_data, &measurement_data,1);

/*
	//Calling closest point.
	column_vector ={ rt(4), rt_lower(4), rt_upper(4)};

	rt = -theta, -cos(theta)*point_x - sin(theta)*point_y, sin(theta)*point_x - cos(theta)*point_y, point_z;
	std::cout<<"rt: "<<rt<<endl;

	rt = -0.0,0.0,0.0,0.0;
	rt_lower = -1.0, -1.0,-1.0,-1.0;
	rt_upper = 1.0, 1.0, 1.0, 1.0;
	

	double final_error = 0;
	// time measurement variables 


	double cpu_starttime , cpu_endtime;
	for(int i = 0; i<20; i++)
	{
		cout<<"iteration #: "<<i<<endl;
		cpu_starttime = clock();
		cal_closest_points(rt);
		cpu_endtime = clock();
		cout<<"The time taken for calculation = "<<((cpu_endtime - cpu_starttime)/CLOCKS_PER_SEC)<<endl;

		final_error = find_optimal_parameters(0.01, 0.000000001,100000, rt, rt_lower, rt_upper,findTotalErrorInCloud);
		cout<<"Rt parameters "<<rt<<endl;
		cout<<"current error: "<<final_error<<endl;
		
	}
	//cout<<"Error after optimization "<<final_error<<endl;
*/
	
	


	return 0;
}




	
