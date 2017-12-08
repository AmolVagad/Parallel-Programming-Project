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

extern "C"

using namespace std;


__constant__ double R_constant[9];

#include "icp_kernel.cu"


// Creating variables to store the measurement and model data

	point_cloud_data measurement_data;
	point_cloud_data model_data;

	int bin_size = 4;

	double min_x =0;
	double min_y = 0;
	double min_z = 0;
	double range_x = 0;
	double range_y = 0;
	double range_z = 0;



	

///// For initial testing purposes carrying out rotation and translation operation on cuda//////////////////



// Define Octree 
Octree<std::vector<double>> octree_icp(bin_size); 

// Define the column vector 
typedef dlib::matrix<double,0,1> column_vector;





// Function to carry out Rotation of given point on the device 

double* PerformRotationOnDevice(const Matrix R_h, const Matrix t_h, const Matrix Point_h, Matrix Rotated_Point_h)
{
	
	int size_R = R_h.width*R_h.height*sizeof(double);
	int size_T = t_h.width*t_h.height*sizeof(double);
	int size_Point = Point_h.width*Point_h.height*sizeof(double);
	// Declare the device variables 
	
	Matrix  t_d,Point_d, Rotated_Point_d;
	
	// Allocate memory on the device
	
	
	cudaMalloc((void **)&t_d.elements,size_T);
	cudaMalloc((void **)&Point_d.elements,size_Point);
	
	// Copy from host to device 
	
	
	cudaMemcpy(t_d.elements,t_h.elements,size_T, cudaMemcpyHostToDevice);
	cudaMemcpy(Point_d.elements,Point_h.elements,size_Point, cudaMemcpyHostToDevice);
	
	
	// Allocate device memory for result 
	cudaMalloc((void **)&Rotated_Point_d.elements,size_T);
	
	
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

	cudaMemcpyToSymbol(R_constant,R_h.elements,3 * 3*sizeof(double));

    // Launch the device computation threads!

     PerformRotationKernel<<<dimGrid,dimBlock>>>(t_d, Point_d, Rotated_Point_d);
		
	// Transfer Rotated Point from device to host
     cudaMemcpy(Rotated_Point_h.elements, Rotated_Point_d.elements, size_T, cudaMemcpyDeviceToHost);
       // Free device memory for all
     cudaFree(t_d.elements); cudaFree (Point_d.elements);cudaFree (Rotated_Point_d.elements);
     
     return Rotated_Point_h.elements;	
	
}






// Function that calls transformation function and stores the transformed values 

void PerformTransformationToAllPoints(const Matrix R,const Matrix t, point_cloud_data * data, point_cloud_data * transformed_data, int skips)
{
	Matrix point, rotated_point;
	rotated_point.height =3;
	rotated_point.width = 1;  
	point.elements = (double*)malloc(3*data->size*sizeof(double));
	rotated_point.elements = (double*)malloc(rotated_point.width*rotated_point.height*sizeof(double));
	for(int i  = 0; i < data->size; i++)
	{
		
		
		point.elements[i+0] = data->x_coord.at(i);
		point.elements[i+1] = data->y_coord.at(i);
		point.elements[i+2] = data->z_coord.at(i);
		
		transformed_data->x_coord.push_back(rotated_point.elements[0]);
		transformed_data->y_coord.push_back(rotated_point.elements[1]);
		transformed_data->z_coord.push_back(rotated_point.elements[2]);		 
	}
	transformed_data->size = transformed_data->x_coord.size();
	rotated_point.elements = PerformRotationOnDevice(R, t, point, rotated_point);
	
}





// The main function 

int main()
{
	
	
	


	ifstream infile1;
  	infile1.open ("icp_model.csv");
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
	// Allocating memory to the matrices 

	R.elements = (double*)malloc(R.width*R.height*sizeof(double));
	t.elements = (double*)malloc(t.width*t.height*sizeof(double));

	R.elements[0] = cos(theta);R.elements[1]= -sin(theta); R.elements[2]= 0;
	R.elements[3] =sin(theta);  R.elements[4]=cos(theta); R.elements[5]= 0;
	R.elements[6] = 0; R.elements[7]= 0; R.elements[8]= 1;
	
	t.elements[0] = point_x;
	t.elements[1]= point_y;
	t.elements[2] = point_z;
	
	
	// Generate mesasurement datra by rorating the model data
	PerformTransformationToAllPoints(R, t, &model_data, &measurement_data,1);

	
/*

	//Calling closest point.
	column_vector  rt(4), rt_lower(4), rt_upper(4);

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




	
