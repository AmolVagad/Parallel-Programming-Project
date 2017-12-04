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


using namespace std;
using namespace dlib;



///////////Globals///////////////////////
double min_x =0;
double min_y = 0;
double min_z = 0;
double range_x = 0;
double range_y = 0;
double range_z = 0;
int bin_size = 16;


// Define Octree 
Octree<std::vector<double>> octree_icp(bin_size); 




//////////////////////////////////////////////////////////
typedef dlib::matrix<double,0,1> column_vector;

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





point_cloud_data measurement_data;
point_cloud_data model_data;
//point_cloud_data transformed_data;
int iterations = 0;
int skips = 5000;


//dlib::matrix<double> point(3,1);
//dlib::matrix<double> Rot_point(3,1);


void cal_closest_points(const column_vector &rt);

double findTotalErrorInCloud(const column_vector &rt);



dlib::matrix<double> PerformRotation(dlib::matrix<double> R,dlib::matrix<double> t, dlib::matrix<double> point)
{
	dlib::matrix<double> point_new(3,1);
	point_new = R*point + t;
	return point_new;
}

void PerformTransformationToAllPoints(dlib::matrix<double> R, dlib::matrix<double> t, point_cloud_data * data, point_cloud_data * transformed_data, int skips)
{
	for(int i  = 0; i < data->size; i++)
	{
		//cout<<"Print x "<<data->x_coord.at(i)<<endl;
		dlib::matrix<double,3,1> point, point_new;
		point = data->x_coord.at(i), data->y_coord.at(i), data->z_coord.at(i);
		//cout<<"Error here above "<<endl;
		point_new = PerformRotation(R, t, point);
		transformed_data->x_coord.push_back(point_new(0));
		transformed_data->y_coord.push_back(point_new(1));
		transformed_data->z_coord.push_back(point_new(2));		 
	}
	transformed_data->size = transformed_data->x_coord.size();
	//cout<<"Transformed data size "<<transformed_data->size<<endl;
}

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
	/*	
	//Checking the octree bin sizes
	for(int i = 0; i < 16; i++)
	{
		for(int j = 0; j < 16; j++)
		{
			for(int k = 0; k < 16; k++)
			{
				cout<<"Bin: "<<i<<"\t"<<j<<"\t"<<k<<"\t"<<octree_icp(i, j, k).size()<<endl;
			}
		}
	}
	*/

	//cout<<"Octree at 0 "<<octree_icp(0,0,0)[0]<<endl;		
		
	

	//Rotational function test
	double theta = 0.03;
	double point_x = 0.003;
	double point_y = 0.005;
	double point_z = 0.0;
	dlib::matrix<double> R(3,3);
	dlib::matrix<double> t(3,1);

	R = cos(theta), -sin(theta), 0,
	    sin(theta), cos(theta), 0,
	    0, 0, 1;

	t = point_x, point_y, point_z;
	
	//R theta = -0.785398
	//x = -0.014142135
	//y = 0
	//z = 0

	// Generate mesasurement datra by rorating the model data
	PerformTransformationToAllPoints(R, t, &model_data, &measurement_data,1);


	//Calling closest point. Currently for testing purpose I am matching i and j
	column_vector rt(4), rt_lower(4), rt_upper(4);

	rt = -theta, -cos(theta)*point_x - sin(theta)*point_y, sin(theta)*point_x - cos(theta)*point_y, point_z;
	cout<<"rt: "<<rt<<endl;

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

	
	


	return 0;
}

void cal_closest_points(const column_vector &rt)
{

	point_cloud_data transformed_data;
	dlib::matrix<double> R(3,3);
	dlib::matrix<double> t(3,1);

	R = cos(rt(0)), -sin(rt(0)), 0,
	    sin(rt(0)), cos(rt(0)), 0,
	    0, 0, 1;

	t = rt(1), rt(2), rt(3);

	float distance = 0.0;
	int best_index ;
	float closest_distance;
	int best_bin_index_x;
	int best_bin_index_y;
	int best_bin_index_z;

	PerformTransformationToAllPoints(R, t, &measurement_data, &transformed_data,1);

	int numPointsUpdated = 0;

	
	for(int i = 0; i < transformed_data.size; i++)
	{
		best_index = 0;
		closest_distance = 65535;
		best_bin_index_x = 0;
		best_bin_index_y = 0;
		best_bin_index_z = 0;
		int index_x_t= 0;
		int index_y_t = 0;
		int index_z_t = 0;
		
		//cout<<"Range x value "<<range_x<<endl;
		index_x_t = floor(((transformed_data.x_coord[i]  - min_x)/range_x)*bin_size);
		index_y_t = floor(((transformed_data.y_coord[i]  - min_y)/range_y)*bin_size);
		index_z_t = floor(((transformed_data.z_coord[i]  - min_z)/range_z)*bin_size);
		
		// Boundary conditon 
		index_x_t = max(min(index_x_t, bin_size - 1),0);
		index_y_t = max(min(index_y_t, bin_size - 1),0);
		index_z_t = max(min(index_z_t, bin_size - 1),0);
		
		
		int p_q_r[3] = {0,0,0};
		int non_empty_bin = 0;
		int max_bin_offset = 2;
//*****************Mandatory search in adjacent bins*******************************************************************		
		for(p_q_r[0] = -1; p_q_r[0] < max_bin_offset; p_q_r[0]++)
		{
			for(p_q_r[1] = -1; p_q_r[1] < max_bin_offset; p_q_r[1]++)
			{
				for(p_q_r[2] = -1; p_q_r[2] < max_bin_offset; p_q_r[2]++)
				{
					bool p_flag = ((p_q_r[0] + index_x_t) >= 0) && ((p_q_r[0] + index_x_t) <= bin_size - 1);
					bool q_flag = ((p_q_r[1] + index_y_t) >= 0) && ((p_q_r[1] + index_y_t) <= bin_size - 1);
					bool r_flag = ((p_q_r[2] + index_z_t) >= 0) && ((p_q_r[2] + index_z_t) <= bin_size - 1);
					if(p_flag && q_flag && r_flag)
					{
						if(octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2]).size()/3 > 0)
							non_empty_bin++;
						for(int l = 0; l < octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2]).size()/3;l++)
						{
					
							distance = sqrt(pow((transformed_data.x_coord[i] - octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2])[3*l]),2) + pow((transformed_data.y_coord[i] - octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2])[3*l+1]),2) + pow((transformed_data.z_coord[i] - octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2])[3*l + 2]),2));
	
							if(distance < closest_distance)
							{
								closest_distance = distance;
								best_index = l;
								best_bin_index_x = index_x_t + p_q_r[0];
								best_bin_index_y = index_y_t + p_q_r[1];
								best_bin_index_z = index_z_t + p_q_r[2];	
							}
						}
					}
				}
			}
			
		}
//*************************************************************************************************************
		
		while(non_empty_bin == 0)
		{
			//cout<<"Entered this if"<<endl;
			max_bin_offset = max_bin_offset + 1;
			
			for(int a = 0; a < 6; a++)
			{
				p_q_r[a%3] = pow(-1,a)*(max_bin_offset - 1);
				for(int b = -max_bin_offset + 1; b < max_bin_offset; b++)
				{
					p_q_r[(a+1)%3] = b;
					for(int c = -max_bin_offset + 1; c < max_bin_offset; c++)
					{
						p_q_r[(a+2)%3] = c;
						bool p_flag = ((p_q_r[0] + index_x_t) >= 0) && ((p_q_r[0] + index_x_t) <= bin_size - 1);
						bool q_flag = ((p_q_r[1] + index_y_t) >= 0) && ((p_q_r[1] + index_y_t) <= bin_size - 1);
						bool r_flag = ((p_q_r[2] + index_z_t) >= 0) && ((p_q_r[2] + index_z_t) <= bin_size - 1);
						if(p_flag && q_flag && r_flag)
						{
							if(octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2]).size()/3 > 0)
								non_empty_bin++;
							for(int l = 0; l < octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2]).size()/3;l++)
							{
					
								distance = sqrt(pow((transformed_data.x_coord[i] - octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2])[3*l]),2) + pow((transformed_data.y_coord[i] - octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2])[3*l+1]),2) + pow((transformed_data.z_coord[i] - octree_icp(index_x_t + p_q_r[0], index_y_t + p_q_r[1], index_z_t + p_q_r[2])[3*l + 2]),2));
	
								if(distance < closest_distance)
								{
									closest_distance = distance;
									best_index = l;
									best_bin_index_x = index_x_t + p_q_r[0];
									best_bin_index_y = index_y_t + p_q_r[1];
									best_bin_index_z = index_z_t + p_q_r[2];	
								}
							}
						}
					}
				}
			}		





		}	


		if(non_empty_bin > 0)
		{						
			measurement_data.index[i] = best_index;
			measurement_data.bin_index_x[i] = best_bin_index_x;
			measurement_data.bin_index_y[i] = best_bin_index_y;
			measurement_data.bin_index_z[i] = best_bin_index_z;
		}						
		

		
	}
	std::cout<<"numPoints Updated "<<numPointsUpdated<<endl;
	//cout<<"Octree at 0 "<<octree_icp(0,0,0)[0]<<endl;
}



double findTotalErrorInCloud(const column_vector &rt)
{
	iterations++;
	double icp_error = 0.0;
	point_cloud_data transformed_data;
	
	dlib::matrix<double> R(3,3);
	dlib::matrix<double> t(3,1);

	R = cos(rt(0)), -sin(rt(0)), 0,
	    sin(rt(0)), cos(rt(0)), 0,
	    0, 0, 1;

	t = rt(1), rt(2), rt(3);
	//cout<<"Check measurement data element "<<measurement_data.x_coord.at(0)<<endl;
	PerformTransformationToAllPoints(R, t, &measurement_data, &transformed_data,1);
	//cout<<"Measurement data size "<<measurement_data.size<<endl;
	double true_map_error = 0.0;

	for(int i = 0; i < measurement_data.size; i++)
	{
		//cout<<"Check 2"<<endl;
		int j = measurement_data.index[i];
		int x_Idx = measurement_data.bin_index_x[i];
		int y_Idx = measurement_data.bin_index_y[i];
		int z_Idx = measurement_data.bin_index_z[i];
		//cout<<"x IDX "<<x_Idx<<endl;
		//cout<<"y IDX "<<y_Idx<<endl;
		//cout<<"z IDX "<<z_Idx<<endl;
		//cout<<"Value of j  "<<j<<endl;


		//cout<<"Octree check "<<octree_icp(x_Idx, y_Idx, z_Idx)[3*j]<<endl;
		icp_error +=sqrt(pow((transformed_data.x_coord[i] - octree_icp(x_Idx, y_Idx, z_Idx)[3*j]),2) + pow((transformed_data.y_coord[i] - octree_icp(x_Idx, y_Idx, z_Idx)[3*j + 1]),2) + pow((transformed_data.z_coord[i] - octree_icp(x_Idx, y_Idx, z_Idx)[3*j + 2]),2)); 

		//i//cp_error += transformed_data.x_coord[i] + transformed_data.y_coord[i] + transformed_data.z_coord[i];
		//cout<<"Check 3"<<endl;
		//true_map_error +=sqrt(pow((transformed_data.x_coord.at(i) - model_data.x_coord.at(i)),2) + pow((transformed_data.y_coord.at(i) - model_data.y_coord.at(i)),2) + pow((transformed_data.z_coord.at(i) - model_data.z_coord.at(i)),2));
		//cout<<"Value of i "<<i<<endl;
	}
	

	return icp_error;
}


	


