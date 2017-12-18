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

//////////////////////////////////////////////////////////



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
}

void cal_closest_points_cpu(const column_vector &rt)
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


	PerformTransformationToAllPoints(R, t, &measurement_data, &transformed_data,1);

	int numPointsUpdated = 0;

	
	for(int i = 0; i < transformed_data.size; i++)
	{
		best_index = 0;
		closest_distance = 65535;

		for(int j = 0; j < model_data.size; j++)
		{
			
			distance = sqrt(pow((transformed_data.x_coord.at(i) - model_data.x_coord.at(j)),2) + pow((transformed_data.y_coord.at(i) - model_data.y_coord.at(j)),2) + pow((transformed_data.z_coord.at(i) - model_data.z_coord.at(j)),2));
			
			if(distance < closest_distance)
			{
				closest_distance = distance;
				best_index = j;
			}
		

		}
		measurement_data.index[i] = best_index;
		
		
	}
	

}

double findTotalErrorInCloud_cpu(const column_vector &rt)
{
	
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
	//cout<<"Check 1"<<endl;
	double true_map_error = 0.0;
	for(int i = 0; i < measurement_data.size; i++)
	{
		//cout<<"Value of i "<<i<<endl;
		int j = measurement_data.index.at(i);

		icp_error +=sqrt(pow((transformed_data.x_coord.at(i) - model_data.x_coord.at(j)),2) + pow((transformed_data.y_coord.at(i) - model_data.y_coord.at(j)),2) + pow((transformed_data.z_coord.at(i) - model_data.z_coord.at(j)),2)); 

		//true_map_error +=sqrt(pow((transformed_data.x_coord.at(i) - model_data.x_coord.at(i)),2) + pow((transformed_data.y_coord.at(i) - model_data.y_coord.at(i)),2) + pow((transformed_data.z_coord.at(i) - model_data.z_coord.at(i)),2));

	}


	//cout<<"current error: "<<total_error<<endl;
	//cout<<"final map error: "<<true_map_error<<endl;
	return icp_error;
}

dlib::matrix<double> compute_gold()
{
	
	

	//Rotational function test
	double theta = 0.03;
	double point_x = 0.003;
	double point_y = 0.005;
	double point_z = 0.0;

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
	cpu_starttime = clock();
	for(int i = 0; i<20; i++)
	{
		cout<<"iteration cpu #: "<<i<<endl;
		cal_closest_points_cpu(rt);
		final_error = find_optimal_parameters(0.01, 0.000000001,100000, rt, rt_lower, rt_upper,findTotalErrorInCloud_cpu);		
		cout<<"Rt parameters "<<rt<<endl;
		cout<<"current error: "<<final_error<<endl;
		
	}
	cpu_endtime = clock();
	//cout<<"Error after optimization "<<final_error<<endl;

	column_vector rt_vec(5);
	rt_vec = rt(0), rt(1), rt(2), rt(3), (cpu_endtime - cpu_starttime)/CLOCKS_PER_SEC ;
	
	return rt_vec;
}

