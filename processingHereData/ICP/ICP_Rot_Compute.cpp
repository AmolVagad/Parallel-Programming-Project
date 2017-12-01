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



using namespace std;
using namespace dlib;



///////////Globals///////////////////////
double min_x =0;
double min_y = 0;
double min_z = 0;
double range_x = 0;
double range_y = 0;
double range_z = 0;

// Define Octree 
Octree<std::vector<double>> octree_icp(16); 




//////////////////////////////////////////////////////////
typedef dlib::matrix<double,0,1> column_vector;

struct point_cloud_data{

	std::vector <double> x_coord;
	std::vector <double> y_coord;
	std::vector <double> z_coord;
	std::vector <int> index;
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
	for(int i  = 0; i < data->size; i+= skips)
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
	
	}
	
	// Calculating the min and max values of x,y,z
	double max_x =  *max_element(model_data.x_coord.begin(),model_data.x_coord.end()) ;
	double max_y =  *max_element(model_data.y_coord.begin(),model_data.y_coord.end()) ;
	double max_z =  *max_element(model_data.z_coord.begin(),model_data.z_coord.end()) ;
	double min_x =  *min_element(model_data.x_coord.begin(),model_data.x_coord.end()) ;
	double min_y =  *min_element(model_data.y_coord.begin(),model_data.y_coord.end()) ;
	double min_z =  *min_element(model_data.z_coord.begin(),model_data.z_coord.end()) ;
	
	// Calculating the range
	range_x = max_x - min_x; 
	range_y = max_y - min_y; 
	range_z = max_z - min_z;
	
	model_data.size = model_data.x_coord.size();
	
	
	// Storing the data into Octrees 
	
	for(int i= 0; i < model_data.size; i++)
	{
		int index_x = 0;
		int index_y = 0;
		int index_z = 0;
		
		index_x = floor(((model_data.x_coord.at(i)  - min_x)/range_x)*16);
		index_y= floor(((model_data.y_coord.at(i)  - min_y)/range_y)*16);
		index_z = floor(((model_data.z_coord.at(i)  - min_z)/range_z)*16);
		
		// Boundary conditon 
		index_x = min(index_x, 15);
		index_y = min(index_y, 15);
		index_z = min(index_z, 15);
		
		
		octree_icp(index_x, index_y, index_z).push_back(model_data.x_coord.at(i));
		octree_icp(index_x, index_y, index_z).push_back(model_data.y_coord.at(i));
		octree_icp(index_x, index_y, index_z).push_back(model_data.z_coord.at(i));
		
	}
		
			
		
	

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

		//final_error = find_optimal_parameters(0.01, 0.000000001,100000, rt, rt_lower, rt_upper,findTotalErrorInCloud);
		//cout<<"Rt parameters "<<rt<<endl;
		//cout<<"current error: "<<final_error<<endl;
		
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


	PerformTransformationToAllPoints(R, t, &measurement_data, &transformed_data,1);

	int numPointsUpdated = 0;

	
	for(int i = 0; i < transformed_data.size; i++)
	{
		best_index = 0;
		closest_distance = 65535;
		
		int index_x_t= 0;
		int index_y_t = 0;
		int index_z_t = 0;
		
		
		index_x_t = floor(((transformed_data.x_coord[i]  - min_x)/range_x)*16);
		index_y_t = floor(((transformed_data.y_coord[i]  - min_y)/range_y)*16);
		index_z_t = floor(((transformed_data.z_coord[i]  - min_z)/range_z)*16);
		
		// Boundary conditon 
		index_x_t = max(min(index_x_t, 15),0);
		index_y_t = max(min(index_y_t, 15),0);
		index_z_t = max(min(index_z_t, 15),0);
		
		
		for(int p = -1; p + index_x_t >0 && p + index_x_t < 15 && p < 2; p++)
		{
					for(int q = -1; q + index_y_t >0 && q + index_y_t < 15 && q < 2; q++)
					{
								for(int r = -1; r + index_z_t >0 && r + index_z_t < 15 && r < 2; r++)
								{
									for(int l = 0; l < octree_icp(index_x_t + p, index_y_t + q, index_z_t + r).size()/3;l++)
									{
										
										distance = sqrt(pow((transformed_data.x_coord[i] - octree_icp(index_x_t + p, index_y_t + q, index_z_t + r)[3*l]),2) + pow((transformed_data.y_coord[i] - octree_icp(index_x_t + p, index_y_t + q, index_z_t + r)[3*l+1]),2) + pow((transformed_data.z_coord[i] - octree_icp(index_x_t + p, index_y_t + q, index_z_t + r)[3*l + 2]),2));
								
										if(distance < closest_distance)
										{
											closest_distance = distance;
											best_index = l;
											
											
										}
									}
								}
					}
		}
								
								
		
		

		
	}
	std::cout<<"numPoints Updated "<<numPointsUpdated<<endl;

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
	//cout<<"Check 1"<<endl;
	double true_map_error = 0.0;
	for(int i = 0; i < measurement_data.size; i++)
	{
		//cout<<"Check 2"<<endl;
		int j = measurement_data.index.at(i);

		icp_error +=sqrt(pow((transformed_data.x_coord.at(i) - model_data.x_coord.at(j)),2) + pow((transformed_data.y_coord.at(i) - model_data.y_coord.at(j)),2) + pow((transformed_data.z_coord.at(i) - model_data.z_coord.at(j)),2)); 

		true_map_error +=sqrt(pow((transformed_data.x_coord.at(i) - model_data.x_coord.at(i)),2) + pow((transformed_data.y_coord.at(i) - model_data.y_coord.at(i)),2) + pow((transformed_data.z_coord.at(i) - model_data.z_coord.at(i)),2));

	}

	double sign_object_localizer = 0.0;

	sign_object_localizer +=  sqrt(pow((transformed_data.x_coord.at(0) - model_data.x_coord.at(0)),2) + pow((transformed_data.y_coord.at(0) - model_data.y_coord.at(0)),2) + pow((transformed_data.z_coord.at(0) - model_data.z_coord.at(0)),2));

	sign_object_localizer += sqrt(pow((transformed_data.x_coord.at(1*skips) - model_data.x_coord.at(1*skips)),2) + pow((transformed_data.y_coord.at(1*skips) - model_data.y_coord.at(1*skips)),2) + pow((transformed_data.z_coord.at(1*skips) - model_data.z_coord.at(1*skips)),2));

	//sign_object_localizer += sqrt(pow((transformed_data.x_coord.at(2*skips) - model_data.x_coord.at(2*skips)),2) + pow((transformed_data.y_coord.at(2*skips) - model_data.y_coord.at(2*skips)),2) + pow((transformed_data.z_coord.at(2*skips) - model_data.z_coord.at(2*skips)),2));

	//sign_object_localizer += sqrt(pow((transformed_data.x_coord.at(3*skips) - model_data.x_coord.at(3*skips)),2) + pow((transformed_data.y_coord.at(3*skips) - model_data.y_coord.at(3*skips)),2) + pow((transformed_data.z_coord.at(3*skips) - model_data.z_coord.at(3*skips)),2));

	//sign_object_localizer += sqrt(pow((transformed_data.x_coord.at(4*skips) - model_data.x_coord.at(4*skips)),2) + pow((transformed_data.y_coord.at(4*skips) - model_data.y_coord.at(4*skips)),2) + pow((transformed_data.z_coord.at(4*skips) - model_data.z_coord.at(4*skips)),2));

	double total_error = (sign_object_localizer*measurement_data.size*0.5) + 1*icp_error;

	//cout<<"current error: "<<total_error<<endl;
	//cout<<"final map error: "<<true_map_error<<endl;
	return total_error;
}


	


