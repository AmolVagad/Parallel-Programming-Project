// File contains all the global variables accessed by various files and functions 

#include "dlib/optimization/optimization.h"
#include "dlib/optimization/find_optimal_parameters_abstract.h"
#include "dlib/optimization/optimization_bobyqa.h"
#include "dlib/optimization/find_optimal_parameters.h" 
#include "octree_code/octree.h"

extern double min_x ;
extern double min_y ;
extern double min_z ;
extern double range_x ;
extern double range_y ;
extern double range_z ;

// Define Octree 
extern Octree<std::vector<double>> octree_icp; 

typedef dlib::matrix<double,0,1> column_vector;

#ifndef _ICP_H_
#define _ICP_H_

// Matrix Structure declaration
typedef struct {
    unsigned int width;
    unsigned int height;
    unsigned int pitch;
    double* elements;
} Matrix;


// Vector structure declaration

typedef struct {

    unsigned int size;
    unsigned int pitch;
    int* elements;
} Vector;


// Point cloud data structure 

struct point_cloud_data{

	std::vector <double> x_coord;
	std::vector <double> y_coord;
	std::vector <double> z_coord;
	std::vector <int> index;
	int size ;	
};

extern point_cloud_data measurement_data;
extern point_cloud_data model_data;

extern dlib::matrix<double> PerformRotation(dlib::matrix<double> R,dlib::matrix<double> t, dlib::matrix<double> point);

extern void PerformTransformationToAllPoints(dlib::matrix<double> R, dlib::matrix<double> t, point_cloud_data * data, point_cloud_data * transformed_data, int skips);

extern void cal_closest_points_cpu(const column_vector &rt);

extern double findTotalErrorInCloud_cpu(const column_vector &rt);

extern dlib::matrix<double> compute_gold();





#endif
 
