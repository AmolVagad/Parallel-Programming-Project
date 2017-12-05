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


#ifndef _ICP_H_
#define _ICP_H_

// Matrix Structure declaration
typedef struct {
    unsigned int width;
    unsigned int height;
    unsigned int pitch;
    double* elements;
} Matrix;





// Point cloud data structure 

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




#endif
 
