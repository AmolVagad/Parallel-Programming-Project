#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "com/here/pb/hdmap/external/alpha/landmark/layer-landmark-obstacles.pb.h"
using namespace std;
using namespace com::here::pb::hdmap::external::alpha::landmark;
using namespace com::here::pb::hdmap::external::v1::geometry;


ofstream ofs;
ofstream ofs2;
std::string filename = "exampleOutput2.csv";
std::string filename_apc = "icp_model_landmarks_intersection.csv";

double car_x = 0.00772 + 0.00005;
double car_y = -0.0033 + 0.00005;
double car_z = 0.0;
double point_x, point_y, point_z;
double dist, theta, phi;


int twos_complement_to_decimal(string input)
{
	
	int n = input.length();
	char a = input[0];
	//cout<<"MSB "<<stol(&a,nullptr,2)<<endl;
	long int sum = -stol(&a,nullptr,2)*pow(2,n-1);
	for(int i = 1; i < n; i++)
	{
	//cout<<"Sum "<<sum<<endl;
	a = input[i];
	sum += stol(&a,nullptr,2)*pow(2, n-1-i);
	}
	//cout<<"Sum value "<<sum<<endl;
	//sum = sum - stoi(&input[n-1],nullptr,2)*pow(2,n);
	return sum;
}

std::pair<double, double> decode_morton_2d(int64_t m)
{
    //Converting number to 64 bits
    string buffer;
    string odd, even;
    int i=0;
    while (i<64)
    {
        buffer.push_back('0' + (m & 1));
	m>>=1;
	i++;
    } 
    std::reverse(buffer.begin(), buffer.end());
    //cout<<"First bit "<<buffer[0]<<endl;
    for(int i = 1; i<64;i++)
    {
	if(i%2 == 0)
	    even.push_back(buffer[i]);//latitude
	else
	    odd.push_back(buffer[i]);//longitude
    }
    //cout<<"Odd is "<<odd<<endl;
    //cout<<"Even is "<<even<<endl;
    long int x = twos_complement_to_decimal(odd);
    long int y = twos_complement_to_decimal(even);
    double x_new = x/pow(2,31)*180;
    double y_new = y/pow(2,31)*180;
   
  return std::make_pair(x_new,y_new);
}


void findPointInCloud(double lat, double lon, double height){

  dist = sqrt(pow((lat),2) + pow((lon),2) + pow(height,2));
  
  theta = atan((lon)/lat);
  if(lat < 0)
	theta += M_PI;

  phi = acos((height)/dist);

}

void findCoordinateInCloud()
{

  theta += 0.06;

  point_x =  dist*cos(theta)*sin(phi);
  point_y =  dist*sin(theta)*sin(phi);
  point_z =  dist*cos(phi);

}



void ListObstacles(const LandmarkObstaclesLayerTile& tile) {
  ofs.open(filename);
  ofs2.open(filename_apc);
 
  int obstacle_id = 0;
  int64_t current_2d_encoding_int = tile.tile_center_here_2d_coordinate();
  int64_t tile_center = tile.tile_center_here_2d_coordinate();
  std::pair<double, double> decoded_tile_center = std::make_pair(0.0,0.0);//decode_morton_2d(tile_center);
  cout<<"Tile center : "<<tile.tile_center_here_2d_coordinate()<<endl;
  cout<<"Tile Id: "<<tile.here_tile_id()<<endl;
  for (int i = 0; i < tile.obstacles_for_lane_groups_size(); i++) {

    for (int j = 0; j < tile.obstacles_for_lane_groups()[i].obstacles_size(); j++) {
      double obstacle_height = tile.obstacles_for_lane_groups()[i].obstacles()[j].height_above_road_cm() /10000000.0;
      LineString2dOffset theGeometry = tile.obstacles_for_lane_groups()[i].obstacles()[j].geometry();
      double lat = decoded_tile_center.second; double lon = decoded_tile_center.first;
      current_2d_encoding_int = 0;
      for (int k = 0; k < theGeometry.here_2d_coordinate_diffs().size(); k++) {
	
	int64_t here_2d_coordinate_diff = theGeometry.here_2d_coordinate_diffs()[k];
	int64_t theXor = (here_2d_coordinate_diff ^ current_2d_encoding_int);
	std::pair<double, double> decoded = decode_morton_2d(theXor);
	//saving for next iteration
	current_2d_encoding_int = theXor;
	lat = decoded.second;
	lon = decoded.first;


	if(lon>(0.02197265625/2))
		lon = -(0.02197265625) +lon;
	if(lat>(0.02197265625/2))
		lat = -(0.02197265625) +lat;
	
	cout<<"Lat: "<<fixed<<lat<<"   Long: "<<fixed<<lon<<" Obstacle height"<<fixed<<obstacle_height<<endl;
	//Writing to CSV
	 ofs.precision(12);
	ofs << fixed<<lat << "," << fixed<<lon << "," << fixed<<obstacle_height<<","<<obstacle_id<<std::endl;

	//Testing functions
	//point_x = lon;
	//point_y = lat;
	//point_z = obstacle_height;
	findPointInCloud(lat,lon,obstacle_height);
	findCoordinateInCloud();
	
	//cout<<"Lat: "<<point_y<<"Long: "<<point_x<<" Height "<<point_z<<endl;

	//something = newfunction(lat,lon,height);
	ofs2<<point_x<<","<<point_y<<","<<fixed<<point_z<<endl;
	//cout<<"here 2d coordinate int: "<<theGeometry.here_2d_coordinate_diffs()[k]<<endl;
	//cout<<"here z level index int: "<<theGeometry.z_level_indexes()[0]<<endl;
      }
      obstacle_id++;
    }
	
  }
  ofs.close(); 
} 

int main(int argc, char* argv[]) {
  // Verify that the version of the library that we linked against is
  // compatible with the version of the headers we compiled against.
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (argc != 2) {
    cerr << "Usage:  " << argv[0] << " LANDMARK_FILE" << endl;
    return -1;
  }

  LandmarkObstaclesLayerTile tile;

  {
    // Read the existing address book.
    fstream input(argv[1], ios::in | ios::binary);
    if (!tile.ParseFromIstream(&input)) {
      cerr << "Failed to parse address book." << endl;
      return -1;
    }
  }

  ListObstacles(tile);

  // Optional:  Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
