#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "com/here/pb/hdmap/external/v1/lanes/layer-lane-geometry-polyline.pb.h"

using namespace std;

using namespace com::here::pb::hdmap::external::v1::geometry;
using namespace com::here::pb::hdmap::external::v1::lanes;

ofstream ofs;
std::string filename = "laneData.csv";


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

    for(int i = 1; i<64;i++)
    {
	if(i%2 == 0)
	    even.push_back(buffer[i]);
	else
	    odd.push_back(buffer[i]);
    }

    int32_t lon = std::stol(odd, nullptr, 2);
    int32_t lat = std::stol(even, nullptr, 2);
   
    int x = (lon - 1);
    int y = (lat - 1);

    double x_new = x/pow(2,31)*180;
    double y_new = y/pow(2,31)*180;
   
  return std::make_pair(x_new,y_new);
}



void ListLanes(const LaneGeometryPolylineLayerTile& tile) {
  ofs.open(filename);
  int lane_id = 0;
  int64_t tile_center = tile.tile_center_here_3d_coordinate().here_2d_coordinate();
  std::pair<double, double> decoded_tile_center = std::make_pair(0.0,0.0);//decode_morton_2d(tile_center);
  //cout<<"Tile center : "<<tile.tile_center_here_2d_coordinate()<<endl;
  //cout<<"Tile Id: "<<tile.here_tile_id()<<endl;
  
  for (int i = 0; i < tile.lane_group_geometries_size(); i++) {
    
    for (int j = 0; j < tile.lane_group_geometries()[i].lane_geometries_size(); j++) {
      double lat = decoded_tile_center.second; double lon = decoded_tile_center.first;
      LineString3dOffset theGeometry = tile.lane_group_geometries()[i].lane_geometries()[j].lane_path_geometry();
      for (int k = 0; k < theGeometry.here_2d_coordinate_diffs().size(); k++) {
        
	int64_t here_2d_coordinate_diff = theGeometry.here_2d_coordinate_diffs()[k];
	std::pair<double, double> decoded = decode_morton_2d(here_2d_coordinate_diff);
	lat += decoded.second;
	lon += decoded.first;
	cout<<"lat: "<<lat<<"long: "<<lon<<endl;
	ofs << lat << "," << lon<<","<<lane_id<<std::endl;
      }
      lane_id++;
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

  LaneGeometryPolylineLayerTile tile;

  {
    // Read the existing address book.
    fstream input(argv[1], ios::in | ios::binary);
    if (!tile.ParseFromIstream(&input)) {
      cerr << "Failed to parse address book." << endl;
      return -1;
    }
  }

  ListLanes(tile);

  // Optional:  Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
