#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "com/here/pb/hdmap/external/alpha/landmark/layer-landmark-poles.pb.h"
using namespace std;
using namespace com::here::pb::hdmap::external::alpha::landmark;
using namespace com::here::pb::hdmap::external::v1::geometry;


ofstream ofs;
ofstream ofs2;
std::string filename = "poles.csv";


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
    //cout<<"Buffer "<<buffer<<endl;
    for(int i = 1; i<64;i++)
    {
	if(i%2 == 0)
	    even.push_back(buffer[i]);
	else
	    odd.push_back(buffer[i]);
    }

    char temp = even[0];
    char temp_odd = odd[0];
    cout<<"even sign"<<even[0]<<endl;
    cout<<"odd sign"<<odd[0]<<endl;
    if(temp == '1')
{  
    even = "";
    buffer = "";
    m-=2;
    i = 0;
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
    }
    //cout<<"even "<<even<<endl;

    for(int i = 0; i<31; i++)
    {
	if(even[i] == '0')
	    even[i] = '1';
	else
	    even[i] = '0';
    }
} 
	if(temp_odd == '1')
	{  

	    odd = "";
	    buffer = "";
	    m-=1;
	    i = 0;
	    while (i<64)
	    {
		buffer.push_back('0' + (m & 1));
		m>>=1;
		i++;
	    } 
	    std::reverse(buffer.begin(), buffer.end());

	    for(int i = 1; i<64;i++)
	    {
		if(i%2 != 0)
		    odd.push_back(buffer[i]);
	    }
		//cout<<"odd "<<odd<<endl;

	    for(int i = 0; i<32; i++)
	    {
		if(odd[i] == '0')
		    odd[i] = '1';
		else
		    odd[i] = '0';
	    }
	} 

    int32_t lon = std::stol(odd/*.substr(1,32)*/, nullptr, 2);
    int32_t lat = std::stoi(even/*.substr(1,31)*/, nullptr, 2);

    if(temp == '1')
	lat = -lat;
    if(temp_odd == '1')
	lon = -lon;
    
    //cout<<"lat: "<<lat<<endl;
    //cout<<"lon: "<<lon<<endl;
   

    double x_new = lon/pow(2,31)*180;
    double y_new = lat/pow(2,31)*180;
   
  return std::make_pair(x_new,y_new);
}




void ListPoles(const LandmarkPolesLayerTile& tile) {
  ofs.open(filename);

  cout.precision(10);
  int pole_id = 0;

  //int64_t tile_center = tile.tile_center_here_2d_coordinate();
  std::pair<double, double> decoded_tile_center = std::make_pair(0.0,0.0);//decode_morton_2d(tile_center);
  //cout<<"Tile center : "<<tile.tile_center_here_2d_coordinate()<<endl;
  cout<<"Tile Id: "<<tile.here_tile_id()<<endl;
  for (int i = 0; i < tile.poles_for_lane_groups_size(); i++) {

    for (int j = 0; j < tile.poles_for_lane_groups()[i].poles_size(); j++) {
	Pole pole = tile.poles_for_lane_groups()[i].poles()[j];
	cout<<"tree: "<<pole.bottom_cross_section_diameter_cm()<<endl;
	double pole_top_diameter = pole.top_cross_section_diameter_cm() /10000000.0;
	double pole_bottom_diameter = pole.bottom_cross_section_diameter_cm() /10000000.0;
	Point3d top_center_point = pole.top_center_point();
	Point3d bottom_center_point = pole.bottom_center_point();
	double pole_height = (top_center_point.cm_from_wgs84_ellipsoid() - bottom_center_point.cm_from_wgs84_ellipsoid())  /10000000.0;
	double average_diameter = (pole_top_diameter+pole_bottom_diameter) /2.0;
	int64_t bottom_center_2d_coordinate = bottom_center_point.here_2d_coordinate();

	int64_t theXor = (bottom_center_2d_coordinate ^ 5960596415578112000); //tile center
	std::pair<double, double> decoded = decode_morton_2d(theXor);

	double lat = decoded.second;
	double lon = decoded.first;


	if(lon>(0.02197265625/2))
		lon = -(0.02197265625) +lon;
	if(lat>(0.02197265625/2))
		lat = -(0.02197265625) +lat;

	cout<<"Lat: "<<fixed<<lat<<"   Long: "<<fixed<<lon<<" pole diameter"<<fixed<<average_diameter<<endl;
	//cout<<"lat "<<decode_morton_2d(4354955124161939766).second<<endl;
	//cout<<"lon "<<decode_morton_2d(4354955124161939766).first<<endl;
	//Writing to CSV
	ofs << lat << "," << lon << "," << fixed<<pole_height<<","<<average_diameter<<std::endl;

      
      pole_id++;
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

  LandmarkPolesLayerTile tile;

  {
    // Read the existing address book.
    fstream input(argv[1], ios::in | ios::binary);
    if (!tile.ParseFromIstream(&input)) {
      cerr << "Failed to parse address book." << endl;
      return -1;
    }
  }

  ListPoles(tile);

  // Optional:  Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
