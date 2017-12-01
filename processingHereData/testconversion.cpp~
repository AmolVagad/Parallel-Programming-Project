#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
//#include "com/here/pb/hdmap/external/v1/lanes/layer-lane-geometry-polyline.pb.h"

using namespace std;

//using namespace com::here::pb::hdmap::external::v1::geometry;
//using namespace com::here::pb::hdmap::external::v1::lanes;

ofstream ofs;
std::string filename = "laneData.csv";

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
    even = even[0] + even;
    string test = "00000000000000110010001011100001";
    long int lon = twos_complement_to_decimal(odd);
    long int lat = twos_complement_to_decimal(even);
   
    long int x = (lon );//- 1);
    long int y = (lat );//- 1);

    double x_new = x/pow(2,31)*180;
    double y_new = y/pow(2,31)*180;
   
  return std::make_pair(x_new,y_new);
}

//0.0784735/2
int main(int argc, char* argv[]) 
{
  int64_t test = 4354955124161939766;

  std::pair<double, double> decoded_tile_center = decode_morton_2d(test);

  cout<<"Longitude "<<decoded_tile_center.first<<endl;
  cout<<"Latitude "<<decoded_tile_center.second<<endl;



  return 0;
}



