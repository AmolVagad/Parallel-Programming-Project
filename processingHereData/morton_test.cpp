
#include <cstdint>
#include <utility>
#include <cassert>
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <algorithm>

using namespace std;
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

    uint32_t lon = std::stol(odd, nullptr, 2);
    uint32_t lat = std::stol(even, nullptr, 2);
   
    int x = (lon - 1);
    int y = (lat - 1);

    double x_new = x/pow(2,31)*180;
    double y_new = y/pow(2,31)*180;
   
  return std::make_pair(x_new,y_new);
}

int main()
{
  std::pair<double, double> decoded = decode_morton_2d((int64_t)5960596415578112000);

  std::cout<<"Decoded Longitude "<<decoded.first<<std::endl;
  std::cout<<"Decoded Latitude "<<decoded.second<<std::endl;

  std::pair<double, double> decoded2 = decode_morton_2d((int64_t)38853464117);

  std::cout<<"Shifted to object Longitude "<<decoded2.first + decoded.first<<std::endl;
  std::cout<<"Shifted to object Decoded Latitude "<<decoded2.second + decoded.second<<std::endl;
   
  return 0;
}

