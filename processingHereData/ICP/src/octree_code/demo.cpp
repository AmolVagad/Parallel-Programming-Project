#include "octree.h"

#ifdef WIN32
#include <tchar.h>
int _tmain(int argc, _TCHAR* argv[])
#else
#include<iostream>
#include<stdio.h>
using namespace std;
int main(int argc, char* argv[])
#endif

{
    Octree<vector<double>> o(4096); /* Create 4096x4096x4096 octree containing doubles. */
    o(1,2,3).push_back(3.1416);      /* Put pi in (1,2,3). */
   //s o.erase(1,2,3);         /* Erase that node. */
   cout<<o(1,2,3).at(0)<<endl;
	return 0;
}


