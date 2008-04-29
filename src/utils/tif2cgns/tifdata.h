#ifndef _tifdata_h_
#define _tifdata_h_
#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

typedef struct
{
	int    index;
	double x;
	double y;
	double h;
} Node_t;


typedef struct
{
	int index;
	int point1;
	int point2;
	int bcode;
	int bc_type;
} Edge_t;


typedef struct
{
	int  index;
	int  node_num;             //the node num of this region
	int  tri_num;		   //the triangle num of this region
	char type[32];
	char name[32];
	vector<int> boundary;   
} Region_t;


typedef struct
{
	int  index;
	char type[32];
	char name[32];
	int  region;
	vector<int> boundary;
}Interface_t;


typedef struct
{
	int  index;
	int  region;
	int  c1,c2,c3;
	int  t1,t2,t3;
} Tri_t;

typedef struct
{
	int  region;
	vector<string> parameter_name_array;
	vector<double> parameter_value_array;
}Parameter_t;

typedef struct
{
	int sol_num;
	vector<string> sol_name_array;
}SolHead_t;


typedef struct
{
	int index;
	string  material;
	vector<double> data_array;
}SolData_t;

typedef struct
{
	char sol_name[32];
	vector<string> data_name_array;
	vector<string> data_unit_array;
	vector<double> data_value_array;
}TranSol_t;

class Component_t
{
  public:
  	int region;
  	double xmin,xmax;
  	double ymin,ymax;
  	int    direction;
  	double mole_begin;
  	double mole_ratio;
  public:
  	void mole(double x, double y, double & mole_x)
  	{
  	  if(x>=xmin-1e-10 && x <=xmax+1e-10 && y>=ymin-1e-10 && y <=ymax+1e-10)
  	  {
  	  	if(direction == 0) 	 mole_x = mole_begin;
  	  	if(direction == 1)	 mole_x = mole_begin + mole_ratio*(x-xmin);
  	  	if(direction == 2)	 mole_x = mole_begin - mole_ratio*(y-ymax);	
  	  	if(mole_x<0) 	         mole_x = 0;	
  	  }
  	}
};

extern vector<Node_t>      node_array;
extern vector<Edge_t>      edge_array;
extern vector<Tri_t>       tri_array;
extern vector<Region_t>    region_array;
extern vector<Interface_t> interface_array;
extern vector<Component_t> component_array;
extern vector<Parameter_t> parameter_array;
extern SolHead_t           solhead;
extern vector<SolData_t>   soldata;
extern TranSol_t           transol;
#endif

