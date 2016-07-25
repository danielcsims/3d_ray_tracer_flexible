// Stretch_Imaging_Ray.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "stretch_imaging_2.h"
#include "stretch_imaging_2.h"
#include <random>

#include <Eigen/Core>
#include <ctime>

using namespace std;

int a;

 //This is the index of the lens under study. 
//double aper_width = 0.1;
int data_index = 0;

int main (int argc, char* argv[])
{
	if (argc != 4) 
	{
		cout << "nope" << endl;
		return -1;
	}
	double delta_u = atof(argv[1]);
	double delta_v = atof(argv[2]);
	double li = atof(argv[3]);
	double aper_width = atof(argv[4]);
	int file_index = atoi(argv[5]);

	//double stretch = atof(argv[4]); //This is the amount the whole lens array is stretched.
	//double samples = atoi(argv[5]); //This is the number of samples we measure as the lens is stretched.
	//double delta = stretch/(samples); //This is the amount each lens is stretched inbetween each sample. 
	test2 (delta_u,delta_v, li,aper_width, file_index);
}




/**
int main ()
{
  double n1 = 1;
  double n2 = 1.41;
  Eigen::Vector2d input_normal_vector = Eigen::Vector2d(0, -1);
  Eigen::Vector2d input_incident_vector = Eigen::Vector2d(-1, -1);
  Eigen::Vector2d result = snell_s(n1,n2,input_normal_vector,input_incident_vector);
  cout << "Surface normal: " << input_normal_vector << endl;
  cout << "Incident Ray: " << input_incident_vector << endl;
  cout << "Incident Ray: " << input_incident_vector << endl;
  getchar();  
  return 0;
}
**/

