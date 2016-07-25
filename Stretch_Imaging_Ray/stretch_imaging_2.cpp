// Stretch_Imaging_Ray.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stretch_imaging_2.h"
#include <stdint.h>
#include <math.h>  
#include <stdio.h>
#include <ctime>

#include <iostream>
#include <random>
#include <string>
#include <sstream>

#include <Eigen/Core>

using namespace std;



int test2 (double delta_u, double delta_v, int32_t lens_index, double aper_width, int32_t number_of_rays_total , int data_index)
{
	//int32_t height = 351; //The number of row in the output array
    //int32_t width = 351; //The number of columns in the output array 
	int32_t number_of_rays  = std::pow(number_of_rays_total, 1.0/4);  //The 4th root of the number of rays along the stretched direction
	int32_t lens_index_s = lens_index;
	//double scene_size = 70; //mm
	double scene_height = 50.0; //mm, just for reference 
	double sensor_width = 0.05; //mm
	double aperture_width = aper_width; //mm
	double aperture_height = 5.0; //mm
	double lens_width = lens_geo_width(delta_u);
/** Examples for using vectors
	Eigen::Vector2d a = Eigen::Vector2d(1.0, 2.0);
	Eigen::Vector2d b = Eigen::Vector2d(2.0, 3.0);

	Eigen::Vector2d x = a+b;

	const double t = a.dot(b);
	const double na = a.norm();

	const double ax = a.x();
	const double ax_alternative = a(0);

	const double ay = a.y();
	const double ay_alternative = a(1);
**/ 
  
  FILE * pFile;
  //double buffer[6250000] = { 0.0 };  

  std::ostringstream oss;
  oss << "rays_" << data_index ;
  std::string file_name_for_rays = oss.str();
  const char *filename = file_name_for_rays.c_str();
  pFile = fopen (filename, "wb");


  //fwrite (&height , sizeof(int32_t), 1, pFile);
  //fwrite (&width , sizeof(int32_t), 1, pFile);

  //int nElems = height*width; //Find the number of elements in the array
  //double *buffer_r = (double*)malloc(nElems*sizeof(double)); //Set the memory needed to store these values
  //double *buffer_g = (double*)malloc(nElems*sizeof(double)); //Set the memory needed to store these values
  //double *buffer_b = (double*)malloc(nElems*sizeof(double)); //Set the memory needed to store these values
  //double *buffer_positions = (double*)malloc(number_of_rays*number_of_rays*8*sizeof(double)+sizeof(double)); //Set the memory needed to store these values
  //for(int i=0; i<nElems; i++) //Fill the arrays with zeros
  //{
	//  buffer_r[i] = 0.0;
	 // buffer_g[i] = 0.0;
	 // buffer_b[i] = 0.0;
  //}
  // for(int i=0; i<number_of_rays*number_of_rays*8+1; i++) //Fill the arrays with zeros
  //{
	//  buffer_positions[i]=0.0;
  //}




  //
  //int32_t *buffer_write_int_info = (int32_t*)malloc(2*sizeof(int32_t));
  //buffer_write_int_info[0] = number_of_rays;
  //buffer_write_int_info[1] = lens_index;
  fwrite (&number_of_rays, sizeof(int32_t),1, pFile);
  fwrite (&lens_index, sizeof(int32_t),1, pFile);
  //cout << "number_of_rays " << buffer_write_int_info[1] << endl;
  
  //
  //
  //double *buffer_write_double_info = (double*)malloc(8*sizeof(double));
  double delta_at_len_index_u = (delta_u/33)*lens_index; 
  double lens_radius = lens_geo_radius_u(delta_u, delta_v);
  double curved_surface_y = lens_geo_lens_center(delta_u, delta_v)+aperture_height;
  //fwrite (&buffer_write_double_info, sizeof(double),8, pFile);
  fwrite (&delta_at_len_index, sizeof(double),1, pFile);
  fwrite (&lens_radius, sizeof(double),1, pFile);
  fwrite (&curved_surface_y, sizeof(double),1, pFile);
  fwrite (&sensor_width, sizeof(double),1, pFile);
  fwrite (&aperture_width, sizeof(double),1, pFile);
  fwrite (&aperture_height, sizeof(double),1, pFile);
  fwrite (&delta_u, sizeof(double),1, pFile);
  fwrite (&delta_v, sizeof(double),1, pFile);
  //free(buffer_write_double_info);
  //later
  Eigen::Vector2d curved_surface_location = Eigen::Vector2d(0,curved_surface_y);
  Eigen::Vector2d lens_array_bottom_surface_normal = Eigen::Vector2d(0, -1); //This defines the first surface normal which is the base of the lens array. The ray is traveling from air to the material of the lens 
  // https://www.guyrutenberg.com/2014/05/03/c-mt19937-example/
  mt19937::result_type seed = time(0);
  cout << "seed: " << seed << endl;
  //std::mt19937 mt_rand(time(0));
  auto real_rand = std::bind(std::uniform_real_distribution<double>(0,1),mt19937(seed));
  //cout << "test: " << real_rand() << endl;
  //std::uniform_real_distribution<> dis(1, 2);
  bool write_rays = true;
  for(int i_u=0;i_u<number_of_rays;i_u++) //For loop covering the change in the top point that defines the starting ray 
  {
	  //cout << "current i: " << i << endl;
	  for(int i_v=0;i_v<number_of_rays;i_v++)
	  {
		  for(int j_u=0;j_u<number_of_rays;j_u++)
		  {
			  for(int j_v=0;j_v<number_of_rays;j_v++)
			  {
				//cout << "--------- " << endl;
				//cout << "current i: " << i << endl;
				//cout << "current j: " << j << endl;
				double shift_amount = aperture_width / (number_of_rays);
				double top_position_start_v = ( (-0.5*aperture_width)) + (shift_amount*i)+delta_at_len_index; //This defines the x position of the top point
				double bottom_position_start_u = ( (-0.5*aperture_width) ) + (shift_amount*j)+delta_at_len_index; //This defines the x position of the bottom point
				double top_position = top_position_start + real_rand()*shift_amount;
				double bottom_position = bottom_position_start + real_rand()*shift_amount;
				//double top_position = top_position_start + 1*shift_amount;
				//double bottom_position = bottom_position_start + 1*shift_amount;
				//cout << "First Ray Start: " << top_position_start << endl;
				//cout << "Second ray Start: " << bottom_position_start << endl;
				//system("pause");
				double y = aperture_height; //This sets the apature height 
				double x = top_position;
				Eigen::Vector2d starting_ray = Eigen::Vector2d(-(top_position-bottom_position), -y); //This defines the first ray leaving the aperture 
				Eigen::Vector2d aperture_ray_normalized = starting_ray.normalized(); 
				//cout << "Ray unnormal: " << aperture_ray << endl;
				//cout << "Ray normal: " << aperture_ray_normalized << endl;
				Eigen::Vector2d internal_ray = snell_s(1.00, 1.41, lens_array_bottom_surface_normal, aperture_ray_normalized); //This is the ray that travels inside the lens. 
				Eigen::Vector2d origin_of_internal = Eigen::Vector2d(top_position, aperture_height);// this defines the point of the first refraction of ray ij 
				Eigen::Vector2d intersection_point = internal_ray_lens_intersection( -1*internal_ray, origin_of_internal, curved_surface_location, lens_radius);//Interaction of internal ray with the curvered surface of the lens 
				//double outside_lens_value = (lens_width/2)+delta_at_len_index;
				double outside_lens_value = (lens_width/2);
				if (intersection_point.x()>outside_lens_value)
				{
					//system("pause");
					fclose (pFile);  
					write_rays = false; 
					throw 20;
				}
				if (intersection_point.x()<(-1*outside_lens_value))
				{
					//system("pause");
					fclose (pFile);  
					write_rays = false;
					throw 30;
				}

				//cout << "Lens Width: " << lens_width/2 << endl;
				//cout << "Ray X: " << intersection_point.x() << endl;
				//cout << "Ray Y: " << intersection_point.y() << endl;
				//system("pause");

				Eigen::Vector2d normal_at_intersection_point = surface_normal_of_a_sphere( intersection_point, curved_surface_location);
				Eigen::Vector2d final_ray = snell_s(1.41, 1.00, normal_at_intersection_point, -1*internal_ray); //This is the ray that is refracted by the spherical part of the lens
				Eigen::Vector2d final_ray_normalized = final_ray.normalized();
				Eigen::Vector2d point_defining_scene = Eigen::Vector2d(0,(scene_height+23+aperture_height)); //A point on the scene plane 
				Eigen::Vector2d scene_ray_intersection = final_ray_scene_intersection( point_defining_scene, final_ray_normalized, intersection_point); //The point on the scene where the ray leaving the lens intersects 
				//double conversion_factor = height/scene_size;  //This conversation factor allows us to convert the ray's intersection with the scene from mm to matrix index 
				//int x_contact = floor(scene_ray_intersection.x()*conversion_factor + height/2); //This is the X position in a matrix for where the ray lands
				//int y_contact = floor(scene_ray_intersection.y()*conversion_factor + width/2); //This is the Y position in a matrix for where the ray lands.
				//int column = floor(nElems/ double(y_contact)); 
				//int array_index = column*height+x_contact;
				//double previous_value = buffer_r[array_index];
				//Now we calculate the energy of this ray
				Eigen::Vector2d d_i = Eigen::Vector2d(top_position, y) - Eigen::Vector2d(bottom_position, 0);
				Eigen::Vector2d w_i = d_i/d_i.norm();
				double c_i = w_i.dot(Eigen::Vector2d(0,1));
				double L_i = d_i.norm();
		
				//buffer_r[array_index] = previous_value+energy_of_ray;
				//cout << "energy_of_ray: " << energy_of_ray << endl;
				//cout << "array index " << array_index << endl;
				//x_contact = t[0]*conversion_factor + float(1024)
				//y_contact = t[1]*conversion_factor + float(1024)
		//		double *buffer_positions = (double*)malloc(9); //Set the memory needed to store these values
		//		buffer_positions[0] = bottom_position;
		//		buffer_positions[1] = 0;
		//		buffer_positions[2] = top_position;
		//		buffer_positions[3] = aperture_height;
		//		buffer_positions[4] = intersection_point.x();
		//		buffer_positions[5] = intersection_point.y();
		//		buffer_positions[6] = scene_ray_intersection.x();
		//		buffer_positions[7] = scene_ray_intersection.y();
				//double *buffer_write_ray_info = (double*)malloc(5*sizeof(double));
				double energy_of_ray = (c_i*c_i)/(L_i*L_i);
				//buffer_write_ray_info[0] = bottom_position;
				//buffer_write_ray_info[1] = 0;
				//buffer_write_ray_info[2] = top_position;
				//buffer_write_ray_info[3] = aperture_height;
				//buffer_write_ray_info[0] = intersection_point.x();
				//buffer_write_ray_info[1] = intersection_point.y();
				//buffer_write_ray_info[2] = scene_ray_intersection.x();
				//buffer_write_ray_info[3] = scene_ray_intersection.y();
				//buffer_write_ray_info[4] = energy_of_ray;
				double bottom_position_y = 0.0;
				if (write_rays)
				{
					fwrite (&bottom_position, sizeof(double),1, pFile);
					fwrite (&bottom_position_y, sizeof(double),1, pFile);
					fwrite (&top_position, sizeof(double),1, pFile);
					fwrite (&aperture_height, sizeof(double),1, pFile);
					fwrite (&intersection_point.x(), sizeof(double),1, pFile);
					fwrite (&intersection_point.y(), sizeof(double),1, pFile);
					fwrite (&scene_ray_intersection.x(), sizeof(double),1, pFile);
					fwrite (&scene_ray_intersection.y(), sizeof(double),1, pFile);
					fwrite (&energy_of_ray, sizeof(double),1, pFile);

				}
				//free(buffer_write_ray_info);
				//cout << "x2: " << top_position << endl;
				//cout << "y1: " << 0 << endl;
				//cout << "x2: " << top_position << endl;
				//cout << "y2: " << aperture_height << endl;
				//cout << "x3: " << intersection_point.x() << endl;
				//cout << "y3: " << intersection_point.y() << endl;
				//cout << "x4: " << scene_ray_intersection.x() << endl;
				//cout << "y4: " << scene_ray_intersection.y() << endl;
				//cout << "Ray count" << ray_count << endl;
				//if (abs(top_position)>(aperture_width/2)||abs(bottom_position)>(aperture_width/2))
				//{
					//cout << "x1: " << bottom_position << endl;
					//cout << "x2: " << top_position << endl;
				//}
			  }
		  }
		}
	}
  //double lens_radius = lens_geo_radius_u(delta_u, delta_v);
  //double curved_surface_y = lens_geo_lens_center(delta_u, delta_v)+aperture_height;
  //fwrite (buffer_r , sizeof(double), nElems, pFile);

 // fwrite (buffer_positions , sizeof(double), number_of_rays*number_of_rays*8+1, pFile);
  //fwrite (buffer_g , sizeof(double), nElems, pFile);
  //fwrite (buffer_b , sizeof(double), nElems, pFile);
  fclose (pFile);  
  //free(buffer_r);
  //free(buffer_g);
  //free(buffer_b);
  //free(buffer_positions);
  //free(filename);
  //free(buffer_write_int_info);

  return 0;
}




Eigen::Vector2d final_ray_scene_intersection( const Eigen::Vector2d &scene_point, const Eigen::Vector2d &final_ray_normalized, const Eigen::Vector2d &intersection_point_on_curved_surface)
{
	double d = (scene_point-intersection_point_on_curved_surface).dot(Eigen::Vector2d(0,-1))/(final_ray_normalized.dot(Eigen::Vector2d(0,-1)));
	return final_ray_normalized*d+intersection_point_on_curved_surface;
}

/**
Eigen::Vector2d internal_ray_lens_intersection( const Eigen::Vector2d &direction_of_internal, const Eigen::Vector2d &origin_of_internal, const Eigen::Vector2d &centerpoint_of_sphere, double radius_of_sphere)
{
	//https://people.cs.clemson.edu/~dhouse/courses/405/notes/raycast.pdf
	Eigen::Vector2d return_intersection;
	Eigen::Vector2d direction_of_internal_normalized = direction_of_internal.normalized();
	Eigen::Vector2d o_minus_c = origin_of_internal-centerpoint_of_sphere;
	double discriminant = (direction_of_internal_normalized.dot(o_minus_c))*(direction_of_internal_normalized.dot(o_minus_c))-o_minus_c.norm()*o_minus_c.norm()+radius_of_sphere*radius_of_sphere;
	if (discriminant < 0)
	{
		throw 30;
	}
	else
	{
		double d_1 = -1*(direction_of_internal_normalized.dot(o_minus_c))+sqrt((direction_of_internal_normalized.dot(o_minus_c))*(direction_of_internal_normalized.dot(o_minus_c))-o_minus_c.norm()*o_minus_c.norm()+radius_of_sphere*radius_of_sphere);
		double d_2 = -1*(direction_of_internal_normalized.dot(o_minus_c))-sqrt((direction_of_internal_normalized.dot(o_minus_c))*(direction_of_internal_normalized.dot(o_minus_c))-o_minus_c.norm()*o_minus_c.norm()+radius_of_sphere*radius_of_sphere);
		if (d_2=d_1)
		{
			return_intersection = origin_of_internal+d_2*direction_of_internal_normalized;
		}
		else
		{
			Eigen::Vector2d return_intersection_1 = origin_of_internal+d_1*direction_of_internal_normalized;
			Eigen::Vector2d return_intersection_2 = origin_of_internal+d_2*direction_of_internal_normalized;
			if (return_intersection_1.y()<return_intersection_2.y())
			{
				return_intersection = return_intersection_1;
			}
			else
			{
				return_intersection = return_intersection_2;
			}

		}
		cout << "lens choice 1: " << d_1 << endl;
		cout << "lens choice 2: " << d_2 << endl;
	}
	return return_intersection;
}
**/
Eigen::Vector2d internal_ray_lens_intersection( const Eigen::Vector2d &direction_of_internal, const Eigen::Vector2d &origin_of_internal, const Eigen::Vector2d &centerpoint_of_sphere, double radius_of_sphere)
{
	// https://people.cs.clemson.edu/~dhouse/courses/405/notes/raycast.pdf
	Eigen::Vector2d return_intersection;
	Eigen::Vector2d direction_of_internal_normalized = direction_of_internal.normalized(); 
	Eigen::Vector2d o_minus_c = origin_of_internal-centerpoint_of_sphere;
	double t_close = direction_of_internal.dot(centerpoint_of_sphere-origin_of_internal);
	Eigen::Vector2d x_close = origin_of_internal + t_close*direction_of_internal;
	double discriminant = (x_close-centerpoint_of_sphere).norm();
	if (discriminant < radius_of_sphere)
	{
		double a = sqrt(radius_of_sphere*radius_of_sphere-discriminant*discriminant);
		return_intersection = origin_of_internal + (t_close-a)*direction_of_internal;
		//cout << "Centerpoint of sphere : " << centerpoint_of_sphere << endl;
		//cout << "radius of sphere : " << radius_of_sphere << endl;
	}
	else
	{
		throw 30;
	}
	return return_intersection;
}

Eigen::Vector2d surface_normal_of_a_sphere(const Eigen::Vector2d &position_on_sphere, const Eigen::Vector2d &centerpoint_of_sphere)
{
	Eigen::Vector2d surface_vector = position_on_sphere - centerpoint_of_sphere;
	Eigen::Vector2d surface_vector_normalized = surface_vector.normalized();
	if (surface_vector_normalized.y()>0)
	{
		return -1*surface_vector_normalized;
	}
	else
	{
		return surface_vector_normalized;
	}
}
double lens_geo_radius_u(double delta_u, double delta_v)//This returns the radius of the lens along a plane through the optical axis in the u direction
{
	double G2 = 0.8667504192892; //initial chord height 
	double A1 = 231.00; //lens width 
	double v = 231.00; //lens depth 
	double n = 33.0;// number of lens along the stretched direction 
	double B2 = delta_u+A1;
	double J2 = (A1+delta_u)/n;
	//double radius = ((2*u*v*a)/((u+delta_u)*(v+delta_v))) + (((u+delta_u)/n)*((u+delta_u)/n)*(u+delta_u)*(v+delta_v))/(8*u*v*a);
	double radius = ((A1*G2)/(2*B2))+(J2*J2*B2)/(8*G2*A1);
	return radius;
}

double lens_geo_width(double delta_u)//This returns the radius of the lens along a plane through the optical axis in the u direction
{ 
	double A1 = 231.00; //lens width 
	double v = 231.00; //lens depth 
	double n = 33.0;// number of lens along the stretched direction 
	double B2 = delta_u+A1;
	double J2 = (A1+delta_u)/n;
	//double radius = ((2*u*v*a)/((u+delta_u)*(v+delta_v))) + (((u+delta_u)/n)*((u+delta_u)/n)*(u+delta_u)*(v+delta_v))/(8*u*v*a);
	return J2;
}

double lens_geo_thickness(double delta_u, double delta_v)//This returns the radius of the lens along a plane through the optical axis in the u direction
{
	double u = 231.00; //lens width initial 
	double v = 231.00; //lens depth initial
	double T = 23; //lens thickness initial 
	int n = 33;// number of lens along the stretched direction 
	double deformed_thickness = (u*v*T)/((u+delta_u)*(v+delta_v));
	return deformed_thickness;
}

double lens_geo_chord_height(double delta_u, double delta_v)
{
	double T = 23; //lens thickness initial 
	double a = 0.8667504192892; //initial chord height 
	double deformed_thickness = lens_geo_thickness(delta_u,  delta_v);
	double deformed_chord_height = ((T+deformed_thickness)/T)*a;
	return deformed_chord_height;
}

double lens_geo_lens_center(double delta_u, double delta_v) //This returns the y position of that circle that makes up the curved geometry of a lens in the flexible lens array 
{
	double y_position = lens_geo_thickness( delta_u,  delta_v) - lens_geo_radius_u( delta_u,  delta_v);
	return y_position;
}


/**
snell_s returns the refracted vector, with the incident and normal
vectors defined in the following system

     ^
     |
     |^
     |r
     |
_____|_________
    /| 
   / |
^ /  |
I/   |
v    | ^
     V n
**/


Eigen::Vector2d snell_s (double n1, double n2, const Eigen::Vector2d &input_normal_vector, const Eigen::Vector2d &input_incident_vector)
{
	Eigen::Vector2d normal_vector = input_normal_vector.normalized();
	Eigen::Vector2d incident_vector = input_incident_vector.normalized();

	if (n1 > n2) //check for critcal angle 
	{
		double critical_angle = asin(n1/n2);
		double angle_between_n_I = acos(normal_vector.dot(incident_vector)/(normal_vector.norm()*incident_vector.norm()));
		if (angle_between_n_I>critical_angle)
		{
			throw 20;
		}
	}
	double c = normal_vector.dot(incident_vector);
	double r = n1/n2;
	Eigen::Vector2d v_refract = (r*c - sqrt(1-r*r*(1-c*c)))*normal_vector - r*incident_vector;
	return v_refract;
}


// def three_d(n1,n2,numpy_normal_vector,numpy_incident_vector):
//    r = float(n1)/float(n2)
//    normal_vector = numpy_normal_vector/np.linalg.norm(numpy_normal_vector)
//    incident_vector = numpy_incident_vector/np.linalg.norm(numpy_incident_vector)
//    if float(n1) > float(n2):
//        critical_angle = math.asin(float(n2)/float(n1))
//        angle_between_n_I = math.acos(np.dot(normal_vector,incident_vector)/(np.linalg.norm(normal_vector)*np.linalg.norm(incident_vector)))
//        if angle_between_n_I > critical_angle:
//            raise Exception("Total Interal Reflection Detected")
//    c = np.dot(normal_vector,incident_vector)
//    v_refract = (r*c - math.sqrt(1-r*r*(1-c*c)))*normal_vector - r*incident_vector
//    return v_refract

int test1 ()
{
	
  
  FILE * pFile;
  double buffer[] = { 1.1 , 1.2 , 1.3 , 1.5 , 1.6, 1.7 };  
  pFile = fopen ("raydata.bin", "wb");
  fwrite (buffer , sizeof(double), 6, pFile);
  
  cout << "sizeof(double): " << sizeof(double) << endl;
  cout << "sizeof(buffer): " << sizeof(buffer) << endl;

  fclose (pFile);
  return 0;
}