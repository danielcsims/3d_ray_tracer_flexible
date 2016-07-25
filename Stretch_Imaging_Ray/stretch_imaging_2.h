#ifndef __STRETCH_IMAGING_2_H__
#define __STRETCH_IMAGING_2_H__

#include "stdafx.h"
#include <iostream>
#include <Eigen/Core>
using namespace std;
typedef __int32 int32_t;

int test1 ();
int test2 (double delta_u, double delta_v, int32_t lens_index, double aper_width, int data_index);
Eigen::Vector2d snell_s (double n1, double n2, const Eigen::Vector2d &input_normal_vector, const Eigen::Vector2d &input_incident_vector);
double lens_geo_radius_u(double delta_u, double delta_v);
double lens_geo_thickness(double delta_u, double delta_v);
double lens_geo_chord_height(double delta_u, double delta_v);
double lens_geo_lens_center(double delta_u, double delta_v); //This returns the y position of that circle that makes up the curved geometry of a lens in the flexible lens array 
double lens_geo_width(double delta_u);
Eigen::Vector2d final_ray_scene_intersection(const Eigen::Vector2d &scene_point, const Eigen::Vector2d &final_ray_normalized, const Eigen::Vector2d &intersection_point_on_curved_surface);
Eigen::Vector2d internal_ray_lens_intersection( const Eigen::Vector2d &direction_of_internal, const Eigen::Vector2d &origin_of_internal, const Eigen::Vector2d &centerpoint_of_sphere, double radius_of_sphere);
Eigen::Vector2d surface_normal_of_a_sphere(const Eigen::Vector2d &position_on_sphere, const Eigen::Vector2d &centerpoint_of_sphere);

extern int a;

#endif
