//
//  Framework for a raytracer
//  File: sphere.cpp
//
//  Created for the Computer Science course "Introduction Computer Graphics"
//  taught at the University of Groningen by Tobias Isenberg.
//
//  Authors:
//    Maarten Everts
//    Jasper van de Gronde
//
//  This framework is inspired by and uses code of the raytracer framework of
//  Bert Freudenberg that can be found at
//  http://isgwww.cs.uni-magdeburg.de/graphik/lehre/cg2/projekt/rtprojekt.html
//
#define PI (3.141592653589793)

#include "sphere.h"
#include <iostream>
#include <math.h>

/************************** Sphere **********************************/

Hit Sphere::intersect(const Ray &ray){
    /****************************************************
    * RT1.1: INTERSECTION CALCULATION
    *
    * Given: ray, position, r
    * Sought: intersects? if true: *t
    *
    * Insert calculation of ray/sphere intersection here.
    *
    * You have the sphere's center (C) and radius (r) as well as
    * the ray's origin (ray.O) and direction (ray.D).
    *
    * If the ray does not intersect the sphere, return false.
    * Otherwise, return true and place the distance of the
    * intersection point from the ray origin in *t (see example).
    ****************************************************/

	double t, tp, tm, D;
    Vector V(ray.O - position);

    D = ((ray.D.dot(V))*(ray.D.dot(V))) - (ray.D.dot(ray.D))*(V.dot(V) - r*r);  // Discriminant

    if (D<0) {
      return Hit::NO_HIT();
    } else {
      tp = (-ray.D.dot(V) + sqrt(D))/ray.D.dot(ray.D);
      tm = (-ray.D.dot(V) - sqrt(D))/ray.D.dot(ray.D);
      if (tp<0 && tm<0){
			return Hit::NO_HIT();
      }
      if (tp<0){
			t = tm;
      } else if (tm<0){
			t = tp;
      } else {
			t = min(tp,tm);
      }
    }

    /****************************************************
    * RT1.2: NORMAL CALCULATION
    *
    * Given: t, C, r
    * Sought: N
    *
    * Insert calculation of the sphere's normal at the intersection point.
    ****************************************************/

    Vector N(ray.at(t) - position);
    N.normalize();

    return Hit(t,N);
}

Color Sphere::calcTexture(Point hit){
	double theta = acos((hit.z-position.z)/r);
	double phi = atan2(hit.y-position.y,hit.x-position.x);

	if (phi<0)
		phi = phi+2*PI;

	double rotX = 0.3*(90.0/360.0);
	double rotY = 0.5*(90.0/360.0);
	double rotZ = 0.7*(90.0/360.0);
	//cout << rotX << " " << rotY << endl;

	double u = rotY + (phi/(2*PI));
	double v = ((PI-theta)/PI) + rotX;

	if (u>0.0)
		u = u-1.0;

	if (v>1.0)
		v = v-1.0;


	return material->texture->colorAt(u,v);
}
