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

#include "sphere.h"
#include <iostream>
#include <math.h>

/************************** Sphere **********************************/

Hit Sphere::intersect(const Ray &ray)
{
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

    Vector diff = ray.O - position;
    //Discriminant
    double discr = ((ray.D.dot(diff)) * (ray.D.dot(diff))) - (ray.D.dot(ray.D)*(diff.dot(diff)) - r*r);
    double t_min, t_plus;

    //Solve the formula to t
    if (discr < 0){
      return Hit::NO_HIT();
    } else {
      t_min = (-ray.D.dot(ray.O-position) - sqrt(discr))/ray.D.dot(ray.D);
      t_plus = (-ray.D.dot(ray.O-position) + sqrt(discr))/ray.D.dot(ray.D);
    }
    double t = min(t_min,t_plus);

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
