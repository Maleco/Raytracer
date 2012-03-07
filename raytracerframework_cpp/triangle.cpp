//
//  Framework for a raytracer
//  File: triangle.cpp
//
//  Created for the Computer Science course "Introduction Computer Graphics"
//  taught at the University of Groningen by Tobias Isenberg.
//
//  Authors:
//    Tom Eijkelenkamp
//    Hans-Rudolf Woldring
//
//  This framework is inspired by and uses code of the raytracer framework of 
//  Bert Freudenberg that can be found at
//  http://isgwww.cs.uni-magdeburg.de/graphik/lehre/cg2/projekt/rtprojekt.html 
//

#include "triangle.h"
#include <iostream>
#include <math.h>

/************************** Sphere **********************************/

Hit Triangle::intersect(const Ray &ray)
{
    /****************************************************
    * RT1.1: INTERSECTION CALCULATION
    *
    * Given: ray, A, B, C
    * Sought: intersects? if true: *t
    * 
    * Insert calculation of ray/triangle intersection here. 
    *
    * If the ray does not intersect the triangle, return false.
    * Otherwise, return true and place the distance of the
    * intersection point from the ray origin in *t (see example).
    ****************************************************/
    
    double a = A.x-B.x,     d = A.x-C.x,     g = ray.D.x,
           b = A.y-B.y,     e = A.y-C.y,     h = ray.D.y,
           c = A.z-B.z,     f = A.z-C.z,     i = ray.D.z,
           
           j = A.x-ray.O.x, k = A.y-ray.O.y, l = A.z-ray.O.z;
   
    double M, B, Y, t;
   
    M = a*(e*i - h*f)   +  b*(g*f - d*i)  +  c*(d*h - e*g);
    B = (j*(e*i - h*f)  +  k*(g*f - d*i)  +  l*(d*h - e*g))  / M;    
    Y = (i*(a*k - j*b)  +  h*(j*c - a*l)  +  g*(b*l - k*c))  / M;
    t = (f*(a*k - j*b)  +  e*(j*c - a*l)  +  d*(b*l - k*c))  / M;

    
    if (Y<0 || Y>1) {
      return Hit::NO_HIT();
    }
    
    if (B<0 || B>1-Y) {
      return Hit::NO_HIT();
    }
        
    return Hit(t,N);
}
