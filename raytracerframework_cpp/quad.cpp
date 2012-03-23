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

#include "quad.h"
#include "triangle.h"
#include <iostream>
#include <math.h>

/************************** Quad **********************************/

Hit Quad::intersect(const Ray &ray)
{
    /****************************************************
    * RT1.1: INTERSECTION CALCULATION
    *
    * Given: ray, A, B, C, D
    * Sought: intersects? if true: *t
    * 
    * Insert calculation of ray/quad intersection here. 
    *
    * If the ray does not intersect the quad, return false.
    * Otherwise, return true and place the distance of the
    * intersection point from the ray origin in *t (see example).
    ****************************************************/
    
    Hit min_hit(std::numeric_limits<double>::infinity(),Vector());
    
    Triangle t1(A,B,C);
    Triangle t2(C,D,A);
    
    Hit hit(t1.intersect(ray));
    if (hit.t<min_hit.t) {
	  hit.N  = N;
	  return hit;
    }  
    Hit hit2(t2.intersect(ray));
    if (hit2.t<min_hit.t) {
	  hit2.N  = N;
	  return hit2;
    }
    return Hit::NO_HIT();

}

Color Quad::calcTexture(Point hit){
	Color n;
	return n;
}

