//
//  Framework for a raytracer
//  File: sphere.h
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

#ifndef TRIANGLE_H_115209AE
#define TRIANGLE_H_115209AE

#include "object.h"

class Triangle : public Object
{
public:
    Point A;
    Point B;
    Point C;
    Vector N;
  
    Triangle(const Point A, const Point B, const Point C) 
	: A(A), B(B), C(C), N((B-A).cross(C-A))
    { }

    virtual Hit intersect(const Ray &ray);
};

#endif /* end of include guard: TRIANGLE_H_115209AE */
