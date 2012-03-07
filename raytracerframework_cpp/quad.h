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

#ifndef QUAD_H_115209AE
#define QUAD_H_115209AE

#include "object.h"

class Quad : public Object
{
public:
    Point A;
    Point B;
    Point C;
    Point D;
    Vector N;
  
    Quad(const Point A, const Point B, const Point C, const Point D) 
	: A(A), B(B), C(C), D(D), N((B-A).cross(C-A))
    { }

    virtual Hit intersect(const Ray &ray);
};

#endif /* end of include guard: QUAD_H_115209AE */
