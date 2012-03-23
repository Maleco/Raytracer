//
//  Framework for a raytracer
//  File: scene.h
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

#ifndef SCENE_H_KNBLQLP6
#define SCENE_H_KNBLQLP6

#include <vector>
#include <string.h>
#include "triple.h"
#include "light.h"
#include "object.h"
#include "image.h"

enum RenderMode {
	phong,
	zbuffer,
	normal,
	gooch
};

class Scene
{
private:
   std::vector<Object*> objects;
   std::vector<Light*> lights;
   Triple eye;
   Triple center;
   Triple up;
   double kBlue;
   double kYellow;
   double alpha;
   double beta;
   RenderMode rendermode;
   bool Shadow;
   int maxRecursionDepth;
   int superSamplingFactor;
   int apertureSamples;
   int apertureRadius;

public:
   Hit findMinHit(const Ray &ray);
   Color traceNormal(const Ray &ray);
	double traceZBuffer(const Ray &ray);
	Color superSample(Point p, Vector right, double pSize);
	Color depthOfField(Point p, Vector right, double pSize);
   Color tracePhongGooch(const Ray &ray, int recursionDepth);
   Color trace(const Ray &ray, Image &img);
   void render(Image &img);
   void addObject(Object *o);
   void addLight(Light *l);
   void setEye(Triple e);
   void setBlue(double b);
   void setYellow(double y);
   void setBeta(double b);
   void setAlpha(double a);
   void setCenter(Triple e);
   void setUp(Triple e);
   void setRenderMode(string rm);
   void setShadowMode(bool shadowMode);
   void setMaxRecursionDepth(int max);
   void setSuperSamplingFactor(int superSampling);
   void setApertureSamples(int as);
   void setApertureRadius(int ar);
   unsigned int getNumObjects() { return objects.size(); }
   unsigned int getNumLights() { return lights.size(); }
};

#endif /* end of include guard: SCENE_H_KNBLQLP6 */
