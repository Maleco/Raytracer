//
//  Framework for a raytracer
//  File: scene.cpp
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

#include "scene.h"
#include "material.h"
#include <iostream>

#include <math.h>
#define CONST_C apertureRadius/(up.length()*sqrt(apertureSamples))
#define GOLDEN_ANGLE 137.508

double minZB = std::numeric_limits<double>::max(); // for zBuffer
double maxZB = std::numeric_limits<double>::min(); //


// ------------------------------------------------------- findMinHit ---------------------------------------------------------------- //
Hit Scene::findMinHit(const Ray &ray){
	// Find nearest object and distance
    Hit min_hit(std::numeric_limits<double>::infinity(),Vector());
    for (unsigned int i = 0; i < objects.size(); ++i) {
        Hit hit(objects[i]->intersect(ray));
        if (hit.t<min_hit.t) {
            min_hit = hit;
            }
    }
    return min_hit;
}
// ------------------------------------------------------- traceNormal --------------------------------------------------------------- //
Color Scene::traceNormal(const Ray &ray){

	//Find the nearest object (hit)
	Hit min_hit = findMinHit(ray);
	// No hit? Return background color.
	if (min_hit.t >= std::numeric_limits<double>::infinity())
		return Color(0.0, 0.0, 0.0);

	//Transform the range <-1,1> to <0, 1>
	double x = (min_hit.N.x+1)/2;
	double y = (min_hit.N.y+1)/2;
	double z = (min_hit.N.z+1)/2;

	//Set the color
	return Color(x, y, z);
}
// ------------------------------------------------------ traceZBuffer --------------------------------------------------------------- //
double Scene::traceZBuffer(const Ray &ray){

	Hit min_hit = findMinHit(ray);
	if (min_hit.t == std::numeric_limits<double>::infinity()){
		return std::numeric_limits<double>::max();
	}
	if (min_hit.t < minZB){
		minZB = min_hit.t;
	}
	if(min_hit.t > maxZB){
		maxZB = min_hit.t;
	}
	return min_hit.t;
}
// ------------------------------------------------------- depthOfField -------------------------------------------------------------- //
Color Scene::depthOfField(Point p, Vector right, double pSize){
	Color color;

	for (int n=0; n<apertureSamples; n++){
		double th = n*GOLDEN_ANGLE;
		double r = (CONST_C*sqrt(n))*cos(th)*pSize;
		p += r*right;
		p += r*up;
		Ray ray(eye, (p-eye).normalized());
		color += tracePhongGooch(ray, 1);
	}
	return color /= apertureSamples;
}
// ------------------------------------------------------ superSample -----------------------------------------------------------------//
Color Scene::superSample(Point p, Vector right, double pSize){	// p = left upper bound of pixel
	Color color;

	for (int i=1; i<=superSamplingFactor; i++){
		p -= up*pSize/(superSamplingFactor+1);					// correct y position for super sampling
		Point buffer = p;
		for (int j=1; j<=superSamplingFactor; j++){
			buffer += right*pSize/(superSamplingFactor+1);	// correct x position for super sampling
			color += depthOfField(buffer, right, pSize);
		}
	}
	color /= superSamplingFactor*superSamplingFactor;
	color.clamp();
	return color;
}
// ----------------------------------------------------- tracePhong ------------------------------------------------------------------ //
Color Scene::tracePhongGooch(const Ray &ray, int recursionDepth){
    // Find hit object and distance
    Hit min_hit(std::numeric_limits<double>::infinity(),Vector());
    Object *obj = NULL;
    for (unsigned int i = 0; i < objects.size(); ++i) {
        Hit hit(objects[i]->intersect(ray));
        if (hit.t<min_hit.t) {
            min_hit = hit;
            obj = objects[i];
        }
    }

    // No hit? Return background color.
    if (!obj) return Color(0.0, 0.0, 0.0);

    Material *material = obj->material;            //the hit objects material
    Point hit = ray.at(min_hit.t);                 //the hit point
    Vector N = min_hit.N;                          //the normal at hit point
    Vector V = -ray.D;                             //the view vector


    /****************************************************
    * This is where you should insert the color
    * calculation (Phong model).
    *
    * Given: material, hit, N, V, lights[]
    * Sought: color
    *
    * Hints: (see triple.h)
    *        Triple.dot(Vector) dot product
    *        Vector+Vector      vector sum
    *        Vector-Vector      vector difference
    *        Point-Point        yields vector
    *        Vector.normalize() normalizes vector, returns length
    *        double*Color        scales each color component (r,g,b)
    *        Color*Color        dito
    *        pow(a,b)           a to the power of b
    ****************************************************/

    Color color;

    //For each light
	for(unsigned int i = 0; i < lights.size(); i++){
		bool isShadow = false;

		if (rendermode == phong) //Ambiant
			color += lights[i]->color * material->color * material->ka;

		//Shadow
		if(Shadow){
			//Ray from light to object
			Ray light(lights[i]->position, hit - lights[i]->position);

			//Hit from light to nearest object
			Hit minHit = findMinHit(light);

			//Hit from current object to light
			Hit currentHit(obj->intersect(light));

			//
			if (minHit.t < currentHit.t) isShadow = true;
		}
		//No Shadow -> Calculate colors
		if(!isShadow){

			//The light vector
			Vector L = lights[i]->position - hit;L.normalize();
			//The R vector (required for specular shading)
			Vector R = (-1*L) + 2 * (L.dot(N)) * N;

			if (rendermode == phong) //Diffuse
				if (material->texture){
					color += (max(0.0,N.dot(L)) * lights[i]->color) * obj->calcTexture(hit) * material->kd;
				} else { 
					color += (max(0.0,N.dot(L)) * lights[i]->color) * material->color * material->kd;
				}
			if (rendermode == gooch && ray.D.dot(N)<-0.2){
				Color kCool = Color(0.0,0.0,kBlue) + alpha * (lights[i]->color * material->color * material->kd);
				Color kWarm = Color(kYellow,kYellow,0.0) + beta * (lights[i]->color * material->color * material->kd);
				color += (kCool *(1 - N.dot(L))/2) + (kWarm * (1 + N.dot(L))/2);
			}

			//Specular		+ dot(R,V)^n LightColor ks
			color += pow(max(0.0,R.dot(V)), material->n) * lights[i]->color * material->ks;
		}
	}

	//Reflection
	if(recursionDepth <= maxRecursionDepth && ray.D.dot(N)<-0.2){
		Color buffer(0.0, 0.0, 0.0);

		//The reflected direction
		Vector reflectV((V + 2*(-N.dot(V))*N));

		for (unsigned int i=0; i<4; i++){

			//The reflected ray
			Ray reflect(hit + 0.001 * N, -reflectV);

			//Calculate reflection through recursion
			buffer += material->ks*tracePhongGooch(reflect, recursionDepth+1);
			reflectV.x += 0.01 * cos(0.5);
			reflectV.y += 0.01 * sin(0.5);
		}

		//Color = average buffer value
		color += buffer/4;
	}

    return color;
}
// ------------------------------------------------------ render -------------------------------------------------------------------- //
void Scene::render(Image &img){

	int w = img.width();
	int h = img.height();

	Vector right = ((center-eye).cross(up)).normalized();			// right vector
	double pSize = up.length();											// pixel size as length of up vector
	Point leftUp(center - (w/2)*pSize*right + (h/2)*pSize*up);	// pixel left upper bound

	double zBuffer[w][h];

	cout << "eye: " << eye.x << "," << eye.y << "," << eye.z << endl;
	cout << "center: " << center.x << "," << center.y << "," << center.z << endl;
	cout << "up: " << up.x << "," << up.y << "," << up.z << endl;
	cout << "right: " << right.x << "," << right.y << "," << right.z << endl;
	cout << "viewSize: " << w << "," << h << endl;
	cout << "leftUp: " << leftUp.x << "," << leftUp.y << "," << leftUp.z << endl;
	cout << "pSize: " << pSize << endl;
	cout << "maxRecursionDepth: " << maxRecursionDepth << endl;
	cout << "superSamplingFactor: " << superSamplingFactor << endl;
	cout << "apertureSamples: " << apertureSamples << endl;
	cout << "apertureRadius: " << apertureRadius << endl;

	Color color;
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			Point pixel = leftUp + (x*pSize*right) - (y*pSize*up);		// go to correct x,y position
			if (rendermode == phong || rendermode == gooch){
				img(x,y) = superSample(pixel, right, pSize);
			}
			if (rendermode == normal){
				Ray ray(eye, ((pixel + right*pSize/2 - up*pSize/2)-eye).normalized());
				img(x,y) = traceNormal(ray);
			}
			if (rendermode == zbuffer){
				Ray ray(eye, ((pixel + right*pSize/2 - up*pSize/2)-eye).normalized());
				zBuffer[x][y] = traceZBuffer(ray);
			}
		}
	}

	if (rendermode == zbuffer){
		for (int y = 0; y < h; y++){
			for (int x = 0; x < w; x++){
				color.set((-1*(zBuffer[x][y]-minZB)/(maxZB-minZB))+1);
				color.clamp();
				img(x,y) = color;
			}
		}
	}
}

// -------------------------------------------------------- setters ------------------------------------------------------------------ //

void Scene::setRenderMode(string rm){
	if(!rm.compare("phong")){
		rendermode = phong;
	}
	if(!rm.compare("zbuffer")){
		rendermode = zbuffer;
	}
	if(!rm.compare("normal")){
		rendermode = normal;
	}
	if(!rm.compare("gooch")){
		rendermode = gooch;
	}
}

void Scene::addObject(Object *o){
    objects.push_back(o);
}

void Scene::addLight(Light *l){
    lights.push_back(l);
}

void Scene::setEye(Triple e){
    eye = e;
}

void Scene::setBlue(double b){
	kBlue = b;
}

void Scene::setYellow(double y){
	kYellow = y;
}

void Scene::setBeta(double b){
	beta = b;
}

void Scene::setAlpha(double a){
	alpha = a;
}

void Scene::setCenter(Triple c){
	center = c;
}

void Scene::setUp(Triple u){
	up = u;
}

void Scene::setShadowMode(bool shadowMode){
	Shadow = shadowMode;
}

void Scene::setMaxRecursionDepth(int max){
	maxRecursionDepth = max;
}

void Scene::setSuperSamplingFactor(int superSampling){
	superSamplingFactor = superSampling;
}

void Scene::setApertureSamples(int as){
	apertureSamples = as;
}

void Scene::setApertureRadius(int ar){
	apertureRadius = ar;
}

