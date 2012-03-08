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

Hit Scene::findMinHit(const Ray &ray){
	// Find nearest object and distance
    Hit min_hit(std::numeric_limits<double>::infinity(),Vector());
    for (unsigned int i = 0; i < objects.size(); ++i) {
        Hit hit(objects[i]->intersect(ray));
        if (hit.t<min_hit.t) {
            min_hit.t = hit.t;
            min_hit.N = hit.N;
            }
    }
    return min_hit;
}

Hit Scene::findMinMaxDistance(const Ray &ray, Image &img){
	//The image dimensions
	int w = img.width();
    int h = img.height();

    double minMaxDistance[] = {std::numeric_limits<double>::infinity(), 0.0}

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
        	Ray ray(eye, (pixel-eye).normalized());
        	Hit min_hit = findMinHit(ray);

        	if(min_hit.z < minMaxDistance[0])
				minMaxDistance[0] = min_hit.z
        }
}

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

Color Scene::traceZBuffer(const Ray &ray, Image &img){

	Color color;

	//Find the nearest object (hit)
	Hit min_hit = findMinHit(ray);

	// No hit? Return background color.
	if (min_hit.t >= std::numeric_limits<double>::infinity())
		return Color(0.0, 0.0, 0.0);

	//The hit point
	Point hit = ray.at(min_hit.t);

	//The distance = z-coordinate divided by it's maximum value (normalized)
	color.set((hit.z/min_hit.t));
	return color;
}

Color Scene::tracePhong(const Ray &ray, int recursionDepth){
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

		//Ambiant		+ LightColor MaterialColor ka
		color += lights[i]->color * material->color * material->ka;

		//Shadow
		if(Shadow){
			//Outgoing light ray (L)
			Ray light(hit + 0.001 * N, lights[i]->position - hit);

			//For each object
			for(unsigned int j = 0; j < objects.size(); ++j){

				//Check to see if there's a hit between the object and the light
				Hit shadow(objects[j]->intersect(light));

				//If hit -> There's a shadow
				if(shadow.t < std::numeric_limits<double>::infinity())
					isShadow = true;
			}
		}

		//No Shadow -> Calculate colors
		if(!isShadow){

			//The light vector
			Vector L = lights[i]->position - hit;L.normalize();
			//The R vector (required for specular shading)
			Vector R = (-1*L) + 2 * (L.dot(N)) * N;

			//Diffuse		+ dot(L,N) LightColor MaterialColor kd
			color += (max(0.0,N.dot(L)) * lights[i]->color) * material->color * material->kd;

			//Specular		+ dot(R,V)^n LightColor ks
			color += pow(max(0.0,R.dot(V)), material->n) * lights[i]->color * material->ks;
		}
	}

	//Reflection
	if(recursionDepth <= maxRecursionDepth){
		Color buffer(0.0, 0.0, 0.0);

		//The reflected direction
		Vector reflectV((V + 2*(-N.dot(V))*N));

		for (unsigned int i=0; i<10; i++){

			//The reflected ray
			Ray reflect(hit + 0.001 * N, -reflectV);

			//Calculate reflection through recursion
			buffer += material->ks*tracePhong(reflect, recursionDepth+1);
			reflectV.x += 0.01 * cos(0.5);
			reflectV.y += 0.01 * sin(0.5);
		}

		//Color = average buffer value
		color += buffer/10;
	}

    return color;
}

Color Scene::trace(const Ray &ray, Image &img){
	Color color;

	//Calculate color values according to the rendermode
	switch(rendermode){
    	case phong:
			color = tracePhong(ray, 1);
			break;
		case zbuffer:
			color = traceZBuffer(ray, img);
			break;
		case normal:
			color = traceNormal(ray);
			break;
	}

	//Return the color
	return color;
}

void Scene::render(Image &img)
{
    int w = img.width();
    int h = img.height();

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            Point pixel(x+0.5, h-1-y+0.5, 0);
            Ray ray(eye, (pixel-eye).normalized());
            Color col = trace(ray, img);
            col.clamp();
            img(x,y) = col;
        }
    }
}

void Scene::addObject(Object *o)
{
    objects.push_back(o);
}

void Scene::addLight(Light *l)
{
    lights.push_back(l);
}

void Scene::setEye(Triple e)
{
    eye = e;
}

void Scene::setRenderMode(string rm){
	if(!rm.compare("phong")){
		rendermode = phong;
	}else if(!rm.compare("zbuffer")){
		rendermode = zbuffer;
	}else if(!rm.compare("normal")){
		rendermode = normal;
	}
}

void Scene::setShadowMode(bool shadowMode){
	Shadow = shadowMode;
	cout << Shadow << endl;
}

void Scene::setMaxRecursionDepth(int max){
	maxRecursionDepth = max;
}
