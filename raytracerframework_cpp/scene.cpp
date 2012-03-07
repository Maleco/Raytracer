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

Color Scene::trace(const Ray &ray, int recursionDepth)
{
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

    switch(rendermode){
    	case phong:{
			color += material->color*material->ka;				//With multiple lights, ambiant only depends on the material values

			Hit shadow_min(std::numeric_limits<double>::infinity(),Vector());
			for(unsigned int i = 0; i < lights.size(); i++){	//For each light
				bool isShadow = false;

				if(Shadow){
					Ray ray(hit + 0.001 * N, lights[i]->position - hit);	//Outgoing light ray (L)
					for(unsigned int j = 0; j < objects.size(); ++j){		//For each object
						Hit shadow(objects[j]->intersect(ray));
						if(!isnan(shadow.t)) {
								cout << "shadow: " << shadow.t << " shadow_min: " << shadow_min.t << endl;
						}
						if(shadow.t < shadow_min.t){
							isShadow = true;
						}
					}
				}
				if(!isShadow){
					//The light vector
					Vector L = lights[i]->position - hit;
					L.normalize();
					//The R vector (required for specular shading)
					Vector R = (-1*L) + 2 * (L.dot(N)) * N;

					//Diffuse		+ dot(L,N) LightColor MaterialColor kd
					color += (max(0.0,N.dot(L)) * lights[i]->color) * material->color * material->kd;
					//Specular		+ dot(R,V)^n LightColor ks
					color += pow(max(0.0,R.dot(V)), material->n) * lights[i]->color * material->ks;
				}
			}

			if(recursionDepth <= maxRecursionDepth){
				Color buffer(0.0, 0.0, 0.0);
				Vector reflectV((V + 2*(-N.dot(V))*N));
					for (unsigned int i=0; i<10; i++){
						Ray reflect(hit + 0.001 * N, -reflectV);
						buffer += material->ks*trace(reflect, recursionDepth+1); // Reflect
						reflectV.x += 0.01 * cos(0.5);
						reflectV.y += 0.01 * sin(0.5);
					}
					color += buffer/10;
			}
    	}
    	break;

		case zbuffer:{
			cout << "henk" << endl;
			double distance = hit.z/min_hit.t; //The distance = z-coordinate divided by it's maximum value (normalized)
			color.set(distance);
		}
		break;

		case normal:{
			double x = (min_hit.N.x+1)/2;
			double y = (min_hit.N.y+1)/2;
			double z = (min_hit.N.z+1)/2;
			color.set(x,y,z);
		}
    }

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
            Color col = trace(ray, 0 /* RecursionDepth*/);
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
	}else {
		rendermode = phong;
	}
}

void Scene::setShadowMode(bool shadowMode){
	Shadow = shadowMode;
	cout << Shadow << endl;
}

void Scene::setMaxRecursionDepth(int max){
	maxRecursionDepth = max;
}
