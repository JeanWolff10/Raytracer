//
//  Classes.hpp
//  raytracer2
//
//  Created by Jean Wolff on 26/03/2020.
//  Copyright © 2020 Jean Wolff. All rights reserved.
//
#pragma once

#ifndef Classes_hpp
#define Classes_hpp
#define M_pi 3.1416
#endif /* Classes_hpp */

#include <stdio.h>
#include <math.h>
#include "Vector.h"


class Ray {
public:
    // Constructeur
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};


class Object {
public :
    Object(){};
    virtual bool intersection (const Ray& d,  Vector& P, Vector& N, double &t) const=0 ;
    Vector albedo;
    bool is_mirror;
    bool is_transparent;
    bool is_bulle;
};

class Sphere : public Object{
public :
    // Attributs de la classe
    Vector O;
    double R;
    // Constructeur
    Sphere(const Vector& origin, const double& rayon, const Vector &couleur, bool mirror = false, bool transparent = false, bool bulle=false) : O(origin), R(rayon) {
        albedo = couleur;
        is_mirror = mirror;
        is_transparent = transparent;
        is_bulle = bulle;
    };

    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const;
};

class Triangle : public Object {
public :
    //const Vector &A, &B, &C; entraine un bug!!
    Vector A, B, C;
    Triangle(const Vector& A, const Vector &B, const Vector& C, const Vector &couleur, bool mirror=false, bool transparent=false, bool bulle=false) : A(A), B(B), C(C) {
        albedo = couleur;
        is_mirror = mirror;
        is_transparent = transparent;
        is_bulle = bulle;
    }
    
    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const;
};

class BBox {
public:
    BBox() {};
    BBox(const Vector& bmin, const Vector& bmax): bmin(bmin), bmax(bmax) {};
    
    // Intersection avec la boite englobante
    bool intersection(const Ray& d) const {
        double t_1_x = (bmin[0] - d.origin[0]) / d.direction[0];
        double t_2_x = (bmax[0] - d.origin[0]) / d.direction[0];
        double t_min_x = std::min(t_1_x, t_2_x);
        double t_max_x = std::max(t_1_x, t_2_x);
        
        double t_1_y = (bmin[1] - d.origin[1]) / d.direction[1];
        double t_2_y = (bmax[1] - d.origin[1]) / d.direction[1];
        double t_min_y = std::min(t_1_y, t_2_y);
        double t_max_y = std::max(t_1_y, t_2_y);
        
        double t_1_z = (bmin[2] - d.origin[2]) / d.direction[2];
        double t_2_z = (bmax[2] - d.origin[2]) / d.direction[2];
        double t_min_z = std::min(t_1_z, t_2_z);
        double t_max_z = std::max(t_1_z, t_2_z);
        
        if(std::min(std::min(t_max_x, t_max_y), t_max_z) > std::max(std::max(t_min_x, t_min_y), t_min_z)) return true;
        return false;
    }
    
    Vector bmin, bmax;
};

class Geometry : public Object {
public:
    Geometry(const char* obj, double scaling, const Vector& offset, const Vector& couleur, bool mirror=false, bool transp=false);
    
    std::vector<int> faceGroup;
    std::vector<int> faces;
    std::vector<int> normalIds;
    std::vector<int> uvIds;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    
    bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const;

private:
    BBox bb;
};

class Scene{
public:
    Scene() {};
    void addSpheres(const Sphere& s) {objects.push_back(&s);}
    void addTriangle(const Triangle& s) {objects.push_back(&s);}
    void addGeometry(const Geometry& s) {objects.push_back(&s);}
    
    bool intersection(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t) const;
    
    // Ensemble des spheres de la scene
    std::vector<const Object*> objects;
};
