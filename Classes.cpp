//
//  Classes.cpp
//  raytracer2
//
//  Created by Jean Wolff on 26/03/2020.
//  Copyright © 2020 Jean Wolff. All rights reserved.
//

#include "Vector.h"
#include "Classes.hpp"
#include <stdio.h>
#include <vector>

#include <map>
#include <list>

Geometry::Geometry(const char* obj, double scaling, const Vector& offset, const Vector& couleur, bool mirror, bool transp) {
    albedo = couleur;
    is_mirror = mirror;
    is_transparent = transp;
    
    FILE* f;
    f = fopen(obj, "r"); // modif : fopen_s(&f,...)
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        fgets(line, 255, f);
        
        if (line[0] == 'u' && line[1] == 's') {
            curGroup++;
        }
        
        // Chargement sommets
        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;
            sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); // modif : sscanf_s
            vertices.push_back(scaling*vec + offset);
        }
    
        // Chargement normales
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);
            normals.push_back(vec);
        }
        
        // Chargement coordonnes de texture
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        
        // Chargement triangles
        if (line[0] == 'f') {
            int i0, i1, i2;
            int j0, j1, j2;
            int k0, k1, k2;
            faceGroup.push_back(curGroup);
            int nn = sscanf(line, "f %u/%u/%u %u/%u/%u %u/%u/%u\n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2); // modif : sscanf_s
            if (nn == 9) {
                faces.push_back(i0-1);
                faces.push_back(i1-1);
                faces.push_back(i2-1);
                uvIds.push_back(j0-1);
                uvIds.push_back(j1-1);
                uvIds.push_back(j2-1);
                normalIds.push_back(k0-1);
                normalIds.push_back(k1-1);
                normalIds.push_back(k2-1);
            }
            else {
                // Cas des quadrilateres a diviser en 2 triangles (que pour le dragon)
            }
        }
    }
    fclose(f);
    
    // Calcul de la boite englobante
    bb.bmax = vertices[0];
    bb.bmin = vertices[0];
    for (int i=1; i<vertices.size(); i++) {
        for (int j=0; j<3; j++) {
            bb.bmin[j] = std::min(bb.bmin[j], vertices[i][j]);
            bb.bmax[j] = std::max(bb.bmax[j], vertices[i][j]);
        }
    }
}

bool Geometry::intersection(const Ray& d, Vector& P, Vector& N, double& t) const {
    if (!bb.intersection(d)) return false;
    
    t = 1E99;
    bool has_inter = false;
    for (int i=0; i<faces.size()/3; i++) {
        int i0 = faces[3*i];
        int i1 = faces[3*i+1];
        int i2  = faces[3*i+2];
        Triangle tri(vertices[i0], vertices[i1], vertices[i2], albedo, is_mirror, is_transparent);
        Vector localP, localN;
        double localt;
        if (tri.intersection(d, localP, localN, localt)) {
            has_inter = true;
            if(localt < t) {
                t = localt;
                P = localP;
                N = localN;
            }
        }
    }
    return has_inter;
}

bool Sphere::intersection(const Ray& d, Vector& P, Vector& N, double &t) const {
    // resolution de l equation du 2nd degre
    double a = 1.;
    double b = 2. * dot(d.direction, d.origin - O);
    double c = (d.origin - O).getNorm2() - R*R;
    double delta =  b*b - 4. * a*c;
    if (delta < 0) return false;
    double t1 = (-b - sqrt(delta)) / (2. * a);
    double t2 = (-b + sqrt(delta)) / (2. * a);
    if (t2 < 0) return false;
    if (t1 > 0)
        t = t1;
    else
        t = t2;
    // P est le point d'intersection entre le rayon incident et la sphere
    P = d.origin + t * d.direction;
    // N est la normale a la sphere au point P
    N = (P - O) / R;
    return true;
}

bool Triangle::intersection(const Ray& d, Vector& P, Vector& N, double &t) const {
    N = -cross(B-A, C-A).getNormalized(); // - ajoute apres l ajout d un maillage
    t = dot(C - d.origin, N) / dot(d.direction, N);
    if (t<0) return false;

    P = d.origin + t*d.direction;
    Vector u = B-A;
    Vector v = C-A;
    Vector w = P-A;
    double m11 = u.getNorm2();
    double m12 = dot(u,v);
    double m22 = v.getNorm2();
    double detm = m11*m22 - m12*m12;
    
    double b11 = dot(w,u);
    double b21 = dot(w, v);
    double detb = b11*m22 - b21 * m12;
    double beta = detb/detm;
    
    double g12 = b11;
    double g22 = b21;
    double detg = m11*g22 - m12*g12;
    double gamma = detg / detm;

    return (gamma > 0 && beta > 0 && (gamma + beta) < 1);
}

bool Scene::intersection(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t) const {
    bool has_inter = false;
    min_t = 1E99;
    sphere_id = -1;
    
    for (int i=0; i<objects.size(); i++) {
        if (i==0) continue;
        Vector localP, localN;
        double t;
        bool local_has_inter = objects[i]->intersection(d, localP, localN, t);
        if (local_has_inter) {
            has_inter = true;
            // Teste si la sphere est plus proche que la plus proche actuelle
            if (t < min_t) {
                min_t = t;
                P = localP;
                N = localN;
                // sphere_id vaut des valeurs genre 3786239801 ??
                sphere_id = i;
            }
        }
    }
    return has_inter;
}
