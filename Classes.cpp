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
    
    char matfile[255];
    
    FILE* f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        fgets(line, 255, f);
        
        if (line[0] == 'u' && line[1] == 's') {
            curGroup++;
        }
        
        // Fichier materiau - pas necessaire
        if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
            sscanf(line, "mtllib %100s", matfile);
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
    
    // Chargement des textures
    f = fopen("BeautifulGirl.mtl","r");
    while (!feof(f)) {
        char line[255];
        fgets(line, 255, f);
        if (line[0] == 'm' && line[4] == 'K' && line[5] == 'd'){
            char texturefile[255];
            sscanf(line, "map_Kd %100s\n", texturefile);
            add_texture((std::string("girlOK/") + std::string(texturefile)).c_str());
        }
    }
    fclose(f);
    
    // Construction du BVH
    build_bvh(&bvh, 0, faces.size()/3);
}

void Geometry::add_texture(const char* filename) {
    
    textures.resize(textures.size() + 1);
    w.resize(w.size() + 1);
    h.resize(h.size() + 1);
    
    FILE* f;
    f = fopen(filename, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    w[w.size() - 1] = *(int*)&info[18]; // extract image height and width from header
    h[h.size() - 1] = *(int*)&info[22];
    
    int size = 3 * w[w.size() - 1] * h[h.size() - 1];
    textures[textures.size() - 1].resize(size); // allocate 3 bytes per pixel
    fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);
    
    for (int i = 0; i < size; i += 3) {
        std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
    }
}


BBox Geometry::build_bbox(int i0, int i1) {
    // Calcul de la boite englobante
    BBox result;
    result.bmax = vertices[faces[i0*3]];
    result.bmin = vertices[faces[i0*3]];
    for (int i=i0; i<i1; i++) { // Indice du triangle
        for (int j=0; j<3; j++) { // Indice du sommet
            for (int k=0; k<3; k++) { // Indice de dimension
                result.bmin[k] = std::min(result.bmin[k], vertices[faces[i*3+j]][k]);
                result.bmax[k] = std::max(result.bmax[k], vertices[faces[i*3+j]][k]);
            }
        }
    }
    return result;
}

void Geometry::build_bvh(BVH* node, int i0, int i1) {
    node->bbox = build_bbox(i0, i1);
    node->i0 = i0;
    node->i1 = i1;
    node->fg = NULL;
    node->fd = NULL;
    
    Vector diag = node->bbox.bmax - node->bbox.bmin;
    
    // Recherche de l axe selon lequel l objet est le plus etendu
    int split_dim;
    if ((diag[0] > diag[1]) && (diag[0] > diag[2])) {
        split_dim = 0;
    }
    else {
        if ((diag[1] > diag[0]) && (diag[1] > diag[2])) {
            split_dim = 1;
        }
        else {
            split_dim = 2;
        }
    }
    
    // Point de separation
    double split_val = node->bbox.bmin[split_dim] + diag[split_dim] * 0.5;
    
    // Partition
    int pivot = i0-1;
    for (int i=i0; i<i1; i++) {
        double center_split_dim = (vertices[faces[i*3]][split_dim] + vertices[faces[i*3+1]][split_dim] + vertices[faces[i*3+2]][split_dim]) / 3.;
        if (center_split_dim < split_val) {
            pivot++;
            std::swap(faces[i*3+0], faces[pivot*3+0]);
            std::swap(faces[i*3+1], faces[pivot*3+1]);
            std::swap(faces[i*3+2], faces[pivot*3+2]);
            
            std::swap(normalIds[i*3+0], normalIds[pivot*3+0]);
            std::swap(normalIds[i*3+1], normalIds[pivot*3+1]);
            std::swap(normalIds[i*3+2], normalIds[pivot*3+2]);
            
            std::swap(uvIds[i*3+0], uvIds[pivot*3+0]);
            std::swap(uvIds[i*3+1], uvIds[pivot*3+1]);
            std::swap(uvIds[i*3+2], uvIds[pivot*3+2]);
            
            std::swap(faceGroup[i], faceGroup[pivot]);
        }
    }
    // A la fin, toutes les valeurs de centre plus petites que split_val sont avant pivot, et les autres sont apres
    
    // Critere d arret
    if (pivot<=i0 || pivot>=i1 || i1==i0+1) return;
    node->fg = new BVH();
    build_bvh(node->fg, i0, pivot); // Dans la description yt il y a des + ou - 1
    node->fd = new BVH();
    build_bvh(node->fd, pivot, i1);
}


bool Geometry::intersection(const Ray& d, Vector& P, Vector& N, double& t, Vector& color) const {
    t = 1E99;
    bool has_inter = false;
    
    // Test d intersection avec la racine de l arbre
    if (!bvh.bbox.intersection(d)) return false;
    
    
    std::list<const BVH*> l;
    l.push_front(&bvh);
    
    while (!l.empty()) {
        const BVH* current = l.front();
        l.pop_front();
        if (current->fg && current->fg->bbox.intersection(d)) {
            l.push_back(current->fg);
        }
        if (current->fd && current->fd->bbox.intersection(d)) {
            l.push_back(current->fd);
        }
        
        // Cas d une feuille de l arbre
        if (!current->fg) {
            
            for (int i=current->i0; i<current->i1; i++) {
                int a = faces[i*3];
                int b = faces[i*3+1];
                int c = faces[i*3+2];
                
                Triangle tri(vertices[a], vertices[b], vertices[c], albedo, is_mirror, is_transparent);
                Vector localP, localN;
                double localt;
                double alpha, beta, gamma;
                if (tri.intersection(d, localP, localN, localt, alpha, beta, gamma)) {
                    has_inter = true;
                    if(localt < t) {
                        t = localt;
                        P = localP;
                        
                        // Interpolation de Phong
                        N = normals[normalIds[i*3+0]] * alpha + normals[normalIds[i*3+1]] * beta + normals[normalIds[i*3+2]] * gamma;
                        N.normalize();
                        
                        // La couleur depend de la texture
                        int textureId = faceGroup[i];
                        // On refait l interpolation comme pour les normales, mais sur les uv
                        int x = (uvs[uvIds[i*3+0]][0] * alpha + uvs[uvIds[i*3+1]][0] * beta + uvs[uvIds[i*3+2]][0] * gamma) * (w[textureId]-1);
                        int y = (uvs[uvIds[i*3+0]][1] * alpha + uvs[uvIds[i*3+1]][1] * beta + uvs[uvIds[i*3+2]][1] * gamma) * (h[textureId]-1);
                        
                        double cr = (textures[textureId][(y * w[textureId] + x)*3]) / 255.; // Composante rouge
                        double cg = (textures[textureId][(y * w[textureId] + x)*3+1]) / 255.; // Composante verte
                        double cb = (textures[textureId][(y * w[textureId] + x)*3+2]) / 255.; // Composante bleue
                        color = Vector(cr, cg, cb);
                    }
                }
            }
        }
    }
    return has_inter;
}

bool Sphere::intersection(const Ray& d, Vector& P, Vector& N, double &t, Vector& color) const {
    // Resolution de l equation du 2nd degre
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
    color = albedo;
    return true;
}

bool Triangle::intersection(const Ray& d, Vector& P, Vector& N, double &t, double &alpha, double &beta, double &gamma) const {
    N = -cross(B-A, C-A).getNormalized(); // Le - ajoute apres l ajout d un maillage peut faire bugger un triangle
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
    beta = detb/detm;
    double g12 = b11;
    double g22 = b21;
    double detg = m11*g22 - m12*g12;
    gamma = detg / detm;
    alpha = 1 - beta - gamma;
    return (gamma > 0 && beta > 0 && (gamma + beta) < 1);
}

bool Scene::intersection(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t, Vector& color) const {
    bool has_inter = false;
    min_t = 1E99;
    sphere_id = -1;
    
    for (int i=0; i<objects.size(); i++) {
        if (i==0) continue;
        Vector localP, localN, localColor;
        double t;
        bool local_has_inter = objects[i]->intersection(d, localP, localN, t, localColor);
        if (local_has_inter) {
            has_inter = true;
            // Teste si la sphere est plus proche que la plus proche actuelle
            if (t < min_t) {
                min_t = t;
                P = localP;
                N = localN;
                sphere_id = i;
                color = localColor;
            }
        }
    }
    return has_inter;
}
  
