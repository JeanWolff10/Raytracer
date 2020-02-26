//
//  Vector.h
//  raytracer2
//
//  Created by Jean Wolff on 11/01/2020.
//  Copyright © 2020 Jean Wolff. All rights reserved.
//
#pragma once
#ifndef Vector_h
#define Vector_h
#define M_pi 3.1416

#endif /* Vector_h */

#include <math.h>
#include <vector>

// Pour afficher des variables
#include <iostream>
using namespace std;


class Vector {
public:
    // Constructeur
    Vector(double x=0, double y=0, double z=0){
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    // Accesseur pour acceder au i eme element de la coordonnee du vecteur
    const double& operator[](int i) const { return coord[i]; }
    double& operator[](int i) { return coord[i]; }
    
    // On evitera le calcul couteux de racines carres
    double getNorm2() {
        return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
    }
    
    // Normaliser un vecteur
    void normalize() {
        double norm = sqrt(getNorm2());
        coord[0] /= norm;
        coord[1] /= norm;
        coord[2] /= norm;
    }
    
    // Renvoie le vecteur normalise, mais sans modifier le vecteur initial
    Vector getNormalized() {
        // Copie du vecteur
        Vector result(*this);
        result.normalize();
        return result;
    }
    
private:
    double coord[3];
};

Vector operator+(const Vector& a, const Vector &b);
Vector operator-(const Vector& a, const Vector &b);
Vector operator*(double a, const Vector &b);
Vector operator*(const Vector &b, double a);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector& a, double b);
Vector operator-(const Vector& a); // ligne oubliee
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);

// Direction aleatoire
Vector random_cos(const Vector &N);
 
class Ray {
public:
    // Constructeur
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};

class Sphere{
public:
    // Attributs de la classe
    Vector O;
    double R;
    // Albedo rend compte du reflet pour chaque couleur
    Vector albedo;
    bool is_mirror;
    bool is_transparent;
    bool is_bulle;
    // Constructeur
    Sphere(const Vector &origin, double rayon, const Vector &couleur, bool is_mirror=false, bool is_transparent=false, bool is_bulle=false) : O(origin), R(rayon), albedo(couleur), is_mirror(is_mirror), is_transparent(is_transparent), is_bulle(is_bulle) {};
    
    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) {
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
        N = (P - O).getNormalized();
        return true;
    }
};
    
class Scene{
public:
    Scene() {};
    void addSpheres(const Sphere& s) {spheres.push_back(s);}

    bool intersection(const Ray& d, Vector& P, Vector& N, int &sphere_id, double& t) {
        bool has_inter = false;
        double min_t = 1E99;
        
        for (int i=0; i<spheres.size(); i++) {
            Vector localP, localN;
            bool local_has_inter = spheres[i].intersection(d, localP, localN, t);
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
    
    bool shadow(Vector& P, Vector& N, int &sphere_id, double &epsilon, Vector& position_lumiere) {
        // Rayon qui part de P et va vers la lumière
        Vector direction_ombre = position_lumiere - P;
        double distance_source_carre = direction_ombre.getNorm2();
        direction_ombre.normalize();
        // On prend un point legerement decolle de la surface
        P = P + epsilon * N;
        // Teste si un objet est rencontre
        Ray r_ombre(P, direction_ombre);
        Vector P1, N1;
        int sphere_id1;
        double t;
        // On laisse N1 et sphere_id1 bien qu'ils soient modifies de maniere inutile
        bool has_inter_ombre = this->intersection(r_ombre, P1, N1, sphere_id1, t);
        if (has_inter_ombre) {
            // Teste si l objet rencontre est plus proche que la lumiere
            double distance_intersection_carre = (P - P1).getNorm2();
            if (distance_intersection_carre < distance_source_carre) {
                return true;
            }
        }
        return false;
    }
    
    Vector getColor2(const Ray& d, int &numero_rebond, int &numero_rebond_transp) {
        double epsilon = 0.00001;
        double intensite_lumiere = 5000000;
        Vector position_lumiere(15, 30, -20);
        Vector intensite_pixel(0,0,0);
        Vector P, N;
        int sphere_id;
        double valeur_dot;
        
        double n_air = 1.03;
        double n_sphere = 1.5;
        double t;
        
        bool has_inter = this->intersection(d, P, N, sphere_id, t);
        // On prend un point legerement decolle de la surface
        P = P + epsilon * N;
        if (has_inter) {
            
            bool is_mirror = this->spheres[sphere_id].is_mirror;
            bool is_transparent = this->spheres[sphere_id].is_transparent;
            
            if (is_mirror and numero_rebond > 0) {
                numero_rebond = numero_rebond - 1;
                Vector new_direction = d.direction - 2 * dot(d.direction, N) * N;
                new_direction.normalize(); // a priori inutile
                Ray new_r(P, new_direction);
                // Produit scalaire entre le rayon vers l objet dont on voit le reflet, et N
                valeur_dot = dot(new_direction, N); // Le prof ne prend pas ca !!
                // Calcul de la couleur du point dont on voit le reflet
                intensite_pixel = valeur_dot * this->getColor2(new_r, numero_rebond, numero_rebond_transp);
                
                // Pour donner une couleur au miroir
                //Vector modulation_couleur(1.5,1,1);
                //intensite_pixel = modulation_couleur * intensite_pixel_reflet;
            }
            
            else if (is_mirror and numero_rebond == 0)
                return Vector(0,0,0);
            
            else if (is_transparent and numero_rebond_transp > 0) {
                double n1 = n_air;
                double n2 = n_sphere;
                // Cas d'une bulle
                bool is_bulle = this->spheres[sphere_id].is_bulle;
                if (is_bulle) {
                    n1 = n_sphere;
                    n2 = n_air;
                }
                // Si on sort de la sphere
                if (dot(d.direction, N) > 0) {
                    n1 = n_sphere;
                    n2 = n_air;
                    N = -N;
                    if (is_bulle) {
                        n1 = n_air;
                        n2 = n_sphere;
                    }
                }
                double temp = 1 - (n1/n2)*(n1/n2) * (1 - dot(d.direction, N)*dot(d.direction, N));
                if (temp>0) {
                    numero_rebond_transp = numero_rebond_transp - 1;
                    Vector new_direction = (n1/n2) * d.direction - ((n1/n2) * dot(d.direction, N) + sqrt(temp)) * N;
                    new_direction.normalize();
                    // On prend un point legerement decolle de la surface
                    P = P - 2 * epsilon * N;
                    Ray new_r(P, new_direction);
                    // Produit scalaire entre le rayon vers le point vers lequel le rayon refracte va, et N
                    //valeur_dot = dot(new_direction, N);
                    intensite_pixel = this->getColor2(new_r, numero_rebond, numero_rebond_transp);
                    if (is_bulle) {
                        //Vector violet2(0.3, 0.04, 0.7);
                        //Vector couleur(1.15, 1.02, 1.3);
                        Vector couleur(1, 1, 1.7);
                        intensite_pixel = couleur * intensite_pixel;
                    }

                }
            }
            
            else if (is_transparent and numero_rebond_transp == 0)
                return Vector(0,0,0);
            
            else {
                // LUMIERE PONCTUELLE
                // OMBRE
                // bool is_shadow = this->shadow(P, N, sphere_id, epsilon, position_lumiere);
                // if (is_shadow) {
                //    intensite_pixel = Vector(0,0,0);
                //}
                //else {
                // CONTRIBUTION DIRECTE
                    // Produit scalaire entre le rayon qui va vers la lumiere, et la normale N a l objet
                    //valeur_dot = dot((position_lumiere - P).getNormalized(), N);
                    // Calcul de la couleur du point dont on voit le reflet
                    //intensite_pixel = this->spheres[sphere_id].albedo * intensite_lumiere * std::max(0., valeur_dot) / (position_lumiere - P).getNorm2();
                    // diviser par pi la ligne ci-dessus normalement (video 3 a 35:00)                    
                //}
                
                // CONTRIBUTION DIRECTE (lumiere spherique)
                // Direction aleatoire entre le centre de la lumiere et son hemisphere du cote de P
                Vector axe_PO = (P - spheres[0].O).getNormalized();
                Vector dir_aleatoire = random_cos(axe_PO);
                // Point a la surface de la sphere lumineuse
                Vector point_aleatoire = dir_aleatoire * spheres[0].R + spheres[0].O;
                Vector wi = (point_aleatoire - P).getNormalized();
                double d_light2 = (point_aleatoire - P).getNorm2();
                // Normale au point aleatoire
                Vector Np = dir_aleatoire;
                double proba = dot(axe_PO, dir_aleatoire);
                
                Ray ray_light(P, wi);
                Vector P_light, N_light;
                int sphere_id_light;
                double t_light;
                bool has_inter_light = this->intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
                // 0.99 pour ne pas prendre la sphere lumineuse avec elle meme...
                if (has_inter_light && t_light*t_light < 0.99*d_light2)
                {
                    intensite_pixel = Vector(0,0,0);
                }
                else {
                    intensite_pixel = intensite_lumiere / (4 * M_pi * d_light2) * std::max(0., dot(N, wi)) * std::max(0., dot(Np, -wi)) / proba * this->spheres[sphere_id].albedo;
                }
                
                // CONTRIBUTION INDIRECTE
                if (numero_rebond > 0) {
                    // LUMIERE SPHERIQUE
                    numero_rebond = numero_rebond - 1;
                    Vector direction_aleatoire = random_cos(N);
                    Ray rayon_aleatoire(P, direction_aleatoire);
                    // Calcul de la couleur du point dont on voit le reflet
                    intensite_pixel = intensite_pixel + this->spheres[sphere_id].albedo * this->getColor2(rayon_aleatoire, numero_rebond, numero_rebond_transp) / 2;
                    // (diviser par pi)
                }
            }
            
            return intensite_pixel;
        }
        else
            return Vector(0,0,0);
    }
    
    
    // Ensemble des spheres de la scene
    std::vector<Sphere> spheres;
};
