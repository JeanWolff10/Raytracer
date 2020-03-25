//
//  main.cpp
//  raytracer2
//
//  Created by Jean Wolff on 08/01/2020.
//  Copyright © 2020 Jean Wolff. All rights reserved.
//
// clang++ -Xpreprocessor -fopenmp main.cpp Vector.cpp -std=c++11 -lomp -O3
// ./a.out 5 image11.png 1024 1024
//

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
// #include "Vector2.h"
#include <math.h>

#define M_pi 3.1416
// Pour afficher des variables
#include <iostream>
using namespace std;

#include <chrono>
// Record start time
auto start = std::chrono::high_resolution_clock::now();

#include <random>
std::default_random_engine engine;
std::uniform_real_distribution<double> uniform(0,1);

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

// & permet de ne pas charger un tableau complet, seulement de prendre la case memoire du bon element
Vector operator+(const Vector& a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(double a, const Vector &b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector &a, const Vector &b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator*(const Vector &b, double a) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator/(const Vector& a, double b) {
    return Vector(a[0]/b, a[1]/b, a[2]/b);
}
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}
double dot(const Vector& a, const Vector& b) {
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}
Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
};
Vector random_cos(const Vector &N) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    Vector direction_aleatoire_repere_local(cos(2 * M_pi * r1) * sqrt(1 - r2), sin(2 * M_pi * r1) * sqrt(1 - r2), sqrt(r2));
    Vector aleatoire(uniform(engine)-0.5, uniform(engine)-0.5, uniform(engine)-0.5);
    Vector tangent1 = cross(N, aleatoire); tangent1.normalize();
    Vector tangent2 = cross(tangent1, N);
    return direction_aleatoire_repere_local[2]*N + direction_aleatoire_repere_local[0]*tangent1 + direction_aleatoire_repere_local[1]*tangent2;
}


class Ray {
public:
    // Constructeur
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};


class Object {
public:
    // Albedo rend compte du reflet pour chaque couleur
    Vector albedo;
    bool is_mirror;
    bool is_transparent;
    bool is_bulle;
    Object(const Vector &couleur, bool is_mirror = false, bool is_transparent = false, bool is_bulle = false) :
    is_mirror(is_mirror), is_transparent(is_transparent), is_bulle(is_bulle) {
        albedo = couleur;
    };
    virtual bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const = 0;
    
};

class Sphere : public Object {
public:
    // Attributs de la classe
    Vector O;
    double R;

    // Constructeur
    Sphere(const Vector &origin, double rayon, const Vector &couleur, bool is_mirror=false, bool is_transparent=false, bool is_bulle=false) : O(origin), R(rayon), Object(couleur, is_mirror, is_transparent, is_bulle) {
    };
    
    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const {
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
        //if ( fabs(N[0]) > 0.2 &&
          //  cout << N[0] << "  " << N[1] << "  " << N[2] << "  " << endl;
        return true;
    }
};

class Triangle : public Object {
public:
    Vector A, B, C;
    
    Triangle(const Vector& A, const Vector &B, const Vector& C, const Vector &couleur, bool is_mirror=false, bool is_transparent=false, bool is_bulle=false) : A(A), B(B), C(C), Object(couleur, is_mirror, is_transparent, is_bulle) {
    };

    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const {
        N = cross(B-A, C-A).getNormalized();
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
};

class Scene{
public:
    Scene() {};
    void addSpheres(const Sphere& s) {objects.push_back(&s);}
    void addTriangle(const Triangle& s) {objects.push_back(&s);}

    bool intersection(const Ray& d, Vector& P, Vector& N, int& sphere_id, double& min_t) const {
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
    
    // Ensemble des spheres de la scene
    std::vector<const Object*> objects;
};

Vector position_lumiere(25, 25, 30);
double rayon_lumiere = 3;

Vector getColor2(const Ray& d, const Scene& s, int& numero_rebond) {
    if(numero_rebond == 0) return Vector(0,0,0);

    double epsilon = 0.00001;
    double intensite_lumiere = 50000000;
    //const Sphere *light = dynamic_cast<Sphere*>(const_cast<Object*>(objects[0]));
    //Vector position_lumiere =  dynamic_cast<Sphere*>(objects[0])->O;
    //double rayon_lumiere = 5;
    Vector intensite_pixel(0,0,0);
    Vector P, N;
    int sphere_id;
    double valeur_dot;
    double t;

    //return (Vector(N[0], N[1], N[2])*0.5 + 0.5)*45500.0 ;
    //if (dot(N, d.direction) > 0) N= N*(-1.);

    // On prend un point legerement decolle de la surface
    P = P + epsilon * N;
    if (s.intersection(d, P, N, sphere_id, t)) {
        
        // --- SURFACE SPECULAIRE ---
        if (s.objects[sphere_id]->is_mirror) {
            numero_rebond = numero_rebond - 1;
            Vector new_direction = d.direction - 2 * dot(d.direction, N) * N;
            new_direction.normalize(); // a priori inutile
            Ray new_r(P, new_direction);
            // Produit scalaire entre le rayon vers l objet dont on voit le reflet, et N
            valeur_dot = dot(new_direction, N); // Le prof ne prend pas ca !!
            // Calcul de la couleur du point dont on voit le reflet
            intensite_pixel = valeur_dot * getColor2(new_r, s, numero_rebond);
            
            // Pour donner une couleur au miroir
            //Vector modulation_couleur(1.5,1,1);
            //intensite_pixel = modulation_couleur * intensite_pixel_reflet;
        }
        
        // --- SURFACE TRANSPARENTE ---
        else if (s.objects[sphere_id]->is_transparent) {
            double n_air = 1.03;
            double n_sphere = 1.5;
            double n1 = n_air;
            double n2 = n_sphere;
            // Cas d'une bulle
            bool is_bulle = s.objects[sphere_id]->is_bulle;
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
                numero_rebond = numero_rebond - 1;
                Vector new_direction = (n1/n2) * d.direction - ((n1/n2) * dot(d.direction, N) + sqrt(temp)) * N;
                new_direction.normalize();
                // On prend un point legerement decolle de la surface
                P = P - 2 * epsilon * N;
                Ray new_r(P, new_direction);
                // Produit scalaire entre le rayon vers le point vers lequel le rayon refracte va, et N
                //valeur_dot = dot(new_direction, N);
                intensite_pixel = getColor2(new_r, s, numero_rebond);
                if (is_bulle) {
                    //Vector violet2(0.3, 0.04, 0.7);
                    //Vector couleur(1.15, 1.02, 1.3);
                    Vector couleur(1, 1, 1.7);
                    intensite_pixel = couleur * intensite_pixel;
                }

            }
        }
        
        // --- SURFACE DIFFUSE ---
        else {
            // --- CONTRIBUTION DIRECTE (lumiere spherique) ---
            // Direction aleatoire entre le centre de la lumiere et son hemisphere du cote de P
            Vector axe_PO = (P - position_lumiere).getNormalized();
            Vector dir_aleatoire = random_cos(axe_PO);
            
            // Point a la surface de la sphere lumineuse
            Vector point_aleatoire = dir_aleatoire * rayon_lumiere + position_lumiere;
            Vector wi = (point_aleatoire - P).getNormalized();
            double d_light2 = (point_aleatoire - P).getNorm2();
            // Normale au point aleatoire
            
            Vector Np = dir_aleatoire;
            // Vector Np = axe_PO;// TOTOTOTOTOT
            double proba = 1;// dot(axe_PO, dir_aleatoire);
            
            
            P = P + 0.001 * N; // TOTOTOTO
            Ray ray_light(P, wi);
            Vector P_light, N_light;
            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
            // 0.99 pour ne pas prendre la sphere lumineuse avec elle meme...
            if (1 && has_inter_light && (t_light*t_light < 0.99*d_light2))
            {
                intensite_pixel = Vector(0,0,0);
                //intensite_pixel = 250*intensite_lumiere / (4 * M_pi * d_light2) * std::max(0., dot(N, wi)) * std::max(0., dot(Np, -wi)) / proba * s.objects[sphere_id]->albedo / 3;
            }
            else {
                //if (sphere_id == 6 ) cout << this->objects[sphere_id]->albedo[0] << endl;
                intensite_pixel = 250*intensite_lumiere / (4 * M_pi * d_light2) * std::max(0., dot(N, wi)) * std::max(0., dot(Np, -wi)) / proba * s.objects[sphere_id]->albedo;
            }
            
            // --- CONTRIBUTION INDIRECTE ---
            // Lumiere spherique
            numero_rebond = numero_rebond - 1;
            Vector direction_aleatoire = random_cos(N);
            Ray rayon_aleatoire(P, direction_aleatoire);
            // Calcul de la couleur du point dont on voit le reflet
            // TOTOTOOTO
            intensite_pixel = intensite_pixel + s.objects[sphere_id]->albedo * getColor2(rayon_aleatoire, s, numero_rebond) / 15;
            // (diviser par pi)
        }
    }
    return intensite_pixel;
}

int main(int argc, char** argv) {
    int nrays = atoi(argv[1]);
    char* image_name = argv[2];
    int W = atoi(argv[3]);
    int H = atoi(argv[4]);
    double fov = 60 * M_pi / 180;
    //cout<<fov<<endl;
    
    Vector beige(255, 162, 89);
    beige.normalize();
    Vector orange(254, 104, 69);
    orange.normalize();
    Vector rose(250, 66, 82);
    rose.normalize();
    Vector vert(145, 189, 58);
    vert.normalize();
    Vector bleu_ciel(210, 235, 233);
    bleu_ciel.normalize();
    
    Sphere s1(Vector(0,-20000-20,0), 20000, vert); // Sol
    Sphere s2(Vector(0,20060,0), 20000, rose); // Plafond
    Sphere s3(Vector(-20050,0,0), 20000, beige); // Mur gauche
    Sphere s4(Vector(20050,0,0), 20000, orange); // Mur droit
    Sphere s5(Vector(0,0,-20080), 20000, bleu_ciel, false); // Mur fond
    Sphere s9(Vector(0,0,20040), 20000, Vector(1,1,0), false); // Mur invisible avant la camera
    Sphere s6(Vector(-13,-15,-40), 5, Vector(0,1,0), false, true); // Sphere de gauche
    Sphere s8(Vector(13,-15,-40), 5, Vector(1,0,0), true); // Sphere de droite
    //Sphere s7(Vector(-10,12,-10), 3.4, Vector(1,1,1), false, true, true); // Bulle intérieure
    //Sphere s10(Vector(-10,12,-10), 3.5, Vector(1,1,1), false, true); // Bulle extérieure
    Sphere s11(Vector(0,-10,-28), 10, Vector(1,1,1)); // Sphere blanche pour voir l eclairage indirect
    
    Triangle tri(Vector(-20,0,-20), Vector(0,0,0), Vector(-10,20,-10), orange);
    
    Sphere slum(position_lumiere, rayon_lumiere, Vector(1,1,1));
    //s.lumiere = &slum;
    
    Vector position_camera(0,7,36);
    double focus_distance = 36 + 40; // tout ce qui est avant ou apres est flou
    
    Scene s;
    s.addSpheres(slum); // A mettre en premiere position
    s.addSpheres(s1);
    s.addSpheres(s2);
    s.addSpheres(s3);
    s.addSpheres(s4);
    s.addSpheres(s5);
    //s.addSpheres(s6);
    //s.addSpheres(s7);
    //s.addSpheres(s8);
    //s.addSpheres(s9);
    //s.addSpheres(s10);
    s.addSpheres(s11);
    s.addTriangle(tri);
    
    int nb_rebonds_max = 5;
    int rebounds_max = nb_rebonds_max;
    //double coef_reflection = 1;
    //double coef_diff = 0;
    //double gamma = 2.2;
    
    std::vector<unsigned char> image(W*H * 3, 0);
    #pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector intensite_pixel(0,0,0);
            
            for (int k=0; k < nrays; k++) {
                // Anti-aliasing - methode Box Muller
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double R = sqrt(-2 * log(r1));
                double dx = R * cos(2 * M_pi * r2);
                double dy = R * sin(2 * M_pi * r2);

                // Variation de la profondeur de champ
                // Ouverture carree
                double dx_aperture = (uniform(engine) - 0.5) * 5.;
                double dy_aperture = (uniform(engine) - 0.5) * 5.;
                
                // Rayon qui part de la camera et passe par le pixel (i,j)
                Vector direction(j - W / 2.0 + 0.5 + dx, i - H / 2.0 + 0.5 + dy, -W / (2 * tan(fov / 2)));
                direction.normalize();
                
                // Mise au point, flou
                // Vector destination = position_camera + focus_distance * direction;
                // Vector new_origin = position_camera + Vector(dx_aperture, dy_aperture, 0);
                // Ray r(new_origin, (destination - new_origin).getNormalized());

                // Sans mise au point
                //Ray r(position_camera, direction);
                Ray r = Ray(position_camera, direction);
                
                int numero_rebond = rebounds_max;
                intensite_pixel = intensite_pixel + getColor2(r, s, numero_rebond) / nrays;
            }

            //cout << intensite_pixel[0] << "  " << intensite_pixel[1] << "  " << intensite_pixel[2] << "  " << endl;
            // Enregistrement des couleurs des pixels
            // Avec la correction gamma
            image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0., (std::pow(intensite_pixel[0], 0.45))));
            image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0., (std::pow(intensite_pixel[1], 0.45))));
            image[((H-i-1)*W + j) * 3 + 2] = std::min(255., std::max(0., (std::pow(intensite_pixel[2], 0.45))));
        }
    }
    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    
    stbi_write_png(image_name, W, H, 3, &image[0], 0); // &image[0] est l adresse du premier pixel
    return 0;
}
