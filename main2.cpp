//
//  main.cpp
//  raytracer2
//
//  Created by Jean Wolff on 08/01/2020.
//  Copyright © 2020 Jean Wolff. All rights reserved.
//
// clang++ -Xpreprocessor -fopenmp main.cpp Vector.cpp -std=c++11 -lomp -O3
// clang++ -Xpreprocessor -fopenmp main2.cpp Vector.cpp Classes.cpp -std=c++11 -lomp -O3
// ./a.out 5 image11.png 1024 1024
//

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <math.h>

#define M_pi 3.1416
#include "Vector.h"
#include "Classes.hpp"

// Pour afficher des variables
#include <iostream>
using namespace std;
//cout<<1<<endl;

// Pour mesurer la duree d execution
#include <chrono>
auto start = std::chrono::high_resolution_clock::now();

// Generer des nombres aleatoires
#include <random>
std::default_random_engine engine;
std::uniform_real_distribution<double> uniform(0,1);

// Champ visuel
double fov = 60 * M_pi / 180;
Vector position_camera(0,7,47);
double focus_distance = 47 + 40; // Tout ce qui est avant ou apres est flou

// Lumiere
Vector position_lumiere(25, 25, 30);
double rayon_lumiere = 3;
double intensite_lumiere = 6000000000;

// P = P + epsilon * N pour eviter les bugs du aux imprecisions numeriques
double epsilon = 0.001;

int rebounds_max = 5;


Vector random_cos(const Vector &N) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    Vector direction_aleatoire_repere_local(cos(2 * M_pi * r1) * sqrt(1 - r2), sin(2 * M_pi * r1) * sqrt(1 - r2), sqrt(r2));
    Vector aleatoire(uniform(engine)-0.5, uniform(engine)-0.5, uniform(engine)-0.5);
    Vector tangent1 = cross(N, aleatoire);
    tangent1.normalize();
    Vector tangent2 = cross(tangent1, N);
    Vector rep = direction_aleatoire_repere_local[2]*N + direction_aleatoire_repere_local[0]*tangent1 + direction_aleatoire_repere_local[1]*tangent2;
    return rep;
}


Vector getColor2(const Ray& d, const Scene& s, int& numero_rebond) {
    if(numero_rebond == 0) return Vector(0,0,0);
    
    Vector intensite_pixel(0,0,0);
    Vector P, N, albedo;
    int sphere_id;
    double t;
    
    bool has_inter = s.intersection(d, P, N, sphere_id, t, albedo);
    
    if (has_inter) {
        // On prend un point legerement decolle de la surface
        P = P + epsilon * N;

        // --- SURFACE SPECULAIRE ---
        if (s.objects[sphere_id]->is_mirror) {
            numero_rebond = numero_rebond - 1;
            Vector dir_miroir = d.direction - 2 * dot(d.direction, N) * N;
            Ray r_miroir(P, dir_miroir);
            // Calcul de la couleur du point dont on voit le reflet
            intensite_pixel = getColor2(r_miroir, s, numero_rebond);
            
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
            Ray new_ray;
            bool entering = true;
            
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
                entering = false;
                if (is_bulle) {
                    n1 = n_air;
                    n2 = n_sphere;
                }
            }
            double temp = 1 - (n1/n2)*(n1/n2) * (1 - dot(d.direction, N)*dot(d.direction, N));
            if (temp>0) {
                numero_rebond = numero_rebond - 1;
                Vector direction_refracte = (n1/n2) * d.direction - ((n1/n2) * dot(d.direction, N) + sqrt(temp)) * N;
                // On prend un point legerement decolle de la surface
                P = P - 2 * epsilon * N;
                
                // Coef de reflexion
                double R0 = ((n1-n2)/(n1+n2)) * ((n1-n2)/(n1+n2));
                double R;
                if (entering) {
                    R = R0 + (1-R0) * std::pow(1 + dot(d.direction, N), 5.);
                }
                else {
                    R = R0 + (1-R0) * std::pow(1 + dot(direction_refracte, N), 5.);
                }
                
                if (uniform(engine) < R) {
                    // Rayon reflechi
                    P = P +  2 * epsilon * N;
                    Vector dir_miroir = d.direction - 2 * dot(d.direction, N) * N;
                    new_ray = Ray(P, dir_miroir);
                }
                else {
                    // Rayon refracte
                    new_ray = Ray(P, direction_refracte);
                }
                
                if (is_bulle) {
                    Vector couleur(1, 1, 1.7); // Donne une artificiellement une couleur violette a la bulle
                    intensite_pixel = couleur * intensite_pixel;
                }
            }
            
            // Reflexion totale
            else {
                Vector dir_miroir = d.direction - 2 * dot(d.direction, N) * N;
                new_ray = Ray(P, dir_miroir);
                
            }
            
            intensite_pixel = getColor2(new_ray, s, numero_rebond);
        }
        
        // --- SURFACE DIFFUSE ---
        else {
            // --- CONTRIBUTION DIRECTE ---
            // Direction aleatoire entre le centre de la lumiere et son hemisphere du cote de P
            Vector axe_PO = (P - position_lumiere).getNormalized();
            Vector dir_aleatoire = random_cos(axe_PO);
            
            // Point a la surface de la sphere lumineuse
            Vector point_aleatoire = dir_aleatoire * rayon_lumiere + position_lumiere;
            Vector wi = (point_aleatoire - P).getNormalized();
            double d_light2 = (point_aleatoire - P).getNorm2();
            
            // Normale au point aleatoire
            Vector Np = dir_aleatoire;
            
            Ray ray_light(P, wi);
            Vector P_light, N_light, albedo_light;
            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light, albedo_light);
            
            // --- OMBRE ---
            if (has_inter_light && (t_light*t_light < 0.99*d_light2)) // 0.99 pour ne pas prendre la lumiere avec elle meme
            {
                intensite_pixel = Vector(0,0,0);
            }
            
            // --- ECLAIRAGE DIRECT ---
            else {
                double proba = dot(axe_PO, dir_aleatoire);
                intensite_pixel = intensite_lumiere / (4 * M_pi * d_light2) * std::max(0., dot(N, wi)) * dot(Np, -wi) / proba * albedo;
            }
            
            // --- CONTRIBUTION INDIRECTE ---
            numero_rebond = numero_rebond - 1;
            Vector direction_aleatoire = random_cos(N);
            Ray rayon_aleatoire(P, direction_aleatoire);
            // Calcul de la couleur du point dont on voit le reflet
            intensite_pixel = intensite_pixel + albedo * getColor2(rayon_aleatoire, s, numero_rebond);
        }
    }
    return intensite_pixel;
}

int main(int argc, char** argv) {
    int nrays = atoi(argv[1]);
    char* image_name = argv[2];
    int W = atoi(argv[3]);
    int H = atoi(argv[4]);
    
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
    
    // Murs de la piece
    Sphere s1(Vector(0,-20000-20,0), 20000, vert, true); // Sol
    Sphere s2(Vector(0,20060,0), 20000, rose); // Plafond
    Sphere s3(Vector(-20050,0,0), 20000, Vector(0,0,1)); // Mur gauche
    Sphere s4(Vector(20050,0,0), 20000, vert); // Mur droit
    Sphere s5(Vector(0,0,-20080), 20000, bleu_ciel); // Mur fond
    //Sphere s9(Vector(0,0,20040), 20000, Vector(1,1,0), false); // Mur invisible avant la camera
    
    // Spheres
    Sphere s6(Vector(-16,-10,-28), 10, Vector(0,1,0), false, true); // Sphere transparente à gauche
    Sphere s8(Vector(0,27,-40), 10, Vector(1,0,0), true); // Sphere miroir en l air
    //Sphere s7(Vector(-10,12,-10), 3.4, Vector(1,1,1), false, true, true); // Bulle intérieure
    //Sphere s10(Vector(-10,12,-10), 3.5, Vector(1,1,1), false, true); // Bulle extérieure
    Sphere s11(Vector(-28,28,-32), 6, Vector(1,1,1)); // Sphere blanche pour voir l eclairage indirect
    
    //Triangle tri(Vector(-20,0,-20), Vector(0,0,0), Vector(-10,20,-10), orange);
    
    Sphere slum(position_lumiere, rayon_lumiere, Vector(1,1,1));
    
    Geometry g1("girl.obj", 26, Vector(20,-20 ,-28), Vector(1.,1.,1.));
    
    Scene s;
    s.addSpheres(slum); // A mettre en premiere position
    s.addSpheres(s1);
    s.addSpheres(s2);
    s.addSpheres(s3);
    s.addSpheres(s4);
    s.addSpheres(s5);
    s.addSpheres(s6);
    //s.addSpheres(s7);
    s.addSpheres(s8);
    //s.addSpheres(s9);
    //s.addSpheres(s10);
    s.addSpheres(s11);
    //s.addTriangle(tri);
    s.addGeometry(g1);
    
    
    std::vector<unsigned char> image(W*H * 3, 0);
    #pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector intensite_pixel(0,0,0);
            
            for (int k=0; k < nrays; k++) {
                // Anti-aliasing - methode Box Muller
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double R = 0.2 * sqrt(-2 * log(r1));
                double dx = R * cos(2 * M_pi * r2);
                double dy = R * sin(2 * M_pi * r2);
                
                // Variation de la profondeur de champ
                // Ouverture carree
                //double dx_aperture = (uniform(engine) - 0.5) * 5.;
                //double dy_aperture = (uniform(engine) - 0.5) * 5.;
                
                // Rayon qui part de la camera et passe par le pixel (i,j)
                Vector direction(j - W / 2.0 + 0.5 + dx, i - H / 2.0 + 0.5 + dy, -W / (2 * tan(fov / 2)));
                direction.normalize();
                
                // Mise au point, flou
                // Vector destination = position_camera + focus_distance * direction;
                // Vector new_origin = position_camera + Vector(dx_aperture, dy_aperture, 0);
                // Ray r(new_origin, (destination - new_origin).getNormalized());

                // Sans mise au point
                Ray r(position_camera, direction);
                
                int numero_rebond = rebounds_max;
                intensite_pixel = intensite_pixel + getColor2(r, s, numero_rebond) / nrays;
            }

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
