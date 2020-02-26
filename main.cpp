//
//  main.cpp
//  raytracer2
//
//  Created by Jean Wolff on 08/01/2020.
//  Copyright © 2020 Jean Wolff. All rights reserved.
//
// clang++ -Xpreprocessor -fopenmp main.cpp Vector.cpp -std=c++11 -lomp -O3
// ./a.out 5 image11.png 1024 1024

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Vector.h"
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
    
    Vector violet(35, 4, 68);
    violet.normalize();
    Vector mauve(144, 48, 61);
    mauve.normalize();
    Vector orange2(237, 130, 64);
    orange2.normalize();
    Vector bleu_ciel(210, 235, 233);
    bleu_ciel.normalize();
    
    // Avec mur fond courbe
    //Sphere s10(Vector(0,30,-160), 70, Vector(0,0,0), true); // Mur fond courbe
    //Sphere s9(Vector(0,0,20090), 20000, Vector(1,1,0), true); // Mur invisible avant la camera
    //Sphere s6(Vector(-12,-12,-40), 8, Vector(0,1,0), true);
    //Sphere s8(Vector(12,-12,-40), 8, Vector(1,0,0), true);
    
    // Avec differentes couleurs
    //Sphere s1(Vector(0,-20000-20,0), 20000, violet); // Sol
    //Sphere s2(Vector(0,20060,0), 20000, bleu_ciel); // Plafond
    //Sphere s3(Vector(-20050,0,0), 20000, orange2); // Mur gauche
    //Sphere s4(Vector(20050,0,0), 20000, mauve); // Mur droit
    
    Sphere s1(Vector(0,-20000-20,0), 20000, vert); // Sol
    Sphere s2(Vector(0,20060,0), 20000, rose); // Plafond
    Sphere s3(Vector(-20050,0,0), 20000, beige); // Mur gauche
    Sphere s4(Vector(20050,0,0), 20000, orange); // Mur droit
    Sphere s5(Vector(0,0,-20080), 20000, bleu_ciel, true); // Mur fond
    Sphere s9(Vector(0,0,20040), 20000, Vector(1,1,0), true); // Mur invisible avant la camera

    Sphere s6(Vector(-13,-15,-40), 5, Vector(0,1,0), false, true);
    Sphere s8(Vector(13,-15,-40), 5, Vector(1,0,0), true);
    // Bulles
    Sphere s7(Vector(-10,12,-10), 3.4, Vector(1,1,1), false, true, true);
    Sphere s10(Vector(-10,12,-10), 3.5, Vector(1,1,1), false, true);
    
    // Immense bulle qui passe juste devant la camera
    //Sphere s11(Vector(0, 7,-2000-5), 2000, Vector(1,1,1), false, true, true);
    //Sphere s12(Vector(0,7,-2000-5), 2000.1, Vector(1,1,1), false, true);
    
    // Bulle juste devant la camera
    //Sphere s11(Vector(0, 7,33), 2, Vector(1,1,1), false, true, true);
    //Sphere s12(Vector(0,7,33), 1.9, Vector(1,1,1), false, true);
    
    // Sphere transparente juste devant la camera
    //Sphere s12(Vector(0,7,33), 1.6, Vector(1,1,1), false, true);
    
    Sphere s11(Vector(0,-10,-28), 10, Vector(1,1,1));
    
    Vector position_lumiere(15, 30, -20);
    Sphere slum(Vector(15, 30, -20), 5, Vector(1,1,1));
    //s.lumiere = &slum;
    
    Vector position_camera(0,7,36);

    Scene s;
    s.addSpheres(slum); // A mettre en premiere position
    s.addSpheres(s1);
    s.addSpheres(s2);
    s.addSpheres(s3);
    s.addSpheres(s4);
    s.addSpheres(s5);
    s.addSpheres(s6);
    s.addSpheres(s7);
    s.addSpheres(s8);
    s.addSpheres(s9);
    s.addSpheres(s10);
    s.addSpheres(s11);
    //s.addSpheres(s12);
    
    
    

    int nb_rebonds_max = 5;
    int rebounds_max = nb_rebonds_max;
    int numero_rebond_transp_max = nb_rebonds_max;
    //double coef_reflection = 1;
    //double coef_diff = 0;
    //double gamma = 2.2;
    
    std::vector<unsigned char> image(W*H * 3, 0);
    #pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector intensite_pixel(0,0,0);
            
            for (int k=0; k < nrays; k++) {
                // Anti-aliasing
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double R = sqrt(-2 * log(r1));
                double dx = R * cos(2 * M_pi * r2);
                double dy = R * sin(2 * M_pi * r2);
                // Rayon qui part de la camera et passe par le pixel (i,j)
                Vector direction(j - W / 2.0 + 0.5 + dx, i - H / 2.0 + 0.5 + dy, -W / (2 * tan(fov / 2)));
                direction.normalize();
                Ray r(position_camera, direction);
                
                int numero_rebond = rebounds_max;
                int numero_rebond_transp = numero_rebond_transp_max;
                intensite_pixel = intensite_pixel + s.getColor2(r, numero_rebond, numero_rebond_transp) / nrays;
            }
            
            // Enregistrement des couleurs des pixels
            image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0]));
            image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1]));
            image[((H-i-1)*W + j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2]));
            
            // Avec la correction gamma
            //image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0., pow(intensite_pixel[0], 1/gamma)));
            //image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0., pow(intensite_pixel[1], 1/gamma)));
            //image[((H-i-1)*W + j) * 3 + 2] = std::min(255., std::max(0., pow(intensite_pixel[2], 1/gamma)));
        }
    }
    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    
    stbi_write_png(image_name, W, H, 3, &image[0], 0); // &image[0] est l adresse du premier pixel
    return 0;
}
