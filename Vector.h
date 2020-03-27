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
Vector operator*(const Vector &a, const Vector &b);
Vector operator*(const Vector &b, double a);
Vector operator/(const Vector& a, double b);
Vector operator-(const Vector& a);
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);
