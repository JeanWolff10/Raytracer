//
//  Vector.cpp
//  raytracer2
//
//  Created by Jean Wolff on 11/01/2020.
//  Copyright © 2020 Jean Wolff. All rights reserved.
//

#include <stdio.h>
#include "Vector.h"
#include <vector>

#include <random>
static std::default_random_engine engine;
static  std::uniform_real_distribution<double> uniform(0, 1);
#define M_pi 3.1416
 

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

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

