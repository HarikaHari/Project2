/*
  File Name    : vector.c
  Assignment   : Raycasting
  Created by   : Harika Hari (hh453). */

#include "custom.h"

double sqr(double v){
    return v*v;
}

double vectorLength(Vector a){
    return sqrt(sqr(a[0])+sqr(a[1])+sqr(a[2]));
}

void vectorUnit(Vector in, Vector out){
    double len = vectorLength(in);
    out[0] = in[0]/len;
    out[1] = in[1]/len;
    out[2] = in[2]/len;
    
}

void VectorAddition(Vector a, Vector b, Vector c){
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

void VectorSubstraction(Vector a, Vector b, Vector c){
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

double VectorDotProduct(Vector a, Vector b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void normalize(double *v) {
    double len = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
    len = sqrt(len);
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
}