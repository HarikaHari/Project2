/*
  File Name    : custom.h
  Assignment   : Raycasting
  Created by   : Harika Hari (hh453). */

#ifndef custom_headers
#define custom_headers


#include <ctype.h>
#include <errno.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//Pre Defined Values 
#define ObjectsCount 20  //maximum number of JSON objects that can be parsed
#define MAX_COLOR_VAL 255
#define CAM 1
#define SPH 2
#define PLN 3
#define QUAD 4

//structures based on the JSON sample
typedef struct CAMERA {
    double width;
    double height;
} CAMERA;

typedef struct SPHERE {
    double *color;
    double *position;
    double radius;
} SPHERE;

typedef struct PLANE {
    double *color;
    double *position;
    double *normal;
} PLANE;

typedef struct QUADRIC {
    double *color;
    double *coefficient;
} QUADRIC;

//structure to store data as an object
typedef struct OBJECT{
    int type;
    union {
        CAMERA camera;
        SPHERE sphere;
        PLANE plane;
        QUADRIC quadric;
    } data;
} OBJECT;

//structure of image 
typedef struct Image
{
    int width;
    int height;
    char *data; 
    int magic_number;
    int depth;
    char *tupltype;
    int maxval;
}Image;

typedef struct ray{
    double origin[3];
    double direction[3];
} ray;

typedef double Vector[3]; 
extern int line;  
extern OBJECT objects[ObjectsCount]; 
void read_scene(const char* filename);

void ImageWrite(Image *image, const char *filename,int format);

int getCameraPosition(OBJECT *objects);
void colorPixel(double *color, int row, int col,Image *image);
double planeIntersection(double *Ro, double *Rd, double *Pos, double *Norm);
double sphereIntersection(double *Ro, double *Rd, double *C, double r);
void raycasting(Image *image, double cam_width, double cam_height, OBJECT *object);


double sqr(double v);
double vectorLength(Vector a);
void VectorAddition(Vector a, Vector b, Vector c);
void VectorSubstraction(Vector a, Vector b, Vector c);
double VectorDotProduct(Vector a, Vector b);
void normalize(double *v);

#endif /* custom_h */