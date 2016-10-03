/*
  File Name    : raycast.c
  Assignment   : Raycasting
  Created by   : Harika Hari (hh453). */

#include "custom.h"
#include "vector.c"
#include "ppm.c"

//Returns the position of camera object from the JSON file
int getCameraPosition(OBJECT *objects){
    int i=0;
    while(i < ObjectsCount && objects[i].type !=0){
        if (objects[i].type == CAM) {
            return i;
        }
        i++;
    }
    return -1;
}


// fill in pixel color to our image - need to flip the y axis due to the placement of the viewplane

void colorPixel(double *color, int row, int col,Image *image){
 //store results in image->data which store pixels
    image->data[row * image->width*4 + col*4] = (char)color[0];
    image->data[row * image->width*4 + col*4+1] = (char)color[1];
    image->data[row * image->width*4 + col*4+2]= (char)color[2];
    image->data[row * image->width*4 + col*4+3]= 255;

}


//Ro means Ray origin - Rd means Ray direction (other end of the ray)
 
double planeIntersection(double *Ro, double *Rd, double *Pos, double *Norm){

    double alph,delta;
    normalize(Norm);
    /*
     *if the plane is parallel to the ray
     * return -1;
     */
    alph = VectorDotProduct(Norm, Rd);
    
    // the plane is parallel to the ray
    if (fabs(alph) <0.0001) {
        return -1;
    }
    
    Vector incident_vector;
    VectorSubstraction(Pos, Ro, incident_vector);
    delta = VectorDotProduct(incident_vector, Norm);
    
    
    double t = delta/alph; // whcih means we check thea1 and thea2
    
    if (t<0.0) { // reflection, no intersection
        return -1;
    }
    
    //printf("planeIntersection finished works\n");
    return t; // return something, but not t , need to figure out it

    
    }

double sphereIntersection(double *Ro, double *Rd, double *C, double r){

    double a, b, c;
    double s0,s1; // we have two solutions if delta >0
    //calculate quadratic formula
    //First find a, b, c
    a = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
    b = 2 * (Rd[0]*(Ro[0]-C[0]) + Rd[1]*(Ro[1]-C[1]) + Rd[2]*(Ro[2]-C[2]));
    c = sqr(Ro[0]-C[0]) + sqr(Ro[1]-C[1]) + sqr(Ro[2]-C[2]) - sqr(r);
     // make sure a is 1 (means the ray direction was normalized)
    if (a > 1.0001 || a < .9999) {
        printf("a = %lf\n", a);
        fprintf(stderr, "Ray direction was not normalized\n");
        exit(1);
    }
    // check that discriminant is <, =, or > 0
    double disc = sqr(b) - 4*a*c;
    
    if (disc < 0) {
        //printf("disc was < 0\n");
        return -1; // no solution
    }
    else if (disc == 0) {
        s0 = -1*(b / (2*a)); // single solution
        s1 = -1*(b / (2*a));
        //printf("t0 = %lf\n", t0);
    }
    else {  // 2 solutions: find the smaller
        s0 = (-1*b - sqrt(sqr(b) - 4*c))/2;
        s1 = (-1*b + sqrt(sqr(b) - 4*c))/2;

    }
    
    
    if (s0 < 0 && s1 < 0) {
        // no intersection
        return -1;
    }
    else if (s0 < 0 && s1 > 0) {
        return s1;
    }
    else if (s0 > 0 && s1 < 0) {
        return s0;
    }
    else { // they were both positive
        if (s0 <= s1)
            return s0;
        else
            return s1;
    }
    
}


void raycasting(Image *image, double cam_width, double cam_height, OBJECT *objects){
    // loop over all pixels and test for intersections with objects
   
    int x, y, counter; 
    
    Vector vp_pos= {0,0,1}; //view plane position
    Vector Ro = {0,0,0}; // Ray origin position or Camera position
    Vector point = {0,0,0}; //point on viewplane where intersection happens
    
    double pixheight = cam_height / image->height;
    double pixwidth = cam_width / image->width;
    
    Vector Rd = {0,0,0}; // direction of Ray
    point[2] = vp_pos[2]; // set intersecting point Z to viewplane Z
    
    for(x=0;x<image->height;x++){
        point[1] = -(vp_pos[1] - cam_height/2.0 + pixheight*(x + 0.5)); 
        for (y=0; y<image->width; y++) {
             point[0] = vp_pos[0] - cam_width/2.0 + pixwidth*(y + 0.5);
            normalize(point);   //normalize the point
            //store normalized point as our direction
            Rd[0] = point[0];
            Rd[1] = point[1];
            Rd[2] = point[2];
            
            int best_counter =0;
            double best_t = INFINITY;
            for (counter=0; objects[counter].type!=0; counter++) {
                // we need to return intersection test on each object
                double t =0;
                switch (objects[counter].type) {
                    case 0:
                        printf("no object found\n");
                        break;
                        
                    case CAM:
                        break;
                        
                    case SPH:
                        t = sphereIntersection(Ro, Rd, objects[counter].data.sphere.position, objects[counter].data.sphere.radius);
                        break;
                        
                    case PLN:
                        t= planeIntersection(Ro, Rd, objects[counter].data.plane.position, objects[counter].data.plane.normal);
                        break;
                        
                    default:
                        exit(1);
                }
                if (t > 0 && t < best_t) {
                    best_t = t;
                    best_counter = counter;
                }
                
            }
            if (best_t > 0 && best_t != INFINITY) {
                // there was an intersection
                if (objects[best_counter].type == PLN) {
                    colorPixel(objects[best_counter].data.plane.color, x, y, image);
                }
                else if (objects[best_counter].type == SPH){
                    colorPixel(objects[best_counter].data.sphere.color, x, y, image);
                }
                
            }
            else{
                Vector white ={1,1,1};
                colorPixel(white,x,y,image);
            }
            
        }
    }
    
}
