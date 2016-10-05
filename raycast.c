/*
  File Name    : raycast.c
  Assignment   : Raycasting
  Created by   : Harika Hari (hh453). */


#include "custom.h"
#include "intersection.c"
#include "json.c"


int main(int argc, const char * argv[]) {

   //verify the length of command line arguments passed		
   if (argc != 5) {
        fprintf(stderr, "Error: Expected 4 arguments Image width, Image height, Input JSON file and output file name \n");
        exit(1);
    }
	//verify the image width and height values	
    if (atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0) {
        fprintf(stderr, "Error: width and height parameters must be > 0\n");
        exit(1);
    }
    
    //Parsing JSON file start
    const char *file = argv[3];
    if (file == NULL) {
        fprintf(stderr, "Error: Failed to open input file '%s'\n", argv[3]);
        exit(1);
    }
    read_scene(file);
	//Parsing JSON file End
	
	//creating Image size for rendering Image
    Image *image =(Image *)malloc(sizeof(Image));
    printf("Creating the Image\n");
    //Reading Image height, width, maxval and data
    image->width = atoi(argv[1]);
    printf("image width is %d \n",image->width);
    image->height = atoi(argv[2]);
    printf("image height is %d \n",image->height);
    image->maxval = 255;
    image->data =(unsigned char*) malloc(sizeof(unsigned char) * image->width * image->height*4);
    //Reading the camera position from parsed JSON data
	int pos = getCameraPosition(objects);
	
    raycast(image, objects[pos].data.camera.width, objects[pos].data.camera.height, objects);
    //writing Image to a PPM file
	const char *output = argv[4];
	ImageWrite(image, output,6);
    return 0;
}
