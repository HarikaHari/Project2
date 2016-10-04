/*
  File Name    : ppm.c
  Assignment   : Raycasting
  Created by   : Harika Hari (hh453). */


#include "custom.h"

 void ImageWrite(Image *buffer, const char *filename,int format) {
    size_t num2;
    int size = buffer->width * buffer->height * 4;
    FILE *fp = fopen(filename, "w");
    if (!fp) { 
        fprintf(stderr,"cannot open file for writing");
		exit(0);
	}
    
	//write the header file
    //image format
    fprintf(fp, "P%d\n",format);

    //comments
    fprintf(fp, "# Created by %s\n","Harika");

    //image size
    fprintf(fp, "%d %d\n",buffer->width, buffer->height);

    // rgb component depth
    fprintf(fp, "%d\n",buffer->maxval);

    // pixel data
	 for(int i=1; i<size+1;i++){    // for each for slots we skip it,
            char ch=buffer->data[i-1];
            if (i%4 !=0) {
               fwrite(&ch, 1, 1, fp);  // which means we skip A at here
            }
   	}
	fclose(fp);
}