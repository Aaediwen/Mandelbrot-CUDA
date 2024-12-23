#include <stdio.h>
#include <X11/Xlib.h>
#include <unistd.h>
#include <stdlib.h>
#include <tiffio.h>
#include <pthread.h>

// =============== STRUCT DEFINITIONS ========================
// the X context
struct drawgc {
     Display *Display;       //     X display connection
     Window window;        //     dialog box
     Pixmap pixmap;    // a Pixmap for later use
     int winx;         // X size
     int winy;         // Y size
     int pixmapx;
     int pixmapy;
     GC gc;            // window Graphics Context
     GC pixmapgc;      // PixMap Graphics Context
     int depth;
} fractwin;

// definition for a specific row in the mandelbrot set
struct row_data{
     double winx;
     double winy;
     double xstart;
     double xscale;
     double xend;
     int max_itterations;
     int complete;
     pthread_t thread;
     double sy;
     unsigned long *colors;
};
// data about the current viewport
struct viewport {
  double viewx;
  double viewy;
  double zoom; 
  double factor;
} fractview;

struct fract_cords {
 int x;
 int y;
} ;


// ===================== HOST X INITIALIZATION =================
int Xinit (){
     char *display_name = NULL;
     XWindowAttributes windowattribs;
     fractwin.Display = XOpenDisplay(display_name);
     if (!fractwin.Display) {
          return 1;
     }
     fractwin.window = XCreateSimpleWindow(fractwin.Display, DefaultRootWindow(fractwin.Display), 10, 10, fractwin.winx,fractwin.winy, 10, 0, 1234);

     fractwin.pixmap = XCreatePixmap(fractwin.Display, DefaultRootWindow(fractwin.Display), fractwin.winx, fractwin.winy, 32);
     XMapWindow(fractwin.Display, fractwin.window);
     XGetWindowAttributes(fractwin.Display, fractwin.window, &windowattribs);
     fractwin.winx=windowattribs.width;
     fractwin.winy=windowattribs.height;
     fractwin.depth=2;
     for (int c=1 ; c < windowattribs.depth ; c++) {
          fractwin.depth *=2;
     }
     fractwin.gc=XCreateGC(fractwin.Display, fractwin.window, 0,0);
     fractwin.pixmapgc=XCreateGC(fractwin.Display, fractwin.pixmap, 0,0);
     XClearWindow(fractwin.Display, fractwin.window);
     XSelectInput(fractwin.Display, fractwin.window, ButtonPressMask|ButtonReleaseMask|KeyPressMask|KeyReleaseMask);
     return 0;  
}




void tiff_write(struct row_data *bitmap, int renderx, int rendery) {
          TIFF *tif=TIFFOpen("test.tif", "w");
          int sampleperpixel=3;
          char * image= (char*)malloc (renderx*rendery*sampleperpixel);
          int index=0;
          for (int y=0 ; y < rendery ; y++) {     
               for (int x=0 ; x < renderx ; x++) {
                    unsigned long i=bitmap[y].colors[x];
                    unsigned long r, g, b;
                    b = i & 255;
                    i = i >> 8;
                    g = i & 255;
                    i = i >> 8;
                    r = i & 255;
                    image[index]=(int)r;
                    index ++;
                    image[index]=(int)g;
                    index ++;
                    image[index]=(int)b;
                    index ++;
               }
          }
          TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, renderx);  // set the width of the image
          TIFFSetField(tif, TIFFTAG_IMAGELENGTH, rendery);    // set the height of the image
          TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
          TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
          TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
          //   Some other essential fields to set that you do not have to understand for now.
          TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
          TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
          tsize_t linebytes = sampleperpixel*renderx;
          unsigned char * buf=(unsigned char*)malloc(linebytes);
          TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, sampleperpixel*renderx));
          index=0;
          for (int y=0 ; y < rendery ; y++) {     
               TIFFWriteScanline(tif, &image[index], y, 0);
               index += renderx*sampleperpixel;
          }
          
          TIFFClose(tif);
          free(image);
          free(buf);
}


__global__ void nv_pixel(unsigned long* color, struct row_data *data) {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;
//    int d_itterations = *max_itterations;
    float sx;
    sx = (index/ data->xend) * data->xscale + data->xstart;
    double ax=0;
    double ay=0;
 
    unsigned long i=0;
    
    while (ax*ax + ay*ay < 4 && i < data->max_itterations) {
         double xtemp = ax*ax - ay*ay + sx;
         ay = 2*ax*ay + data->sy;
         ax = xtemp;
         i++;
    }
    
                 i=i*2;
               unsigned long blue= i & 255;
               i = i >> 8;
               unsigned long green= i & 255;
               i = i >> 8;
               unsigned long red = i & 255;
               i=0;
      //         printf ("%i, %i, %i\n", red, green, blue);
               i=(green << 16) | (blue << 8) | red;

    
    color[index] = i;
}


// test nv kernel ===============================================================================
__global__ void nv_pixel_test(unsigned long* color, struct row_data *data) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
//    printf ("In NV Pixel Test for row: %i pixel %i\n", (int) data->winy, index);
    float sx;
    sx = (index/ data->xend) * data->xscale + data->xstart;
    double ax=0;
    double ay=0;
 
    unsigned long i=0;
    while (ax*ax + ay*ay < 4 && i < data->max_itterations) {
         double xtemp = ax*ax - ay*ay + sx;
         ay = 2*ax*ay + data->sy;
         ax = xtemp;
         i++;
    }
    i=i*2;
    unsigned long blue= i & 255;
    i = i >> 8;
    unsigned long green= i & 255;
    i = i >> 8;
    unsigned long red = i & 255;
    i=0;
    i=(green << 16) | (blue << 8) | red;
    __syncthreads();
    color[(int)data->winx * (int)data->winy + index] = i;
//    printf ("In NV Pixel Test for %i, %i i is %i -- %i\n", (int) data->winy, index, (int)i, (int)color[index]);
}

// end test nv kernel ==========================================================================
void *row_thread (void *args) {
     struct row_data *data=(struct row_data*)args;
     for (int x=0 ; x < data->xend ; x++){
          double sx=((x/data->xend)*data->xscale)+data->xstart;
          double ax=0;
          double ay=0;
          unsigned long i=0;
          // per pixel loop
          while (ax*ax + ay*ay < 4 && i < data->max_itterations) {
               double xtemp = ax*ax - ay*ay + sx;
               ay = 2*ax*ay + data->sy;
               ax=xtemp;
               i++;
          }
          data->colors[x]=i;
//          if (i > 255) {
  //             data->colors[x]=0;
  //        } else {
               i=i*2;
               unsigned long blue= i & 255;
               i = i >> 8;
               unsigned long green= i & 255;
               i = i >> 8;
               unsigned long red = i & 255;
               i=0;
      //         printf ("%i, %i, %i\n", red, green, blue);
               i=(green << 16) | (blue << 8) | red;
               data->colors[x]=i;
               
    //      }
//          printf("debug: %f %i %i %i\n", temp, i, data->max_itterations, data->colors[x] );
     }
     data->complete=1;
     return data;

}


//--------------------------
void mandelbrot_threaded(int full, int pixmap) {

     struct timespec starttime;
     struct timespec stoptime;
     int renderx, rendery;
     GC rendergc;
     if (pixmap ==1 ) {
          renderx=fractwin.pixmapx;
          rendery=fractwin.pixmapy;
          rendergc=fractwin.pixmapgc;
     } else {
          renderx=fractwin.winx;
          rendery=fractwin.winy;
          rendergc=fractwin.gc;
     }
     clock_gettime(CLOCK_MONOTONIC, &starttime);
     int max_threads=5;
     // CUDA override
     double max_itterations=100000/fractview.zoom;
     if (!full) {
          max_itterations=100+fractview.factor;
          max_threads=50;
     }
printf("============================================\nCPU Pthread Render with %i Threads\n============================================\n",max_threads);
     printf("Dims: %i, %i Framebuffer %i Bytes\n", renderx, rendery, ((sizeof(struct row_data)*rendery)+(sizeof(unsigned long)*renderx*rendery)));
     int active_threads=0;
     struct row_data *bitmap=(struct row_data*)malloc(sizeof(struct row_data)*rendery);
     double scalex=3.5*fractview.zoom;
     double scaley=2*fractview.zoom;
     printf("Scale: %f, %f Viewport: %f, %f Zoom: %la\n", scalex, scaley, fractview.viewx, fractview.viewy, fractview.zoom);     
     printf("Max Itterations: %f Factor %la \n", max_itterations, fractview.factor);
     for (double y=0 ; y < rendery ; y++){
          double sy=((y/rendery)*scaley)+fractview.viewy;
          bitmap[(int)y].winx=renderx;
          bitmap[(int)y].xstart=fractview.viewx;
          bitmap[(int)y].xscale=scalex;
          bitmap[(int)y].xend=renderx;
          bitmap[(int)y].max_itterations=max_itterations;
          bitmap[(int)y].sy=sy;
          bitmap[(int)y].complete=0;
          bitmap[(int)y].winy=y;
          bitmap[(int)y].colors=(unsigned long*)malloc(sizeof(long)*(bitmap[(int)y].winx+2));
//          printf ("\rspawning thread: %i", (int)y);
          pthread_create(&(bitmap[(int)y].thread),NULL,row_thread,&bitmap[(int)y]);
          active_threads++;
//          printf ("Drawing\n");
          while (active_threads > max_threads) {
                for (int j=0 ; j < (int)y ; j++) {
                     if (bitmap[j].complete==1) {
                          pthread_t thisthread = bitmap[j].thread;
                          struct row_data *returnval=(struct row_data*)malloc(sizeof(struct row_data));
                          pthread_join(thisthread, (void**)returnval);
                          bitmap[j].complete=2;
                          active_threads--;
                          for (int k=0 ; k < renderx ; k++) {
                               XSetForeground (fractwin.Display, rendergc, bitmap[j].colors[k]);
                               if (pixmap ==0) {
                                    XDrawPoint(fractwin.Display, fractwin.window, rendergc, k, bitmap[j].winy);
                               } else {
                                    XDrawPoint(fractwin.Display, fractwin.pixmap, rendergc, k, bitmap[j].winy);
                               }
                          } // draw rows

                     } // row complete
                     
                } // itterate through launched rows
                XFlush(fractwin.Display);
          } // if we have exceeded max_threads
     } // Y loop
    for (int y=0 ; y < rendery ; y++) {     
          if (bitmap[y].complete==1) {
               pthread_t thisthread = bitmap[y].thread;
               struct row_data *returnval=(struct row_data*)malloc(sizeof(struct row_data));
               pthread_join(thisthread, (void**)returnval);
               bitmap[y].complete=2;
          }
          for (int x=0 ; x < renderx ; x++) {
               XSetForeground (fractwin.Display, fractwin.gc, bitmap[y].colors[x]);
               if (pixmap==0) {
                    XDrawPoint(fractwin.Display, fractwin.window, fractwin.gc, x, bitmap[y].winy);
               } else {
                    XDrawPoint(fractwin.Display, fractwin.pixmap, fractwin.pixmapgc, x, bitmap[y].winy);
               }
               } // draw rows
     } // ensure everything is drawn
     XFlush(fractwin.Display);
     clock_gettime(CLOCK_MONOTONIC, &stoptime);
     printf("\ncompleted in %ld seconds.\n", (stoptime.tv_sec - starttime.tv_sec));
          if (pixmap==1) {
//          XWriteBitmapFile(fractwin.Display, "test.xbm", fractwin.pixmap, renderx, rendery, -1, -1);
          tiff_write(bitmap, renderx, rendery);
     }
     for (int y=0 ; y < (int)rendery ; y++) {
	free(bitmap[y].colors);
	}
     free (bitmap);    
     
printf("============================================\nThread Render Complete\n============================================\n");
}

//---------------------------
void mandelbrot_cuda(int full, int pixmap) {
// The Mandlebrot set is interesting in the real region x = -2 to +1 and y= -1 to +1.

     struct timespec starttime;
     struct timespec stoptime;
     int renderx, rendery;
     double max_itterations;
     GC rendergc;
     if (pixmap ==1 ) {
          renderx=fractwin.pixmapx;
          rendery=fractwin.pixmapy;
          rendergc=fractwin.pixmapgc;
     } else {
          renderx=fractwin.winx;
          rendery=fractwin.winy;
          rendergc=fractwin.gc;
     }
     clock_gettime(CLOCK_MONOTONIC, &starttime);
     max_itterations=100000/fractview.zoom;
     if (!full) {
          max_itterations=100+fractview.factor;
     }

     
     struct row_data *bitmap= (struct row_data*)malloc(sizeof(struct row_data)*rendery);
     double scalex=3.5*fractview.zoom;
     double scaley=2*fractview.zoom;

     struct row_data * d_row;
     unsigned long *d_color;
     
     int threadcount = 512; //multiples of 32 max 1024 per Nvidia documentation
          //  smallest executable unit of parallelism on a CUDA device comprises 32 threads
          // http://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html
     int blockcount=((renderx + threadcount - 1)/threadcount);
     printf("============================================\nCUDA Render with %i blocks of %i threads\n============================================\n", blockcount, threadcount);
     printf("Dims: %i, %i Framebuffer %i Bytes\n", renderx, rendery, ((sizeof(struct row_data)*rendery)+(sizeof(unsigned long)*renderx*rendery)));
     printf("Scale: %f, %f Viewport: %f, %f Zoom: %la\n", scalex, scaley, fractview.viewx, fractview.viewy, fractview.zoom);
     printf("Max Itterations: %f Factor %la \n", max_itterations, fractview.factor);
     size_t mysize = sizeof(long) * renderx;
     cudaMalloc ((void**)&d_row, sizeof(struct row_data));
     cudaMalloc ((void**)&d_color, mysize);
     for (double y=0 ; y < rendery ; y++) {
          double sy=((y/rendery)*scaley)+fractview.viewy;
          bitmap[(int)y].winx=renderx;
          bitmap[(int)y].xstart=fractview.viewx;
          bitmap[(int)y].xscale=scalex;
          bitmap[(int)y].xend=renderx;
          bitmap[(int)y].max_itterations=max_itterations;
          bitmap[(int)y].sy=sy;
          bitmap[(int)y].winy=y;
          bitmap[(int)y].colors=(unsigned long*)malloc(sizeof(unsigned long)*(bitmap[(int)y].winx+2));
          
          //----------------- CUDA_ROW -----------
          
          cudaMemcpy(d_row, &bitmap[(int)y], sizeof(struct row_data), cudaMemcpyHostToDevice);
//          cudaMemcpy(d_color, bitmap[(int)y].colors, mysize, cudaMemcpyHostToDevice);

          nv_pixel<<<blockcount,threadcount>>>(d_color, d_row);
          
          // draw previous row
          if (y > 0) {
               int k;
               for (k=0 ; k < fractwin.winx ; k++) {
                    XSetForeground (fractwin.Display, fractwin.gc, bitmap[(int)y-1].colors[k]);
                    if (pixmap == 0) {
                         XDrawPoint(fractwin.Display, fractwin.window, rendergc, k, bitmap[(int)y-1].winy);
                    } else {
                         XDrawPoint(fractwin.Display, fractwin.pixmap, rendergc, k, bitmap[(int)y-1].winy);
                    }
                }
                XFlush(fractwin.Display);
           }
           // wait for Cuda
           
//          cudaDeviceSynchronize();
          cudaMemcpy(bitmap[(int)y].colors,d_color, mysize, cudaMemcpyDeviceToHost);
          
          //------------------------------------
          
          
     } // main render loop
     cudaFree(d_row); 
     cudaFree(d_color);
     cudaDeviceSynchronize();
     for (int y=0 ; y < rendery ; y++) {
          for (int x=0 ; x < renderx ; x++) {
               XSetForeground (fractwin.Display, fractwin.gc, bitmap[y].colors[x]);
               if (pixmap==0) {
                    XDrawPoint(fractwin.Display, fractwin.window, fractwin.gc, x, bitmap[y].winy);
               } else {
                    XDrawPoint(fractwin.Display, fractwin.pixmap, fractwin.pixmapgc, x, bitmap[y].winy);
               }
          } // draw rows
     } // ensure everything is drawn
     XFlush(fractwin.Display);
     clock_gettime(CLOCK_MONOTONIC, &stoptime);
     printf("\ncompleted in %ld seconds.\n", (stoptime.tv_sec - starttime.tv_sec));
//     printf("factor is: %lf, %lf\n", fractview.factor, fractview.zoom);
     if (pixmap==1) {
          tiff_write(bitmap, renderx, rendery);
     }
     for (int y=0 ; y < (int)rendery ; y++) {
	free(bitmap[y].colors);
	}
     free (bitmap); 
     printf("============================================\nCUDA Render Complete\n============================================\n");  
}


// test routine =====================================================================================================================

void mandelbrot_test(int full, int pixmap) {
// The Mandlebrot set is interesting in the real region x = -2 to +1 and y= -1 to +1.
     struct timespec starttime;
     struct timespec stoptime;
     int renderx, rendery;
     double max_itterations;
     GC rendergc;
     if (pixmap ==1 ) {
          renderx=fractwin.pixmapx;
          rendery=fractwin.pixmapy;
          rendergc=fractwin.pixmapgc;
     } else {
          renderx=fractwin.winx;
          rendery=fractwin.winy;
          rendergc=fractwin.gc;
     }
     clock_gettime(CLOCK_MONOTONIC, &starttime);
     max_itterations=100000/fractview.zoom;
     if (!full) {
          max_itterations=1000+fractview.factor;
     }
     
     struct row_data *bitmap= (struct row_data*)malloc(sizeof(struct row_data)*rendery);
     double scalex=3.5*fractview.zoom;
     double scaley=2*fractview.zoom;

     struct row_data *d_row;
     unsigned long *d_color;
  //   unsigned long  *h_color;
     
     int threadcount = 64; //multiples of 32 max 1024 per Nvidia documentation
          //  smallest executable unit of parallelism on a CUDA device comprises 32 threads
          // http://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html
     int blockcount=((renderx + threadcount - 1)/threadcount);
     printf("============================================\nCUDA Test Render with %i blocks of %i threads\n============================================\n", blockcount, threadcount);
     
     printf("Dims: %i, %i Framebuffer %i Bytes\n", renderx, rendery, (int)(rendery * renderx * sizeof(unsigned long)));
     printf("Scale: %f, %f Viewport: %f, %f Zoom: %l\n", scalex, scaley, fractview.viewx, fractview.viewy, fractview.zoom);
     printf("Max Itterations: %f Factor %l \n", max_itterations, fractview.factor);
     
     // initialization 
     
     cudaStream_t *streams = (cudaStream_t *)malloc(sizeof(cudaStream_t) * rendery);
     cudaMalloc ((void**)&d_row, (sizeof(struct row_data)*rendery));
 //    cudaMalloc ((void**)&d_color, renderx * sizeof(unsigned long));
//     cudaMallocHost ((void**)&h_color, (sizeof(unsigned long)*rendery*renderx));
     
     // collect input
	printf ("memory definition\n");
     for (double y=0 ; y < rendery ; y++) {
          double sy=((y/rendery)*scaley)+fractview.viewy;
          bitmap[(int)y].winx=renderx;
          bitmap[(int)y].xstart=fractview.viewx;
          bitmap[(int)y].xscale=scalex;
          bitmap[(int)y].xend=renderx;
          bitmap[(int)y].max_itterations=max_itterations;
          bitmap[(int)y].sy=sy;
          bitmap[(int)y].winy=y;
	  bitmap[(int)y].complete=0;
          bitmap[(int)y].colors=(unsigned long*)malloc(sizeof(unsigned long)*(bitmap[(int)y].winx));
     }
     
     // send to CUDA
     cudaMemcpy((void*)d_row, (void*)bitmap, (sizeof(struct row_data)*rendery), cudaMemcpyHostToDevice);
printf ("Spawning Kernels \n");
int running = 0;
     cudaMalloc (&d_color, rendery * renderx * sizeof(unsigned long));
     for (int y=0 ; y < (int)rendery ; y++) {
          printf ("\rSpawning Kernel %i ", y);
          if (running < 100 ) {
          	cudaStreamCreate(&streams[y]);
	          nv_pixel_test<<<blockcount,threadcount,0, streams[y]>>>(d_color, &d_row[y]);
		running++;
          } else {
		for (int t=0 ; t < y ; t++) {
			if (!bitmap[t].complete && (cudaStreamQuery(streams[t])==cudaSuccess)) {
				bitmap[t].complete=1;
				cudaMemcpyAsync(bitmap[t].colors, &d_color[t*(int)renderx], sizeof (unsigned long)*renderx, cudaMemcpyDeviceToHost, streams[t]);
				cudaStreamDestroy(streams[t]);
				running--;
                                printf ("\rDrawing Row %i",t);
				for (int x=0 ; x < renderx ; x++) {
			               XSetForeground (fractwin.Display, fractwin.gc, bitmap[t].colors[x]); // working
			               if (pixmap == 0) {
                        			 XDrawPoint(fractwin.Display, fractwin.window, rendergc, x, t);
			               } else {
                        			 XDrawPoint(fractwin.Display, fractwin.pixmap, rendergc, x, t);
			               }
			          }
			}
                }
		y--;
          }
          
     }
     
     printf ("\rcleaning up       \n");
     // get returns and render
     for (int y=0 ; y < (int)rendery ; y++) {
          
//          printf ("\rCompleting Cuda Row %i.",y);
//          cudaStreamSynchronize(streams[y]);   
          if (bitmap[y].complete==0) {     
	        cudaMemcpyAsync(bitmap[y].colors, &d_color[y*(int)renderx], sizeof (unsigned long)*renderx, cudaMemcpyDeviceToHost, streams[y]); // working
          	cudaStreamDestroy(streams[y]);
		printf ("\rDrawing Row %i",y);
		for (int x=0 ; x < renderx ; x++) {
               XSetForeground (fractwin.Display, fractwin.gc, bitmap[y].colors[x]); // working
               if (pixmap == 0) {
                         XDrawPoint(fractwin.Display, fractwin.window, rendergc, x, y);
                    } else {
                         XDrawPoint(fractwin.Display, fractwin.pixmap, rendergc, x, y);
                    }
          }
	  }
     }

	cudaFree(d_color);
	cudaFree(d_row);
        free (streams);
// refreshing screen

     XFlush(fractwin.Display);
     clock_gettime(CLOCK_MONOTONIC, &stoptime);
     printf("\ncompleted in %ld seconds.\n", (stoptime.tv_sec - starttime.tv_sec));
     if (pixmap==1) {
          tiff_write(bitmap, renderx, rendery);
     }
     for (int y=0 ; y < (int)rendery ; y++) {
	free(bitmap[y].colors);
	}
     free (bitmap); 
     printf("============================================\nCUDA Test Render Complete\n============================================\n");  
}

// end test routine =================================================================================================================




void zoom (int direction) {
//double zoomnew;
     if (direction == 0) { 	// zoom out
          fractview.zoom += fractview.factor/1000;
          if (fractview.zoom > (fractview.factor/1000)) {
               fractview.factor *= 10;
          }
          if (fractview.factor > 10) {
               fractview.factor = 10;
          }
//         fractview.factor += .1/fractview.factor;
//        fractview.factor += 1;
     } else {			// zoom in
//          zoomnew = fractview.zoom - (fractview.factor/1000);
//          if (zoomnew < 0) {
//          fractview.factor -= 1/fractview.factor;
//          }
          if (fractview.zoom - (fractview.factor/1000) <=0) {
               fractview.factor /= 10;
          }
          fractview.zoom -= fractview.factor/1000;
//          fractview.factor -= 1/fractview.factor;
//          fractview.factor -= 1;
          if (fractview.zoom < 0 ) {
//               fractview.zoom=1;
//               fractview.factor=10;
          }
     }
}

struct fract_cords start;
     struct fract_cords end;

int main (int argc, char* argv[]) {
// initialization 
     fractwin.winx=1024;
     fractwin.winy=768;
     if (argc == 3) {
          fractwin.winx=atoi(argv[1]);
          fractwin.winy=atoi(argv[2]);
          
     }
     fractwin.pixmapx=fractwin.winx;
     fractwin.pixmapy=fractwin.winy;
      if (fractwin.winx > 1920) {
          fractwin.winx=1920;
     }
     if (fractwin.winy > 1080) {
          fractwin.winy=1080;
     }

     fractview.viewx = -2;
     fractview.viewy = -1;
     fractview.zoom=1;
     fractview.factor=10;

// initialize X
     Xinit();

     

     int quit=1;
     float scalex;
     float scaley;
     while (quit) {
          Window root_return, child_return;
          int root_x_return, root_y_return;
          unsigned int mask_return;
          int mouse_x, mouse_y;
          XEvent event;
          XNextEvent(fractwin.Display, &event);
/*
q -- Exit
c -- zoom in
d -- zoom out
r -- cuda single xwin
f -- cuda single tiff
y -- test xwin
h -- test tiff
t -- cpu xwin
g -- cpu tiff
' ' -- preview render
right click -- cuda xwin
middle click -- cuda tiff
left click -- preview render
scroll -- zoom
*/
          switch (event.type) {
               case (ButtonPress): {
//                    printf("Buttonpress %i\n", event.xbutton.button);
                    switch(event.xbutton.button) {
                         case (1): { // left click
                              scalex=3.5*fractview.zoom;
                              scaley=2*fractview.zoom;
                              XQueryPointer(fractwin.Display, fractwin.window, &child_return, &root_return, &root_x_return, &root_y_return, &mouse_x, &mouse_y, &mask_return);
                              start.y=mouse_y;
                              start.x=mouse_x;
 //                             printf("\n--------------------\n\nStart is: %i, %i\n\n\n", start.x, start.y);
                              break;
                         }
                         
                    }
                    break;
               } // ButtonPress
               case (ButtonRelease): {
 //                   printf("Buttonrelease %i\n", event.xbutton.button);
                    switch(event.xbutton.button) {
                        case(1): // left click
                        {
                              scalex=3.5*fractview.zoom;
                              scaley=2*fractview.zoom;
                              XQueryPointer(fractwin.Display, fractwin.window, &child_return, &root_return, &root_x_return, &root_y_return, &mouse_x, &mouse_y, &mask_return);
                              end.y=mouse_y;
                              end.x=mouse_x;
                              int deltax = start.x-end.x;
                              int deltay = start.y-end.y;
                              float xtemp=((float)deltax/(float)fractwin.winx)*scalex;
                              float ytemp=((float)deltay/(float)fractwin.winy)*scaley;
//                              printf ("------\ndebug output:  \n%i, %i, %f\n-------\n", deltay, fractwin.winy, ((float)deltay/(float)fractwin.winy));
//                              printf("translation: %f, %f, %i \n", deltay, ytemp, fractwin.winy);
                              fractview.viewx += xtemp;
                              fractview.viewy += ytemp;
//                              printf("\n\n\nDelta: %i, %i to %i, %i\n", start.x, start.y, end.x, end.y);
//                              printf("Delta: %i, %i \n", deltax, deltay);
                              
//                              printf("Viewport moving to: %f, %f \n", fractview.viewx, fractview.viewy);
                              if (fractview.viewy > 1) {
                                   fractview.viewy = 1;
                              }
                              if (fractview.viewy < -1) {
                                   fractview.viewy = -1;
                              }
                              if (fractview.viewx > 1 ) {
                                   fractview.viewx = 1;
                              }
                              if (fractview.viewx < -2) {
                                   fractview.viewx = -2;
                              }
                              if ((start.x != end.x) || (start.y != end.y)) {
                              start.x=0;
                              start.y=0;
                              }
                              end.x=0;
                              end.y=0;
                              if (XPending(fractwin.Display)==0) {
                                   mandelbrot_cuda(0,0);
                              }
                              break;           
                         }
                         case (2): // middle click
                         {
                              mandelbrot_cuda(1,1);
                         }
                         case(3): // right click
                         {
                              mandelbrot_cuda(1,0);
                              break;
                         }
                         case(4): // scroll up
                         {
                              zoom(0);
                              if (XPending(fractwin.Display)==0) {
                                   mandelbrot_cuda(0,0);
                              }
                              break;
                         }
                         case(5): // scroll down
                         {
                              zoom(1);
                              if (XPending(fractwin.Display)==0) {
                                   mandelbrot_cuda(0,0);
                              }
                              break;
                         } 
                    } // switch button
                    break;
               } // ButtonRelease
          
               case (KeyPress): {
//                    printf("Keypress %i\n", event.xbutton.button);
                    switch(event.xkey.keycode) {
                         case (40):
                         {
                              zoom(0);
                              if (XPending(fractwin.Display)==0) {
                                   mandelbrot_cuda(0,0);
                              }
                              break;
                         }
                         case (54):
                         {
                              zoom(1);
                              if (XPending(fractwin.Display)==0) {
                                   mandelbrot_cuda(0,0);
                              }
                              break;
                         }

                    } // keycode
                    break;
               } // KeyPress
                    
               case (KeyRelease): {
                    printf("Keyrelease %i\n", event.xbutton.button);
                    switch (event.xkey.keycode) {
                         case(24):
                         {
                              quit=0;
                              break;
                         }
			case(27):
                         {
                              mandelbrot_cuda(1,0);
                              break;
                         }
			case(41):
                         {
                              mandelbrot_cuda(1,1);
				break;
                         }
                         case(28):
                         {
                              mandelbrot_threaded(1,0);
                              break;
                         }
			case(42):
                         {
                              mandelbrot_threaded(1,1);
				break;
                         }
                         case(29):
                         {
                              mandelbrot_test(1,0);
				break;
                         }
			case(43):
                         {
                              mandelbrot_test(1,1);
				break;
                         }
			 case(65):
                         {
                              mandelbrot_test(0,0);
				break;
                         }

                    } // keycode
                    break;
               } // KeyRelease
          } // switch event type      
     } // main event loop

// normal exit
     return 0;
}