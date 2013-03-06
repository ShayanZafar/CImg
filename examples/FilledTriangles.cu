/*
 #
 #  File        : CImg_demo.cpp
 #                ( C++ source file )
 #
 #  Description : A multi-part demo demonstrating some of the CImg capabilities.
 #                This file is a part of the CImg Library project.
 #                ( http://cimg.sourceforge.net )
 #
 #  Copyright   : David Tschumperle
 #                ( http://tschumperle.users.greyc.fr/ )
 #
 #  License     : CeCILL v2.0
 #                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed by the CeCILL  license under French law and
 #  abiding by the rules of distribution of free software.  You can  use,
 #  modify and/ or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
 #
*/

// Include static image data, so that the exe does not depend on external image files.
#include "img/CImg_demo.h"
#include <iostream>
#include <iomanip>
//Include the nVidia CUDA runtime for Parallel programming

#include <curand_kernel.h>

// Include CImg library header.
#include "CImg.h"
using namespace cimg_library;
#undef min
#undef max

/*
 * Setup and initialize curand with a seed
 */
__global__ void initCurand(curandState* state){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(100, idx, 0, &state[idx]);
	__syncthreads();
}

/*
 * CUDA kernel that will execute 100 threads in parallel
*/

__global__ void initializeArrays(float* posx, float* posy,float* rayon, float* veloc, float* opacity
								,float* angle, unsigned char* color, int height, int width, curandState* state, size_t pitch){

	int idx =  blockIdx.x * blockDim.x + threadIdx.x;
	curandState localState = state[idx];

	posx[idx] = (float)(curand_normal(&localState)*width);
    posy[idx] = (float)(curand_normal(&localState)*height);
    rayon[idx] = (float)(10 + curand_normal(&localState)*50);
    angle[idx] = (float)(curand_normal(&localState)*360);
    veloc[idx] = (float)(curand_uniform(&localState)*20 - 10);
    color[idx*pitch] = (unsigned char)(curand_normal(&localState)*255);
    color[(idx*pitch)+1] = (unsigned char)(curand_normal(&localState)*255);
    color[(idx*pitch)+2] = (unsigned char)(curand_normal(&localState)*255);
    opacity[idx] = (float)(0.3f + 1.5f *curand_normal(&localState));
	
	__syncthreads();
}

void errCheck(cudaError_t err, const char* msg){
	 if (err != cudaSuccess)
        std::cout<< msg << ": " << cudaGetErrorString(err) << std::endl;
}

/*---------------------------

  Main procedure

  --------------------------*/
int main() {

    // Create a colored 640x480 background image which consists of different color shades.
    CImg<float> background(640,480,1,3);
    cimg_forXY(background,x,y) background.fillC(x,y,0,
                                                x*std::cos(6.0*y/background.height()) + y*std::sin(9.0*x/background.width()),
                                                x*std::sin(8.0*y/background.height()) - y*std::cos(11.0*x/background.width()),
                                                x*std::cos(13.0*y/background.height()) - y*std::sin(8.0*x/background.width()));
    background.normalize(0,180);
    
    // Init images and create display window.
    CImg<unsigned char> img0(background), img;
    unsigned char white[] = { 255, 255, 255 }, color[100][3];
    CImgDisplay disp(img0,"[#6] - Filled Triangles (Click to shrink)");

    // Define random properties (pos, size, colors, ..) for all triangles that will be displayed.
    float posx[100];
	float posy[100];
	float rayon[100];
	float angle[100];
	float veloc[100];
	float opacity[100];
	// Define the same properties but for the device
	float* d_posx;
	float* d_posy;
	float* d_rayon;
	float* d_angle;
	float* d_veloc;
	float* d_opacity;
	unsigned char* d_color;

	// CURAND state
	curandState* devState;
	// error handling
	cudaError_t err;

	// allocate memory on the device for the device arrays
	err = cudaMalloc((void**)&d_posx, 100 * sizeof(float));
	errCheck(err, "cudaMalloc((void**)&d_posx, 100 * sizeof(float))");
	err = cudaMalloc((void**)&d_posy, 100 * sizeof(float));
	errCheck(err,"cudaMalloc((void**)&d_posy, 100 * sizeof(float))");
	err = cudaMalloc((void**)&d_rayon, 100 * sizeof(float));
	errCheck(err,"cudaMalloc((void**)&d_rayon, 100 * sizeof(float))");
    err = cudaMalloc((void**)&d_angle, 100 * sizeof(float));
	errCheck(err,"cudaMalloc((void**)&d_angle, 100 * sizeof(float))");
	err = cudaMalloc((void**)&d_veloc, 100 * sizeof(float));
	errCheck(err,"cudaMalloc((void**)&d_veloc, 100 * sizeof(float))");
	err = cudaMalloc((void**)&d_opacity, 100 * sizeof(float));
	errCheck(err,"cudaMalloc((void**)&d_opacity, 100 * sizeof(float))");
	err = cudaMalloc((void**)&devState, 100*sizeof(curandState));
	errCheck(err,"cudaMalloc((void**)&devState, 100*sizeof(curandState))");
	size_t pitch;
	//allocated the device memory for source array  
	err = cudaMallocPitch(&d_color, &pitch, 3 * sizeof(unsigned char),100);
	errCheck(err,"cudaMallocPitch(&d_color, &pitch, 3 * sizeof(unsigned char),100)");
	// launch grid of threads
	dim3 dimBlock(100);
	dim3 dimGrid(1);
	  
	/* Kernel for initializing CURAND */
	initCurand<<<1,100>>>(devState);

	// synchronize the device and the host
    cudaDeviceSynchronize();
     
	/*Kernel for initializing Arrays */
	initializeArrays<<<1, 100>>>(d_posx, d_posy, d_rayon, d_veloc, d_opacity, d_angle,
										d_color, img0.height(), img0.width(), devState, pitch);
	// synchronize the device and the host
    cudaDeviceSynchronize();
	
	// get the populated arrays back to the host for use
	err = cudaMemcpy(posx,d_posx, 100 * sizeof(float), cudaMemcpyDeviceToHost);
	errCheck(err,"cudaMemcpy(posx,d_posx, 100 * sizeof(float), cudaMemcpyDeviceToHost)");
	err = cudaMemcpy(posy,d_posy, 100 * sizeof(float), cudaMemcpyDeviceToHost);
	errCheck(err,"cudaMemcpy(posy,d_posy, 100 * sizeof(float), cudaMemcpyDeviceToHost)");
	err = cudaMemcpy(rayon,d_rayon, 100 * sizeof(float), cudaMemcpyDeviceToHost);
	errCheck(err,"cudaMemcpy(rayon,d_rayon, 100 * sizeof(float), cudaMemcpyDeviceToHost)");
	err = cudaMemcpy(veloc,d_veloc, 100 * sizeof(float), cudaMemcpyDeviceToHost);
	errCheck(err,"cudaMemcpy(veloc,d_veloc, 100 * sizeof(float), cudaMemcpyDeviceToHost)");
	err = cudaMemcpy(opacity,d_opacity, 100 * sizeof(float), cudaMemcpyDeviceToHost);
	errCheck(err,"cudaMemcpy(opacity,d_opacity, 100 * sizeof(float), cudaMemcpyDeviceToHost)");
	err = cudaMemcpy(angle,d_angle, 100 * sizeof(float), cudaMemcpyDeviceToHost);
	errCheck(err,"cudaMemcpy(angle,d_angle, 100 * sizeof(float), cudaMemcpyDeviceToHost)");
	// pitch of color array is 3+1 padded
	err = cudaMemcpy2D(color,4,d_color,pitch,3 *sizeof(unsigned char),3, cudaMemcpyDeviceToHost);
	errCheck(err,"cudaMemcpy2D(color,pitch,d_color,100*3,3 *sizeof(unsigned char),100* sizeof(unsigned char), cudaMemcpyDeviceToHost)");
    // measuring time it takes for triangle animations in 1000 iterations
    int i = 0, num = 1;
    
    // Start animation loop.
    while (!disp.is_closed() && !disp.is_keyQ() && !disp.is_keyESC() && i < 1000) {
        img = img0;
        
        i++;
        // Draw each triangle on the background image.
        for (int k = 0; k<num; ++k) {
            const int
            x0 = (int)(posx[k] + rayon[k]*std::cos(angle[k]*cimg::PI/180)),
            y0 = (int)(posy[k] + rayon[k]*std::sin(angle[k]*cimg::PI/180)),
            x1 = (int)(posx[k] + rayon[k]*std::cos((angle[k] + 120)*cimg::PI/180)),
            y1 = (int)(posy[k] + rayon[k]*std::sin((angle[k] + 120)*cimg::PI/180)),
            x2 = (int)(posx[k] + rayon[k]*std::cos((angle[k] + 240)*cimg::PI/180)),
            y2 = (int)(posy[k] + rayon[k]*std::sin((angle[k] + 240)*cimg::PI/180));
            if (k%10) img.draw_triangle(x0,y0,x1,y1,x2,y2,color[k],opacity[k]);
            else img.draw_triangle(x0,y0,x1,y1,x2,y2,img0,0,0,img0.width()-1,0,0,img.height()-1,opacity[k]);
            img.draw_triangle(x0,y0,x1,y1,x2,y2,white,opacity[k],~0U);
            
            // Make the triangles rotate, and check for mouse click event.
            // (to make triangles collapse or join).
            angle[k]+=veloc[k];
            if (disp.mouse_x()>0 && disp.mouse_y()>0) {
                float u = disp.mouse_x() - posx[k], v = disp.mouse_y() - posy[k];
                if (disp.button()) { u = -u; v = -v; }
                posx[k]-=0.03f*u, posy[k]-=0.03f*v;
                if (posx[k]<0 || posx[k]>=img.width()) posx[k] = (float)(cimg::rand()*img.width());
                if (posy[k]<0 || posy[k]>=img.height()) posy[k] = (float)(cimg::rand()*img.height());
            }
        }
        
        // Display current animation framerate, and refresh display window.
        img.draw_text(5,5,"%u frames/s",white,0,0.5f,13,(unsigned int)disp.frames_per_second());
        img0.resize(disp.display(img).resize(false).wait(20));
        if (++num>100) num = 100;
        
        // Allow the user to toggle fullscreen mode, by pressing CTRL+F.
        if (disp.is_keyCTRLLEFT() && disp.is_keyF()) disp.resize(640,480,false).toggle_fullscreen(false);
    }

	// free allocated device memory
	cudaFree(d_posy);
	cudaFree(d_posx);
	cudaFree(d_rayon);
	cudaFree(d_veloc);
	cudaFree(d_opacity);
	cudaFree(d_color);
	cudaFree(d_angle);
	cudaFree(devState);
  return 0;
}