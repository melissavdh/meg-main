#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <complex>
#include <vector>
#include <stdio.h>
#include <fftw3.h>
#include "Array.h"
#include "fftw++.h"
#include "megio.h"
#include "FFT.h"
#include "globals.h"
#include "vars.h"

    using namespace std;
    using namespace Array;
    using namespace fftwpp;
    
    typedef vector<vector<double> > pixelArray;
    
    // pad kernel to new dimensions and implement column and row wrapping
    pixelArray padKernel(pixelArray actual_valuesPSF,int pady,int padx) {
        
        int yrange = actual_valuesPSF.size();
        int xrange = actual_valuesPSF[0].size();
        
        pixelArray kernel;
        kernel.resize(pady);
        for (int i = 0; i < pady; i++) kernel[i].resize(padx);
        for (int i = 0; i < pady; i++)
            for (int j = 0; j < padx; j++)
                kernel[i][j] = 0;
        
        int x0 = ( (xrange % 2 != 0) ? (xrange+1)/2 : xrange/2 );
        int y0 = ( (yrange % 2 != 0) ? (yrange+1)/2 : yrange/2 );
        int xx = xrange - x0;
        int yy = yrange - y0;
               
        for (int i = 0; i < y0; i++) { // bottom left
            for (int j = 0; j < x0; j++) {
                kernel[i][j] = actual_valuesPSF[yy+i][xx+j];
            }
        }

        for (int i = 0; i < y0; i++) { // bottom right
            for (int j = x0; j < xrange; j++) {
                kernel[i][padx-xrange+j] = actual_valuesPSF[yy+i][(x0-j)*-1];
            }
        }

        for (int i = y0; i < yrange; i++) { // top left
            for (int j = 0; j < x0; j++) {
                kernel[pady-yrange+i][j] = actual_valuesPSF[(y0-i)*-1][xx+j];
            }
        }
            
        for (int i = y0; i < yrange; i++) { // top right
            for (int j = x0; j < xrange; j++) {
                kernel[pady-yrange+i][padx-xrange+j] = actual_valuesPSF[(y0-i)*-1][(x0-j)*-1];
            }
        }

        return kernel;
        
    }
            
    // create Gaussian kernel (only used if PSF not input as image)
    pixelArray gaussianPSF(int dim, double sigma) {
        
        pixelArray kernel;
        kernel.resize(dim);
        for (int i = 0; i < dim; i++) kernel[i].resize(dim);
        double c = ((dim-1)/2.0);
        
        double norm = 0;
        double power = 0;
        for (int i = -c; i <= c; i++) {
            for (int j = -c; j <= c; j++) {
                power = (pow((i),2) + pow((j),2))/(2*sigma*sigma);
                norm = 2.0*M_PI*sigma*sigma;
                kernel[i+c][j+c] = (exp(-power))/norm;
            }
        }
        
        return kernel;
        
    }
    
    // pad image to new dimensions on each side so it remains central
    pixelArray padImage(pixelArray sersic,int pady,int padx) {
    
        int yrange = sersic.size();
        int xrange = sersic[0].size();
        
        pixelArray image;
        image.resize(pady);
        for (int i = 0; i < pady; i++) image[i].resize(padx);
        for (int i = 0; i < pady; i++) {
            for (int j = 0; j < padx; j++) {
                image[i][j] = 0;
            }
        }
        
        int left_pad = ( (xrange % 2 != 0) ? (padx - xrange + 1)/2 : (padx - xrange)/2 );
        int right_pad = ( (xrange % 2 != 0) ? (padx - left_pad + 1) : (padx - left_pad) );
        int top_pad = ( (yrange % 2 != 0) ? (pady - yrange + 1)/2 : (pady - yrange)/2 );
        int bottom_pad = ( (yrange % 2 != 0) ? (pady - top_pad + 1) : (pady - top_pad) );

        for (int i = top_pad; i < bottom_pad; i++) {
            for (int j = left_pad; j < right_pad; j++) {
                image[i][j] = sersic[i - top_pad][j - left_pad];
            }
        }

        return image;
        
    }
    
    // reverse the previous padding operation to extract convolved model
    pixelArray unpadImage(pixelArray conv_values,int xrange,int yrange) {
    
        int pady = conv_values.size();
        int padx = conv_values[0].size();
        
        pixelArray image;
        image.resize(yrange);
        for (int i = 0; i < yrange; i++) image[i].resize(xrange);
        
        int left_pad = ( (xrange % 2 != 0) ? (padx - xrange + 1)/2 : (padx - xrange)/2 );
        int top_pad = ( (yrange % 2 != 0) ? (pady - yrange + 1)/2 : (pady - yrange)/2 );
        
        for (int i = 0; i < yrange; i++) {
            for (int j = 0; j < xrange; j++) {
                image[i][j] = conv_values[i + top_pad][j + left_pad];
            }
        }
        
        return image;
        
    }
    
    // call the above functions and perform convolution
    pixelArray convolve(pixelArray sersic, int xrange, int yrange, double sigma) {
        
        // FIND NEW DIMENSIONS (POWER OF 2)
        double newdim = 2 * max(xrange,yrange);
        double x = ceil(log(newdim)/log(2));
        int X = pow(2,x);
        int Y = X;
        
        pixelArray kernel;
        
        // SET UP PSF KERNEL
        if (psfimage == 1) { // create Gaussian
            double width = sqrt(2*ktol*log(10));
            int kdim = (int)(width*sigma + 0.5) * 2 + 1;
            kernel = gaussianPSF(kdim,sigma);
            int psfxrange = kernel.size();
            int psfyrange = kernel[0].size();
//            string out_kernel = "kernel.fits";
//            sersicImage(out_kernel,kernel,psfyrange,psfxrange);
        } else { // get values of psf image
            int k_x = psfkernel[0].size();
            int k_y = psfkernel.size();
            kernel.resize(k_y);
            for (int i = 0; i < k_y; i++) kernel[i].resize(k_x);
            for (int i = 0; i < k_y; i++) {
                for (int j = 0; j < k_x; j++) {
                    kernel[i][j] = psfkernel[i][j];
                }
            }
        }
        
//        string out_img_orig = "image.fits";
//        sersicImage(out_img_orig,sersic,xrange,yrange);
       
        // PAD IMAGE TO NEW DIMENSIONS
        pixelArray padded_image = padImage(sersic,X,Y);
//        string out_img_pad = "imagePAD.fits";
//        sersicImage(out_img_pad,padded_image,X,Y);
        
        // PAD KERNEL TO NEW DIMENSIONS
        pixelArray padded_kernel = padKernel(kernel,X,Y);
//        string out_ker_pad = "kernelPAD.fits";
//        sersicImage(out_ker_pad,padded_kernel,X,Y);

        // FFT IMAGE
        unsigned int pad = X/2 + 1;
        size_t align = sizeof(Complex);
        array2<double> f(Y,X,align);
        array2<Complex> g(Y,pad,align);
        rcfft2d Forward(X,f,g);
        crfft2d Backward(X,g,f);
        
        for (int i = 0; i < Y; i++)
            for (int j = 0; j < X; j++)
                f(i,j) = padded_image[j][i];
            
        Forward.fft(f,g);
        
        // FFT KERNEL
        array2<double> psff(Y,X,align);
        array2<Complex> psfg(Y,pad,align);
        rcfft2d psfForward(X,psff,psfg);
        
        for (int i = 0; i < Y; i++)
            for (int j = 0; j < X; j++)
                psff(i,j) = padded_kernel[j][i];
            
        psfForward.fft(psff,psfg);
        
        // CONVOLVE
        array2<Complex> decon(Y,pad,align);
        for (int i = 0; i < Y; i++) {
            for (int j = 0; j < pad; j++) {
                real(decon(i,j)) = real(g(i,j))*real(psfg(i,j)) - imag(g(i,j))*imag(psfg(i,j));
                imag(decon(i,j)) = imag(g(i,j))*real(psfg(i,j)) + real(g(i,j))*imag(psfg(i,j));
            }
        }
        
        // INVERSE FFT
        Backward.fftNormalized(decon,f);
        
        // UNPAD RESULT
        pixelArray pixels;
        pixels.resize(X);
        for (int i = 0; i < X; i++) pixels[i].resize(Y);
        for (int i = 0; i < Y; i++) {
            for (int j = 0; j < X; j++) {
                 pixels[j][i] = f(i,j);
            }
        }
        pixelArray done = unpadImage(pixels, xrange, yrange);
        
//        string out_img = "convdmodel.fits";
//        sersicImage(out_img,done,xrange,yrange); 
        
        // RETURN CONVOLVED IMAGE
        return done;
        
    }

