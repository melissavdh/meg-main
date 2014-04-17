#ifndef __FFT_H_INCLUDED__
#define __FFT_H_INCLUDED__
    
    typedef std::vector<std::vector<double> > pixelArray;

    pixelArray padKernel(pixelArray actual_valuesPSF,int pady,int padx);
    
    pixelArray removePadding(pixelArray image, int xrange, int yrange);        
    
    pixelArray gaussianPSF(int dim,double sigma);
    
    pixelArray padImage(pixelArray actual_values,int pady,int padx);
    
    pixelArray unpadImage(pixelArray actual_values,int xrange,int yrange);
    
    pixelArray convolve(pixelArray sersic,int xrange,int yrange,double sigma);

#endif
