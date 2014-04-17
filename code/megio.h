#ifndef __MEGIO_H_INCLUDED__
#define __MEGIO_H_INCLUDED__

#include <vector>

    typedef std::vector<std::vector<double> > pixelArray;

    // read in co-ordinates of ROI from text file
    pixelArray readin(char file[]);    

    // read in pixel values from region in FITS file
    pixelArray read_values(char file[],pixelArray coords,int xrange,int yrange);
    
    // read in all pixels values in whole image
    pixelArray read_all_values(char file[]);

    // copy original image and replace region with new values.  Keeps all headers and extensions of FITS file.
    void newImage(std::string in_name,pixelArray new_region,int xrange,int yrange,pixelArray coords);

    // create a FITS image showing only the simulated Sersic profile.  Note no header file.
    void sersicImage(std::string out_name,pixelArray new_region,int xrange,int yrange);
    
#endif
