#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include "fitsio.h"
#include "globals.h"
#include "megio.h"
	
    using namespace std;


    typedef vector<vector<double> > pixelArray;   
    
    pixelArray readin(char file[]) {

        // read in coordinates and put in 2d array coords to be passed to other
        // functions called from here.
        pixelArray coordinates;
        coordinates.resize(2);
        for (int i = 0; i < 2; i++) {
            coordinates[i].resize(2);
        }
        int count1 = 0;
        int count2 = 0;

        ifstream infile(file);
        char line[10];

        // check the file exists
        if (file == NULL) {
            cerr << "ERROR: No co-ordinate file specified." << endl;
            cerr << "       See --help." << endl;
            exit(1);
        }
        if (!infile) {
            cerr << "ERROR: Can't find co-ordinate file " << file << "." << endl;
            exit(1);
        }
        
        // read in pixel co-ordinates to coords array
        while ( infile.getline(line,10) ) {
            coordinates[count1][count2] = atof(line);
            if ( count1 < count2 ) {
                count1++;
                count2--;
            } else count2++;
        }
        //cout << "bottom left: (" << coordinates[0][0] << "," << coordinates[0][1] << ")" << endl;
        //cout << "top right  : (" << coordinates[1][0] << "," << coordinates[1][1] << ")" << endl;
        
        return coordinates;
        
    }
    
    
    // read in pixel values from region in FITS file
    pixelArray read_values(char file[],pixelArray coords,int xrange,int yrange) {

        // check the fits file exists
        if (file == NULL) {
            cerr << "ERROR: No FITS file specified." << endl;
            cerr << "       See --help." << endl;
            exit(1);
        }
            
        // variables needed for FITSIO functions
        double values[yrange][xrange];
        fitsfile *fptr;
        int status = 0;
        long int fpixel[2] = {coords[0][0],coords[0][1]};
        long int lpixel[2] = {coords[1][0],coords[1][1]};
        long int inc[2] = {1,1};
        
        // vector array to pass back to main function
        pixelArray pixels;
        pixels.resize(yrange);
        for (int i = 0; i < yrange; i++) {
            pixels[i].resize(xrange);
        }
        
        // open FITS file and read in pixel values in correct co-ords
        if (fits_open_file(&fptr, file, READONLY, &status)) {
            cerr << "ERROR: Can't find FITS file " << file << "." << endl;
            exit(1);
        } else {
            fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,NULL,values,NULL,&status);
            fits_close_file(fptr, &status);
        }
        if (status) {
            fits_report_error(stderr,status);
            exit(1);
        }

        // copy pixel values into vector array
        for (int i = 0; i < yrange; i++) {
            for (int j = 0; j < xrange; j++) {
                pixels[i][j] = values[i][j];
            }
        }
        
        return pixels;
        
    }
    
    
    // read in pixel values from entire file
    pixelArray read_all_values(char file[]) {
        
        // check the fits file exists
        if (file == NULL) {
            cerr << "ERROR: No FITS file specified." << endl;
            cerr << "       See --help." << endl;
            exit(1);
        }
            
        // variables needed for FITSIO functions
        //double values[yrange][xrange];
        fitsfile *fptr;
        int status = 0;
        long fpixel[2] = {1,1};
        long naxes[2];
        
        // open FITS file, get dimensions, and read in pixel values in correct co-ords
        if (fits_open_file(&fptr, file, READONLY, &status)) {
            cerr << "ERROR: Can't find FITS file " << file << "." << endl;
            exit(1);
        } else {
            fits_get_img_size(fptr,2,naxes,&status);
            //fits_read_pix(fptr,TDOUBLE,fpixel,(naxes[0]*naxes[1]),NULL,values,NULL,&status);
            //fits_close_file(fptr, &status);
        }
        
        double values[naxes[1]][naxes[0]];
        long nelements = naxes[0]*naxes[1];
        
        // vector array to pass back to main function
        pixelArray pixels;
        pixels.resize(naxes[1]);
        for (int i = 0; i < naxes[1]; i++) {
            pixels[i].resize(naxes[0]);
        }
        
        if (!fits_read_pix(fptr,TDOUBLE,fpixel,nelements,NULL,values,NULL,&status)) {
            fits_close_file(fptr, &status);
        }
        if (status) {
            cerr << "ERROR: problem reading FITS file " << file << "." << endl;
            fits_report_error(stderr,status);
            exit(1);
        }
        
        // copy pixel values into vector array
        for (int i = 0; i < naxes[1]; i++) {
            for (int j = 0; j < naxes[0]; j++) {
                pixels[i][j] = values[i][j];
            }
        }
        
        return pixels;
        
    }
    
    // copy original image and replace region with new values.  Keeps all extensions of FITS file.
    void newImage(string in_name,pixelArray new_region,int xrange,int yrange,pixelArray coords) {
        
        double region[yrange][xrange];
        
        // copy region values from 2d vector to 2d array
        for (int i = 0; i < yrange; i++) {
            for (int j = 0; j < xrange; j++) {
                region[i][j] = new_region[i][j];
            }
        }
        
        string out_name = in_name;             
        size_t position;
        position = out_name.find(".fit");
        out_name.insert(position,"_residual");
        string rewrite_name = out_name;        
        out_name.insert(0,"!");

        char *in_file = new char[in_name.size() + 1];
        in_file[in_name.size()] = 0;
        memcpy(in_file, in_name.c_str(), in_name.size());
        
        char *out_file = new char[out_name.size() + 1];
        out_file[out_name.size()] = 0;
        memcpy(out_file, out_name.c_str(), out_name.size());
        
        char *rewrite_file = new char[rewrite_name.size() + 1];
        rewrite_file[rewrite_name.size()] = 0;
        memcpy(rewrite_file, rewrite_name.c_str(), rewrite_name.size());
        
        fitsfile *infptr;
        fitsfile *outfptr;
        fitsfile *rewritefptr;
        int status = 0;
        long int fpixel[2] = {coords[0][0],coords[0][1]};
        long int lpixel[2] = {coords[1][0],coords[1][1]};
        
        if (!fits_open_file(&infptr, in_file, READONLY, &status)) {
            if (!fits_create_file(&outfptr, out_file, &status)) {
                fits_copy_file(infptr, outfptr, 1, 1, 1, &status);
                fits_close_file(outfptr, &status);
                fits_close_file(infptr, &status);
            }
        }
        if (!fits_open_file(&rewritefptr, rewrite_file, READWRITE, &status)) {
            fits_write_subset(rewritefptr, TDOUBLE, fpixel, lpixel, region, &status);
            fits_close_file(rewritefptr, &status);
        }
        if (status) {
            cerr << "ERROR: Couldn't write residual to FITS file " << out_name.substr(1) << endl;
            fits_report_error(stderr,status);
            exit(1);
        }
        
        cout << "Residual saved as new image: " << out_name.substr(1) << endl;
        
    }
    
    
    // create a FITS image showing only the simulated Sersic profile
    void sersicImage(string out_name,pixelArray new_region,int xrange,int yrange) {
        
        double region[yrange][xrange];

        // copy region values from 2d vector to 2d array
        for (int i = 0; i < yrange; i++) {
            for (int j = 0; j < xrange; j++) {
                region[i][j] = new_region[i][j];
            }
        }
                   
        size_t position;
        position = out_name.find(".fit");
        out_name.insert(position,"_model");
        out_name.insert(0,"!");
        
        char *out_file = new char[out_name.size() + 1];
        out_file[out_name.size()] = 0;
        memcpy(out_file, out_name.c_str(), out_name.size());
        
        fitsfile *outfptr;
        int status = 0;
        long int fpixel[2] = {1,1};
        long naxes[2] = {xrange,yrange};
        long nelements = xrange*yrange;
        long naxis = 2;
        
        if (!fits_create_file(&outfptr, out_file, &status)) {
            if (!fits_create_img(outfptr, DOUBLE_IMG, naxis, naxes, &status)) {
                fits_write_pix(outfptr, TDOUBLE, fpixel, nelements, region, &status);
                fits_close_file(outfptr, &status);
            }
        }
        if (status) {
            cerr << "ERROR: Couldn't write model to FITS file " << out_name.substr(1) << endl;
            fits_report_error(stderr,status);
            exit(1);
        }
        
        cout << "Model saved as new image: " << out_name.substr(1) << endl;
        
    }

