#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <CCfits/CCfits>
#include "lookup.cpp"

    using namespace std;

    typedef vector<vector<double> > pixelArray;
    
    using namespace CCfits;
    
    // read in pixel values from region in FITS file
    pixelArray read_values(char argv[],pixelArray coords,int xrange,int yrange) {
            
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
        if (!fits_open_file(&fptr, argv, READONLY, &status)) {
            cout << "file open: status " << status << endl;
            fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,NULL,values,NULL,&status);
            cout << "file read: status " << status << endl;
            fits_close_file(fptr, &status);
            cout << "file closed: status = " << status << endl;
        }
        if (status) fits_report_error(stderr,status);
        
        // copy pixel values into vector array
        for (int i = 0; i < yrange; i++) {
            for (int j = 0; j < xrange; j++) {
                pixels[i][j] = values[i][j];
            }
        }
        //std::cout << "file " << argv << " opened and values read successfully" << std::endl;
        //std::cout << "(read values) values = " << values[67][73] << std::endl;
        //std::cout << "(read values) pixels = " << pixels[67][73] << std::endl;
        
        return pixels;
        
    }
    
    // estimate total light
    double total_light(int b,pixelArray actual_values,int xrange,int yrange) {
        
        double total_light = 0;
        
        for(int i = 0; i < yrange; i++) {
            for(int j = 0; j < xrange; j++) {
                total_light += (actual_values[i][j] - b);
            }
        }
        
        return total_light;
        
    }
    
    // calculate light at certain distance from centre
    double half_light(int distance,double array_x,double array_y,pixelArray actual_values) {
        
        double half_light = 0;
        for (int i = (array_y - distance); i < (array_y + distance + 1); i++) {
            for (int j = (array_x - distance); j < (array_x + distance + 1); j++) {
                half_light += actual_values[i][j];
            }
        }
       
        return half_light;
        
    }
    
    // estimate half light radius by stepping 1 pixel at a time.  Estimate is just
    // below actual half light radius rather than just above.
    double estimate_rs(double array_x,double array_y,pixelArray actual_values,double Lt) {
    
        double rs = 0;
        
        for (int i = 1; i < 100; i++) {
            double light = half_light(i,array_x,array_y,actual_values);
            if (light > (Lt/2)) {
                rs = i-1;
                break;
            }
        }
        
        return rs;
        
    }

    
    // fill lookup table with gamma function values
    void createLookUp() {
               
        for (int i = 0; i < 1000001; i++) {
            double n = i*0.0001;
            value[i] = GAMMA(n);
        }
                
    }
    
    // estimate start light for sersic calculation
    double estimate_l0(double total_light, double rs, double n) {
        
        double a = 2*M_PI*(rs*rs)*n*gamma(2*n);
        double l0 = total_light/a;
        
        return l0;
        
    }
       
    
    // calculate sersic profile based on distance from centre
    double contsersic(double l0, double rs, double n, double distance) {
        
        double contvalue = 0;

        contvalue = l0 * exp( (pow( (distance/rs),(1/n) ))*(-1) );
        
        return contvalue;
        
    }
    
    
    // calculate fitness statistic
    double fitnessStatistic(pixelArray actual_values,pixelArray sersic_values,int xrange,int yrange,double b) {
        
        double statistic = 0;
        
        for (int y = 0; y < yrange; y++) {
            for (int x = 0; x < xrange; x++) {
                statistic += pow( (actual_values[y][x] - sersic_values[y][x] - b), 2 );
            }
        }
        
        return statistic;
        
    }
    
    // calculate sersic profile in 360 degrees
    pixelArray newProfile(int xrange, int yrange, vector<double> estimates) {
        
        pixelArray sersic_values;
        sersic_values.resize(yrange);
        for (int i = 0; i < yrange; i++) {
            sersic_values[i].resize(xrange);
        }
        int array_x = estimates[4];
        int array_y = estimates[5];

        double distance = 0;
        for (int y = 0; y < yrange; y++) {
            for (int x = 0; x < xrange; x++) {
                distance = sqrt(pow((array_x - x),2) + pow((array_y - y),2) );
                sersic_values[y][x] = contsersic(estimates[3],estimates[2],estimates[0],distance);
            }
        }
        
        return sersic_values;
        
    }
    
    // calculate new image region by subtracting sersic values from actual values
    pixelArray createRegion(pixelArray actual_values, pixelArray sersic_values, int xrange, int yrange) {
        
        pixelArray new_region;
        new_region.resize(yrange);
        for (int i = 0; i < yrange; i++) {
            new_region[i].resize(xrange);
        }
        
        for (int y = 0; y < yrange; y++) {
            for (int x = 0; x < xrange; x++) {
                new_region[y][x] = actual_values[y][x] - sersic_values[y][x];
            }
        }
        
        return new_region;
        
    }
    
    // copy original image and replace region with new values
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
        out_name.insert(position,"_sersic");
        string rewrite_name = out_name;
        out_name.insert(0,"!");
        cout << "outname = " << out_name << endl;

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
            cout << "file opened, status = " << status << endl;
            if (!fits_create_file(&outfptr, out_file, &status)) {
                cout << "file created, status = " << status << endl;
                fits_copy_file(infptr, outfptr, 1, 1, 1, &status);
                cout << "file copied, status = " << status << endl;
                fits_close_file(outfptr, &status);
                cout << "outfile closed, status = " << status << endl;
                fits_close_file(infptr, &status);
                cout << "infile closed, status = " << status << endl;
            }
        }
        if (!fits_open_file(&rewritefptr, rewrite_file, READWRITE, &status)) {
            cout << "outfile opened, status = " << status << endl;
            fits_write_subset(rewritefptr, TDOUBLE, fpixel, lpixel, region, &status);
            cout << "region written: status " << status << endl;
            fits_close_file(rewritefptr, &status);
            cout << "outfile closed, status = " << status << endl;
        }
        if (status) fits_report_error(stderr,status);
        
    }
    

