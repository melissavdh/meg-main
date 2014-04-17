#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

    using namespace std;
    
    typedef vector<vector<double> > pixelArray;
    
    pixelArray readin(char file[]) {

        // read in coordinates and put in 2d array coords to be passed to other
        // functions called from here.
        pixelArray coords;
        coords.resize(2);
        for (int i = 0; i < 2; i++) {
            coords[i].resize(2);
        }
        int count1 = 0;
        int count2 = 0;
        cout << "file = " << file << endl;
        ifstream infile(file);
        char line[4];

        // check the file exists
        if (!infile) {
            cout << "Usage: [fits file] [co-ordinate file]" << endl;
            return coords;
        }
        //cout << "Opened file to obtain region co-ordinates." << endl;
        
        // read in pixel co-ordinates to coords array
        while ( infile.getline(line,4) ) {
            coords[count1][count2] = atof(line);
            if ( count1 < count2 ) {
                count1++;
                count2--;
            }
            else count2++;
        }
        
        return coords;
        
    }
