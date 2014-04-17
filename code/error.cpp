#include <iostream>
#include <stdlib.h>
#include <string>

    using namespace std;
    
    struct Error {
        
        string error;
        string file;
        int line;
        Error(string e, string f, int l) : error(e), file(f), line(l) {}
        
    };
    
    void catchError(Error err) {
        
        cout << "Exception: " << err.error << " in file " << err.file << " at line " << err.line;
        exit(1);
        
    }
    
    #define throw(error) throw(Error(error,__FILE__,__LINE__));
