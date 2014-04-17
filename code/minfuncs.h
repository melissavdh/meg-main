#ifndef __MIN_FUNCS_H_INCLUDED__
#define __MIN_FUNCS_H_INCLUDED__
   
    template<class T>
    inline T SIGN(const T &a, const T &b) {
        return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
    }
    
    inline void shft3(double &a, double &b, double &c, const double d) {
        a = b;
        b = c;
        c = d;
    } 
    
    template<class T>
    inline void SWAP(T &a, T &b) {
        T dum=a; a=b; b=dum;
    }
    
    template<class T>
    inline const T &MAX(const T &a, const T &b) {
        return b > a ? (b) : (a);
    }

    inline float MAX(const double &a, const float &b) {
        return b > a ? (b) : float(a);
    }

    inline float MAX(const float &a, const double &b)  {
        return b > a ? float(b) : (a);
    }
    
    template<class T>
    inline const T &MIN(const T &a, const T &b) {
        return b < a ? (b) : (a);
    }

    inline float MIN(const double &a, const float &b) {
        return b < a ? (b) : float(a);
    }

    inline float MIN(const float &a, const double &b) {
        return b < a ? float(b) : (a);
    }
    
#endif
