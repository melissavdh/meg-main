#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <CCfits/CCfits>
#include "matrix.h"
#include "minfuncs.h"
#include "callfitness.h"
#include "findiff.h"

    using namespace std;
       
    template <class T>
    void lnsrch(vector<double> &xold, const double fold, vector<double> &g, vector<double> &p, vector<double> &x, double &f, const double stpmax, bool &check, T &func) {
        
        cout << "n = " << xold[0] << ", b = " << xold[1] << ", rs = " << xold[2] << ", x = " << xold[3] << endl;
        cout << "y = " << xold[4] << ", q = " << xold[5] << ", th = " << xold[6] << ", l = " << xold[7] << endl;
        
        const double ALF = 1.0e-4, TOLX = numeric_limits<double>::epsilon();
        double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
        double rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
        int i;
        
        int n = xold.size();
        check = false;
        for (i = 0; i < n; i++) sum += p[i]*p[i];
        sum = sqrt(sum);
        if ( sum > stpmax ) for (i = 0; i < n; i++) p[i] *= stpmax/sum;
        for (i = 0; i < n; i++) slope += g[i]*p[i];
        
        if ( slope >= 0.0 ) throw ("roundoff problem in lnsrch");
        test = 0.0;
        for (i = 0; i < n; i++) {
            temp = abs(p[i])/max(abs(xold[i]),1.0);
            if ( temp > test ) test = temp;
        }
        
        alamin = TOLX/test;
        alam = 1.0;
        for (;;) {
            for (i = 0; i < n; i++) x[i] = xold[i]+alam*p[i];
            f = func(x);
            if ( alam < alamin ) {
                for (i = 0; i < n; i++) x[i] = xold[i];
                check = true;
                return;
            } else if ( f <= fold+ALF*alam*slope ) return;
            else {
                if ( alam == 1.0 ) tmplam = -slope/(2.0*(f-fold-slope));
                else {
                    rhs1 = f-fold-alam*slope;
                    rhs2 = f2-fold-alam2*slope;
                    a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                    b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                    if ( a == 0.0 ) tmplam = -slope/(2.0*b);
                    else {
                        disc = b*b-3.0*a*slope;
                        if ( disc < 0.0 ) tmplam = 0.5*alam;
                        else if ( b <= 0.0 ) tmplam = (-b+sqrt(disc));
                    }
                    if ( tmplam > 0.5*alam ) tmplam = 0.5*alam;
                }
            }
            alam2 = alam;
            f2 = f;
            alam = max(tmplam,0.1*alam);
        }
    }
    
    
    template <class T>
    void dfpmin(vector<double> &p, const double gtol, int &iter, double &fret, T &funcd) {
        cout <<"!!!!!!!!!!!!!"<< endl;
        const int ITMAX = 200;
        const double EPS = numeric_limits<double>::epsilon();
        const double TOLX = 4*EPS, STPMX = 100.0;
        bool check;
        double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
        int n = p.size();
        vector<double> dg(n),g(n),hdg(n),pnew(n),xi(n);
        MatDoub hessin(n,n);
        
        fp = funcd(p);
        funcd.df(p,g);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) hessin[i][j] = 0.0;
            hessin[i][i] = 1.0;
            xi[i] = -g[i];
            sum += p[i]*p[i];
        }
        
        stpmax = STPMX*max(sqrt(sum),double(n));
        
        for (int its = 0; its < ITMAX; its++) {
            iter = its;
            lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,funcd);
            fp = fret;
            
            for (int i = 0; i < n; i++) {
                xi[i] = pnew[i] - p[i];
                p[i] = pnew[i];
            }
            test = 0.0;
            
            for(int i = 0; i < n; i++) {
                temp = abs(xi[i])/max(abs(p[i]),1.0);
                if ( temp > test ) test = temp;
            }
            if ( test < TOLX ) return;
            
            for (int i = 0; i < n; i++) dg[i] = g[i];
            funcd.df(p,g);
            test = 0.0;
            den = max(abs(fret),1.0);
            
            for (int i = 0; i < n; i++) {
                temp = abs(g[i])*max(abs(p[i]),1.0)/den;
                if ( temp > test ) test = temp;
            }
            if ( test < gtol ) return;
            
            for (int i = 0; i < n; i++) dg[i] = g[i]-dg[i];
            for (int i = 0; i < n; i++) {
                hdg[i] = 0.0;
                for (int j = 0; j < n; j++) hdg[i] += hessin[i][j]*dg[j];
            }
            fac = fae = sumdg = sumxi = 0.0;
            
            for (int i = 0; i < n; i++) {
                fac += dg[i]*xi[i];
                fae += dg[i]*hdg[i];
                sumdg += pow(dg[i],2);
                sumxi += pow(xi[i],2);
            }
            
            if ( fac > sqrt(EPS*sumdg*sumxi) ) {
                fac = 1.0/fac;
                fad = 1.0/fae;
                for (int i = 0; i < n; i++) dg[i] = fac*xi[i] - fad*hdg[i];
                for (int i = 0; i < n; i++) {
                    for (int j = i; j < n; j++) {
                        hessin[i][j] += fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
                        hessin[j][i] = hessin[i][j];
                    }
                }
            }
            
            for (int i = 0; i < n; i++) {
                xi[i] = 0.0;
                for (int j = 0; j < n; j++) xi[i] -= hessin[i][j]*g[j];
            }
        }
        throw("too many iterations in dfpmin");
    }
            
            
            
        
        
        
        
        
