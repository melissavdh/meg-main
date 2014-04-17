#include <stdlib.h>
#include <limits>
#include <iostream>
#include "callfitness.h"
#include "powell.cpp"
#include "amoeba.cpp"
#include "quasi.cpp"
#include "callminimise.h"
#include "dlib/optimization.h"
#include "levmar.h"


    using namespace std;
    
    typedef dlib::matrix<double,0,1> column_vector;
    int itmax = 2;
    
    // start Nelder-Mead minimisation
    vector<double> callAmoeba(vector<double> estimates) {
          
        // step size for each variable
        vector<double> del(10);
        del[0] = 0.5;
        del[1] = 10;
        del[2] = 0.5;
        del[3] = 0.5;
        del[4] = 0.5;
        del[5] = -0.01;
        del[6] = 0.1;
        del[7] = 1000; // try 1000 upwards was 500
        del[8] = 0.01;
        del[9] = -0.00001; // was 50
      
        Amoeba am(0.0001);
        int its = 0;
        
        while ( (its < itmax) ) {
            estimates = am.minimise(estimates,del,fitness);
            its++;
        }
        
        return estimates;
        
    }
    
    // Start minimsation using Powell's method
    vector<double> callPowell(vector<double> estimates) {
        
        
        Powell<double(vector<double>&)> powell(fitness,0.1);
        int its = 0;

        while ( (its < itmax) ) {
            estimates = powell.minimise(estimates);
            its++;
        }
        
        return estimates;
            
    }
    
    // set up functor for BOBYQA
    struct testFunction {
        
        double operator ()(const column_vector & starting_point) const {

            vector<double> estimates(10);
            estimates[0] = starting_point(0)/10; estimates[1] = starting_point(1)*100; estimates[2] = starting_point(2)/8; 
            estimates[3] = starting_point(3)/2; estimates[4] = starting_point(4)/2; estimates[5] = starting_point(5)/100; 
            estimates[6] = (starting_point(6)/15)*M_PI; estimates[7] = starting_point(7)*1000; estimates[8] = starting_point(8)/100;
            estimates[9] = starting_point(9)/50000;
            
            return fitness(estimates);
        }
        
    };
    
    // Minimise using BOBYQA
    vector<double> callBobyqa(vector<double> estimates) {
        
        testFunction f;
        column_vector starting_point;
        starting_point.set_size(10);
        
        // starting point defined by scaled initial estimates
        starting_point = estimates[0]*10,estimates[1]/100,estimates[2]*8,estimates[3]*2,
            estimates[4]*2,estimates[5]*100,(estimates[6]/M_PI)*15,estimates[7]/1000,estimates[8]*100,
            estimates[9]/100;

            
        double l0lower = estimates[7]/1000;
        double siglow = (estimates[8]-0.5)*100;
        double sighigh = (estimates[8]+0.5)*100;
        // lower and upper bounds.  Must be equal ranges for each parameter.
        dlib::matrix<double> lower_bound(10,1), upper_bound(10,1);
        lower_bound = 0,0,0,0,0,0,0,l0lower,siglow,0;
        upper_bound = 100,100,500,1000,1000,200,100,10000,sighigh,11;
        
        dlib::find_min_bobyqa(f, 
                        starting_point, 
                        20,    // number of interpolation points
                        lower_bound,  // lower bound constraint
                        upper_bound,   // upper bound constraint
                        5,    // initial trust region radius rho_begin
                        0.001,  // stopping trust region radius rho_end
                        10000    // max number of objective function evaluations
        );
        
        estimates[0] = starting_point(0)/10; estimates[1] = starting_point(1)*100; estimates[2] = starting_point(2)/8; 
        estimates[3] = starting_point(3)/2; estimates[4] = starting_point(4)/2; estimates[5] = starting_point(5)/100; 
        estimates[6] = (starting_point(6)/15)*M_PI; estimates[7] = starting_point(7)*1000; estimates[8] = starting_point(8)/100;
        estimates[9] = starting_point(9)/50000;
        
        vector<double> second(10);
        dlib::find_min_bobyqa(f, 
                        starting_point, 
                        20,    // number of interpolation points
                        lower_bound,  // lower bound constraint
                        upper_bound,   // upper bound constraint
                        3,    // initial trust region radius rho_begin
                        0.0001,  // stopping trust region radius rho_end
                        10000    // max number of objective function evaluations
        );

        second[0] = starting_point(0)/10; second[1] = starting_point(1)*100; second[2] = starting_point(2)/8; 
        second[3] = starting_point(3)/2; second[4] = starting_point(4)/2; second[5] = starting_point(5)/100; 
        second[6] = (starting_point(6)/15)*M_PI; second[7] = starting_point(7)*1000; second[8] = starting_point(8)/100;
        second[9] = starting_point(9)/50000;

        if ( fitness(estimates) < fitness(second) ) return estimates;
        else return second;
        //return estimates;
        
    }
    
    // set up functor for L-M
    void fitnessFunc(double *p, double *hx, int m, int n, void *adata) {
        
        vector<double> estimates(10);
        /*estimates[0] = p[0]; estimates[1] = p[1]*100; estimates[2] = p[2]; estimates[3] = p[3]*5;
        estimates[4] = p[4]*5; estimates[5] = p[5]/100; estimates[6] = (p[6]/100)*M_PI; estimates[7] = p[7]*1000;
        estimates[8] = p[8]/10; estimates[9] = p[9]/1000;*/
        
        estimates[0] = p[0]; estimates[1] = p[1]; estimates[2] = p[2]*5; estimates[3] = p[3];
        estimates[4] = p[4]; estimates[5] = p[5]; estimates[6] = p[6]*50; estimates[7] = p[7]*25000;
        estimates[8] = p[8]; estimates[9] = p[9]*100;
        
        /*cout << "n = " << estimates[0] << ", b = " << estimates[1] << ", rs = " << estimates[2] << endl;
        cout << "x0 = " << estimates[3] << ", y0 = " << estimates[4] << ", q = " << estimates[5] << endl;
        cout << "theta = " << estimates[6] << ", l0 = " << estimates[7] << ", sigma = " << estimates[8] << endl;
        cout << "delta = " << estimates[9] << endl;*/
        
        register double stat = fitness(estimates);
        for (int i = 0; i < 10; i++) hx[i] = 0.0;
        hx[0] = stat;
        
    }
    
    // Minimise using L-M
    vector<double> callLM(vector<double> estimates) {
        
        // starting point defined by scaled initial estimates.
        double p[10];
        
        p[0] = estimates[0]; p[1] = estimates[1]; p[2] = estimates[2]/5; p[3] = estimates[3];
        p[4] = estimates[4]; p[5] = estimates[5]; p[6] = estimates[6]/50; p[7] = estimates[7]/25000;
        p[8] = estimates[8]; p[9] = estimates[9]/100;
        
        // lower bounds
        double l0lower = estimates[7]/25000;
        double lb[10];
        lb[0] = 0; lb[1] = 0; lb[2] = 0; lb[3] = 0; lb[4] = 0; lb[5] = 0;
        lb[6] = 0; lb[7] = l0lower; lb[8] = 0; lb[9] = 0.0000000000001;
        
        // upper bounds
        double ub[10];
        ub[0] = 100; ub[1] = 10000; ub[2] = 25; ub[3] = 500; ub[4] = 500; ub[5] = 1;
        ub[6] = M_PI*50; ub[7] = 250; ub[8] = 50; ub[9] = 0.001;
        
        // finite difference approximation (will be calculated by the minimisation)
        double x[10];
        x[0] = x[1] = x[2] = x[3] = x[4] = x[5] = x[6] = x[7] = x[8] = x[9] = 0.0;
        
        // set tolerance criterion
        cout << "mu = " << LM_INIT_MU << endl;
        double opts[5];
        opts[0] = LM_INIT_MU; opts[1] = 0.000000000001; opts[2] = 0.000000000001; opts[3] = 0.000001; opts[4] = LM_DIFF_DELTA;
        
        // minimise
        int result = dlevmar_bc_dif(fitnessFunc, p, x, 10, 10, lb, ub, 1000, opts, NULL, NULL, NULL, NULL);
        result = dlevmar_bc_dif(fitnessFunc, p, x, 10, 10, lb, ub, 1000, opts, NULL, NULL, NULL, NULL);
               
        estimates[0] = p[0]; estimates[1] = p[1]; estimates[2] = p[2]*5; estimates[3] = p[3];
        estimates[4] = p[4]; estimates[5] = p[5]; estimates[6] = p[6]*50; estimates[7] = p[7]*25000;
        estimates[8] = p[8]; estimates[9] = p[9]*100;
        
        while ( (estimates[6] >= M_PI) || (estimates[6] < 0) ) {
            if (estimates[6] < 0) estimates[6] = M_PI - estimates[6];
            if (estimates[6] >= M_PI) estimates[6] -= M_PI;
        }
        
        /*cout << "n = " << estimates[0] << ", b = " << estimates[1] << ", rs = " << estimates[2] << endl;
        cout << "x0 = " << estimates[3] << ", y0 = " << estimates[4] << ", q = " << estimates[5] << endl;
        cout << "theta = " << estimates[6] << ", l0 = " << estimates[7] << ", sigma = " << estimates[8] << endl;
        cout << "delta = " << estimates[9] << endl;*/
        
        return estimates;
    }
               
