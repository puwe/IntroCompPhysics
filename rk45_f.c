#include "rk45_f.h"

int rk45_s = 7;
float rk45_ab[]={
            0.,
            1./5.,  1./5.,
            3./10., 3./40.,         9./40.,
            4./5.,  44./45.,        -56./15.,       32./9.,
            8./9.,  19372./6561.,   -25360./2187.,  64448./6561.,   -212./729.,
            1.,     9017./3168.,    -355./33.,      46732./5247.,   49./176.,   -5103./18656.,
            1.,     35./384.,        0.,             500./1113.,     125./192.,  -2187./6784.,   11./84.
            };

float rk45_c5[]={
            35./384.,
            0.,
            500./1113.,
            125./192.,
            -2187./6784.,
            11./84.,
            0.
            };

float rk45_c4[]={
            1951./21600.,
            0.,
            22642./50085.,
            451./720.,
            -12231./42400.,
            649./6300.,
            1./60.,
};

float rk45_h(float t, float h, float *y, int N){ 

    float *yh = (float *) malloc(sizeof(float)*N);
    float *dydh = (float *) malloc(sizeof(float)*rk45_s*N);

    rk45_f(t, dydh, y, N);
    for(int j=1; j<rk45_s; j++){
            int a = j*(j+1)/2;
            int b = a + 1; 

            for(int i=0; i<N; i++){
                float sum = 0;
                for(int k=0; k<j; k++){
                    sum += h*rk45_ab[b+k]*dydh[k*N+i];
                }
                yh[i] = y[i] + sum;
            }
            rk45_f(t+rk45_ab[a]*h, dydh+N*j, yh, N); 
    }
        
    for(int i=0; i<N; i++){
        float dy = 0;
        for(int j=0; j<rk45_s; j++){;
            // c4, c5
            dy += h*rk45_c4[j]*dydh[N*j+i];
        }
        y[i] += dy;
    }
    
    // time update
    t = t+h;
   
    free(yh);
    free(dydh);
    return t; 
}

void rk45_p(float t0, float t, int n, float *y0, int N){
    float h = (t-t0)/(n-1);
    FILE *pf = fopen("ode", "w");
    
    for(int i=0; i<n; i++){
        fprintf(pf, "%f ", t0);
        for(int j=0; j<N; j++){
            fprintf(pf, "%f ", y0[j]);
        }
        fprintf(pf, "\n");
        t0 = rk45_h(t0, h, y0, N); 
    }
   fclose(pf); 
} 
