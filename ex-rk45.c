#include "rk45_f.h"

void rk45_f (float t, float *k, float *y, int N){
    
    float f = t*t - t - 3;
    float g = 3*sin(t - 0.25);
    k[0] = -f*y[0] + g;
    
    /*
    float mu = 1;
    k[0] = y[1];
    k[1] = mu*(1 - y[0]*y[0])*y[1] - y[0];
    */
}

int main(){

    int n = 101;
    float t0 = 1;
    float t = 5;
    int N = 1;
    float *y0 = (float*) malloc(sizeof(float)*N);
    y0[0] = 1;
    //y0[1] = 0;

    rk45_p(t0, t, n, y0, N);

    return 0;
}
