#ifndef RK45_H
#define RK45_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// y' = f(t, y)
// y: N by 1 array
void rk45_f(float t, float *dydt, float *y, int N);

// step: t->t+h 
// return t+h
float rk45_h(float t, float h, float *y, int N);

// time: [t0, t] with n points
// ini_value: y0 with size N
// output file: "ode"  
void rk45_p(float t0, float t, int n, float *y0, int N);

#endif

