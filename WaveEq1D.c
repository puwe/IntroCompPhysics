#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double Wave1DGaussian(double x, double t, double c){
    double u;
    //u = exp(-pow(x+c*t-10,2));
    u = exp(-pow(x+c*t-10,2))+exp(-pow(x-c*t-10,2));

    return u;
}

int main(int argc, char* argv[]){
    double c=1;
    double b=1.0; 

    double t=0;
    double xsize=10;
    int nx = 10001; //20+1;
    double dx = 2*xsize/(nx-1);
    //int dx = 1;
    //int xsize = 100;
    //int nx = 2*xsize/dx+1;
    double dt = sqrt(b)*dx/c;

    double* u = (double*) malloc(nx*sizeof(double)); // u(x,t)
    double* u0 = (double*) malloc(nx*sizeof(double)); // u(x,t-dt)
    double u_tmp;
    double tmax=2;
    int i;
    double x;
    double t0=0;
    FILE* pFile;
    //initial
    for(i=0; i<nx; i++){
        x = 10-xsize+i*dx;
        u0[i] = Wave1DGaussian(x,t0-dt,c);
        u[i] = Wave1DGaussian(x,t0,c);
        //printf("%f\t%f\n",u0[i],u[i]);
    }
    for(t=0;t<tmax;t+=dt){
        //char filename[128];
        //sprintf(filename,"./data/Wave_t%0.3f.txt",t);
        //pFile = fopen(filename,"w"); 
        for(i=0; i<nx; i++){
            u0[i] = 2*(1-b)*u[i] + b*(u[(i+1)%nx]+u[(i-1)%nx])-u0[i];
        }
        for(i=0; i<nx; i++){
            x = 10-xsize+i*dx;
            //fprintf(pFile,"%f\t%f\n",x,u[i]);
            u_tmp = u[i];
            u[i] = u0[i];
            u0[i] = u_tmp;
        }
        //fclose(pFile);
    }
    
    printf("t:\t%f\n",t);
    pFile = fopen("WaveGaussian1.txt","w");
    for(i=0; i<nx; i++){
        x = 10-xsize+i*dx;
        u_tmp = Wave1DGaussian(x,t0,c);
        fprintf(pFile,"%f\t%f\t%f\n",x,u[i],u_tmp);
        printf("%f\t%f\t%f\n",x,u[i],u_tmp);
    }
    fclose(pFile);

    return 0;
}

