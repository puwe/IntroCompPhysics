#define _GNU_SOURCE
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<omp.h>

double RandomWalk2D(int N, double dR)
{
    double theta;
    double dx, dy;
    double x0=0;
    double y0=0;
    double x=x0;
    double y=y0;

    for(int i=0; i<N; i++){
        theta = drand48()*2*M_PI;
        dx = dR*cos(theta);
        dy = dR*sin(theta);
        x += dx;
        y += dy;
    }
    
    double R2 = pow(x-x0,2)+pow(y-y0,2);

    return R2;
}

int main(int argc, char** argv)
{
    int M_begin=1;
    int M_end=1e4;
    int N;
    double dR=1;
    double sum=0;
    double sum2=0;
    double mean=0;
    double error=0;
    double R2;
    int M;
    FILE* pFile;
    char filename[128];
    srand48((unsigned)time(NULL));
    
    M=1e4;
    sprintf(filename, "RandomWalk2D_M%d",M);
    pFile = fopen(filename,"w");
    for(N=1; N<500;N++){
        sum = 0;
        sum2 = 0;
        for(int i=0; i<M; i++){
            R2 = RandomWalk2D(N, dR);
            sum += R2;
            sum2 += pow(R2,2);
        }

        mean = sum/M;
        error = sqrt(1.0/M*(sum2/M-pow(mean,2)));
        fprintf(pFile,"%d\t%f\n",N,mean);
        printf("%d\t%f\t%f\t%f\n",N,mean,error,error/mean);
    }

    /**N=10;
    sprintf(filename, "RandomWalk2D_N%d",N);
    pFile = fopen(filename,"w");
    for(M=M_begin; M<M_end;M++){
        sum = 0;
        sum2 = 0;
        for(int i=0; i<M; i++){
            R2 = RandomWalk2D(N, dR);
            sum += R2;
            sum2 += pow(R2,2);
        }

        mean = sum/M;
        error = sqrt(1.0/M*(sum2/M-pow(mean,2)));
        fprintf(pFile,"%d\t%f\n",M,error/mean);
        printf("%d\t%f\t%f\t%f\n",M,mean,error,error/mean);
    }**/

    return 0;
}


    
