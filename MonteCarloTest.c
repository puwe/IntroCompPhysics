#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <omp.h>

double ParticlePopulate(double* pParticles, int n, double L, double R)
{

    double d_sum=0;
    double d_mean;
    double d=0;
    unsigned int count=0;
    unsigned int count_max=1e5;

    LOOP1:for(int i=0;i<n;i++){
        count = 0;
        LOOP2:for(int j=0;j<3;j++){
            pParticles[3*i+j] = rand()*L/RAND_MAX; 
            count++;
        }
        for(int k=0;k<i;k++){
            d = pow(pParticles[3*i+0]-pParticles[3*k+0],2)+\
                pow(pParticles[3*i+1]-pParticles[3*k+1],2)+\
                pow(pParticles[3*i+2]-pParticles[3*k+2],2);
            d = sqrt(d);
            if(d<2*R&&count<count_max){
                //printf("%f\n",d);
                goto LOOP2;
            }
        }
        if(count>=count_max){
            //printf("New sphere allocate failure\n");
            goto LOOP1;
        }
    }

        
    d_sum = 0;
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            d = pow(pParticles[3*i+0]-pParticles[3*j+0],2)+\
                pow(pParticles[3*i+1]-pParticles[3*j+1],2)+\
                pow(pParticles[3*i+2]-pParticles[3*j+2],2);
            d = sqrt(d);
            d_sum += d;
        }
    }
    
    
    d_mean = 2.0/(n*(n-1))*d_sum;

    return d_mean;
}

double test(int M, int n, double L, double R)
{
    double res[M];
    clock_t t_start, t_end;
    double cpu_time_used;
    FILE* pFile; 
    t_start=clock();

    for(int k=0;k<M;k++){
        int seed = k+1;
        double* pParticles = (double*)malloc(3*n*sizeof(double)); 
        if(pParticles==NULL)printf("Allocation failure\n");
        res[k] = ParticlePopulate(pParticles, n, L, R);
        free(pParticles);
        printf("%d\t%f\n",k,res[k]);
    }

    t_end=clock();
    cpu_time_used = ((double)(t_end-t_start))/CLOCKS_PER_SEC;
    double dmean;
    pFile = fopen("MonteCarlo.txt","w"); 
    for(int k=0;k<M;k++){
        double sum=0.0;
        for(int j=0;j<=k;j++){
            sum += res[k];
        }
        dmean = sum/(k+1);
        //printf("%d\t%f\n",k+1,dmean); 
        fprintf(pFile,"%d\t%f\n",k+1,dmean);
    }

    printf("Time = %f,\t M = %d,\t dmean = %f \n",cpu_time_used,M,dmean);
    
    return dmean; 

}


int main(int argc, char** argv)
{
    int M=500;
    int n=16;
    double L=1;
    double R=0.1;
    double nu = (n*4/3*M_PI*pow(R,3))/pow(L,3);
    double dmean;
    
    assert(nu < M_PI/(3*sqrt(2))); // fcc close packing
    printf("nu = %f, nu(close packing) = %f\n",nu,M_PI/(3*sqrt(2)));

    
    //dmean = test(M,n,L,R);
    double test[M];
    clock_t t_start, t_end;
    double cpu_time_used;
    FILE* pFile;

    t_start=clock();
    for(int m=1;m<=M;m++){
        double res[m];
        for(int k=0;k<m;k++){
            double* pParticles = (double*)malloc(3*n*sizeof(double)); 
            if(pParticles==NULL)printf("Allocation failure\n");
            res[k] = ParticlePopulate(pParticles, n, L, R);
            free(pParticles); 
        }
        double sum=0.0;
        for(int k=0;k<m;k++){
            sum += res[k];
        }
        double dmean = sum/m;
        test[m-1] = dmean;
        printf("M = %d,\t dmean = %f \n", m, dmean);
    }
    t_end=clock();
    cpu_time_used = ((double)(t_end-t_start))/CLOCKS_PER_SEC; 
    printf("Time = %f,\t max(M) = %d \n",cpu_time_used,M);

    pFile = fopen("MonteCarlo.txt","w"); 
    for(int m=0;m<M;m++){
        fprintf(pFile,"%d\t%f\n",m+1,test[m]);
    }
    return 0;
}
