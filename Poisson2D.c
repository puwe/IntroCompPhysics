#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

int main(int argc, char* argv[])
{
    int N = 2e2;
    double phiJacob[N*N];
    double phiGauss[N*N];
    double rho[N*N];
    double h = 1.0/(N-1);
    double phitmp[N*N];
    // initialize
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            phiGauss[i+j*N]=0;
            phiJacob[i+j*N]=0;
            rho[i+j*N]=0;
            phitmp[i+j*N]=0;
            
        }
    }
    rho[3*N/4+3*N/4*N] = 1.0/pow(h,2);
    rho[N/2+N/2*N] = 1.0/pow(h,2);

    for(int step=0; step<1e4; step++){
        for(int i=1; i<N-1; i++){
            for(int j=1; j<N-1; j++){
                phiGauss[i+j*N] = 1.0/4.0*(phiGauss[i+1+j*N]+phiGauss[i-1+j*N]+phiGauss[i+(j+1)*N]+phiGauss[i+(j-1)*N])+pow(h,2)/4*rho[i+j*N]; 
                phiJacob[i+j*N] = 1.0/4.0*(phitmp[i+1+j*N]+phitmp[i-1+j*N]+phitmp[i+(j+1)*N]+phitmp[i+(j-1)*N])+pow(h,2)/4*rho[i+j*N];
            }
        }

        for(int i=1; i<N-1; i++){
            for(int j=1; j<N-1; j++){
                phitmp[i+j*N] = phiJacob[i+j*N];
            }
        }
    }


    FILE* pFile;
    char filename[128];

    sprintf(filename, "PhiJacob.dat");
    pFile = fopen(filename,"w");
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            fprintf(pFile, "%f\t", phiJacob[i+j*N]);
        }
        fprintf(pFile,"\n");
    }

    memset(filename,0,sizeof(filename));
    
    sprintf(filename, "PhiGauss.dat");
    pFile = fopen(filename,"w");
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            fprintf(pFile, "%f\t", phiGauss[i+j*N]);
        }
        fprintf(pFile,"\n");
    }



    return 0;
}
