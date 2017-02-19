#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<assert.h>
//#include<omp.h>
//#include<mpi.h>

void Metropolis(int* S, int L, double J, double kB, double T, double* M, double* E)
{
    double E1=0;
    double M1=0; 
    for(int k=0; k<L; k++){
        for(int l=0; l<L; l++){
            double e0=0;
            e0 = -J*(S[k+l*L]*S[(k-1+L)%L+l*L]+S[k+l*L]*S[(k+1+L)%L+l*L]+ \
                      S[k+l*L]*S[k+((l-1+L)%L)*L]+S[k+l*L]*S[k+((l+1+L)%L)*L]);
        
            // single spin slip
            S[k+l*L] = -S[k+l*L];
            double e1=0;
            e1 = -J*(S[k+l*L]*S[(k-1+L)%L+l*L]+S[k+l*L]*S[(k+1+L)%L+l*L]+ \
                      S[k+l*L]*S[k+((l-1+L)%L)*L]+S[k+l*L]*S[k+((l+1+L)%L)*L]);
        
            if(e1<e0){
                continue;
            }else{
                double r = rand()*1.0/RAND_MAX;
                double delta = exp(-(e1-e0)/(kB*T));
                double c = (delta<1) ? delta:1;
                if(r<c){
                    continue;
                }else{
                    S[k+l*L] = -S[k+l*L];
                }
            }

            M1 += S[k+l*L];
            E1 += -J*(S[k+l*L]*S[(k-1+L)%L+l*L]+S[k+l*L]*S[(k+1+L)%L+l*L]+ \
                     S[k+l*L]*S[k+((l-1+L)%L)*L]+S[k+l*L]*S[k+((l+1+L)%L)*L]);

        }
    }
                
    //printf("M1: %f\t E1: %f\n", M1, E1);
    *M = M1;
    *E = E1;
}



int main(int argc, char** argv)
{
        
    int L=32;
    int* S = (int*) malloc(L*L*sizeof(int));
    int N=1e2;
    double E=0;
    double M=0;
    double T_start=0;
    double T_end=5;
    unsigned int nT = 1e2;
    double dT = (T_end-T_start)/nT;
    FILE* pFile;
    char filename[128];
    sprintf(filename,"IsingModel2D%dx%d.txt",L,L);
    pFile = fopen(filename,"w");
    // Initialize
    //printf("=== Initialization ===\n");
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++){
            S[i+j*L]=1;
        }
    }
    // Calculation
    for(double T=T_start; T<T_end; T+=0.01){
        // Relaxation
        //printf("=== Relaxatoin ===\n");
        for(int n=0;n<N;n++){
            double M_relax;
            double E_relax;
            Metropolis(S, L, 1.0, 1.0, T,&M_relax,&E_relax);
        }
        //printf("%f\t %f\n",M_relax,E_relax);

        // Sweep
        //printf("=== Calcualtion ===\n");
        M=0;
        E=0;
        for(int n=0; n<N; n++){
            double M0;
            double E0;
            Metropolis(S, L, 1.0, 1.0, T,&M0,&E0);
            M += M0;
            E += E0;        
        }
        M = M/N/(L*L);
        E = E/N/(L*L*2);
        printf("%f\t%f\t%f\n",T,M,E);
        fprintf(pFile,"%f\t%f\t%f\n",T,M,E);
    }

    fclose(pFile);
    return 0;
}

