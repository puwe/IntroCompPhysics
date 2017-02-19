#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<stdbool.h>
#include<math.h>
#include<assert.h>
double Chain(int N, double r, int M)
{
    double* Spheres = (double*) malloc(2*N*sizeof(double));
    //FILE* pFile;
    //char filename[128];

    if(Spheres==NULL){
        printf("Allocation Failure\n");
        exit(EXIT_FAILURE);
    }

    double theta = drand48()*2*M_PI;
    //double phi = acos(-2*drand48()+1);
    double x0 = 0;
    double y0 = 0;
    //double z0 = 0;
    double x = 0;
    double y = 0;
    //double z = 0;
    //Spheres[0] = x0;
    //Spheres[1] = y0;
    //Spheres[2] = z0;
    //x = x0 + 2*r*sin(phi)*cos(theta);
    //y = y0 + 2*r*sin(phi)*sin(theta);
    //z = z0 + 2*r*cos(phi);
    //Spheres[3] = x;
    //Spheres[4] = y;
    //Spheres[5] = z;
    Spheres[0] = x0;
    Spheres[1] = y0;
    x = x0 + 2*r*cos(theta);
    y = y0 + 2*r*sin(theta);
    int maxsteps = N;
    bool overlap = true;
    int num=N*10;
    int n=0;
    //sprintf(filename, "./test/SphereChain_N%d_M%d", N, M);
    //pFile = fopen(filename, "w");
    
LOOP:for(int i=2; i<N&&n<num; i++){
    
        overlap = true;
        for(int step=0; step<maxsteps&&overlap==true; step++){
            theta = drand48()*2*M_PI;
            //phi = acos(-2*drand48()+1);
            
            //x = Spheres[3*(i-1)+0] + 2*r*sin(phi)*cos(theta);
            //y = Spheres[3*(i-1)+1] + 2*r*sin(phi)*sin(theta);
            //z = Spheres[3*(i-1)+2] + 2*r*cos(phi);
            x = Spheres[2*(i-1)+0] + 2*r*cos(theta);
            y = Spheres[2*(i-1)+1] + 2*r*sin(theta);
            overlap=false;
            for(int j=0; j<i-1; j++){
                //if(pow(x-Spheres[3*j+0],2)+pow(y-Spheres[3*j+1],2)+pow(z-Spheres[3*j+2],2)<pow(2*r,2)){
                if(pow(x-Spheres[2*j+0],2)+pow(y-Spheres[2*j+1],2)<pow(2*r,2)){    
                    overlap=true;
                }
            }
        }
        if(overlap==false){
            //Spheres[3*i+0] = x;
            //Spheres[3*i+1] = y;
            //Spheres[3*i+2] = z;
            //printf("%d\t%d\t%f\t%f\t%f\n",N,i,x,y,z);
            Spheres[2*i+0] = x;
            Spheres[2*i+1] = y;
        }
        else{
            //printf("Reach max steps %d\n", maxsteps);
            n++;
            goto LOOP;
        }

    }
    if(n==num){
        printf("Chain failure\n");
    }
    //x = Spheres[3*(N-1)+0];
    //y = Spheres[3*(N-1)+1];
    //z = Spheres[3*(N-1)+2];
    //double R2 = pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2);
    //printf("%d\t%f\t%f\t%f\n",N,x,y,z);
    x = Spheres[2*(N-1)+0];
    y = Spheres[2*(N-1)+1];
    double R2 = pow(x-x0,2)+pow(y-y0,2);
    free(Spheres);
    return R2;
}


int main(int argc, char* argv[])
{
    int M = 1e4;
    int N = 10;
    double r=1;
    double R2;
    double sum=0;
    double sum2=0;
    double mean=0;
    double error=0;
    char filename[128];
    FILE* pFile;
    
    sprintf(filename,"SphereChain2DM%d",M);
    pFile = fopen(filename,"w");
    srand48((unsigned)time(NULL));
    
    for(N=3; N<1e2; N++){
        sum = 0;
        sum2 = 0;
        for(int i=0; i<M; i++){
            R2 = Chain(N,r,i);
            sum += R2;
            sum2 += pow(R2,2);
            //printf("%d\t%f\n", i, R2);
        }
        mean = sum/M;
        error = sqrt(1.0/M*(sum2/M-pow(mean,2)));
        fprintf(pFile,"%d\t%f\n",N,mean);
        printf("%d\t%f\t%f\t%f\n",N,mean,error,error/mean);
    }
    fclose(pFile);
    return 0;
}
