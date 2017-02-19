#include<cstdio>
#include<cstdlib>
#include<cassert>
#include"latticeview.h"

int main(int argc, char** argv)
{
    int N=10;
    int* plat;
    int i,j;
    int ImageWidth=1e3;
    int ImageHeight=1e3;
    double p=0.6;

    plat = (int*) calloc(N*N,sizeof(int));
    if (plat==NULL)
    {
        printf("allocation failure\n"); 
        exit(1);
    }

    //srand48(3);
    assert(p>=0&&p<1);

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            if(drand48() < p){
                plat[i+j*N] = 1;
            }
        }
    }


    Print_lattice(plat,N,N,ImageWidth,ImageHeight,"percolation.ppm");


    return 0;
}



