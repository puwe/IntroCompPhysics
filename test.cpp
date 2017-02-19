#include<cstdio>
#include<cstdlib>
#include<cassert>
#include"latticeview.h"

#define ImageWidth (1000)
#define ImageHeight (1000)
int percolation(int* plat, int N, double p )
{
    //int N=10;
    //int* plat;
    int i,j;
    //int ImageWidth=1e3;
    //int ImageHeight=1e3;
    //double p=0.6;
    int count=0;

    //plat = (int*) calloc(N*N,sizeof(int));
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
                count++;
            }
            else{
                plat[i+j*N] = 0;
            }
        }
    }

    printf("Number of forests: %d\n",count);
    Print_lattice(plat,N,N,ImageWidth,ImageHeight,"percolation.ppm");


    return count;
}

void forestfire(int* plat, int N, int num, int* pspan, int* plife)
{
    int* plat_g;
    int gi,gj;
    int i,j;
    int span;
    int life;
    int t;
    int count=0;
    int count_tmp=0;
    // ghost lattice
    plat_g = (int*) calloc((N+2)*(N+2),sizeof(int));
     
    for(gi=1;gi<N+1;gi++){
        i = gi - 1;
        for(gj=1;gj<N+1;gj++){
            j = gj - 1;
            plat_g[gi+gj*(N+2)] = plat[i+j*N];
        }
    }
    // burn fire at first row
    gj = 1;
    for(gi=1;gi<N+1;gi++){
        if(1==plat_g[gi+gj*(N+2)]){ // occupied
            plat_g[gi+gj*(N+2)] = 2;
        }
    }
    // burning
    for(t=0;t<num;t++){
        count = 0;
        for(gj=1;gj<N+1;gj++){
            for(gi=1;gi<N+1;gi++){
                // fire
                if(plat_g[gi+gj*(N+2)]>1){
                    // check surrounding place whether there is forest occupied and burn it
                    if(1==plat_g[gi+(gj-1)*(N+2)]){
                        plat_g[gi+(gj-1)*(N+2)] = plat_g[gi+gj*(N+2)] + 1;
                    }
                    if(1==plat_g[gi+(gj+1)*(N+2)]){
                        plat_g[gi+(gj+1)*(N+2)] = plat_g[gi+gj*(N+2)] + 1;
                    }
                    if(1==plat_g[(gi+1)+gj*(N+2)]){
                        plat_g[(gi+1)+gj*(N+2)] = plat_g[gi+gj*(N+2)] + 1;
                    }
                    if(1==plat_g[(gi-1)+gj*(N+2)]){
                        plat_g[(gi-1)+gj*(N+2)] = plat_g[gi+gj*(N+2)] + 1;
                    }
                    // unburned forest number reduced
                    count++;
                }
            }
        }
        if(count_tmp!=count){
            count_tmp=count;
        }
        else{
            // fire stops
            printf("fire stops: %d\n",t);
            break;
        }

    }
    // check if spann
    gj = N;
    span = 1;
    for(gi=1;gi<N+1;gi++){
        if(plat_g[gi+gj*(N+2)]>span){
            span = plat_g[gi+gj*(N+2)];
        }
    }
    // life time
    life = 1;
    for(gi=1;gi<N+1;gi++){
        for(gj=1;gj<N+1;gj++){
            if(plat_g[gi+gj*(N+2)]>life){
                life = plat_g[gi+gj*(N+2)];
            }
        }
    }
    
    printf("Span: %d\n",span);
    printf("Life: %d\n",life);

    // return to lattice
    for(gi=1;gi<N+1;gi++){
        i = gi - 1;
        for(gj=1;gj<N+1;gj++){
            j = gj - 1;
            plat[i+j*N] = plat_g[gi+gj*(N+2)];
        }
    }

    *pspan = span;
    *plife = life;

    free(plat_g);

}

double test(int N, double p, int ntest)
{
    //int N=10;
    //double p=0.7;
    int* plat;
    int span;
    int life;
    int num;
    //int ntest=50;
    char filename[128];
    int* vecspan;
    int* veclife;
    FILE* pFile;
    char outputfile[128];

    plat = (int*) calloc(N*N,sizeof(int));
    vecspan = (int*) calloc(ntest,sizeof(int));
    veclife = (int*) calloc(ntest,sizeof(int));
    
    assert(ntest>0);

    sprintf(outputfile,"./test/%dtests_N%d_p%.2f.txt",ntest,N,p);
    pFile = fopen(outputfile,"w");

    for(int t=0;t<ntest;t++){
        num = percolation(plat,N,p);
        forestfire(plat,N,num,&span,&life);
        vecspan[t] = span;
        veclife[t] = life;

        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                if(plat[i+j*N]>2&&plat[i+j*N]!=life){
                    plat[i+j*N]=3; //black
                }
                else if(plat[i+j*N]==life){
                    plat[i+j*N]=4; //blue
                }
            }
        }
        fprintf(pFile,"%d\t%d\t%d\n",t,span,life);
        sprintf(filename,"./test/test%03d.ppm",t);
        Print_lattice(plat,N,N,ImageWidth,ImageHeight,filename);
    
    }

    double span_avg=0;
    double life_avg=0;
    int span_count=0;
    for(int t=0;t<ntest;t++){
        span_avg += vecspan[t];
        life_avg += veclife[t];
        if(vecspan[t]>1){
            span_count++;
        }
    }
    span_avg = span_avg/ntest;
    life_avg = life_avg/ntest;

    printf("ntest: %d\t size: %d\t probability: %f\n",ntest,N,p);
    fprintf(pFile,"%d\t%f\t%f\n",span_count,span_avg,life_avg);
    printf("span_num: %d\t span_avg: %f\t life_avg: %f\t\n",span_count,span_avg,life_avg);
    

    free(plat);
    free(vecspan);
    free(veclife);

    return span_count*1.0/ntest;

}

int main(int argc, char** argv)
{
    int N;
    double p;
    int ntest;
    double res;
    FILE *pFile;
    
    pFile = fopen("test.txt","w");
    N=100;
    ntest=50;
    for(p=0.1;p<1;p+=0.1){
        res=test(N,p,ntest);
        fprintf(pFile,"%f\t%f\n",p,res);
    }

    return 0;
}    
