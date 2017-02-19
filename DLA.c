#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<assert.h>
#include"latticeview.h"
#define ImageWidth (1000)
#define ImageHeight (1000)
int sandbox(int* plat, int N, unsigned int center, unsigned int R, int maxsizelabel)
{
    int i,j;
    int ic,jc;
    int max;
    int count=0;

    if(plat[center+center*N]>0){
        center=center;
    }else if(plat[(center-1)+(center-1)*N]>0){
        center=center-1;
    }else if(plat[(center+1)+(center+1)*N]>0){
        center=center+1;
    }else{
        printf("center:%d\toccupation:%d\n",center,plat[center+center*N]);
        assert(plat[center+center*N]>0);
    }
            
    assert(R<=N/2);

    for(i=0;i<N;i++){
        ic = i-center;
        for(j=0;j<N;j++){
            jc = j-center;
            max = (abs(ic)>abs(jc)) ? abs(ic):abs(jc);
            if(max<=R&&plat[i+j*N]==maxsizelabel){
                count++;
            }
        }
    }

    return count;
}
int boxcount(int* plat, int N, unsigned int R, int maxsizelabel)
{
    int* box;
    int length=N/R;
    int i,j;
    int i_b,j_b;
    int i_r,j_r;
    int count;

    box = (int*) calloc(R*R,sizeof(int));
    if(box==NULL){printf("allocate box failure\n");exit(EXIT_FAILURE);}
    
    for(i=0;i<N;i++){
        i_b = i/length;
        for(j=0;j<N;j++){
            j_b = j/length;
            if(i_b<R&&j_b<R&&plat[i+j*N]==maxsizelabel){
                box[i_b+R*j_b]++;
            }
        }
    }

    for(i_r=0;i_r<R;i_r++){
        for(j_r=0;j_r<R;j_r++){
            if(box[i_r+j_r*R]>0){
                count++;
            }
        }
    }

    free(box);
    return count;
            
}

int DLA(int* plat, int N, unsigned int m)
{
    //int N=100;
    //int* plat = (int*) calloc(N*N,sizeof(int));
    int x_c = N/2;
    int y_c = N/2;
    //unsigned int m=1;
    //center seed
    plat[x_c+y_c*N] = m;

    int dR=3;
    int R=0;
    R+=dR;
    int Rmax=2*R<N/2?2*R:N/2;
    while(R<Rmax){
        double theta = rand();///RAND_MAX*2*M_PI;
        double x = x_c+R*cos(theta);
        double y = y_c+R*sin(theta);
        int x_int = (int) x;
        int y_int = (int) y;
        //printf("%d\t%d\n",x_int,y_int); 
        while(pow(x_int-x_c,2)+pow(y_int-y_c,2)<pow(Rmax,2)){
            int random_walk = rand()%4;
            switch (random_walk){
            case 0:
                y_int++;
                break;
            case 1:
                x_int--;
                break;
            case 2:
                y_int--;
                break;
            case 3:
                x_int++;
                break;
            }
            
            if(plat[(x_int-1)+y_int*N]==m||plat[x_int+(y_int-1)*N]==m||plat[(x_int+1)+y_int*N]==m||plat[x_int+(y_int+1)*N]==m){
                plat[x_int+y_int*N]++;
                assert(plat[x_int+y_int*N]<=m);
                //printf("%d\t%d\n",x_int,y_int);
                double dist = sqrt(pow(x_int-x_c,2)+pow(y_int-y_c,2));
                if(dist+dR>=R){
                    R += dR;
                    Rmax = 2*R<N/2 ? 2*R: N/2;
                }

                break;
            }
        }
        

    }
    
    FILE* pFile;
    pFile = fopen("DLA.txt","w");
    for(int i=0;i<N;++i){
        for(int j=0;j<N;++j){
            //printf("%d\t",plat[i+j*N]);
            fprintf(pFile,"%d\t",plat[i+j*N]);
        }
        //printf("\n");
        fprintf(pFile,"\n");
    }
    fclose(pFile);

    Print_lattice(plat,N,N,ImageWidth,ImageHeight,"DLA.ppm");
    
    return 0;
}

int main(int argc, char** argv)
{
    int N=600;
    //int* plat = (int*) malloc(N*N*sizeof(int));
    int plat[N*N];
    assert(plat!=NULL);
    unsigned int m=2;
    int res;
    int x_c = N/2;
    int y_c = N/2;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            plat[i+j*N]=0;
        }
    }
    //center seed
    plat[x_c+y_c*N] = m;

    int dR=3;
    int R=0;
    R+=dR;
    int Rmax=2*R<N/2?2*R:N/2;
    while(R<Rmax){
        double theta = rand();///RAND_MAX*2*M_PI;
        double x = x_c+R*cos(theta);
        double y = y_c+R*sin(theta);
        int x_int = (int) x;
        int y_int = (int) y;
        //printf("%d\t%d\n",x_int,y_int); 
        while(pow(x_int-x_c,2)+pow(y_int-y_c,2)<pow(Rmax,2)){
            int random_walk = rand()%4;
            switch (random_walk){
            case 0:
                y_int++;
                break;
            case 1:
                x_int--;
                break;
            case 2:
                y_int--;
                break;
            case 3:
                x_int++;
                break;
            }
            
            if(plat[(x_int-1)+y_int*N]==m||plat[x_int+(y_int-1)*N]==m||plat[(x_int+1)+y_int*N]==m||plat[x_int+(y_int+1)*N]==m){
                plat[x_int+y_int*N]++;
                assert(plat[x_int+y_int*N]<=m);
                //printf("%d\t%d\n",x_int,y_int);
                double dist = sqrt(pow(x_int-x_c,2)+pow(y_int-y_c,2));
                if(dist+dR>=R){
                    R += dR;
                    Rmax = 2*R<N/2 ? 2*R: N/2;
                }

                break;
            }
        }
    }
    
    FILE* pFile;
    pFile = fopen("DLA.txt","w");
    for(int i=0;i<N;++i){
        for(int j=0;j<N;++j){
            //printf("%d\t",plat[i+j*N]);
            fprintf(pFile,"%d\t",plat[i+j*N]);
        }
        //printf("\n");
        fprintf(pFile,"\n");
    }
    fclose(pFile);

    Print_lattice(plat,N,N,ImageWidth,ImageHeight,"DLA.ppm");
    
    pFile = fopen("DLA_Sandbox.txt","w");
    
    for(int R=1;R<=N/2;R*=2){
        res = boxcount(plat, N, R, m);
        fprintf(pFile,"%d\t%d\n",R,res);
        printf("%d\t%d\n",R,res);
    }
    fclose(pFile);

    return 0;
}
        
                
            
            


