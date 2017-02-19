#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<algorithm>
#include<iterator>
//#include"latticeview.h"

//#define ImageWidth (1000)
//#define ImageHeight (1000)
int percolation(int* plat, int N, double p )
{
    int i,j;
    int count=0;

    if (plat==NULL)
    {
        printf("allocation failure\n"); 
        exit(1);
    }

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
    //Print_lattice(plat,N,N,ImageWidth,ImageHeight,"percolation.ppm");


    return count;
}

int Hoshen_Kopelman(int* plat, int N, int* M, int kmax, double p)
{
    int i,j;
    int* plat_g;
    int k,k0,k1,k2,ktmp;
    int* n;
    FILE* pFile0;
    FILE* pFile1;

    plat_g = (int*) calloc((N+1)*(N+1),sizeof(int));
    if(!plat_g){printf("allocate plat_g failure!\n");}

    pFile0 = fopen("./test/percolation.txt","w");
    for(i=1;i<N+1;i++){
        for(j=1;j<N+1;j++){
            plat_g[i+j*(N+1)] = plat[(i-1)+(j-1)*(N)];
            fprintf(pFile0,"%d\t",plat[(i-1)+(j-1)*N]);
        }
        fprintf(pFile0,"\n");
    }
    fclose(pFile0);
    k=1;
    for(i=1;i<N+1;i++){
        for(j=1;j<N+1;j++){
            if(plat_g[i+j*(N+1)]==1){
                if(0==plat_g[(i-1)+j*(N+1)]&&0==plat_g[i+(j-1)*(N+1)]){
                    //while(M[k]<0){
                    //    printf("i,j:\t%d\t%d\tk:\t%dM[k]:\t%d\n",i,j,k,M[k]);
                    //    assert(k!=-M[k]);
                    //    k=-M[k];
                    //}
                    k = k+1;
                    plat_g[i+j*(N+1)]=k;
                    M[k]=1;
                }
                if(plat_g[(i-1)+j*(N+1)]>0&&plat_g[i+(j-1)*(N+1)]==0){
                    k0 = plat_g[(i-1)+j*(N+1)];
                    while(M[k0]<0){
                        //printf("i,j:\t%d\t%d\tk0:\t%d\tM[k0]:\t\%d\n",i,j,k0,M[k0]);
                        assert(k0!=-M[k0]);
                        k0=-M[k0];
                    }
                    plat_g[i+j*(N+1)]=k0;
                    M[k0]++;
                    
                }
                if(plat_g[(i-1)+j*(N+1)]==0&&plat_g[i+(j-1)*(N+1)]>0){
                    k0 = plat_g[i+(j-1)*(N+1)];
                    while(M[k0]<0){
                        //printf("i,j:\t%d\t%d\tk0:\t%d\tM[k0]:\t%d\n",i,j,k0,M[k0]);
                        assert(k0!=-M[k0]);
                        k0=-M[k0];
                    }
                    plat_g[i+j*(N+1)]=k0;
                    M[k0]++;
                }
                if(plat_g[(i-1)+j*(N+1)]>0&&plat_g[i+(j-1)*(N+1)]>0){
                    k1 = plat_g[(i-1)+j*(N+1)];
                    k2 = plat_g[i+(j-1)*(N+1)];
                    while(M[k1]<0){
                        //printf("i,j:\t%d\t%d\tk1:\t%d\tM[k1]:\t%d\n",i,j,k1,M[k1]);
                        assert(k1!=-M[k1]);
                        k1=-M[k1];
                    }
                    while(M[k2]<0){
                        //printf("i,j:\t%d\t%d\tk2:\t%d\tM[k2]:\t%d\n",i,j,k2,M[k2]);
                        assert(k2!=-M[k2]);
                        k2=-M[k2];
                    }

                    if(k1!=k2){
                        ktmp = std::min(k1,k2);
                        k2 = std::max(k1,k2);
                        k1 = ktmp;
                        plat_g[i+j*(N+1)]=k1;
                        M[k1] += M[k2]+1;
                        M[k2] = -k1;
                    }
                    else{
                        k0 = k1;
                        //while(M[k0]<0){
                        //    //printf("i,j:\t%d\t%d\tk0:\t%d\tM[k0]:\t%d\n",i,j,k0,M[k0]);
                        //    k0=-M[k0];
                        //}
                        plat_g[i+j*(N+1)]=k0; 
                        M[k0]++;
                    }
                }

            }
            //printf("%d\t",plat_g[i+j*(N+1)]);
        }
        //printf("\n");
    }

    int maxsize = 0;
    int maxlabel = 0;
    for(k=0;k<kmax;k++){
        if(M[k]>0&&M[k]>maxsize){
            maxsize=M[k];
            maxlabel=k;
        }
    }
    n = (int*) calloc(maxsize+1,sizeof(int));
    if(!n){printf("allocate n failure!\n");}
    //printf("label\tsize\n");
    for(k=2;k<kmax;k++){
        if(M[k]>0){
            n[M[k]]++;
            //printf("%d\t%d\n",k,M[k]);
        }
    }
    pFile1 = fopen("./test/cluster0.txt","w");
    for(i=1;i<N+1;i++){
        for(j=1;j<N+1;j++){
            k = plat_g[i+j*(N+1)];
            fprintf(pFile1,"%d\t",k);
        }
        fprintf(pFile1,"\n");
    }
    fclose(pFile1);
    pFile1 = fopen("./test/cluster1.txt","w");
    for(i=1;i<N+1;i++){
        for(j=1;j<N+1;j++){
            k = plat_g[i+j*(N+1)];
            while(M[k]<0){
                k=-M[k];
            }
            plat[(i-1)+(j-1)*N] = k; 
            fprintf(pFile1,"%d\t",plat[(i-1)+(j-1)*N]);
        }
        fprintf(pFile1,"\n");
    }
    fclose(pFile1);
    char sizedist[128];
    sprintf(sizedist,"./test/sizedist_N%d_p%.2f.txt",N,p);
    pFile0 = fopen(sizedist,"w");
    printf("size\tnum\n");
    for(i=0;i<=maxsize;i++){
        if(n[i]>0){
            fprintf(pFile0,"%d\t%d\n",i,n[i]);
            printf("%d\t%d\n",i,n[i]);
        }
    }
    fclose(pFile0);
    //Print_lattice(plat,N,N,ImageWidth,ImageHeight,"cluster.ppm");
    free(n);
    free(plat_g);

    printf("label:%d\tmaxsize:%d\n",maxlabel,M[maxlabel]);
    return maxlabel;
}

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

int sandboxtest(int N, double p)
{
    //int N=1e4;
    //double p=0.6;
    int* plat;
    int kmax;
    int* M;
    int k;
    int res;
    FILE* pFile;

    plat = (int*) calloc(N*N,sizeof(int));
    
        
    kmax = percolation( plat, N,  p );

    M = (int*) calloc(kmax+1,sizeof(int));
    if(M==NULL){
        printf("allocate M failure!\n");
        exit(EXIT_FAILURE);
    }

    k=Hoshen_Kopelman( plat, N,  M, kmax+1, p);
    printf("Sandbox\n");
    printf("R\tM\n");
    pFile = fopen("./test/Sandbox_small.txt","w");
    for(int R=1;R<=N/2;R*=2){
        res = sandbox(plat, N, N/2, R, k);
        fprintf(pFile,"%d\t%d\n",R,res);
        printf("%d\t%d\n",R,res);
    }
    fclose(pFile);


    free(M);
    free(plat);
    return 0;
} 

int boxcounttest(int N, double p)
{
    //int N=1e4;
    //double p=0.6;
    int* plat;
    int kmax;
    int* M;
    int k;
    int res;
    FILE* pFile;

    plat = (int*) calloc(N*N,sizeof(int));
    
        
    kmax = percolation( plat, N,  p );

    M = (int*) calloc(kmax+1,sizeof(int));
    if(M==NULL){
        printf("allocate M failure!\n");
        exit(EXIT_FAILURE);
    }

    k=Hoshen_Kopelman( plat, N,  M, kmax+1, p);
    printf("Boxcount\n");
    printf("R\tM\n");
    pFile = fopen("./test/Boxcount_small.txt","w");
    for(int R=1;R<=N;R*=2){
        res = boxcount(plat, N, R, k);
        fprintf(pFile,"%d\t%d\n",R,res);
        printf("%d\t%d\n",R,res);
    }
    fclose(pFile);


    free(M);
    free(plat);
    return 0;
} 

int main(int argc, char** argv)
{
    int N=1e4;
    double p=0.6;

    boxcounttest(N,p);
    //sandboxtest(N,p);

    return 0;
}
