#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cassert>

// congruential random number generator
// p is prime number
double random(int c, int p, int x0)
{
    int x;
    
    assert(x0!=0);

    x = (c*x0)%p;

    return x;
}


int main(int argc, char** argv)
{
    int k=5;
    int c=3;
    int p=pow(2,17)-1;
    int x0=1;
    int* pN;
    int x=x0;
    int N=1e4;
    int l = (p-1)/k;
    int n;
    FILE *pFile;
    
    double chi2=0;

    pN = (int*) calloc(k,sizeof(int));
    pFile = fopen("chi2.txt","w");

    for(int i=0;i<N;i++)
    {
        n = x/l;
        //printf("%d\t%d\n",n,x);
        pN[n]++;
        x = random(c,p,x);
    }

    for(int i=0;i<k;i++)
    {
        chi2 += pow(pN[i]-N*1.0/k,2)*1.0/(N*1.0/k);
        printf("%d\t%d\n",i,pN[i]);
        fprintf(pFile,"%d\t%f\n",i,pN[i]*1.0/N);
    }
    printf("c: %d p: %d x0: %d binwidth: %d chi2: %f\n",c,p,x0,l,chi2);

    return 0;
}






