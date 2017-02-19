#include<cstdio>
#include<cassert>
#include<cstdlib>
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
    int c=3;
    int p=101;
    int x0=1;
    int N=p-1;
    int x=x0;
    FILE *pFile;
    double r;
    double Rand;

    printf("c: ");
    scanf("%d",&c);
    printf("p: ");
    scanf("%d",&p);
    printf("x0: ");
    scanf("%d",&x0);
    printf("N: ");
    scanf("%d",&N);

    printf("RAND_MAX: %d\n",RAND_MAX);

    pFile = fopen("squaretest.txt","w");
    if(pFile==NULL) perror("Error opening file");

    for(int i=0;i<N;i++)
    {
        r = x*1.0/p;
        printf("%d\t%f\t",i,r);
        fprintf(pFile, "%f ",r);

        Rand = rand()*1./RAND_MAX;
        printf("%f\n",Rand);

        x=random(c,p,x);
        fprintf(pFile, "%f ",x*1.0/p);

        fprintf(pFile, "%f ",Rand);
        fprintf(pFile, "%f \n",rand()*1./RAND_MAX);
    }

    return 0;
}


