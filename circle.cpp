#include<cstdio>
#include<random>
#include<cmath>

double rand_r(double r,double R)
{
    return R*sqrt(r);
}

double rand_phi(double r)
{
    return 2*M_PI*r;
}

int main(int argc, char** argv)
{
    double r,phi;
    int N=1e4;
    FILE *pFile;
    double r_n,phi_n;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    pFile = fopen("circle.txt","w");

    for(int i=0;i<N;i++)
    {
        r_n = distribution(generator);
        r = rand_r(r_n,1.0);

        phi_n = distribution(generator);
        phi = rand_phi(phi_n);

        fprintf(pFile,"%f\t%f\n",phi,r);
    }


    return 0;
}


