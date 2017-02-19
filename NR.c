#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

double GaussFunc(double x[2], double mu[2], double sig[2])
{
    double res;

    res = exp(-(pow(x[0]-mu[0],2)/(2*pow(sig[0],2))+pow(x[1]-mu[1],2)/(2*pow(sig[1],2))));

    return res;
}

void dGaussdX(double Dgauss[2], double x[2], double mu[2], double sig[2])
{
    double gauss = GaussFunc(x, mu, sig);

    Dgauss[0] = -(x[0]-mu[0])/pow(sig[0],2)*gauss;
    Dgauss[1] = -(x[1]-mu[1])/pow(sig[1],2)*gauss;

}

void d2GaussdX2(double Hgauss[4], double x[2], double mu[2], double sig[2])
{ 
    double gauss = GaussFunc(x, mu, sig);

    Hgauss[0] = (-1.0/pow(sig[0],2)+pow(x[0]-mu[0],2)/pow(sig[0],4))*gauss;
    Hgauss[1] = (x[0]-mu[0])*(x[1]-mu[1])/pow(sig[0],2)/pow(sig[1],2)*gauss;
    Hgauss[2] = (x[1]-mu[1])*(x[0]-mu[0])/pow(sig[1],2)/pow(sig[0],2)*gauss;
    Hgauss[3] = (-1.0/pow(sig[1],2)+pow(x[1]-mu[1],2)/pow(sig[1],4))*gauss;

}

void d2GaussdX2approx(double HgaussApprox[4], double x[2], double mu[2], double sig[2])
{
    double Dgauss[2];
    dGaussdX(Dgauss,x,mu,sig);
    double eps=pow(2,-53);
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            double Dgauss_dx[2];
            double hj = x[j]*sqrt(eps);
            double x_dx[2] = {x[0],x[1]};
            x_dx[j] += hj;
            dGaussdX(Dgauss_dx,x_dx,mu,sig);
            HgaussApprox[i+j*2] = (Dgauss_dx[i]-Dgauss[i])/hj;
        }
    }


}

void Newton(double x0[2], double mu[2], double sig[2], int maxsteps)
{
    double df[2];
    double d2f[4];
    double d2finv[4];
    double det;
    double x1[2];
    double eps = 1e-6;
    FILE* pFile;
    pFile = fopen("NewtonRaphson.txt","w");
    for(int step=0; step<maxsteps; step++){
        dGaussdX(df, x0, mu, sig);
        //d2GaussdX2(d2f,x0,mu,sig);
        d2GaussdX2approx(d2f,x0,mu,sig);
        det = d2f[0]*d2f[3]-d2f[1]*d2f[2];
        //printf("det: %f\n",det);
        assert(fabs(det)>0);
        d2finv[0] = 1.0/det*d2f[3];
        d2finv[1] = 1.0/det*-d2f[1];
        d2finv[2] = 1.0/det*-d2f[2];
        d2finv[3] = 1.0/det*d2f[0];

        x1[0] = x0[0] - (d2finv[0]*df[0]+d2finv[1]*df[1]);
        x1[1] = x0[1] - (d2finv[2]*df[0]+d2finv[3]*df[1]);
        
        printf("%d\t%f\t%f\n",step,x0[0],x0[1]);

        fprintf(pFile,"%d\t%f\t%f\n",step,x0[0],x0[1]);
        if(sqrt(pow(df[0],2)+pow(df[1],2))<eps){
            break;
        }
        else{
            x0[0] = x1[0];
            x0[1] = x1[1];
        }
    }

}

int main(int argc, char* argv[])
{
    double mu[2] = {0.0,0.0};
    double sig[2] = {sqrt(0.5),sqrt(0.5)};
    int maxsteps;
    double x0[2] = {sqrt(0.5)/2,sqrt(0.5)/2};
    
    scanf("%lf",&x0[0]);
    scanf("%lf",&x0[1]);
    scanf("%d", &maxsteps);
    
    Newton(x0,mu,sig,maxsteps);

    return 0;
}
