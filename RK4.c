#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<omp.h>

void ydot(double t, double* y, double* ydot, int N, double gamma, double g)
{
    // xdot = vx
    // zdot = vz
    // vxdot = -gamma*v*vx
    // vydot = -g-gamma*v*vz
    double v = sqrt(pow(y[2],2)+pow(y[3],2));
    ydot[0] = y[2];
    ydot[1] = y[3];
    ydot[2] = -gamma*v*y[2];
    ydot[3] = -g-gamma*v*y[3];



}


void RK4(int N, double t0, double t1, double h, double y0[N], double gamma, double g, char filename[128])
{
    double k1[N];
    double k2[N];
    double k3[N];
    double k4[N];
    double y[N];
    double y1[N];
    double y2[N];
    double y3[N];

    //FILE* pFile;

    //pFile = fopen(filename,"w");
    for(int i=0; i<N; i++){
        y[i] = y0[i];
        y1[i] = y0[i];
        y2[i] = y0[i];
        y3[i] = y0[i];
        //fprintf(pFile, "%f\t", y[i]);
        //printf("%f\t", y[i]);
    }
    //fprintf(pFile,"\n");
    //printf("\n");
    for(double t = t0; t<t1; t+= h)
    {
        ydot(t, y, k1, N, gamma, g);
        for(int i=0; i<N; i++){
            y1[i] = y[i] + h/2.0*k1[i];
        }
        ydot(t+h/2.0, y1, k2, N, gamma, g);
        for(int i=0; i<N; i++){
            y2[i] = y[i] + h/2.0*k2[i];
        }
        ydot(t+h/2.0, y2, k3, N, gamma, g);
        for(int i=0; i<N; i++){
            y3[i] = y[i] + h*k3[i];
        }
        ydot(t+h, y3, k4, N, gamma, g);
        for(int i=0; i<N; i++){
            y[i] += h/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
            //fprintf(pFile, "%f\t", y[i]);
            //printf("%f\t", y[i]);
        }
        //fprintf(pFile,"\n");
        //printf("\n");
        if(y[1]<y0[1]){
            break;
        }
    }

    for(int i=0; i<N; i++){
        y0[i] = y[i];
    }

    //fclose(pFile);
}

int main(int argc, char* argv[])
{
    int N=4;
    double x0 = 0;
    double z0 = 0;
    double v0 = 40;
    double alpha = M_PI/4;
    double gamma = 0.7; // in [0,5]
    double g = 9.8;
    double t0 = 0;
    double t1 = 100; //1.2;
    double h = 1e-4;
    double y[N];
    int Ngamma=1e2;
    double alpha_gamma[Ngamma];
    double hgamma = 5.0/Ngamma;
    int Nalpha=1e3;
    double halpha = M_PI/2.0/Nalpha;
    char filename[128]; 
    sprintf(filename, "Trajectory.log");
    FILE* pFile;
    pFile = fopen("alpha_gamma.txt","w");

    for(int j=0; j<Ngamma; j++){
        double xrange=0;
        double alphamax=0;
        gamma = j*hgamma;
        alpha_gamma[j] = alphamax;
        //printf("-----------gamma=%f--------------\n",gamma);
#pragma omp parallel for private(y,alpha) shared(xrange, alphamax)
        for(int k=0; k<Nalpha; k++){
            alpha = k*halpha;
            double vx = v0*cos(alpha);
            double vz = v0*sin(alpha);
            y[0] = x0;
            y[1] = z0;
            y[2] = vx;
            y[3] = vz;
            //sprintf(filename, "Trajectory_%0.2f_%0.2f_%0.2f.txt",v0, alpha, gamma);
            
            RK4(N, t0, t1, h, y, gamma, g, filename);
            //printf("%f\t%f\n",alpha,y[0]-x0);
            if(y[0]-x0>xrange){
                xrange=y[0]-x0;
                alphamax=alpha;
            }

        }
        alpha_gamma[j] = alphamax;
        
        printf("%f\t%f\n",gamma,alphamax);
    }

    for(int j=0; j<Ngamma; j++){
        gamma = j*hgamma;
        fprintf(pFile,"%f\t%f\n",gamma, alpha_gamma[j]);
    }
    fclose(pFile);
    return 0;
}
