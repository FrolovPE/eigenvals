#include <iostream>
#include "lib.h"
#include <fenv.h>
int feenableexcept(int excepts);
int fedisableexcept(int excepts);
int fegetexcept(void);


int main(int argc, char **argv)
{
    int n,m,k,its{};
    double eps;
    double t1{},t2{},r1{},r2{};
    char *filename{};
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    if(!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n)==1 && sscanf(argv[2], "%d", &m)==1 && sscanf(argv[3], "%lf", &eps)==1 && sscanf(argv[4], "%d", &k)==1)) 
    {
        cout<<"Usage : "<<argv[0]<<" <n> "<<" <m> "<<" <eps> "<<" <k>"<<endl;
      
        return 0;
    }

    printf("n = %d m = %d eps = %e k = %d\n",n,m,eps,k);
    

    if(n<0 || eps<0 || k<0)
    {
        printf("<n> or <eps> or <k> <= 0, usage n,eps,k > 0");
        return 0;
    }

    double *a = new double[n*n];

    if(argc == 5 && k>=1 && k<=4)
    {
        init(a,n,k);
    }else if(argc == 5 && (k<1 || k>4))
    {
        r1 = -1; r2 = -1;
        printf("bad argument for k, k=1,2,3,4 \n");
        report(argv[0],r1,r2,t1,t2,its,n); 
        delete []a;
        return -1;
    }

    if(argc == 6 && k == 0 ) 
    {
        filename = argv[5];
        if(readarray(a,n,filename)<0)
        {
            return -1;
        }
    }

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            if(i == j) continue;
            if( fabs(a[i * n + j] - a[j * n + i]) > eps )
            {
                printf("Matrix is not symetric\n");
                r1 = -1;r2 = -1;
                delete []a;
                report(argv[0],r1,r2,t1,t2,its,n); 
                return -1;
            }
        }
    }

    


    printf("Matrix A;\n");
    printlxn(a,n,n,m);

    double trA = trace(a,n);

    threediag(a,n,eps*normofmatrix(a,n));

    double trdA = trace(a,n);

    printlxn(a,n,n,m);

    printf("trA = %lf trdA = %lf\n",trA,trdA);


    delete []a;
}