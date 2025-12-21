#include <iostream>
#include "lib.h"
#include <time.h>
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
            r1 = -1; r2 = -1;
            report(argv[0],r1,r2,t1,t2,its,n); 
            delete []a;
            return -1;
        }
    }else if(argc == 6 && k!=0)
    {
        r1 = -1; r2 = -1;
        report(argv[0],r1,r2,t1,t2,its,n); 
        delete []a;
        return -1;
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

    
    double norma = normofmatrix(a,n);
    for(int i = 0; i < n; i++)
    {
        for(int j = 0 ; j < n ;j++)
        {
            r2 += a[i * n + j] * a[j * n + i];
        }
    }

    printf("Matrix A;\n");
    printlxn(a,n,n,m);

    double trA = trace(a,n);
    t1 = clock();
    threediag(a,n,eps*normofmatrix(a,n));
    t1 = (clock() - t1)/CLOCKS_PER_SEC;

    // double trdA = trace(a,n);

    // printlxn(a,n,n,m);


    // printf("\nd: ");
    // for(int i = 0 ; i < n ; i++)
    // {
    //     printf("%lf ",d[i]);
    // }
    // printf("\nund:  ");
    // for(int i = 0 ; i < n - 1 ; i++)
    // {
    //     printf("%lf ",und[i]);
    // }

    // printf("\ntrA = %lf trdA = %lf\n",trA,trdA);

    t2 = clock();
    QR_refl(a,n,eps * norma,its);
    t2 = (clock() - t2)/CLOCKS_PER_SEC;

    r1 = (fabs(trA - trace(a,n)))/norma;

    double ll = 0;
    for(int i = 0 ; i < n; i++)
    {
        ll += a[i * n + i] * a[i * n + i];
    }

    // printf("sqrt(r2) = %lf sqrt(ll) = %lf\n",sqrt(r2),sqrt(ll));
    
    r2 = fabs(sqrt(r2) - sqrt(ll))/norma;

    // printlxn(a,n,n,m);
    printf("Eigen vals:\n");
    for(int i = 0; i < min(n,m); i++)
    {
        printf("%10.3e ",a[i * n + i]);
    }
    printf("\n");

    report(argv[0],r1,r2,t1,t2,its,n);

    delete []a;
}
