#include <iostream>
#include "lib.h"
using namespace std;

void printlxn(const double *a, int l, int n, int r)
{
    int l1 = (l < r ? l : r), n1 = (n < r ? n : r);
    int i,j;

    for(i = 0; i < l1; ++i)
    {
        printf("\n");
        for(j = 0; j < n1; ++j)
        {
            printf("%10.3e ",a[i * n + j]);
        }
    }

    printf("\n");

}



int readarray(double *a, int n, char* filename){
        int _c=0;
        double el;

        FILE *file = fopen(filename,"r");
        if(!file){
            
            printf("File %s doesnt exist or wrong file name!\n",filename);
            // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
            // delete []a;
            // delete []b;
            // delete []x;
            // delete []realx;
            
            return -1;
        }

        while(fscanf(file,"%lf",&el)==1)
        {
            if(_c<n*n) 
            {
                a[_c] = el;
                _c++;
            }else
            {
                printf("Bad scan from file %s\n",filename);
                // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
                // delete []a;
                // delete []b;
                // delete []x;
                // delete []realx;
                fclose(file);
                return -2;
            }

        
        }
        if(!feof(file) || _c!=n*n){
            printf("Bad file %s\n",filename);
            // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
            fclose(file);
            // delete []a;
            // delete []b;
            // delete []x;
            // delete []realx;
            return -3;
        }
        fclose(file);
        return 0;
    }


double f (int k , int n , int i , int j)
{
    if(k == 1) 
        return n - max(i+1,j+1) +1;
    else if (k == 2)
        if(i == j)
            return 2;
        else if(abs(i-j) == 1)
            return -1;
        else
            return 0;
    else if (k == 3)
        if( i == j && i < n)
            return 1;
        else if(j == n)
            return i;
        else if( i == n)
            return j;
        else
            return 0;
    else if (k == 4)
        return static_cast<double>(1) / (static_cast<double>(i+1) + static_cast<double>(j+1) - static_cast<double>(1));
    else
        return -1;
}

void init(double *a, int n, int s )
{
    if(!a)
    {
        printf("a = nullptr\n");
        return;
    }

    for(int i = 0; i < n ; i++)
    {
        for(int j=0 ; j < n; j++)
        {
            a[i*n+j] = f(s,n,i,j);
        }
    }
}

double vectornorm(double *a , int n)
{
    double res = 0;

    if(!a)
    {
        printf("nullptr in vector norm\n");
        return -1;
    }
    for(int i = 0; i<n; i++)
    {
        res += fabs(a[i]);
    }

    return res;
}

void mat_x_vector(double *res,double *a, double *b, int n)
{
    

    if(!a || !b)
    {
        return;
    }

    for(int i = 0; i< n; i++)
    {
        double s = 0;
        for(int j = 0; j < n; j++)
        {
            s += a[i*n +j] * b[j];
        }
        res[i] = s;
    }

  
}

void vectorsub(double *res,double *a,double *b, int n)
{
    // double *res = new double[n];

    if(!a || !b)
    {
        printf("nullptr in vector subtract\n");
        return;
    }

    for(int i = 0 ; i< n ; i++)
    {
        
        res[i] = a[i] - b[i];
    }

    // return res;
}



void residuals(double &r1,double &r2,double *a,double *b,double *x,double *realx,int n,double *Ax, double *Ax_b, double *x_realx)
{
    mat_x_vector(Ax,a,x,n);// double *Ax = mat_x_vector(a,x,n);
    vectorsub(Ax_b,Ax,b,n);// double *Ax_b = vectorsub( Ax , b, n);
    vectorsub(x_realx,x,realx,n);// double *x_realx = vectorsub( x , realx, n);

    // cout<<"\nVector Ax :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",Ax[i]);
    // }

    // cout<<"\nVector b :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",b[i]);
    // }

    // cout<<"\nVector Ax_b :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",Ax_b[i]);
    // }

    r1 = vectornorm(Ax_b ,  n) / vectornorm(b,n);
    r2 = vectornorm(x_realx ,  n) / vectornorm(realx,n);

    // delete []Ax;
    // delete []Ax_b;
    // delete []x_realx;
    
}






void report(char *title, double r1, double r2 ,double t1,  double t2 ,int its, int n  )
{
    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
        title, r1, r2, its, its / n, t1, t2);
}


void matmult(double *res,double *a, double *b, int n, int m,int l) //a^{n*m} * b^{m*l} = res^{n*l}
{
    int i = 0,j=0;
    
   
    if(!a || !b)
    {
        return;
    }

    for(i = 0 ; i < n ; i++)
    {
        for(j = 0; j< l; j++)
        {
            double s = 0;
            for(int k = 0 ; k < m ; k++)
            {
                s += a[i*m+k]*b[k*l+j];
            }
            res[i*l+j] = s;
        }
    }

}

 void get_block(double *a, double *b, int n, int m, int i, int j)
{
    int i1=0, j1=0, k, l, r, h;
    if(m == 0) {
        printf("m == 0 in i = %d , j = %d\n",i,j);
        return;
    }
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    if(j < k) h = m;
    else h = l;


    double *bl = a + i*n*m + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < r ; i1++)
    {
        for (j1=0 ; j1 < h ; j1++)
        {
            b[i1*h + j1] = bl[i1*n + j1];
            

        }
    }
    

}



void set_block(double *a, double *b, int n, int m, int i, int j)
{
    int i1=0, j1=0, k, l, r, h;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    if(j < k) h = m;
    else h = l;


    double *bl = a + i*n*m + j*m; //start of block

    for (i1 =0 ; i1 < r ; i1++)
    {
        for (j1=0 ; j1 < h ; j1++)
        {
            bl[i1*n + j1] = b[i1*h + j1];
            

        }
    }
    

}

double normofmatrix(double *a , int size)
{
    if (!a)
    { 
        printf("nullptr in norm of matrix\n");
        return -1;
    }

    double mm = -1;

    for(int j = 0; j< size ; j++)
    {
        double s = 0;
        for(int i = 0 ;i <size; i++)
        {
            s+=fabs(a[i*size+j]);
        }
        if(mm < s) mm = s;
    }

    return mm;

}

double* inverse(double *result,const double* A, int size,double eps)
{
    if(!A)
    {
        printf("nullptr in inverse\n");
        return nullptr;
    }

    // for(int i = 0; i<size*size; i++)
    // {
    //     if(A[i])
    //     {
    //         printf("cant inverse non square matrix\n");
    //         return nullptr;
    //     }
    // }

    double* a = new double[size*size];
    memcpy(a, A, size*size * sizeof(double));


    double* E = new double[size*size];
    int* colsw = new int[size];
    
    // init E
    for(int i = 0; i<size; i++)
    {
        for(int j = 0; j< size; j++)
        {
            if (i!=j) E[i*size+j] = 0;
            else E[i*size +j] = 1;
        }
    }

    //init colsw
    for(int i = 0; i <size ;i++) colsw[i] = i;


    for(int i = 0 ; i< size ; i++)
    {
        double maxEl = fabs(a[i*size+i]);
        int pivot = i;
        for(int j = i+1; j<size;j++)
        {
            if(fabs(a[i*size +j]) > maxEl)
            {
                maxEl = fabs(a[i*size +j]);
                pivot = j;
            }
        }

        
        //swap columns
        if(pivot!=i)
        {
            for(int k = 0; k<size ; k++)
            {
                swap(a[k*size+i],a[k*size+pivot]);
                swap(E[k*size+i],E[k*size+pivot]);
            }
            swap(colsw[i],colsw[pivot]);
        }

        if(fabs(a[i*size+i]) < eps)
        {
            // printf("matrix has no inverse\n");
            delete []E;
            delete []colsw;
            delete []a;
            return nullptr;
        }

        //devide row i 
        // cout<<"A matr before devideing"<<endl;
        // printlxn(a,size,size,size,size);

        double mainEl = a[i*size+i];
        for(int k = 0 ; k<size; k++)
        {
            a[i*size+k] /= mainEl;
            E[i*size+k] /= mainEl;
        }
        // cout<<"A matr after devideing"<<endl;
        // printlxn(a,size,size,size,size);


        for(int k = 0; k< size;k++)
        {
            if(k!=i)
            {
                double factor = a[k*size+i];
                for(int j = 0; j <size;j++)
                {
                    a[k*size+j] -= factor * a[i*size+j];
                    E[k*size+j] -= factor * E[i*size+j];
                }
            }
        }

    }

// cout<<"Last a presentation before swap"<<endl;
// printlxn(a,size,size,size,size);


    // swap back columns
    

    for (int i = 0; i < size*size; ++i) result[i] = 0.0;

    
    
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            result[ colsw[i]*size + colsw[j] ] = E[i*size + j];
// cout<<"Last a presentation after swap"<<endl;

// printlxn(a,size,size,size,size);

    delete []a;
    delete [] colsw;
    delete []E;
    return result;
    
    
}

void blocksize(int i, int j,int n,int m, int &r, int &h)
{
    int k,l;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
        else r = l;

    if(j < k) h = m;
    else h = l;
}


void swap_block_columns(double *a, int n,int m, int i, int j)
{
    for(int p =0 ; p < m ; p++)
                {
                    for(int c = 0; c<n ; c++)
                    {
                        //должны свапнуть столбцы блоков 
                        swap(a[c*n+i*m+p],a[c*n+j*m+p]);
                        
                    }
                }
                
                // printf("swapped blocks %d and %d\n",i,j);
}

void swap_block_vec(double *a, int n,int m, int i, int j)
{
    n=n;
    for(int p =0 ; p < m ; p++)
                {
                    
                        //должны свапнуть столбцы блоков 
                        swap(a[i*m+p],a[j*m+p]);
                        
                    
                }
                
               
}

void get_vec_block(double *b,double *block,int n, int m, int i)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    for (i1 =0 ; i1 < r ; i1++)
    {
            block[i1] = b[i*m+i1];
    }
    

}


void mat_mult_sub(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_B, const int col_row) {
    int row_l = row_A % 3;
    int col_l = col_B % 3;
    int row_k = (row_A - row_l) / 3;
    int col_k = (col_B - col_l) / 3;
    double res_00 = 0, res_01 = 0, res_02 = 0;
    double res_10 = 0, res_11 = 0, res_12 = 0;
    double res_20 = 0, res_21 = 0, res_22 = 0;

    for (int b_i = 0; b_i < row_k; b_i++) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            res_20 = 0, res_21 = 0, res_22 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];

                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_12 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_22 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];
            }

            Result[b_i * 3 * col_B + b_j * 3] -= res_00; Result[b_i * 3 * col_B + b_j * 3 + 1] -= res_01; Result[b_i * 3 * col_B + b_j * 3 + 2] -= res_02;
            Result[(b_i * 3 + 1) * col_B + b_j * 3] -= res_10; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] -= res_11; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
            Result[(b_i * 3 + 2) * col_B + b_j * 3] -= res_20; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] -= res_21; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] -= res_22;
        }

        if (col_l != 0) {
            res_00 = 0, res_01 = 0, res_10 = 0;
            res_11 = 0, res_20 = 0, res_21 = 0;

            for (int s = 0; s < col_row; s++) {
                if(col_l > 1) {
                    res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }

                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3];
            }

            Result[b_i * 3 * col_B + col_k * 3] -= res_00; 
            Result[(b_i * 3 + 1) * col_B + col_k * 3] -= res_10;
            Result[(b_i * 3 + 2) * col_B + col_k * 3] -= res_20;

            if(col_l > 1) {
                Result[b_i * 3 * col_B + col_k * 3 + 1] -= res_01; 
                Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
                Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] -= res_21;
            }
        }
            
    }

    if(row_l != 0) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                    res_11 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                    res_12 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 
                }
            }

            Result[row_k * 3 * col_B + b_j * 3] -= res_00; Result[row_k * 3 * col_B + b_j * 3 + 1] -= res_01; Result[row_k * 3 * col_B + b_j * 3 + 2] -= res_02;

            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + b_j * 3] -= res_10; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] -= res_11; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
            }
        }

        if(col_l != 0) {
            res_00 = 0, res_01 = 0;
            res_10 = 0, res_11 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                if (col_l > 1) {
                    res_01 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                }
                if (col_l > 1 && row_l > 1) {
                    res_11 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

            Result[row_k * 3 * col_B + col_k * 3] -= res_00;

            if (col_l > 1) {
                Result[row_k * 3 * col_B + col_k * 3 + 1] -= res_01;
            }
            
            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3] -= res_10;
            }

            if (row_l > 1 && col_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
            }
        }     
    }
}

void set_vec_block(double *b,double *block,int n, int m, int i)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    


    //start of block 

    for (i1 =0 ; i1 < r ; i1++)
    {
            b[i*m+i1]=block[i1] ;
    }
    

}

void get_block_ml(double *a, double *b, int n, int m,int l, int i)
{
    
    int i1=0, j1=0;
   

    double *bl = a + i*n*m + (n-l); //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < m ; i1++)
    {
        for (j1=0 ; j1 < l ; j1++)
        {
            b[i1*l + j1] = bl[i1*n + j1];
            

        }
    }
    

    
}

void set_block_ml(double *a, double *b, int n, int m,int l, int i)
{
    
    int i1=0, j1=0;
   

    double *bl = a + i*n*m + (n-l); //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < m ; i1++)
    {
        for (j1=0 ; j1 < l ; j1++)
        {
             bl[i1*n + j1]=b[i1*l + j1] ;
            

        }
    }
    

    
}

void get_block_lm(double *a, double *b, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (n-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < l ; i1++)
    {
        for (j1=0 ; j1 < m ; j1++)
        {
            b[i1*m + j1] = bl[i1*n + j1];
            

        }
    }
    
}

void set_block_lm(double *a, double *b, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (n-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < l ; i1++)
    {
        for (j1=0 ; j1 < m ; j1++)
        {
            bl[i1*n + j1]= b[i1*m + j1];
            

        }
    }
    
}

void matsub(double *res,double *a, double *b, int m,int l)
{
    for(int i = 0 ;i<m;i++)
    {
        for(int j = 0; j<l; j++)
        {
            res[i*l+j] = a[i*l+j] - b[i*l+j];
        }
    }
}

void vec_mult_sub(double* Result, double* A, double* vec, int m) {
    double* temp = new double[m];

    for (int i = 0; i < m; i++)
    {
        
        temp[i] = 0;
        double sum = 0;
        for (int j = 0; j < m; j++) {
            
            if(fabs( A[i*m + j] ) < EPS64 ) 
            {
                A[i*m + j] = 0;
            }

                sum += A[i*m + j]*vec[j];

        }
        temp[i] += sum;
            // printf("temp[%d] = %lf\n",i,temp[i]);
        
    }

    for (int i = 0; i < m; i++)
        Result[i] -= temp[i];

    // cout<<"vec_mult:\n";
    // for(int i = 0 ; i < m ; i++)
    // {
    //     cout<<Result[i]<<" ";
    // }
    delete[] temp;
}

void vec_mult_sub_lm(double* Result, double* A, double* vec, int l,int m) {
    double* temp = new double[l];
    for (int i = 0; i < l; i++)
    {
        temp[i] = 0;
        double sum = 0;
        for (int j = 0; j < m; j++) {
            
            if(fabs(A[i*m + j]) < EPS64 ) 
            {
                A[i*m + j] = 0;
            }
                sum += A[i*m + j]*vec[j];

        }
        temp[i] += sum;
            // printf("temp[%d] = %lf\n",i,temp[i]);
        
    }

    for (int i = 0; i < l; i++)
        Result[i] -= temp[i];

    // cout<<"vec_mult:\n";
    // for(int i = 0 ; i < m ; i++)
    // {
    //     cout<<Result[i]<<" ";
    // }
    delete[] temp;
}







void multiplication(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_row,const int col_B)
{
    int row_l = row_A % 3;
    int col_l = col_B % 3;
    int row_k = (row_A - row_l) / 3;
    int col_k = (col_B - col_l) / 3;
    double res_00 = 0, res_01 = 0, res_02 = 0;
    double res_10 = 0, res_11 = 0, res_12 = 0;
    double res_20 = 0, res_21 = 0, res_22 = 0;

    for (int b_i = 0; b_i < row_k; b_i++) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            res_20 = 0, res_21 = 0, res_22 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];

                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_12 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_22 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];
            }

            Result[b_i * 3 * col_B + b_j * 3] = res_00; Result[b_i * 3 * col_B + b_j * 3 + 1] = res_01; Result[b_i * 3 * col_B + b_j * 3 + 2] = res_02;
            Result[(b_i * 3 + 1) * col_B + b_j * 3] = res_10; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] = res_11; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
            Result[(b_i * 3 + 2) * col_B + b_j * 3] = res_20; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] = res_21; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] = res_22;
        }

        if (col_l != 0) {
            res_00 = 0, res_01 = 0, res_10 = 0;
            res_11 = 0, res_20 = 0, res_21 = 0;

            for (int s = 0; s < col_row; s++) {
                if(col_l > 1) {
                    res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }

                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3];
            }

            Result[b_i * 3 * col_B + col_k * 3] = res_00; 
            Result[(b_i * 3 + 1) * col_B + col_k * 3] = res_10;
            Result[(b_i * 3 + 2) * col_B + col_k * 3] = res_20;

            if(col_l > 1) {
                Result[b_i * 3 * col_B + col_k * 3 + 1] = res_01; 
                Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
                Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] = res_21;
            }
        }
            
    }

    if(row_l != 0) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                    res_11 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                    res_12 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 
                }
            }

            Result[row_k * 3 * col_B + b_j * 3] = res_00; Result[row_k * 3 * col_B + b_j * 3 + 1] = res_01; Result[row_k * 3 * col_B + b_j * 3 + 2] = res_02;

            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + b_j * 3] = res_10; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] = res_11; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
            }
        }

        if(col_l != 0) {
            res_00 = 0, res_01 = 0;
            res_10 = 0, res_11 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                if (col_l > 1) {
                    res_01 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                }
                if (col_l > 1 && row_l > 1) {
                    res_11 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

            Result[row_k * 3 * col_B + col_k * 3] = res_00;

            if (col_l > 1) {
                Result[row_k * 3 * col_B + col_k * 3 + 1] = res_01;
            }
            
            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3] = res_10;
            }

            if (row_l > 1 && col_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
            }
        }     
    }      
}



int solution(int n, int m, double *a, double *b, double *x,
    double *block_mm, double *block_ml, double *block_ll,double *invblock_mm, double *diaginvblock_mm, 
    double *invblock_ll,double *diagblock_mm,
    int *colsw,double *vecb_m,double *vecb_l,
    double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *tmpvecb_m,double *tmpvecb_l,const double eps)
{
    
    // m=m;
    // a=a;
    // b=b;
    // block_ml=block_ml;
    // diaginvblock_mm=diaginvblock_mm;
    // diagblock_mm=diagblock_mm;
    // vecb_l=vecb_l;
    // tmpblock_ml1 = tmpblock_ml1;

    

    

    int  k, l/*, r, h*/;
    // int i = 1,
    // j=0;

    k = n/m; l = n - m*k ;
    
    int is_l = (l == 0) ? 0:1; 
    
    
    for(int c =0; c < k  ;c++) colsw[c]=c;

    for(int i = 0 ; i < k + is_l; i++)
    {   
        double minNorm = 1e64;

        int mainBlock = i;

        if(i != k)
        {
            for(int j = i ; j< k ; j++)
            {

            get_block(a,block_mm,n,m,i,j);

    //             printf("Block[%d,%d]\n",i,j);
    //             printlxn(block_mm,m,m,m,m);

            
            // printlxn(invblock_mm,m,m,m,m);

            if(inverse(invblock_mm,block_mm,m,eps))
                {

                        
//                     cout<<"inverse "<<i<<" "<<j<<" with norm = "<<normofmatrix(invblock_mm,m)<<endl;
// 
//                 printlxn(invblock_mm,m,m,m,m);

                if(normofmatrix(invblock_mm,m) < minNorm) 
                    {
                        
                        minNorm = normofmatrix(invblock_mm,m);
                        mainBlock = j;
                       
                    }
                }

            }

        }else{
            get_block(a,block_ll,n,m,k,k);

            if(!inverse(invblock_ll,block_ll,l,eps))
            {
                printf("Block [%d,%d] (block[l,l] in our matrix)  has no inverse after the transformations\n",k,k);
                return -1;
            }
            
            
//             printlxn(invblock_ll,l,l,l,l);
            minNorm = normofmatrix(invblock_ll,l);
            
        }

        if((fabs(minNorm - 1e64) < eps))
        {   
            if(i!=0)
                printf("No inverse matrix in row %d after the transformations\n",i);
            else
                printf("No inverse matrix in row %d\n",i);
            return -1;
        }

        if(mainBlock != i)
            {
                swap_block_columns(a,n,m,i,mainBlock);
                // printlxn(a,n,n,n,n);
                swap(colsw[i],colsw[mainBlock]);
                // cout<<"swapped "<< i<<" "<<mainBlock<<" in row "<<i<<endl;
            }
        
        
        // printlxn(a,n,n,n,n);
        // cout<<"TEST1"<<endl;
        if(i<k)
        {
            get_block(a,diagblock_mm,n,m,i,i);
            
            if(!(inverse(diaginvblock_mm,diagblock_mm,m,eps)))
            {
                        printf("no blocks in row has inverse block\n");
                        
                        return -1;
                    }

            get_vec_block(b,vecb_m,n,m,i);
            mat_x_vector(tmpvecb_m,diaginvblock_mm,vecb_m,m);// double *resvec = mat_x_vector(diaginvblock_mm,vecb_m,m);
            // cout<<"tmpvecb_m : "<<endl;
            // printlxn(tmpvecb_m,m,1,m,m);    
            set_vec_block(b,tmpvecb_m,n,m,i);

            for(int j = i ; j < k ; j++) //mb try j = i
            {
                get_block(a,block_mm,n,m,i,j);
                
               multiplication(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// matmult(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// double *resmult = matmult(diaginvblock_mm,block_mm,m,m,m)

                set_block(a,tmpblock_mm,n,m,i,j);
                
                            if (!block_mm || !vecb_m || !invblock_mm) {
                fprintf(stderr, "Error: temporary buffers not initialized!\n");
                return -1;
            }
            }
            // printlxn(a,n,n,n,n);
            if(is_l != 0)
            {
                get_block_ml(a,block_ml,n,m,l,i);
                multiplication(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);// matmult(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);
                set_block_ml(a,tmpblock_ml,n,m,l,i);
            }
            
        }else
            {   
                // printlxn(a,n,n,n,n);
                // printlxn(b,n,1,n,n);
                get_block(a,block_ll,n,m,i,i);
                get_vec_block(b,vecb_l,n,m,i);

                // printlxn(block_ll,l,l,l,n);

                // cout<<"vecb_l:"<<endl;
                // printlxn(vecb_l,l,1,l,n);

                if(!(inverse(invblock_ll,block_ll,l,eps)))
                    {
                        printf("ll block has no inverse\n");
                        return -1;
                    }

                
                // printlxn(invblock_ll,l,l,l,n);

                multiplication(tmpblock_ll,invblock_ll,block_ll,l,l,l);// matmult(tmpblock_ll,invblock_ll,block_ll,l,l,l);

                mat_x_vector(tmpvecb_l,invblock_ll,vecb_l,l);
                
                // printlxn(tmpvecb_l,l,1,l,n);

                set_block(a,tmpblock_ll,n,m,i,i);
                set_vec_block(b,tmpvecb_l,n,m,i);
                // cout<<"WE ARE IN i = k"<<endl;
            }

            // printlxn(a,n,n,n,n);
            // printlxn(b,n,1,n,n);
            //начинаем обнулять столбцы
            
        for(int r = i+1 ; r < k + is_l ; r++)
        {
            // cout<<"TEST "<<r<<endl;
            if(r < k)
            {
                get_block(a,block_mm,n,m,r,i);
                get_block(a,tmpblock_mm,n,m,r,i);
                memset(tmpblock_mm,0, m*m*sizeof(double));
                set_block(a,tmpblock_mm,n,m,r,i);

                // not in i for
                get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
                get_vec_block(b,tmpvecb_m,n,m,r);
                vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
                set_vec_block(b,tmpvecb_m,n,m,r);

                // cout<<"tmpvecb_m in subtract i= "<<i<<" r="<<r<<endl;
                // printlxn(tmpvecb_m,m,1,m,m);
                // printlxn(b,n,1,n,n);

                for (int j = i + 1; j < k; j++) {
                    get_block(a,invblock_mm,n,m,i,j);
                    get_block(a,diagblock_mm,n,m,r,j);
                    mat_mult_sub(diagblock_mm,block_mm,invblock_mm,m,m,m);
                    set_block(a,diagblock_mm,n,m,r,j);
                }

                if (is_l!= 0) {
                get_block_ml(a,tmpblock_ml,n,m,l,i);
                get_block_ml(a,tmpblock_ml1,n,m,l,r);
                mat_mult_sub(tmpblock_ml1,block_mm,tmpblock_ml,m,l,m);
                set_block_ml(a,tmpblock_ml1,n,m,l,r);
                }
            }else
            {
                // printlxn(a,n,n,n,n);
               get_block_lm(a, block_ml, n, m, l, i);

            //    printf("block_lm in col %d\n",i);
            //    printlxn(block_ml,m,l,m,n);

               get_block_lm(a, tmpblock_ml, n, m, l, i);
               memset(tmpblock_ml,0,m*l*sizeof(double));
               set_block_lm(a, tmpblock_ml, n, m, l, i);

               get_vec_block(b,vecb_m,n,m,i);// get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
            //    cout<<"vecb_m in subtract i= "<<i<<" r="<<r<<endl;
            //     printlxn(vecb_m,m,1,m,m);
               get_vec_block(b,tmpvecb_l,n,m,r);// get_vec_block(b,tmpvecb_m,n,m,r);
               vec_mult_sub_lm(tmpvecb_l,block_ml,vecb_m,l,m);// vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
               set_vec_block(b,tmpvecb_l,n,m,r);  // set_vec_block(b,tmpvecb_m,n,m,r);


                for(int j = i + 1; j < k; j++) {
                get_block(a,tmpblock_mm,n,m,i,j);
                get_block_lm(a, tmpblock_ml, n, m, l, j);

                // cout<<"tmpblock_ml in col "<<j<<endl;
                // printlxn(tmpblock_ml,m,l,m,m);

                mat_mult_sub(tmpblock_ml,block_ml,tmpblock_mm,l,m,m);

                
                // get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
                // get_vec_block(b,tmpvecb_m,n,m,r);//
                vec_mult_sub_lm(tmpvecb_m,block_ml,vecb_m,l,m);//
                set_block_lm(a, tmpblock_ml, n, m, l, j);
                // set_vec_block(b,tmpvecb_m,n,m,r);//set_vec...
                }

                if (is_l != 0) {
                    get_block_ml(a,tmpblock_ml,n,m,l,i);
                    get_block(a,tmpblock_ll,n,m,k,k);
                    mat_mult_sub(tmpblock_ll,block_ml,tmpblock_ml,l,l,m);
                    set_block(a,tmpblock_ll,n,m,k,k);
                }

            }
            
        }

        // cout<<"LAST PRINT"<<endl;
        // printlxn(a,n,n,n,n);
        // printlxn(b,n,1,n,n);
         
        
    }
       
    
    //начало обратного хода
        

    // cout<<"colsw : "<<endl;
    // for(int i = 0 ; i < k ; i++) cout<<colsw[i]<<" ";


    for(int i = n-1; i >= 0 ; i--)
    {
        if(i == n-1) x[i] = b[i];
        
        else
        {
            x[i] = b[i];
            for(int j = n-1 ; j >i;j--)
            {
                x[i] -= a[i*n + j]*x[j];
            }
        }
    }

    // cout<<"Vector x before swap :"<<endl;
    // printlxn(x,n,1,n,n);

    

    for(int i = 0 ; i < k ; i++)
    {
        

        if(i != colsw[i]){ 
            int t;
            swap_block_vec(x,n,m,i,colsw[i]);
            t = colsw[colsw[i]];
            colsw[colsw[i]] = colsw[i];
            colsw[i] = t; 
        }


    }

    

    return 0;

}










void threediag(double *a,int n,double eps)
{
    if(!a)
    {
        printf("error in threediag\n");
        return;
    }


    for(int j = 0 ; j < n ; j++)
    {
        for(int i = n - 1; i >= j + 2; i--)
        {
            double x = a[(i-1) * n + j];
            double y = a[(i) * n + j];

            // printf("fabs(y) = %lf eps = %lf\n",fabs(y),eps);

            if(fabs(y) < eps) continue;

            double r = sqrt(x*x+y*y);
            double x1 = x/r;
            double x2 = -y/r;

            a[(i-1) * n + j] = x1 * x - y * x2;
            a[(i) * n + j] = x2 * x + y * x1;

            //mult all others columns T*A
            for(int s = j + 1; s < n; s++)
            {
                double s1 = a[(i-1) * n + s];
                double s2 = a[(i) * n + s];
                a[(i-1) * n + s] = x1 * s1 - s2 * x2;
                a[(i) * n + s] = x2 * s1 + s2 * x1;
            }

             // A*T^{t}
            for(int s = 0; s < n; s++)
            {
                double s1 = a[s * n + (i-1)];
                double s2 = a[s * n + i];
                a[s * n + (i-1)] = x1 * s1 - s2 * x2;
                a[s * n + i] = x2 * s1 + s2 * x1;
            }

            a[i * n + j] = 0;
            a[j * n + i] = 0;
        }
    }

}

double trace(double *a, int n)
{
    if(!a)
    {
        printf("error in trace\n");
        return -1;
    }

    double s = 0;

    for(int i = 0 ; i < n ; i++)
    {
        s += a[i * n + i];
    }

    return s;
}


void QR_refl(double *a,int n, double eps,int &itr)
{
    
    int it,j,k;
    double *x0,*x1;
    x0 = new double[n]();
    x1 = new double[n]();

    for(int r = 0 ; r < n - 1; r++)
    {
        int index = n - r;

        if(fabs(a[(index - 1)*n + index - 2]) < eps)
            continue;
        
        double sk = a[(index - 1) * n + (index - 1)] + 0.5 * a[(index - 1) * n + (index - 2)];

        
        // printf("sk = %lf\n",sk);

        for( it = 0;  it < ITER; it++)
        {
            for(int i = 0; i < index; i++)
                a[i*n + i] -= sk;
            
            for(int k = 0 ; k < index -1; k++)
            {
                double c = a[(k+1)*n + k];
                double b = a[k*n + k];
                double nrm = sqrt(c*c + b*b);
                nrm = (nrm < 0 ? -nrm:nrm);
                x0[k] = a[k*n+k] - nrm;
                x1[k] = a[(k + 1)*n + k];
                double sqr = sqrt(x0[k] * x0[k] + x1[k] * x1[k]);

                if(sqr < eps)
                    continue;
                
                x0[k] = x0[k]/sqr;
                x1[k] = x1[k]/sqr;

                a[k * n + k] = nrm;
                a[(k + 1) * n + k] = 0;
                int cond = (k + 3 < index ? k + 3 : index);
                for (j = k + 1; j < cond; j++)
                {
                    double qq = 2 * (a[k * n + j] * x0[k] + a[(k + 1) * n + j] * x1[k]);
                    if (fabs(qq) < 1e-64)
                        continue;
                    
                    a[k * n + j] = a[k * n + j] - qq * x0[k];
                    a[(k + 1) * n + j] = a[(k + 1) * n + j] - qq * x1[k];
                }
            }

            for (k = 0; k < index - 1; k++)
            {
                int cond = (k + 3 < index ? k + 3 : index);
                for (j = k; j < cond; j++)
                {
                    double qq = 2 * (a[j * n + k] * x0[k] + a[j * n + k + 1] * x1[k]);
                    if (fabs(qq) < 1e-64)
                        continue;
                    
                    a[j * n + k] = a[j * n + k] - qq * x0[k];
                    a[j * n + (k + 1)] = a[j * n + (k + 1)] - qq * x1[k];
                }
            }

            for (int i = 0; i < index; i++)
            {
                for (j = i + 1; j < (i + 3 < index ? i + 3 : index); j++)
                {
                    a[i * n + j] = a[j * n + i];
                }
            }

            for (int i = 0; i < index; i++)
            {
                a[i * n + i] += sk;
            }

            if (fabs(a[(index - 1) * n + index - 2]) < eps)
            {
                break;
            }
            sk = a[(index - 1) * n + (index - 1)] + 0.5 * a[(index - 1) * n + (index - 2)];
            // printlxn(a,n,n,n);
        }

        itr += it;

    }

    delete [] x0;
    delete [] x1;

}



<<<<<<< HEAD
void threediag(double *a,int n,double eps)
{
    if(!a)
    {
        printf("error in threediag\n");
        return;
    }


    for(int j = 0 ; j < n ; j++)
    {
        for(int i = n - 1; i >= j + 2; i--)
        {
            double x = a[(i-1) * n + j];
            double y = a[(i) * n + j];

            // printf("fabs(y) = %lf eps = %lf\n",fabs(y),eps);

            if(fabs(y) < eps) continue;

            double r = sqrt(x*x+y*y);
            double x1 = x/r;
            double x2 = -y/r;

            a[(i-1) * n + j] = x1 * x - y * x2;
            a[(i) * n + j] = x2 * x + y * x1;

            //mult all others columns T*A
            for(int s = j + 1; s < n; s++)
            {
                double s1 = a[(i-1) * n + s];
                double s2 = a[(i) * n + s];
                a[(i-1) * n + s] = x1 * s1 - s2 * x2;
                a[(i) * n + s] = x2 * s1 + s2 * x1;
            }

             // A*T^{t}
            for(int s = 0; s < n; s++)
            {
                double s1 = a[s * n + (i-1)];
                double s2 = a[s * n + i];
                a[s * n + (i-1)] = x1 * s1 - s2 * x2;
                a[s * n + i] = x2 * s1 + s2 * x1;
            }

            a[i * n + j] = 0;
            a[j * n + i] = 0;
        }
    }
=======
>>>>>>> 3591345 ('Mon Dec 22 00:44:48 ')




<<<<<<< HEAD
}

double trace(double *a, int n)
{
    if(!a)
    {
        printf("error in trace\n");
        return -1;
    }

    double s = 0;

    for(int i = 0 ; i < n ; i++)
    {
        s += a[i * n + i];
    }

    return s;
}
=======












>>>>>>> 3591345 ('Mon Dec 22 00:44:48 ')
