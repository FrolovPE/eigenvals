#include <iostream>
#include <chrono>
#include <unistd.h>
#include <cstring>
#include <vector>
#include <pthread.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <cmath>
#include "mytime.h"


#define EPS64 1e-64
#define CLEAR     clear(block_mm,block_ml,block_ll,tmpblock_mm,tmpblock_ml,tmpblock_ml1,tmpblock_ll,invblock_mm,invblock_ll,diagblock_mm,diaginvblock_mm,vecb_m,vecb_l,tmpvecb_m, tmpvecb_l,colsw)
using namespace std;


void printlxn(const double *a, int l, int n, int r);
int readarray(double *a, int n, char* filename);
double f (int k , int n , int i , int j);
void init(double *a,  int n, int s);
int solution(int n, int m, double *a, double *b, double *x,
    double *block_mm, double *block_ml, double *block_ll,double *invblock_mm, double *diaginvblock_mm, 
    double *invblock_ll,double *diagblock_mm,
    int *colsw,double *vecb_m,double *vecb_l,
    double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *tmpvecb_m,double *tmpvecb_l,const double eps);
double vectornorm(double *a , int n);
void mat_x_vector(double *res,double *a, double *b, int n);
void vectorsub(double *res,double *a,double *b, int n);
void residuals(double &r1,double &r2,double *a,double *b,double *x,double *realx,int n,double *Ax, double *Ax_b, double *x_realx);
void report(char *title, double r1, double r2 ,double t1,  double t2 ,int its, int n  );
void matmult(double *res,double *a, double *b, int n, int m,int l) ;
void get_block(double *a, double *b, int n, int m, int i, int j);
void set_block(double *a, double *b, int n, int m, int i, int j);
double normofmatrix(double *a , int size);
double* inverse(double *result,const double* A, int size,double eps);
double* inverse1(const double* a_in, int n,const double eps );
void blocksize(int i, int j,int n,int m, int &r, int &h);
void swap_block_columns(double *a, int n,int m, int i, int j);
void swap_block_vec(double *a, int n,int m, int i, int j);
void get_vec_block(double *b,double *block,int n, int m, int i);
void mat_mult_sub(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_B, const int col_row);
void set_vec_block(double *b,double *block,int n, int m, int i);
void get_block_ml(double *a, double *b, int n, int m,int l, int i);
void set_block_ml(double *a, double *b, int n, int m,int l, int i);
void get_block_lm(double *a, double *b, int n, int m,int l, int j);
void set_block_lm(double *a, double *b, int n, int m,int l, int j);
void matsub(double *res,double *a, double *b, int m,int l);
void vec_eq(double *x, double *b, int m);
void vec_sub(double *x, double *b, int m);
void vec_mult_sub(double* Result, double* A, double* vec, int m);
void vec_mult_sub_lm(double* Result, double* A, double* vec, int l,int m);
void multiplication(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_row,const int col_B);
void* parallelSolve(void* ptr);
void* parallelSolve1(void* ptr);
void pllinit_matrix(double *a, int s,int n , int m , int k, int p);
void pllinit_vectorb(double *b,double *a,int n , int m , int k, int p);
void clear(double *block_mm,double *block_ml,double *block_ll,double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *invblock_mm,double *invblock_ll,double *diagblock_mm,double *diaginvblock_mm,double *vecb_m,double *vecb_l,double *tmpvecb_m, double *tmpvecb_l,int *colsw);
void threediag(double *a,int n,double eps);
double trace(double *a, int n);