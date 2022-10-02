#include <stdio.h>
#include <math.h>
#define N 2

void test_matrix(double array[N][N]);
double calculate_det(double array[N][N]);

int main()
{
  int i,j,k,count;
  double* p[N];
  double A[N][N]= {{ 6.0, 3.0},  
                    { 3.0, 4.0}};
  double I[N][N]= {{ 1.0, 0.0},  
                    { 0.0, 1.0}};
  double b[N]={-3.0,-9.0}, v_1[N]={0.0,0.0},v_2[N];
  double u[N],C[N],errors[N],residuals[N];
  double A_inversed[N][N],B[N][N];
  double det,temp,alpha=0.2;
  double error,residual;
  int num_steps=100;

  printf("original matrix is: \n");
  test_matrix(A);
  det= calculate_det(A);
  printf("det = %f\n",det);
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      A_inversed[i][j]=pow(-1,i+j)*A[N-i-1][N-1-j]/det;
    }
  }
  printf("\ninversed matrix is: \n");
  test_matrix(A_inversed);
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      u[i]=u[i]+A_inversed[i][j]*b[j];
    }
  }
  printf("exact solution u=(%f,%f)\n",u[0],u[1]);
  
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      B[i][j]=I[i][j]-alpha*A[i][j];
    }
  }
  test_matrix(B);

  for(i=0;i<N;i++)
  {
    C[i]=b[i]*alpha;
  }
  
  // Richardson iteration
  for(k=0;k<num_steps;k++)
  {
    for(i=0;i<N;i++)
    {
      temp=0.0;
      for(j=0;j<N;j++)
      {
        temp=temp+B[i][j]*v_1[j];
      }
      v_2[i]=temp+C[i];
    }

    // errors
    error=0.0;
    for(i=0;i<N;i++)
    {
      errors[i]=v_2[i]-u[i];
      error=error+pow(errors[i],2);
    }
    error=pow(error,0.5);

    // residuals
    residual=0.0;
    for(i=0;i<N;i++)
    {
      temp=0.0;
      for(j=0;j<N;j++)
      {
        temp=temp+A[i][j]*v_2[j];
        residuals[i]=b[i]-temp;
        residual=residual+pow(residuals[i],2);
      }
    } 
    residual=pow(residual,0.5);

    // update
    for(i=0;i<N;i++)
    {
      v_1[i]=v_2[i];
    }
    printf("numerical solution v=(%f,%f) of No. %d iteration, error=%f, residual=%f \n",v_2[0],v_2[1],k+1,error,residual);
  }
  return 0;
}

// output a matrix's elements
void test_matrix(double array[N][N]) 
{
  int i,j;
  for(i=0;i<2;i++)
  {
    for(j=0;j<2;j++)
    {
          printf("%lf\t",array[i][j]);
    }
    printf("\n");
  }
}

// calculate the value of a determint
double calculate_det(double array[N][N])
{
  return array[0][0]*array[1][1]-array[0][1]*array[1][0];;
}