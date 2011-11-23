#include "matrix.h"

Mat::Mat(int p1, int p2)
{
  int i;
  p_r = p1;
  p_c = p2;
  Data = new double[p_r * p_c];
  for(i = 0; i < p_r * p_c; i++)Data[i] = 0.0;
}

Mat::Mat(int p1, int p2, double* X)
{
  int i,j;
  p_r = p1;
  p_c = p2;
  Data = new double[p_r * p_c];
  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < p_c; j++)
	{
	  Data[i + j * p_r] = X[i + j * p_r];
	}
    }
}

Mat::~Mat()
{
  delete[] Data;
}

double Mat::get(int i,int j)
{
  return(Data[i + j * p_r]);
}

void Mat::set(int i, int j, double s)
{
  Data[i + j * p_r] = s;
  return;
}

void Mat::set_iden()
{
  int i;
  for(i = 0; i < p_r * p_r; i++) Data[i] = 0;
  for(i = 0; i < p_r; i++) Data[i + i * p_r] = 1.0;
  return;
}

Matrix Mat::copy()
{
  int i,j;
  Matrix A = new Mat(p_r,p_c);
  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < p_c; j++)
	{
	  A->set(i,j,get(i,j));
	}
    }
  return(A);
}
//  Inverts a matrix
Matrix Mat::invert()
{
  Matrix A_copy = copy();
  Matrix A_inv = new Mat(p_r,p_c);
  char uplo = 'U';
  int info;
  int p = p_r;
  A_inv->set_iden();
  
  //-----  Use LAPACK  -------
  F77_NAME(dposv)(&uplo, &p, &p, A_copy->Data, &p, A_inv->Data, &p, &info);
  //--------------------------
  delete A_copy;

  return(A_inv);
}

void Mat::print()
{
	/*
  int i,j;
  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < p_c; j++)
	{
	  printf("%8.5f, ",get(i,j));
	}
      printf( "\n" );
    }
*/
}

//  Utility that prints out a matrix full of doubles
void util_print_mat_dbl(double *X,int p_i,int p_j)
{
/*
  int i,j;
  for(i = 0; i < p_i; i++)
    {
      for(j = 0; j < p_j; j++)
	{
	  printf("%8.8f, ",X[i + j * p_i]);
	}
      printf( "\n" );
    }
*/
}//util_print_mat_dbl

//  Utility that prints out a matrix full of ints
void util_print_mat_int(int *X,int p_i,int p_j)
{
/*
  int i,j;

  for(i = 0; i < p_i; i++)
    {
      for(j = 0; j < p_j; j++)
	{
	  printf("%d, ",X[i + j * p_j]);
	}
      printf("\n");
    }
*/
}//util_print_mat_dbl

// Utility that prints out a vector full of ints
void util_print_vec_int(int *X,int n)
{
/*
  int i;

  for(i = 0; i < n; i++)
    {
      printf("%d, ",X[i]);
    }
  printf("\n");
*/
}//util_print_mat_dbl

// Utility that prints out a vector full of doubles
void util_print_vec_dbl(double *X,int n)
{
/*
  int i;

  for(i = 0; i < n; i++)
    {
      printf("%4.3f, ",X[i]);
    }
  printf("\n");
*/
}//util_print_mat_dbl



//----------- Calculates the Log Determinant ------------
// WARNING: THIS ROUTINE IS SPECIFICALLY DESIGNED FOR
//          SYMMETRIC PD MATRICES
//
//  Here I use the fact that for any symmetric
//  PD Matrix D, then:
//    
//    |D| = 2 |T|
//
//  where T is the cholesky decomposition of D
//  note that |T| = \prod_{i = 1}^p T_{ii}
//  which makes this quite easy.
double Mat::log_det()
{
  int i;
  char uplo = 'U';
  int info;
  Matrix A = copy();
  int p = p_r;
 //----Use Lapack-----
  F77_NAME(dpotrf)(&uplo,&p,A->Data,&p,&info);
  //-------------------

  double result = 0;
  for(i = 0; i < p; i++) result += 2 * log(A->get(i,i));
  delete A;
  return(result);
}

//  Takes a matrix and returns its cholesky
//  Note: the matrix you pass this function is overwritten with the result.
Matrix Mat::chol()
{
  int i,j;
  char uplo = 'U';
  int info;
  Matrix T = copy();
  int p = p_r;
  //----Use Lapack-----
  F77_NAME(dpotrf)(&uplo,&p,T->Data,&p,&info);
  //-------------------

  //---- The lapack function above is messy, it doesn't zero out the lower diagonal -----
  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < i; j++)
	{
	  T->set(i,j,0.0);
	}
    }
  //-------------------------------------------------------------------------------
  return(T);
}

//  Multiplies a p_i x p_k matrix by a p_k x p_j matrix to give a p_i x p_j matrix
Matrix mult_mats(Matrix A, Matrix B)
{
  char trans = 'N';
  double a=1.0;
  double c=0.0;
  int p_i = A->p_r;
  int p_k = A->p_c;
  int p_j = B->p_c;
  Matrix C = new Mat(p_i, p_j);
  F77_NAME(dgemm)(&trans, &trans, &p_i, &p_j, &p_k, &a, A->Data, &p_i, B->Data, &p_k, &c, C->Data, &p_i);
  //---------------------
  return(C);
}

//  Returns the transpose of a matrix
Matrix Mat::transpose()
{
  int i, j;
  Matrix At = new Mat(p_c, p_r);

  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < p_c; j++)
	{
	  At->set(j,i,get(i,j));
	}
    }
  return(At);
}

// Takes a square matrix A (of size p x p) 
// And retrieves the square submatrix B of size p_sub x p_sub, dictated by the vector sub
Matrix Mat::sub_mat(int p_sub, int *sub)
{
  int i,j;
  Matrix B = new Mat(p_sub, p_sub);

  for(i = 0; i < p_sub; i++)
    {
      for(j = 0; j < p_sub; j++)
	{
	  B->set(i, j, get(sub[i],sub[j]));
	}
    }
  return(B);
}

void Mat::set_sub_mat(int p_clique, int* clique, Matrix submat)
{
  int i,j;
  for(i = 0; i < p_clique; i++)
    {
      for(j = 0; j < p_clique; j++)
	{
	  set(clique[i],clique[j],submat->get(i,j));
	}
    }


}
//  Returns the complement of a set clique_ID and stores it in V_ID
void get_complementary_set(int p, int p_clique, int *clique_ID, int *V_ID)
{

  int i,j,k,temp_in_clique;

  k = 0;

  for(i = 0; i < p; i++)
    {
      temp_in_clique = 0;
      for(j = 0; j < p_clique; j++)
	{
	  if(i == clique_ID[j])
	    {
	      temp_in_clique = 1;
	    }
	}
      if(temp_in_clique == 0) 
	{
	  V_ID[k] = i;
	  k++;
	}
    }
  return;
}

//  Standard gaussian stuff.  Takes a Cov matrix Full and returns the conditional
//  Cov matrix for a subset clique_ID.
Matrix get_cond_matrix(int p_cl,int *clique_ID,int *V_ID, Matrix Full)
{
  Matrix Result;

  int i, j;
  int p_V = Full->p_r - p_cl;
  Matrix cV = new Mat(p_cl, p_V);
  Matrix Vc;
  Matrix VV = new Mat(p_V,p_V);
  Matrix VVinv;
  Matrix cV_VVinv;

  //----  Start: Form Matrix Objects --------
  //-----  Make the two cross parts  ------
  //  printf("Here Full\n");
  //  Full->print();
  for(i = 0; i < p_cl; i++)
    {
      for(j = 0;j < p_V; j++)
	{
	  cV->set(i,j,Full->get(clique_ID[i],V_ID[j]));
	}
    }
  //  printf("cV\n");
  //  cV->print();
  Vc = cV->transpose();
  //----------------------------------------------
  
  //---Now the inverse of the remaining elements----
  for(i = 0; i < p_V; i++)
    {
      for(j = 0; j < p_V; j++)
	{
	  VV->set(i,j,Full->get(V_ID[i],V_ID[j]));
	}
    }
  VVinv = VV->invert();
  //  printf("VVinv\n");
  //  VVinv->print();
  //-------------------------------------------------

  //------  Form the conditioning matrix --------
  cV_VVinv = mult_mats(cV,VVinv);
  //  printf("cV_VVinv\n");
  //  cV_VVinv->print();
  Result = mult_mats(cV_VVinv,Vc);
  //  printf("Result\n");
  //  Result->print();
  //---------------------------------------------

  delete cV;
  delete Vc;
  delete VV;
  delete VVinv;
  delete cV_VVinv;

  return(Result);

}//get_cond_matrix

//Takes a matrix A and returns A'A
Matrix square_mat(Matrix A)
{
  Matrix B = new Mat(A->p_c, A->p_c);
  char trans = 'N';
  char transY = 'T';

  double a=1.0;
  double c=0.0;
  int p_i = A->p_r;
  int p_k = A->p_c;
  int p_j = A->p_c;

  F77_NAME(dgemm)(&transY, &trans, &p_i, &p_j, &p_k, &a, A->Data, &p_i, A->Data, &p_k, &c, B->Data, &p_i);
  return(B);
}

void Mat::Add(Matrix A)
{
  int i,j;
  for(i= 0; i < p_r; i++)
    {
      for(j = 0; j < p_c; j++)
	{
	  set(i,j,get(i,j) + A->get(i,j));
	}
    }
  return;
}

void Mat::mult_const(double a)
{
  int i,j;
  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < p_c;j++)
	{
	  set(i,j, get(i,j) * a);
	}
    }
}

double quad_form(Matrix A, Matrix B)
{
  double result;
  Matrix At = A->transpose();
  Matrix AB = mult_mats(A,B);
  Matrix ABAt = mult_mats(AB,At);
  result = ABAt->get(0,0);
  delete At;
  delete AB;
  delete ABAt;

  return(result);
}


intMat::intMat(int p1, int p2)
{

  int i;
  p_r = p1;
  p_c = p2;
  Data = new int[p_r * p_c];
  for(i = 0; i < p_r * p_c; i++)Data[i] = 0;

}

intMat::intMat(int p1, int p2, int* X)
{
  int i,j;
  p_r = p1;
  p_c = p2;
  Data = new int[p_r * p_c];

  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < p_c; j++)
	{
	  Data[i + j * p_r] = X[i + j * p_r];
	}
    }
}

intMat::~intMat()
{
  delete[] Data;
}

int intMat::get(int i,int j)
{
  return(Data[i + j * p_r]);
}

void intMat::set(int i, int j, int s)
{
  Data[i + j * p_r] = s;
  return;
}

void intMat::print()
{
/*
  int i,j;
  for(i = 0; i < p_r; i++)
    {
      for(j = 0; j < p_c; j++)
	{
	  printf("%d, ",get(i,j));
	}
      printf( "\n" );
    }
*/
}
