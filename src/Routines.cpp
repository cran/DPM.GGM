#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

//----------- Code that Interfaces with R to run the iHMM ------------
// This takes all the information from R, creates an iHMMState and iHMMResults object
// And then updates them.
//---------------------------------------------------------------------

// Alex Lenkoski  lenkoski@uni-heidelberg.de


//----- Global information about the Stirling numbers -----------------
double *GLOBAL_STIRLING;
int GLOBAL_STIRLING_DIM;
//--------------------------------------------------------------------

#include "matrix.Rcpp.h"
#include "random.Rcpp.h"
#include "newgraph.Rcpp.h"
#include "gwish.Rcpp.h"
#include "GGM.Rcpp.h"
#include "iHMMState.Rcpp.h"
#include "iHMMResults.Rcpp.h"
#include "DPMState.Rcpp.h"
#include "DPMSliceState.Rcpp.h"

extern "C"{
  void run_iHMM_routine(double *X, int *n, int *p, int *burn, int *reps, double *Stirling_numbers, int *Stirling_dim, 
			int *cluster_mat, int* edge_mat, double* alpha, double* alpha0, int* cluster_total, 
			int *predict, double *stepahead_x, double *stepahead_mu, double *stepahead_K, int* print_every, int *all_k, double *Kall)
{
  int i,k;

  GetRNGstate();


  GLOBAL_STIRLING = Stirling_numbers;
  GLOBAL_STIRLING_DIM = *Stirling_dim;

  State state = new iHMMState(X,5,*n,*p);
  Results results = new iHMMResults(cluster_mat, edge_mat, alpha, 
				    alpha0, cluster_total, predict, 
				    stepahead_x, stepahead_mu, stepahead_K,
				    all_k, reps, Kall);

  int ee = *p * (*p - 1) / 2;

  for(i = 0; i < *n * *n; i++) cluster_mat[i] = 0;
  for(i = 0; i < *n * ee; i++) edge_mat[i] = 0;

  for(i = 0; i < *reps; i++) alpha[i] = 0.0;
  for(i = 0; i < *reps; i++) alpha0[i] = 0.0;
  for(i = 0; i < *reps; i++) cluster_total[i] = 0;

  for(i = 0; i < *n * *p * *p; i++) Kall[i] = 0;

  //---------------- Burn In ------------------
  for(k = 0; k < *burn; k++)
    {
      if((k % *print_every) == 0)Rprintf("***************** Burn Step %d *************\n", k);
      state->UpdateState();
      R_CheckUserInterrupt();
    }
  //----------------------------------------------

  //----------------- Repititions ----------------
  for(k = 0; k < *reps; k++)
    {
      if((k % *print_every) == 0)Rprintf("***************** Rep Step %d *************\n", k);
      state->UpdateState();
      results->UpdateResults(k, state);
      R_CheckUserInterrupt();
    }
  //---------------------------------------------------------------------

  delete state;
  delete results;

  PutRNGstate();

}//run_iHMM_routine


void run_dmp_routine(double *X, int *n, int *p, int *burn, int *reps, int *cluster_mat, int* edge_mat, int* print_every)
{
  int i,j,k,q,r,t;
  GetRNGstate();
  int state_i, state_j;
  DPMState state = new DPMState_CLASS(X,5,*n,*p);
  int ee = *p * (*p - 1) / 2;

  for(i = 0; i < *n * *n;i++) cluster_mat[i] = 0;
  for(i = 0; i < *n * ee; i++) edge_mat[i] = 0;

  //---------------- Burn In ------------------
  for(k = 0; k < *burn; k++)
    {
      if((k % *print_every) == 0)Rprintf("***************** Burn Step %d *************\n", k);
      state->UpdateState();
      R_CheckUserInterrupt();
    }
  //----------------------------------------------

  //----------------- Repititions ----------------
  for(k = 0; k < *reps; k++)
    {
      if((k % *print_every) == 0)Rprintf("***************** Rep Step %d *************\n", k);
      state->UpdateState();
      R_CheckUserInterrupt();
      //------------------- Update Results ------------------------------------
      for(i = 0; i < *n; i++)
	{
	  state_i = state->xi[i];
	  for(j = 0; j < *n; j++)
	    {
	      state_j = state->xi[j];
	      cluster_mat[i * *n + j] += (state_i == state_j);
	    }
	}
      for(i = 0; i < *n; i++)
	{
	  t = 0;
	  for(q = 0; q < *p - 1; q++)
	    {
	      for(r = q + 1; r < *p; r++)
		{
		  edge_mat[i * ee + t] += state->GGMlist[ state->xi[i] ]->graph->Edge[q][r];
		  t++;
		}
	    }
	}
      //----------------------------------------------------------------------------
    }
  //---------------------------------------------------------------------
  
  delete state;
  //
  PutRNGstate();
}//run_dmp_routine;

void run_dmp_slice_routine(double *X, int *n, int *p, int *burn, int *reps, int *cluster_mat, int* edge_mat, int* print_every)
{

  int i,j,k,q,r,t;

  GetRNGstate();
  int state_i, state_j;

  DPMSliceState state = new DPMSliceState_CLASS(X,5,*n,*p);
  int ee = *p * (*p - 1) / 2;

  for(i = 0; i < *n * *n;i++) cluster_mat[i] = 0;
  for(i = 0; i < *n * ee; i++) edge_mat[i] = 0;

  //---------------- Burn In ------------------
  for(k = 0; k < *burn; k++)
    {
      if((k % *print_every) == 0)Rprintf("***************** Burn Step %d *************\n", k);
      state->UpdateState();
      //      printf("N:%d\n",state->N);
      R_CheckUserInterrupt();
    }
  //----------------------------------------------

  //----------------- Repititions ----------------
  for(k = 0; k < *reps; k++)
    {
      if((k % *print_every) == 0)Rprintf("***************** Rep Step %d *************\n", k);
      state->UpdateState();
      R_CheckUserInterrupt();
      //------------------- Update Results ------------------------------------
      for(i = 0; i < *n; i++)
	{
	  state_i = state->xi[i];
	  for(j = 0; j < *n; j++)
	    {
	      state_j = state->xi[j];
	      cluster_mat[i * *n + j] += (state_i == state_j);
	    }
	}
      for(i = 0; i < *n; i++)
	{
	  t = 0;
	  for(q = 0; q < *p - 1; q++)
	    {
	      for(r = q + 1; r < *p; r++)
		{
		  edge_mat[i * ee + t] += state->GGMlist[ state->xi[i] ]->graph->Edge[q][r];
		  t++;
		}
	    }
	}
      //----------------------------------------------------------------------------
    }
  //---------------------------------------------------------------------
  
  delete state;
  //
  PutRNGstate();
}//run_dmp_slice_routine;



}
