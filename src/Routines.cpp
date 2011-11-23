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
  int i,j,k,l,q,r,t;
  GetRNGstate();
  int state_i, state_j;

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
  int i,j,k,l,q,r,t;
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

  int i,j,k,l,q,r,t;
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
 /*
extern "C"{
  void test(double *X, int *p, int *n, double *stirling)
  {
    int i,j;

    //    double *x = sample_mult_norm(mu, Xmat);
    //    for(i = 0; i < *p; i++)Rprintf("%f ",x[i]);
    //    Rprintf("\n");

    /*
    printf("GET: And getting 3,4: %f\n",Xmat->get(2,3));
    Xmat->set(0,0,10.0);
    printf("SET: And now we have: %f\n",Xmat->get(0,0));
    Matrix Q = new Mat(*p,*p);
    Q->set_iden();
    printf("IDEN:\n");
    Q->print();
    delete Q;
    //    Xmat->print();
    Q = Xmat->copy();
    printf("Copy:\n");
    Q->print();
    delete Q;

    printf("Invert:\n");
    Q = Xmat->invert();
    Q->print();
    
    printf("Log det: %f\n", Q->log_det());
    delete Q;

    Q = Xmat->chol();
    printf("Chol:\n");
    Q->print();
    delete Q;

    Q = Xmat->transpose();
    printf("transpose: \n");
    Q->print();

    Q->Add(Xmat);
    printf("Add: \n");
    Q->print();
    
    delete Xmat;
    delete Q;
    */
    /*
    LPGraph graph = new Graph();
    graph->InitGraph(*p);
    for(i = 0; i < *p - 1; i++)
      {
	for(j = i + 1; j < *p; j++)
	  {
	    graph->Edge[i][j] = (unif_rand() < 0.5);
	    graph->Edge[j][i] = graph->Edge[i][j];
	  }
      }
    graph->GetMPSubgraphs();
    graph->PrintA();

    double temp = gwish_nc_complete(102,*p,Xmat);
    printf("Norm Constant: %f\n", temp);
    
    Matrix Ksm = rwish(99, Xmat);
    Ksm->print();

    Ksm->set_iden();
    Ksm = gwish_blgibbs_decomposable(graph, 99, Xmat);

    Ksm->print();
    delete graph;
    */
    /*
    GGM G = new GGM_CLASS(*p);
    util_print_vec_dbl(G->xbar, *p);
    //    G->PrintInfo();
    for(i = 0; i < *n; i++)
      {
	G->AddObs(Xmat, i);	
      }
    for(i =0; i < 50; i++)
      {
	G->DropObs(Xmat,i);
      }
    for(i = 0; i < 100; i++)
      {
	printf("i: %d\n",i);
	G->graph->PrintA();
	G->UpdateG();
      }
    delete G;
    
    State S;
    S = new iHMMState(X,5,*n,*p);
    S->UpdateXi();
    //    S->UpdateGamma();
    S->UpdateG();
    S->Print();
    return;


  }
*/



