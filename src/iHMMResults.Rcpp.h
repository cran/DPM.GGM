#include "iHMMresults.h"
// This is the code for the results class of the iHMM routine

// Alex Lenkoski  lenkoski@stat.washington.edu

//Initialize, basically copy pointers over from R
iHMMResults::iHMMResults(int *cluster_mat_R,int *edge_mat_R,double *alpha_R,
			 double *alpha0_R, int *cluster_total_R, int *predict_R,
			 double *stepahead_R_x, double *stepahead_R_mu, double *stepahead_R_K, int *all_k_R, int *tot_reps_R, double *Kall_R)
{
  cluster_mat = cluster_mat_R;
  edge_mat = edge_mat_R;
  alpha = alpha_R;
  alpha0 = alpha0_R;
  cluster_total = cluster_total_R;
  predict = *predict_R;
  stepahead_x = stepahead_R_x;
  stepahead_mu = stepahead_R_mu;
  stepahead_K = stepahead_R_K;
  all_k = *all_k_R;
  tot_reps = *tot_reps_R;
  Kall = Kall_R;
}

//Deconstructor, let go of pointers so R can have them back
iHMMResults::~iHMMResults()
{
  cluster_mat = NULL;
  edge_mat = NULL;
  alpha = NULL;
  alpha0 = NULL;
  cluster_total = NULL;
  stepahead_x = NULL;
  stepahead_mu = NULL;
  stepahead_K = NULL;
  Kall = NULL;
}

//Main update routine.  Just stores some information about the state to send back to R
void iHMMResults::UpdateResults(int rep, State state)
{

  int i,j,n,p,ee;
  n = state->n;
  p = state->p;
  ee = p * (p - 1) / 2;
  int state_i, state_j;
  int t,q,r;
  //------------------- Update Cluster Information -------------------------------
  for(i = 0; i < n; i++)
    {
      state_i = state->xi[i];
      for(j = 0; j < n; j++)
	{
	  state_j = state->xi[j];
	  cluster_mat[i * n + j] += (state_i == state_j);
	}
    }
  //-----------------------------------------------------------------------------

  //------------------ Update Edge Information for each obs --------------------
  for(i = 0; i < n; i++)
    {
      t = 0;
      for(q = 0; q < p - 1; q++)
	{
	  for(r = q + 1; r < p; r++)
	    {
	      edge_mat[i * ee + t] += state->GGMlist[ state->xi[i] ]->graph->Edge[q][r];
	      t++;
	    }
	}
    }
  //----------------------------------------------------------------------------

  //------------- Update some parameter information --------
  alpha[rep] = state->alpha;
  alpha0[rep] = state->alpha0;
  cluster_total[rep] = state->L;
  //--------------------------------------------------------

  //--------- If we're getting step ahead predictions ------
  /*
  if(predict)
    {
      x_ahead = new double[p];
      mu_ahead = new double[p];
      K_ahead = new double[p * p];
      state->SampleNext(x_ahead, mu_ahead, K_ahead);
      for(r = 0; r < p; r++) stepahead_x[rep * p + r] = x_ahead[r];
      for(r = 0; r < p; r++) stepahead_mu[rep * p + r] = mu_ahead[r];
      for(r = 0; r < p * p; r++) stepahead_K[rep * p * p + r] = K_ahead[r];
      delete[] x_ahead;
      delete[] mu_ahead;
      delete[] K_ahead;
    }
  
  //---------------------------------------------------------

  if(all_k)
    {
      Kall_temp = new double[p * p * state->n];
      state->SampleAllK(Kall_temp);
      for(r = 0; r < state->n; r++)
	{
	  for(i = 0; i < p * p; i++)
	    {
	      Kall[r * p * p + i] += Kall_temp[r * p * p + i] / tot_reps;
	    }
	}
      delete[] Kall_temp;
    }
  */
}
