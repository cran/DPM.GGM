#include "iHMMState.h"

//------ The functions for the iHMMState class -----------------------
// These functions take the current state of the iHMM and update it
//--------------------------------------------------------------------

//Alex Lenkoski  lenkoski@stat.washington.edu


//-------------------------- Utility Functions -----------------------
//Count the transitions between two clusters
int transition_count(int *xi, int n, int from, int to)
{
  int count = 0;
  int i;
  for(i = 0; i < n - 1; i++)
    {
      if( (xi[i] == from) & (xi[i + 1] == to) ) count++;
    }
  return(count);
}

//Count the transitions out of a cluster
int transition_out_count(int *xi, int n, int from)
{
  int count = 0;
  int i;
  for(i = 0; i < n - 1; i++)
    {
      if( (xi[i] == from) & (xi[i + 1] != from)) count++;
    }
  return(count);
}

//  Calculates the transition coefficient for obs j being in cluster l
//  Note: this function looks long, but its really just to take care
//  of all the side cases correctly.  There's not too much to it.
double transition_coefficient(int *xi, int n, int j, int l, int L, double *gamma, double alpha)
{
  double a,b,c;
  int *xi_before;
  int *xi_after;
  int i;

  //--------- For the new cluster ------------------
  if(l == L)
    {
      if(j == (n - 1)) return(alpha * gamma[l]);
      return(alpha * gamma[l] * gamma[ xi[j + 1] ]);
    }
  //------------------------------------------------

  //-------- First time point ----------------------
  if(j == 0)
    {
      xi_after = new int[n - 1];
      for(i = 0; i < n - 1; i++)xi_after[i] = xi[i + (j + 1)];
      a = transition_count(xi_after, n - 1, xi[j + 1], l) + alpha * gamma[l];
      b = transition_count(xi_after, n - 1, l, xi[j + 1]) + alpha * gamma[ xi[j + 1] ];
      c = transition_out_count(xi_after, n - 1, l) + alpha + 1;
      delete[] xi_after;
      if(xi[j] != xi[j + 1])return(a * b / c);
      if(xi[j] == xi[j + 1])return(a * (b + 1) / c);
    }
  //-----------------------------------------------

  //-------- Last time point ----------------------
  if(j == (n - 1) )
    {
      xi_before = new int[n - 1];
      for(i =0; i < j; i++)xi_before[i] = xi[i];
      a = transition_count(xi_before, n - 1, xi[j - 1], l) + alpha * gamma[l];
      delete[] xi_before;
      return(a);
    }
  //-----------------------------------------------			 

  //--------- Time point 2 ------------------------
  if( j == 1 )
    {
      xi_after = new int[n - 2];
      for(i = 0; i < n - 2; i++) xi_after[i] = xi[i + j + 1];
      a = transition_count(xi_after, n - 2, xi[j - 1], l) + alpha * gamma[l];
      b = transition_count(xi_after, n - 2, xi[j + 1], l) + alpha * gamma[l];
      c = transition_out_count(xi_after, n - 2, l) + alpha;
      delete[] xi_after;
      if(xi[j - 1] != l) return(a * b / c);
      if(xi[j - 1] == xi[j + 1]) return( a * (b + 1) / (c + 1));
      if(xi[j - 1] != xi[j + 1]) return( a * b / (c  + 1));
    }
  //-----------------------------------------------

  //--------- Time point n - 1 --------------------
  if( j == (n - 2) )
    {
      xi_before = new int[n - 2];
      for(i = 0; i < n - 2; i++)xi_before[i] = xi[i];
      a = transition_count(xi_before, n - 2, xi[j - 1], l) + alpha * gamma[l];
      b = transition_count(xi_before, n - 2, xi[j + 1], l) + alpha * gamma[l];
      c = transition_out_count(xi_before, n - 2, l) + alpha;
      delete[] xi_before;
      if(xi[j - 1] != l) return( a * (b + 1) / (c + 1));
      if(xi[j - 1] == xi[j + 1]) return(a * (b + 1) / (c + 1));
      if(xi[j - 1] != xi[j + 1]) return(a * b / (c + 1) );
    }
  //--------------------------------------------------

  //--------- If we're at a middle time point --------
  if( (j > 1) & (j < n - 2) )
    {
      
      xi_before = new int[j];
      xi_after = new int[n  - (j + 1)];
      for(i = 0; i < j; i++) xi_before[i] = xi[i];
      for(i = 0; i < n - (j + 1); i++) xi_after[i] = xi[j + i + 1];
      a = transition_count(xi_before, j, xi[j - 1], l) + transition_count(xi_after, n - (j + 1), xi[j - 1], l) + alpha * gamma[l];
      b = transition_count(xi_before, j, xi[j + 1], l) + transition_count(xi_after, n - (j + 1), xi[j + 1], l) + alpha * gamma[l];
      c = transition_out_count(xi_before, j, l) + transition_out_count(xi_after, n - (j + 1), l) + alpha;
      delete[] xi_before;
      delete[] xi_after;
      if(xi[j - 1] != l) return(a * b / c);
      if(xi[j - 1] == xi[j + 1]) return(a * (b + 1) / (c + 1));
      if(xi[j - 1] != xi[j + 1]) return( a * b / (c + 1));
    }
  //--------------------------------------------------	
}//transition_coefficient


// Samples the parameter m
int sample_m(int r, double alpha, double gamma)
{
  double *prob_m;
  double max_m;
  double sum_m;
  int sample_m;
  int i;
  //  printf("Sampling m with alpha: %f and gamma: %f\n", alpha, gamma);
  if(r == 0)return(1);

  prob_m = new double[r];

  for(i = 0; i < r; i++)
    {
      //      printf("for %d and %d stirling was: %f\n", r, i, GLOBAL_STIRLING[(r - 1) * GLOBAL_STIRLING_DIM + i]);
      prob_m[i] = GLOBAL_STIRLING[(r - 1) * GLOBAL_STIRLING_DIM + i] + (i + 1) * log(alpha) + (i + 1) * log(gamma);
    }
  //  util_print_mat_dbl(prob_m, 1, r,0);
  max_m = prob_m[0];
  for(i = 1; i < r; i++)if(max_m < prob_m[i])max_m = prob_m[i];
  for(i = 0; i < r; i++) prob_m[i] = exp(prob_m[i] - max_m);
  sum_m = 0;
  for(i = 0; i < r; i++) sum_m += prob_m[i];
  for(i = 0; i < r; i++) prob_m[i] = prob_m[i] / sum_m;
  //  util_print_mat_dbl(prob_m,1,r,0);
  double m_avg = 0.0;
  sample_m = rand_int_weighted(r, prob_m) + 1;
  //  printf("We had %d and %d sampled with %f and %f and the average was %f\n", r, sample_m, alpha, gamma, m_avg);
  delete[] prob_m;
  return(sample_m);
}




//--------------------------------------------------------------------
//  This calculates the predictive distribution score for
//  observation k belonging to cluster l, which has graph thegraph
double iHMMState::PredictiveDistribution(int k, int l, GGM  G)
{

  //-------  Declarations -------------
  int i,j;
  double a;
  double n0 = .01;
  double *mu0 = new double[p];
  int n_group;
  double *mu_bar = new double[p];
  double *mu_tilde = new double[p];
  for(i = 0; i < p; i++)mu0[i] = 0;
  Matrix D_post;
  Matrix D_prior;
  double score;
  double J_G;
  double Norm_terms;
  //------------------------------------

  //---- Form sufficient statistics for observations currently in l -------
  //     Note: this is the "prior" information as it excludes obs k 
  for(i = 0; i < p; i++) mu_bar[i] = (G->n * G->xbar[i] + n0 * mu0[i]) / (G->n + n0);
  D_prior = G->DplusU->copy();
  for(i = 0; i < p; i++)
    {
      for(j = 0; j < p; j++)
	{
	  a = -(G->n + n0) * mu_bar[i] * mu_bar[j] + G->n * G->xbar[i] * G->xbar[j] + n0 * mu0[i] * mu0[j];
	  D_prior->set(i,j,D_prior->get(i,j) + a);
	}
    }
  //------------------------------------------------------------------------

  //----- Now get "posterior" parameters, which include obs k --------------
  D_post = D_prior->copy();
  for(i = 0; i < p; i++) mu_tilde[i] = (X->get(k,i) + (G->n + n0) * mu_bar[i]) / (G->n + n0 + 1);
  for(i = 0; i < p; i++)
    {
      for(j = 0; j < p; j++)
	{
	  a = -(G->n + 1 + n0) * mu_tilde[i] * mu_tilde[j] + X->get(k,i) * X->get(k,j) + (G->n + n0) * mu_bar[i] * mu_bar[j];
	  D_post->set(i, j, D_post->get(i, j) + a);
	}
    }
  //--------------------------------------------------------------------------
  
  //------- Calculate the ratio of the two normalizing constants -------------
  J_G = j_g_decomposable(G->graph, D_prior, D_post, 3 + G->n, 1);
  //---------------------------------------------------------------------------

  //------ The stuff that dangles on this, due to the mean -------------------
  Norm_terms = - p / 2 * log(2 * PI) + p / 2 * log(( (double) G->n + n0) / ( (double) G->n + 1.0 + n0));
  //--------------------------------------------------------------------------

  //------- And return the score ------
  score =  Norm_terms + J_G;
  //-----------------------------------

  //------- Clean up ------------------
  delete D_post;
  delete D_prior;
  delete[] mu_bar;
  delete[] mu_tilde;
  delete[] mu0;
  //------------------------------------

  return(score);
}

// Initialization Routine 		    
iHMMState::iHMMState(double *data, int L_start, int n_obs, int p_model)
{
  int i,k,l;
  int temp;
  n = n_obs;
  p = p_model;

  X = new Mat(n,p,data);
  L = L_start;
  GGMlist = new GGM[L];
  for(i = 0; i < L; i++)
    {
      GGMlist[i] = new GGM_CLASS(p);
    }
  xi = new int[n];
  gamma = new double[L + 1];
  //-------------------------------------------------------

  //-------------- Initalize xi ---------------------------
  for(i = 0; i < n; i++)
    {
      xi[i] = rand_int(L);
      GGMlist[ xi[i] ]->AddObs(X,i);
    } 
  //-------------------------------------------------------

  //------------- Initialize gamma ------------------------
  for(i = 0; i < L + 1; i++) gamma[i] = 1.0 / (L + 1);
  //  printf("************* And our starting value *****\n");
  //  util_print_mat_dbl(gamma, 1, L + 1, 0);
  //-------------------------------------------------------

  //------------- Initialize alpha ------------------------
  alpha = 1;
  alpha0 = 1;
  //-------------------------------------------------------

}

// Deconstructor
iHMMState::~iHMMState()
{
  int i;
  delete X;
  delete[] xi;
  for(i = 0; i < L; i++) delete GGMlist[i];
  delete[] GGMlist;
  delete[] gamma;
}


// Code to update gamma
void iHMMState::UpdateGamma()
{
  int *M, l, lprime;
  int r_llprime;
  int i;
  double *M_dot;

  M_dot = new double[L + 1];

  for(i = 0; i < L + 1; i++)M_dot[i] = 0;

  for(l = 0; l < L; l++)
    {
      for(lprime = 0; lprime < L; lprime++)
	{
	  r_llprime = transition_count(xi, n, l, lprime);
	  M_dot[l] += sample_m(r_llprime, alpha, gamma[lprime]);
	}
    }
  M_dot[L] = alpha0;
  rdirichlet(M_dot, L + 1, gamma);
  delete[] M_dot;
}

//  Code to update Alpha and Alpha0
void iHMMState::UpdateAlpha()
{
  double a = 1;
  double b = 1;
  double a0 = 1;
  double b0 = 1;
  double eta;
  double odds_d_eta;
  double d_eta;

  int m_sum;
  int r_llprime;
  int l, lprime;
  double sum_u, sum_log_zeta;

  int *rldot;
  rldot = new int[L];
  for(l = 0; l < L; l++)rldot[l] = 0;
  //  printf("Here's xi\n");
  //  util_print_mat_int(xi,1,n,0);

  //------------- Get the sum of the matrix M --------------
  m_sum = 0;
  for(l = 0; l < L; l++)
    {
      for(lprime = 0; lprime < L; lprime++)
	{
	  r_llprime = transition_count(xi, n, l, lprime);
	  m_sum += sample_m(r_llprime, alpha, gamma[lprime]);
	  //	  printf("l is: %d lprime is: %d rllprime is: %d m_sum is: %d\n",l,lprime, r_llprime, m_sum);
	}
    }
  //  printf("M sum is: %d\n", m_sum);
  //---------------------------------------------------------

  //---------- Update Alpha ---------------------------------
  sum_log_zeta = 0;
  sum_u = 0;
  for(l = 0; l < L; l++) rldot[l] = transition_out_count(xi, n, l) + transition_count(xi, n, l, l);
  //  printf("Here's rldot:\n");
  //  util_print_mat_int(rldot, 1, L, 0);
  for(l = 0; l < L; l++)
    {
      if(rldot[l] != 0)
	{
	  sum_log_zeta +=  log(rbeta(alpha + 1, rldot[l]));
	  sum_u += rbinom(1, rldot[l] / (alpha + rldot[l]));
	}
    }
  //  printf("Sum u: %f Sum Zeta: %f\n", sum_u, sum_log_zeta);
  //  printf("And parameters are: %f and %f\n", a + m_sum - sum_u, 1.0 / (b - sum_log_zeta));

  alpha = rgamma(a + m_sum -  sum_u, 1.0 / (b - sum_log_zeta));
  //  printf("and now alpha: %f\n", alpha);
  //---------------------------------------------------------

  //---------- Update Alpha0 --------------------------------
  eta = rbeta(alpha0 + 1, m_sum);
  odds_d_eta = (a0 + L - 1) / (m_sum * (b - log(eta)));
  d_eta = odds_d_eta / (1 + odds_d_eta);
  //  printf("eta: %f log(eta): %f odds_eta: %f d_eta: %f\n",eta, log(eta), odds_d_eta, d_eta); 

    if(unif_rand() < d_eta)
      {
        alpha0 = rgamma(a0 + L, 1.0 / (b0 - log(eta)) );
      }
    else
      {
	alpha0 = rgamma(a0 + L - 1, 1.0 / (b0 - log(eta)) );
      }
  //---------------------------------------------------------

  delete[] rldot;

  return;
}

//  Code to update the cluster vector, xi
void iHMMState::UpdateXi()
{
  int i, j, k, l;
  GGM *GGMlist_new;
  GGM newGGM;
  GGM *tempGGMlist;
  double *gamma_temp, *gamma_new;
  double nu;
  double temp;

  LPGraph newgraph;
  int xi_old, xi_new;
  int others;
  double *qs;
  double maxq, sumq;

  //  The big loop through the observations
  for(i = 0; i < n; i++)
    {
      //      printf("Observation %d Clusters %d\n",i, L);      
      //      util_print_mat_int(xi, 1, n, 0);
      xi_old = xi[i];
      xi[i] = -1;
      GGMlist[xi_old]->DropObs(X,i);

      //----------------- Check to see if we need to drop a group -------------
      if(GGMlist[xi_old]->n == 0)
	{
	  //	  printf("Dropping a cluster, current size is: %d\n",L);
	  GGMlist_new = new GGM[L - 1];
	  gamma_new = new double[L];
	  //	  printf("Here's what gamma looks like now\n");
	  //	  util_print_mat_dbl(gamma, 1, L + 1,0);
	  for(l = 0; l < L; l++)
	    {
	      if(l < xi_old)
		{
		  GGMlist_new[l] = GGMlist[l];
		  gamma_new[l] = gamma[l];
		}
	      if(l > xi_old)
		{
		  GGMlist_new[l - 1] = GGMlist[l];
		  gamma_new[l - 1] = gamma[l];
		}
	    }
	  gamma_new[L - 1] = gamma[L];
	  gamma_new[L - 1] += gamma[xi_old];
	  gamma_temp = gamma;
	  gamma = gamma_new;
	  delete[] gamma_temp;
	  //	  printf("And here's gamma now\n");
	  //	  util_print_mat_dbl(gamma, 1, L, 0);
	  tempGGMlist = GGMlist;
	  GGMlist = GGMlist_new;
	  delete tempGGMlist[xi_old];
	  delete[] tempGGMlist;
	  for(j = 0; j < n; j++)
	    {
	      if(xi[j] > xi_old) xi[j] = xi[j] - 1;
	    }
	  L--;
	}//what a hassle.
      //--------------------------------------------------------------------------
      
      //----------- Predictive Distributions, Existing Clusters ------------------ 
      qs = new double[L + 1];
      for(l = 0; l < L; l++)
	{
	  qs[l] = PredictiveDistribution(i,l,GGMlist[l]);
	  qs[l] += log(transition_coefficient(xi, n, i, l, L, gamma, alpha));
	}
     //---------------------------------------------------------------------------

      //----------- New Cluster ---------------------------------------------------
      newGGM = new GGM_CLASS(p);
      qs[L] = PredictiveDistribution(i, L, newGGM);
      qs[L] += log(transition_coefficient(xi, n, i, L, L, gamma, alpha));
      //------------------------------------------------------------------------------

      //------------ Make Proposal ---------------------------------------------------
      logs_to_probs(qs, L + 1);
      //      util_print_mat_dbl(qs, 1, L + 1, 0);
      xi_new = rand_int_weighted(L + 1, qs);
      //      printf("And the new value is: %d\n",xi_new);
      ///--------------------------------------------------------------------------------

      //------------ Clean Up ----------------------------------------------------------
      xi[i] = xi_new;
      if(xi[i] == L)//If we're adding a cluster, we have to do a bit of work, this is all bookkeeping
	{
	  GGMlist_new = new GGM[L + 1];
	  gamma_new = new double[L + 2];

	  for(l = 0; l < L; l++)
	    {
	      GGMlist_new[l] = GGMlist[l];
	      gamma_new[l] = gamma[l];
	    }
	  
	  GGMlist_new[L] = newGGM;
	  tempGGMlist = GGMlist;
	  GGMlist = GGMlist_new;
	  delete[] tempGGMlist;

	  nu = rbeta(alpha0, 1);
	  temp = gamma[L];
	  gamma_new[L] = nu * temp;
	  gamma_new[L + 1] = (1 - nu) * temp;
	  gamma_temp = gamma;
	  gamma = gamma_new;
	  newGGM = NULL;
	  delete[] gamma_temp;
	  L++;
	}
      else //If we're not adding a new cluster, its pretty easy
	{
	delete newGGM;
	}
      delete[] qs;
      //----------------------------------------------------------------------------------
      GGMlist[xi_new]->AddObs(X,i);
    }// End of big loop through the observations
}

//   This code loops through the graphs for each cluster and updates them
void iHMMState::UpdateG()
{
  int l;
  for(l = 0; l < L; l++)
    {
      GGMlist[l]->UpdateG();
    }
}

//  Administrative routine that calls the routines to update each parameter
void iHMMState::UpdateState()
{

  UpdateXi();

  UpdateAlpha();
  UpdateGamma();
  UpdateG();

}

void iHMMState::Print()
{
  int i;
  for(i = 0; i < L;i++)
    {
      Rprintf("Model %d of %d\n",i + 1,L);
      GGMlist[i]->PrintInfo();
    }
}

/*
//  This routine samples one step ahead
void iHMMState::SampleNext(double *x, double *mu, double *K)
{
  //  printf("I got here\n");

  //-----  Set-up -------------
  int i,j,k,l;
  double a;
  double n0 = .01;
  double *mu0 = new double[p];
  double *qs;
  double maxq, sumq;
  int xi_new;
  LPGraph thegraph;
  for(i = 0; i < p; i++) mu0[i] = 0;
  int temp;
  int n_group;
  double *xbar = new double[p];
  double *mu_bar = new double[p];
  double *D = new double[p * p];
  int n_sub = 0;
  double score;
  double J_G;
  double Norm_terms;
  //------------------------------

  //----- Predict the cluster that the next observation will land in -------------
  qs = new double[L + 1];
  for(l = 0; l < L; l++)
    {
      qs[l] = log(transition_coefficient(xi, n, n - 1, l, L, gamma, alpha));
    }
  
  qs[L] = log(transition_coefficient(xi, n, n - 1, L, L, gamma, alpha));

  maxq = qs[0];
  for(l = 0; l < L + 1; l++) if(qs[l] > maxq) maxq = qs[l];
  sumq = 0;
  for(l = 0; l < L + 1;l++) sumq += exp(qs[l] - maxq);
  for(l = 0; l < L + 1;l++) qs[l] = exp(qs[l] - maxq) / sumq;
  //      util_print_mat_dbl(qs, 1, L + 1, 0);
  //  util_print_vec_dbl(qs, L + 1);
  xi_new = rand_int_weighted(L + 1, qs);
  //------------------------------------------------------------------------------
  
  //------- In the off chance that we think we'll be in a new cluster ------------
  if(xi_new == L)
    {
      //      printf("I better not be here\n");
      thegraph = new Graph;
      thegraph->InitGraph(p);
      for(k = 0; k < p - 1; k++)
	{
	  for(l = k + 1; l < p; l++)
	    {
	      temp = (unif_rand() < .5);
	      thegraph->Edge[k][l] = temp;
	      thegraph->Edge[l][k] = temp;
	    }
	}
      TurnFillInGraph(thegraph);
      if(thegraph->IsDecomposable())
	{
	  //	  printf("Worked for this \n");
	}
    }
  else
    {
      //      printf("And I did get here\n");
      thegraph = graphlist[xi_new];
    }
  //-------------------------------------------------------------------------------

  //--------  Now we sample -------------------------------------------------------
  for(i = 0; i < n; i++) if(xi[i] == xi_new) n_sub++;
  make_sub_means(X, xi, xi_new, p, n, n_sub, xbar);
  make_sub_cov(X, xi, xi_new, p, n, n_sub, D);
  for(i = 0; i < p; i++)D[i * p + i] += 1;
  for(i = 0; i < p; i++) mu_bar[i] = (n_sub * xbar[i] + n0 * mu0[i]) / (n_sub + n0);
  if(thegraph->nCliques == 1)
    {
      rwish(p, 3 + n_sub, D, K);
    }
  else
    {
      gwish_blgibbs_decomposable(thegraph, 3 + n_sub, D, K);
    }
  for(i = 0; i < p * p; i++) K[i] = K[i] * (n_sub + n0);
  sample_mult_norm(p, mu_bar, K, mu);
  for(i = 0; i < p * p; i++) K[i] = K[i] / (n_sub + n0);
  sample_mult_norm(p, mu, K, x);
  //--------------------------------------------------------------------------------
  delete[] mu0;
  delete[] qs;
  delete[] xbar;
  delete[] mu_bar;
  delete[] D;
  if(xi_new == L)delete thegraph;
}

//  This routine samples one step ahead
void iHMMState::SampleAllK(double *Kall)
{

  //-----  Set-up -------------
  int i,j,k,l;
  double a;
  double n0 = .01;
  LPGraph thegraph;
  int n_group;
  double *D = new double[p * p];
  int n_sub = 0;
  double *K = new double[p * p];
  //------------------------------
  for(l = 0; l < L; l++)
    {
      //--------  Now we sample ------------------
      
      for(i = 0; i < p * p; i++)K[i] = 0;
      for(i = 0; i < p * p; i++)D[i] = 0;
      thegraph = graphlist[l];
      for(i = 0; i < n; i++) if(xi[i] == l) n_sub++;
      make_sub_cov(X, xi, l, p, n, n_sub, D);
      if(thegraph->nCliques == 1)
	{
	  rwish(p, 3 + n_sub, D, K);
	}
      else
	{
	  gwish_blgibbs_decomposable(thegraph, 3 + n_sub, D, K);
	}
      for(i = 0; i < n; i++)
	{
	  if(xi[i] == l)
	    {
	      for(j = 0; j < p * p; j++)
		{
		  Kall[j + (p * p) * i] = K[j];
		}
	    }
	}
    }
  //--------------------------------------------------------------------------------
  delete[] D;
  delete[] K;
}
*/
