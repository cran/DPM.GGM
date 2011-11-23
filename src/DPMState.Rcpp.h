#include "DPMState.h"

double DPMState_CLASS::PredictiveDistribution(int k, int l, GGM G)
{

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
  double coef = log(alpha);
  double score;
  double J_G;
  double Norm_terms;

  if(G->n > 0) coef = log(G->n);

  D_prior = G->DplusU->copy();

  for(i = 0; i < p; i++) mu_bar[i] = (G->n * G->xbar[i] + n0 * mu0[i]) / (G->n + n0);
  for(i = 0; i < p; i++)
    {
      for(j = 0; j < p; j++)
	{
	  D_prior->set(i,j, D_prior->get(i,j) - (G->n + n0) * mu_bar[i] * mu_bar[j] + G->n * G->xbar[i] * G->xbar[j] + n0 * mu0[i] * mu0[j]);
	}
    }
  
  D_post = D_prior->copy();
  for(i = 0; i < p; i++) mu_tilde[i] = (X->get(k,i) + (G->n + n0) * mu_bar[i]) / (G->n + n0 + 1);
  for(i = 0; i < p; i++)
    {
      for(j = 0; j < p; j++)
	{
	  a = -(G->n + 1 + n0) * mu_tilde[i] * mu_tilde[j] + X->get(k,i) * X->get(k,j) + (G->n + n0) * mu_bar[i] * mu_bar[j];
	  D_post->set(i,j,D_post->get(i,j) + a);
	}
    }

  J_G = j_g_decomposable(G->graph, D_prior, D_post, 3 + G->n, 1);
  Norm_terms = - p / 2 * log(2 * PI) + p / 2 * log((G->n + n0) / (G->n + 1 + n0));
  score =  Norm_terms + J_G;
  
  delete D_post;
  delete D_prior;
  delete[] mu_bar;
  delete[] mu_tilde;
  delete[] mu0;
  
  return(coef + score);
}
		    
DPMState_CLASS::DPMState_CLASS(double *data, int L_start, int n_obs, int p_model)
{
  int i,k,l;
  int temp;

  X = new Mat(n_obs,p_model,data);
  L = L_start;
  n = n_obs;
  p = p_model;
  GGMlist = new GGM[L_start];
  xi = new int[n];

  //--------------- Initialize the graphs ---------------
  for(i = 0; i < L; i++)
    {
      GGMlist[i] = new GGM_CLASS(p);
      TurnFillInGraph(GGMlist[i]->graph);
      if(GGMlist[i]->graph->IsDecomposable())
	{
	  //	  printf("Worked for this \n");
	}
    }
  //-------------------------------------------------------

  //-------------- Initalize xi ---------------------------
  for(i = 0; i < n; i++)
    {
      xi[i] = rand_int(L);
      GGMlist[ xi[i] ]->AddObs(X,i);
    } 
  //-------------------------------------------------------

  //------------- Initialize alpha ------------------------
  alpha = 1;
  //-------------------------------------------------------

}

DPMState_CLASS::~DPMState_CLASS()
{
  int i;
  delete X;
  delete[] xi;
  for(i = 0; i < L; i++) delete GGMlist[i];
  delete[] GGMlist;
}



void DPMState_CLASS::UpdateAlpha()
{
  double a = 1;
  double b = 1;
  double eta = rbeta(alpha + 1, n);
  double odds = (a + L - 1) / (n * b - log(eta));
  double pi_eta = odds / (1 + odds);
  if(unif_rand() < pi_eta)
    {
      alpha = rgamma(a + L, 1 / (b - log(eta)));
    }
  else
    {
      alpha = rgamma(a + L - 1, 1/ (b - log(eta)));
    }
  return;
}

void DPMState_CLASS::UpdateXi()
{
  int i, j, k, l;
  GGM *GGMlist_new;
  GGM *tempGGMlist;

  GGM newGGM;
  int xi_old, xi_new;
  double *qs;
  double temp;

  for(i = 0; i < n; i++)
    {
      //      printf("Observation %d\n",i);      
      //      util_print_mat_int(xi, 1, n, 0);
      xi_old = xi[i];
      xi[i] = -1;
      GGMlist[ xi_old ]->DropObs(X,i);
      //----------------- Check to see if we need to drop a group -------------
      if(GGMlist[ xi_old ]->n == 0)
	{
	  GGMlist_new = new GGM[L - 1];
	  for(l = 0; l < L; l++)
	    {
	      if(l < xi_old)
		{
		  GGMlist_new[l] = GGMlist[l];
		}
	      if(l > xi_old)
		{
		  GGMlist_new[l - 1] = GGMlist[l];
		}
	    }
	  tempGGMlist = GGMlist;
	  GGMlist = GGMlist_new;
	  delete tempGGMlist[xi_old];
	  delete[] tempGGMlist;

	  for(j = 0; j < n; j++)
	    {
	      if(xi[j] > xi_old) xi[j] = xi[j] - 1;
	    }
	  L--;
	}
      //--------------------------------------------------------------------------
      
      //----------- Predictive Distributions, Existing Clusters ------------------
      qs = new double[L + 1];
      for(l = 0; l < L; l++)
	{
	  //	  printf("Here's the graph\n");
	  //	  Print_A(graphlist[l]);
	  qs[l] = PredictiveDistribution(i,l,GGMlist[l]);
	  //	  printf("And the score: %f\n", qs[l]);
	}
      //---------------------------------------------------------------------------

      //----------- New Cluster ---------------------------------------------------
      newGGM = new GGM_CLASS(p);
      TurnFillInGraph(newGGM->graph);
      qs[L] = PredictiveDistribution(i, L, newGGM);
      //------------------------------------------------------------------------------

      //      util_print_mat_dbl(qs, 1, L + 1, 0);

      //------------ Make Proposal ---------------------------------------------------
      logs_to_probs(qs, L + 1);
      xi_new = rand_int_weighted(L + 1, qs);
      //--------------------------------------------------------------------------------

      //------------ Clean Up ----------------------------------------------------------
      xi[i] = xi_new;
      if(xi[i] == L)
	{
	  //	  printf("Adding Cluster\n");
	  GGMlist_new = new GGM[L + 1];
	  for(l = 0; l < L; l++)
	    {
	      GGMlist_new[l] = GGMlist[l];
	    }
	  GGMlist_new[L] = newGGM;

	  tempGGMlist = GGMlist;
	  GGMlist = GGMlist_new;
	  delete[] tempGGMlist;
	  L++;
	}
      else
	{
	  delete newGGM;
	}
      delete[] qs;
      //      util_print_mat_int(xi, 1, n, 0);
      //----------------------------------------------------------------------------------
      GGMlist[xi_new]->AddObs(X,i);
    }
}


void DPMState_CLASS::UpdateG()
{
  int i;
  for(i = 0; i < L; i++)GGMlist[i]->UpdateG();
}



void DPMState_CLASS::UpdateState()
{
  //  printf("Updating Alpha\n");
  UpdateAlpha();
  //  printf("Alpha is %f\n", alpha);
  //  printf("Updating Xi\n");
  UpdateXi();
  //  util_print_mat_int(xi, 1, n,0);
  //  printf("Updating G\n");
  UpdateG();
  //  printf("Number of Clusters: %d\n",L);
}














  
