#include "DPMSliceState.h"

double DPMSliceState_CLASS::PredictiveDistribution(int k, GGM G)
{

  int i;

  double score;
  Matrix xmat = new Mat(1,p);
  for(i = 0; i < p; i++) xmat->set(0,i,X->get(k,i) - G->mu[i]);
  
  score = 0.5 * G->K->log_det();
  //  printf("%f\n", score);
  score -= 0.5 * p * log(2 * M_PI);
  //  printf("%f\n", score);
  score -= 0.5 *  quad_form(xmat, G->K);
  //  G->K->print();
  /*  for(i = 0; i < p; i++)
    {
      printf("%f, ", X->get(k,i));
    }
    printf("\n");
  for(i = 0; i < p; i++)
    {
      printf("%f, ", G->mu[i]);
    }
  printf("\n");
  printf("Score: %f\n",score);exit(1);*/
  delete xmat;
  
  return(score);
}
		    
DPMSliceState_CLASS::DPMSliceState_CLASS(double *data, int N_start, int n_obs, int p_model)
{
  int i;


  X = new Mat(n_obs,p_model,data);
  N = N_start;
  n = n_obs;
  p = p_model;
  GGMlist = new GGM[N];
  xi = new int[n];

  //--------------- Initialize the graphs ---------------
  for(i = 0; i < N; i++)
    {
      //    printf("i: %d\n",i);
      GGMlist[i] = new GGM_CLASS(p);
      TurnFillInGraph(GGMlist[i]->graph);
      if(GGMlist[i]->graph->IsDecomposable())
	{
	  //	  printf("Worked for this \n");
	}
      GGMlist[i]->SampleParams();
      //      printf("Mean\n");
      //      util_print_vec_dbl(GGMlist[i]->mu, p);
      //      printf("Precision\n");
      //      GGMlist[i]->K->print();
    }
  //-------------------------------------------------------

  //-------------- Initalize xi ---------------------------
  r = new int[N];
  for(i = 0; i < N; i++)r[i] = 0;
  for(i = 0; i < n; i++)
    {
      xi[i] = rand_int(N);
      GGMlist[ xi[i] ]->AddObs(X,i);
      r[ xi[i] ]++;
    } 
  //  util_print_vec_int(xi, n);
  //  util_print_vec_int(r, N);
  //-------------------------------------------------------

  //------------- Initialize alpha ------------------------
  alpha = 1;
  //-------------------------------------------------------

  //------------- Initialize W ----------------------------
  w = new double[N];
  w[0] = 0.5;
  double w_sum = 0.5;
  double z;
  for(i = 1; i < N; i++)
    {
      z = rbeta(1,1);
      w[i] = z * (1 - w_sum);
      w_sum += w[i];
    }
  //  util_print_vec_dbl(w, N);  
  //  printf("WSUM: %f\n", w_sum);exit(1);
  UpdateW();

  u = new double[n];
  UpdateU();
  //  util_print_vec_dbl(u,n);
  //-------------------------------------------------------
}

DPMSliceState_CLASS::~DPMSliceState_CLASS()
{
  int i;
  delete X;
  delete[] xi;
  delete[] r;
  delete[] w;
  delete[] u;
  for(i = 0; i < N; i++) delete GGMlist[i];
  delete[] GGMlist;
}



void DPMSliceState_CLASS::UpdateAlpha()
{
  int l;
  int L = 0;
  for(l = 0; l < N;l++)
    {
      L += (r[l] > 0);
    }

  double a = 1.0;
  double b = 1.0;
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

void DPMSliceState_CLASS::UpdateXi()
{
  int i, l;




  int xi_old, xi_new;
  double *qs;
  double temp;
  int check;
  int *ones;
  for(i = 0; i < n; i++)
    {
      //      printf("Observation %d N: %d\n",i,N);      
      //      util_print_mat_int(xi, 1, n, 0);
      check = 1;
      while(check)
	{
	  temp = 1;
	  for(l = 0; l < N;l++)
	    {
	      temp -= w[l];
	    }
	  if(temp > u[i])
	    {
	      ExpandModel();
	    }
	  else
	    {
	      check=0;
	    }
	}
      xi_old = xi[i];
      xi[i] = -1;
      GGMlist[ xi_old ]->DropObs(X,i);
      r[xi_old]--;
      //      printf("Observation %d N: %d\n",i,N);      
      //--------------------------------------------------------------------------
      
      //----------- Predictive Distributions, Existing Clusters ------------------
      qs = new double[N];
      ones = new int[N];
      for(l = 0; l < N; l++)ones[l] = 0;
      //      printf("Here\n");
      for(l = 0; l < N; l++)
	{
	  //	  printf("l: %d\n",l);
	  if(u[i] < w[l])
	    {
	      ones[l] = 1;
	      qs[l] = PredictiveDistribution(i,GGMlist[l]);
	    }
	  else
	    {
	      qs[l] = 0.0;
	    }
	}

      //---------------------------------------------------------------------------

      logs_to_probs_ind(qs, N, ones);
      weights_to_probs(qs,N);
      xi_new = rand_int_weighted(N, qs);
      //--------------------------------------------------------------------------------

      //------------ Clean Up ----------------------------------------------------------
      xi[i] = xi_new;
      GGMlist[xi_new]->AddObs(X,i);
      r[xi_new]++;
      delete[] qs;
    }
}


void DPMSliceState_CLASS::UpdateGGM()
{
  int i;
  for(i = 0; i < N; i++)
    {
      GGMlist[i]->UpdateG();
      GGMlist[i]->SampleParams();
    }
}

void DPMSliceState_CLASS::UpdateW()
{
  int l;
  double r_sum = 0.0;
  double* z = new double[N];
  double z_prod = 1.0;
  for(l=0; l < N;l++) r_sum += (double) r[l];
  
  for(l = 0; l < N; l++)
    {
      r_sum -= (double) r[l];
      z[l] = rbeta(1 + r[l], alpha + r_sum);
    }
  for(l = 0; l < N; l++)
    {
      w[l] = z[l] * z_prod;
      z_prod = (1 - z[l]) * z_prod;
    }
  
}

void DPMSliceState_CLASS::UpdateU()
{
  int i;
  for(i = 0; i < n; i++)
    {
      u[i] = unif_rand() * w[ xi[i] ];
    }
}

void DPMSliceState_CLASS::ExpandModel()
{
  int l;
  double z;
  double w_sum = 0.0;
  double *w_new, *w_temp;
  int *r_new, *r_temp;
  GGM *newGGMlist,*tempGGMlist;

  w_new = new double[N + 1];
  newGGMlist = new GGM[N + 1];
  r_new = new int[N + 1];

  for(l = 0; l < N; l++)
    {
      w_new[l] = w[l];
      w_sum += w[l];
      newGGMlist[l] = GGMlist[l];
      r_new[l] = r[l];
    }

  z = rbeta(1,alpha);
  w_new[N] = z * (1 - w_sum);
  w_temp = w;
  w = w_new;

  newGGMlist[N] = new GGM_CLASS(p);
  newGGMlist[N]->SampleParams();
  tempGGMlist = GGMlist;
  GGMlist = newGGMlist;

  r_new[N] = 0;
  r_temp = r;
  r = r_new;

  N++;
  delete[] w_temp;
  delete[] tempGGMlist;
  delete[] r_temp;

}
void DPMSliceState_CLASS::UpdateState()
{
  //  printf("Updating Alpha\n");
  UpdateAlpha();
  //  printf("Alpha is %f\n", alpha);
  //  printf("Updating Xi\n");
  UpdateXi();
  //  util_print_vec_int(xi,n);
  //  printf("\n");
  //  util_print_vec_int(r,N);
  //  return;
  //  util_print_mat_int(xi, 1, n,0);
  //  printf("Updating G\n");
  UpdateGGM();
  /*  int i;
  for(i=0;i < N;i++)
    {
      printf("Graph number: %d\n",i);
      GGMlist[i]->PrintInfo();
    }
  */
  UpdateW();
  //  util_print_vec_dbl(w,N);
  UpdateU();
  //  printf("\n");
  //  util_print_vec_int(r, N);
  //  printf("\n");
  //  util_print_vec_dbl(u,n);
}














  
