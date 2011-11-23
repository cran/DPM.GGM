#include "random.h"

double* sample_mult_norm(double *mu, Matrix K)
{
  int i;
  int p = K->p_r;
  double* x = new double[p];
  Matrix Sigma = K->invert();
  Matrix x_normed = new Mat(1,p);
  Matrix Phi = Sigma->chol();

  for(i = 0; i < p; i++) x_normed->set(0,i,rnorm(0,1));
  Matrix x_new = mult_mats(x_normed, Phi);
  for(i = 0; i < p; i++) x[i] = x_new->get(0,i) + mu[i];

  delete Sigma;
  delete Phi;
  delete x_normed;
  delete x_new;

  return(x);
}

//Samples an object from a subset dictacted by the vector s (binary) where there are n_sub objects in s that equal 1
int sample_from(int *s, int n, int n_sub)
{
  int i;
  int k = 0;
  int sample = rand_int(n_sub);
  for(i = 0; i < n; i++)
    {
      if(s[i])
	{
	  if(k == sample)
	    {
	      return(i);
	    }
	  k++;
	}
    }
  return(-1);
}

//Your compiler might complain a little about this function but whatever.
int rand_int(int n)
{
  double temp;
  double alpha = unif_rand();
  modf(n * alpha, &temp);
  int value = (int) temp;
  return(value);
}

// Samples from a dirichlet distribution. l is the dimension, a is the parameter vector
void rdirichlet(double *a, int l, double *x)
{
  int i;
  double sum;
  for(i = 0; i < l; i++) x[i] = rgamma(a[i],1);
  sum = 0;
  for(i = 0; i < l; i++)sum+=x[i];
  for(i = 0; i < l; i++) x[i] = x[i] / sum;
}

//  Returns an integer between 0 and n-1 with prob according to weights
int rand_int_weighted(int n, double *weights)
{
  int i;
  double r = unif_rand();
  for(i = 0; i < n; i++)
    {
      if(r < weights[i])return(i);
      r -= weights[i];
    }
  return(-1);
}

//Not really about randomness, but quite related.  Turns log probabilities into real probabilities
void logs_to_probs(double *q, int L)
{
  double maxq, sumq;
  int l;
  maxq = q[0];
  for(l = 0; l < L; l++) if(q[l] > maxq) maxq = q[l];
  sumq = 0;
  for(l = 0; l < L; l++) sumq += exp(q[l] - maxq);
  for(l = 0; l < L; l++) q[l] = exp(q[l] - maxq) / sumq;
  return;
}

void logs_to_probs_ind(double *q, int L, int* ones)
{
  double maxq, sumq;
  int l;
  maxq = R_NegInf;
  for(l = 0; l < L; l++)
    {
      if(ones[l])
	{
	  if(q[l] > maxq) maxq = q[l];
	}
    }
  //  printf("maxq: %f\n",maxq);
  sumq = 0;
  for(l = 0; l < L; l++)
    {
      if(ones[l])
	{
	  sumq += exp(q[l] - maxq);
	}
    }

  for(l = 0; l < L; l++)
    {
      if(ones[l])
	{
	  q[l] = exp(q[l] - maxq) / sumq;
	}
    }
  return;
}

void weights_to_probs(double *q, int L)
{
  int l;
  double sumq = 0.0;
  for(l = 0; l < L;l++)sumq += q[l];
  for(l = 0; l < L; l++) q[l] = q[l] / sumq;
  return;
}
