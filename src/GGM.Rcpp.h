#include "GGM.h"

//The constructor for the GGM Class
// v is the number of variables
// var is the variance associated with updates of K and G.
GGM_CLASS::GGM_CLASS(int v)
{
  
  int i, k, l;
  int temp;
  //  printf("Making new graph, var is: %f\n",var);
  //------ General Model Info ------
  p = v;//Give the dimension
  n = 0;//Right now the model has no data
  delta = 3;//A magic number!
  //---------------------------------

  //---- Propose a graph from the prior ---
  graph = new Graph;
  graph->InitGraph(p);
  for(k = 0; k < p - 1; k++)
    {
      for(l = k + 1; l < p; l++)
	{
	  temp = (unif_rand() < .5);
	  graph->Edge[k][l] = temp;
	  graph->Edge[l][k] = temp;
	}
    }
  TurnFillInGraph(graph);
  graph->GetMPSubgraphs();//Important to remember this!
  //---------------------------------------

  //---- Find sufficient statistics -------
  DplusU = new Mat(p, p);
  DplusU->set_iden();
  //----------------------------------------

  //------ Initialize parameter K ----------
  K = new Mat(p,p);//gwish_blgibbs_decomposable(graph, 3, DplusU);
  K->set_iden();
  //----------------------------------------

  //------- Initialize parameter mu --------
  mu = new double[p];
  for(i = 0; i < p; i++)mu[i] = 0.0;
  //----------------------------------------

  //------- xbar ---------------------------
  xbar = new double[p];
  for(i = 0; i < p; i++) xbar[i] = 0.0;
  //----------------------------------------

}

//Nothing special here, just kill it all
GGM_CLASS::~GGM_CLASS()
{
  delete graph;
  delete K;
  delete DplusU;
  delete[] mu;
  delete[] xbar;
}

// This adds observations to the current GGM. X is the full data observation, k is the row we're adding.
void GGM_CLASS::AddObs(Matrix X, int k)
{
  double a;
  int i,j;
  double* xbar_old = new double[p];
  double ss;
  //-------- Update DplusU, the outer product matrix -----------
  for(i= 0; i < p; i++)
    {
      xbar_old[i] = xbar[i];
      ss = n * xbar[i] + X->get(k,i);
      //      printf("%f\n", ss);
      xbar[i] = ss / (n + 1);
      //      printf("%d %f %f %f\n", i, X->get(k,i), xbar_old[i], xbar[i]);
    }

  for(i = 0; i < p; i++)
    {
      for(j = 0; j < p; j++)
	{
	  a = DplusU->get(i,j);
	  a += n * (xbar_old[i] * xbar_old[j]);
	  a -= (n + 1) * xbar[i] * xbar[j];
	  a += X->get(k,i) * X->get(k,j);
	  DplusU->set(i,j,a);
	}
    }
  //-------------------------------------------------------
  n++;
  delete[] xbar_old;
  return;
}

// This adds observations to the current GGM.
void GGM_CLASS::DropObs(Matrix X, int k)
{
  double a;
  int i,j;
  double* xbar_old = new double[p];
  double ss;
  //-------- Update DplusU, the outer product matrix -----------
  if(n == 1)
    {
      for(i = 0; i < p; i++)
	{
	  xbar_old[i] = xbar[i];
	  xbar[i] = 0.0;
	}
    }
  else
    {
      for(i= 0; i < p; i++)
	{
	  xbar_old[i] = xbar[i];
	  ss = n * xbar[i] - X->get(k,i);
	  xbar[i] = ss / (n - 1);
	}
    }

  for(i = 0; i < p; i++)
    {
      for(j = 0; j < p; j++)
	{
	  a = DplusU->get(i,j);
	  a += n * xbar_old[i] * xbar_old[j];
	  a -= (n - 1) * xbar[i] * xbar[j];
	  a -= X->get(k,i) * X->get(k,j);
	  DplusU->set(i,j,a);
	}
    }
  //-------------------------------------------------------
  n--;
  delete[] xbar_old;
  return;
}

void GGM_CLASS::UpdateG()
{
  //----------  Some set up ----------------------
  int n_dec, n_dec_new;
  int which_change;
  int ee = p * (p - 1) / 2;
  int *which_dec = new int[ee];
  LPGraph newgraph, tempgraph;
  Matrix D_prior = new Mat(p,p);
  D_prior->set_iden();
  double alpha,r;
  double numerator, denominator;
  //----------------------------------------------

  //---------------  Propose a Neighbor ---------------------
  n_dec = FindDecomposableNeighbors(graph, which_dec);
  which_change = sample_from(which_dec,ee,n_dec);
  newgraph = new Graph(graph);
  FlipEdge(newgraph, which_change);
 /* if(!newgraph->IsDecomposable())
    {
       printf("HELP\n");
    }
*/
  n_dec_new = FindDecomposableNeighbors(newgraph, which_dec);//Necessary for the M-H ratio
  //-----------------------------------------------------------

  //-------------- Compute a ratio ----------------------------
  numerator = -log( (double) n_dec_new) +  j_g_decomposable(newgraph, D_prior, DplusU, delta, n);
  denominator = -log( (double) n_dec) + j_g_decomposable(graph, D_prior, DplusU, delta, n);
  alpha = numerator - denominator;
  r = log(unif_rand());
  if(r < alpha)
    {
      //	  printf("r %f was less than %f\n", r, alpha);
      tempgraph = graph;
      graph = newgraph;
      delete tempgraph;
    }
  else
    {
      delete newgraph;
    }
  delete D_prior;
  delete[] which_dec;

}

void GGM_CLASS::SampleParams()
{
  
  //-----  Set-up -------------
  int i;
  //  double a;
  double n0 = .01;
  double *mu0 = new double[p];
  double *mu_bar = new double[p];
  for(i = 0; i < p; i++) mu0[i] = 0;
  Matrix Ktemp;
  //------------------------------

  //--------- Forget me ----------
  delete K;
  delete[] mu;
  //------------------------------

  //----- Predict the cluster that the next observation will land in -------------
  for(i = 0; i < p; i++) mu_bar[i] = (n * xbar[i] + n0 * mu0[i]) / (n + n0);
  if(graph->nCliques == 1)
    {
      K = rwish(3 + n, DplusU);
    }
  else
    {
      K = gwish_blgibbs_decomposable(graph, 3 + n, DplusU);
    }
  Ktemp = K->copy();
  Ktemp->mult_const(n + n0);
  //  printf("Sampling a Mult Norm, Mean\n");
  //  util_print_vec_dbl(mu_bar,p);
  //  printf("Precision\n");
  //  Ktemp->print();
  mu = sample_mult_norm(mu_bar, Ktemp);
  //  printf("Result\n");
  //  util_print_vec_dbl(mu,p);
  //--------------------------------------------------------------------------------

  delete Ktemp;

}

void GGM_CLASS::PrintInfo()
{
  int i,j,k;
  k = 0;
  // printf("GGM with %d Observations\n", n);
  for(i = 0; i < p -1; i++)
    {
      for(j = i + 1; j < p; j++)
	{
	  k += graph->Edge[i][j];
	}
    }
  // printf("Current Edges: %d\n", k);
  if(p < 26)
    {
      // printf("Graph\n");
      graph->PrintA();
    }
  if(p < 11)
    {
      // printf("Precision\n");
      K->print();
      // printf("Mu\n");
      util_print_vec_dbl(mu,p);
    }

}


// Takes a graph and fills the vector "which" with an indication of whether a particular
// edge-away neighbor is decomposable.  Returns total number of decomposable neighbors
int FindDecomposableNeighbors(LPGraph graph, int *which)
{
  int p = graph->nVertices;
  int ee = p * (p - 1) / 2;
  int i;
  LPGraph temp;
  int num = 0;

  temp = new Graph(graph);
  for(i = 0; i < ee; i++)which[i] = 0;

  for(i = 0; i < ee; i++)
    {
      FlipEdge(temp, i);
      if(temp->IsDecomposable())
	{
	  which[i] = 1;
	  num++;
	}
      FlipEdge(temp, i);
    }
  delete temp;
  return(num);
}

//Flips an edge in a graph, where which is a number from 0 to the number of edges in the graph (minus 1)
void FlipEdge(LPGraph graph, int which)
{
  int p = graph->nVertices;
  int i,j,k = 0;

  for(i = 0; i < p - 1; i++)
    {
      for(j = i + 1; j < p; j++)
	{
	  if(k == which)
	    {
	      graph->Edge[i][j] = 1 - graph->Edge[i][j];
	      graph->Edge[j][i] = 1 - graph->Edge[j][i];
	      return;
	    }
	  k++;
	}
    }
}
