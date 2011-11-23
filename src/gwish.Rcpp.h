//--------------------------------------------------------
// gwish.Rcpp:  This is a collection of functions that manipulate graph objects
//   according to the G-Wishart distribution.
//   Note: this depends on the newgraph.Rcpp library, and the matrix.Rcpp library
//
//-----------------------------------------------------------

//Alex Lenkoski  lenkoski@stat.washington.edu

// Samples a multivariate normal vector
#include "gwish.h"

// Calculates the normalizing constant for a G-Wishart distribution
// when the graph is full.
double gwish_nc_complete(int delta, int p, Matrix D)
{
  //  double I;
  double d,c,a,g;
  double dblDelta;//Recasting the inputs makes life easier below
  double dblP;
  int i;

  dblDelta = delta;
  dblP = p;

  //----    Calculate Each Piece --------
  //  util_print_mat_dbl(D, p, p, 0);
  d = D->log_det();
  d = ( (dblDelta + dblP - 1) / 2) * d;

  c = (dblP * (dblDelta + dblP - 1)) / 2 * M_LN2;

  a = (dblDelta + dblP - 1) / 2;
  g = dblP * (dblP - 1) / 4 * log(M_PI);
  for(i = 0; i < dblP; i++) g += lgamma(a - (double)i / 2);
  //  printf("d: %f, c: %f, a: %f, g: %f\n", d,c,a,g);
  return(-d + c + g);
  //------------------------------------
}

/*
//Makes a covariance matrix over a subset of the data
void make_sub_cov(Matrix X, int *sub, int sub_match, int n_sub, Matrix D)
{
  int i,j,k;

  D->set_zero();

  if(n_sub < 2)return;

  double *means = new double[p];
  for(i = 0; i < p; i++) means[i] = 0;

  for(i = 0; i < p; i++)
    {
      for(k = 0; k < n;k++)
	{
	  if(sub[k] == sub_match) means[i] += X->get(k,i) / n_sub;
	}
    }

  for(i = 0; i < p; i ++)
    {
      for(j = 0; j < p; j++)
	{
	  for(k = 0; k < n; k++)
	    {
	      if(sub[k] == sub_match)
		{
		  temp = (X->get(k,i) - means[i]) * (X->get(k,j) - means[j]);
		  D->set(i,j,D->get(i,j) + temp);
		}
	    }
	}
    }

  delete[] means;
  return;
}

//Makes a mean v3ector over a subset of the data
void make_sub_means(Matrix X, int *sub, int sub_match, int n_sub, double *means)
{
  int i,j,k;


  for(i = 0; i < p; i++) means[i] = 0;
  if(n_sub == 0) return;

  for(i = 0; i < p; i++)
    {
      for(k = 0; k < n; k++)
	{
	  if(sub[k] == sub_match) means[i] += X->get(k,i) / n_sub;
	}
    }

}
*/

//Calculates the ratio of normalizing constants when the graph G is decomposable
double j_g_decomposable(LPGraph graph, Matrix D_prior, Matrix D_post, int delta, int n)
{
  double mypost = 0;
  int p = graph->nVertices;
  int i, j;//, p_i, p_j;
  Matrix sub_D_post;
  Matrix sub_D_prior;
  int sub_p;
  //  int sub_ee;
  int *clique = new int[p];

  //  util_print_mat_int(CliqueMat, *p,*p,0);
  for(i = 0; i < graph->nCliques; i++)
    {
      sub_p = graph->CliquesDimens[i];
      for(j = 0; j < p; j++) clique[j] = graph->Cliques[i][j];

      //----------  Form Sub-Matrices  -----------
      sub_D_post = D_post->sub_mat(sub_p, clique);
      sub_D_prior = D_prior->sub_mat(sub_p,clique);
      //-------------------------------------------

      mypost += gwish_nc_complete(delta + n, sub_p, sub_D_post);
      mypost -= gwish_nc_complete(delta, sub_p, sub_D_prior);

      delete sub_D_post;
      delete sub_D_prior;
    }//Loop through the prime components

 for(i = 0; i < graph->nSeparators; i++)
    {

      sub_p = graph->SeparatorsDimens[i];
      for(j = 0; j < p; j++) clique[j] = graph->Separators[i][j];

      //----------  Form Sub-Matrices  -----------
      sub_D_post = D_post->sub_mat(sub_p, clique);
      sub_D_prior = D_prior->sub_mat(sub_p,clique);
      //-------------------------------------------
 
      mypost -= gwish_nc_complete(delta + n, sub_p, sub_D_post);
      mypost += gwish_nc_complete(delta, sub_p, sub_D_prior);

      delete sub_D_post;
      delete sub_D_prior;
      
    }//Loop through the separators

 delete[] clique;
 return(mypost);

}//j_g_decomposable

//-------- Samples a Wishart Variate (Implies Full Graph) ---------
// Note that I use the Barlett decomposition to achieve this
Matrix rwish(int delta, Matrix D)
{
  int i,j;
  int p = D->p_r;

  Matrix D_inv = D->invert();
  Matrix T = D_inv->chol();
  Matrix psi = new Mat(p,p);

  // ***  Sample Values in Psi   ******

  //------Diagaonal Elements----------
  for (i = 0; i < p; i++) psi->set(i, i, sqrt( rchisq(p - i - 1 + delta) ) );
  //----------------------------------
  
  //------Offdiagonal in G------------
  for (i = 0; i < p-1; i++)
    {
      for(j = i + 1; j < p; j++)
	{
	  psi->set(i, j, rnorm(0,1) );
	}
    }

  //----------------------------------
  
  //******** End Sampling  *************

  //----  Now Complete things  ----------
  Matrix psiT = mult_mats(psi,T);
  Matrix K = square_mat(psiT);
  //-------------------------------------

  //----- And Clean up ------------------
  delete D_inv;
  delete T;
  delete psi;
  delete psiT;
  //-------------------------------------

  return(K);
}

//Performs one iteration of the Block Gibbs Sampler
void gwish_blgibbs_iterate(LPGraph graph, int delta, Matrix D, Matrix K)
{
  int p = graph->nVertices;
  int i, p_clique;
  int *clique_ID, *V_ID;
  int clique_num;
  Matrix submatrix_D;
  Matrix submatrix_cc;
  Matrix submatrix_cond;

  //------------ Loop through cliques and update accordingly -------
  for(clique_num = 0; clique_num < graph->nCliques; clique_num++)
    {

      //*******Just recording Indices***********
      p_clique = graph->CliquesDimens[clique_num];
      clique_ID = new int[p_clique];
      V_ID = new int[p - p_clique];
      for(i = 0; i < p_clique; i++) clique_ID[i] = graph->Cliques[clique_num][i];
      get_complementary_set(p, p_clique, clique_ID, V_ID);
      //******************************************
      
      //******  Start: Making Submatrices  ************
	  
      //---First make the part specific to the clique-----
      submatrix_D = D->sub_mat(p_clique, clique_ID);
      submatrix_cc = rwish(delta, submatrix_D);
      //--------------------------------------------------
	  
      //-----   Now form the part that's based ------
      //----- on the entries outside of the clique --
      if(p_clique < p)
	{
	  submatrix_cond = get_cond_matrix(p_clique, clique_ID, V_ID, K);
	  submatrix_cc->Add(submatrix_cond);
	  delete submatrix_cond;
	}
      //---------------------------------------------
	  
      //********  End Forming Matrices ********************
	  
      //---Update the result matrix for the clique---
      K->set_sub_mat(p_clique, clique_ID, submatrix_cc);
      //---------------------------------------------

      //----- Clean up from updating this clique-----
      delete[] clique_ID; delete[] V_ID;
      delete submatrix_cc; delete submatrix_D;
      //----------------------------------------------
      
    }//End Loop through cliques
  //--We've now looped through all cliques in the graph and updated---

  return;  
}

// Samples from a G-Wishart for a decomposable graph. Note that due to Piccioni (2000)
// if you run the block Gibbs sampler twice over a decomposable graph, you get a sample
// from the posterior.
Matrix gwish_blgibbs_decomposable(LPGraph graph, int delta, Matrix D)
{

  //--------------- Set mat to identity and then run twice-------------------
  int p = graph->nVertices;
  Matrix K = new Mat(p,p);
  K->set_iden();
  gwish_blgibbs_iterate(graph, delta, D, K);
  gwish_blgibbs_iterate(graph, delta, D, K);
  //-------------------------------------------------------------------------
  return(K);
}//gwish_blgibbs_decomposable


