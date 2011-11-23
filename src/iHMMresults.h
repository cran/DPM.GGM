typedef class iHMMResults* Results;


class iHMMResults{

 public:
  int *cluster_mat;
  int *edge_mat;
  int predict;
  int all_k;
  int tot_reps;
  double *alpha;
  double *alpha0;
  int *cluster_total;
  double *stepahead_x;
  double *stepahead_mu;
  double *stepahead_K;
  double *Kall;
  iHMMResults(int *cluster_mat_R,int *edge_mat_R,double *alpha_R,
	      double *alpha0_R, int *cluster_total_R, int *predict_R, 
	      double *stepahead_R_x, double *stepahead_R_mu, double *stepahead_R_K, int *all_k_R, int *tot_reps, double *Kall_R);
  ~iHMMResults();
    
  void UpdateResults(int rep, State state);


};
