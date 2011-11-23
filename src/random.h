double* sample_mult_norm(double *mu, Matrix K);
int sample_from(int *s, int n, int n_sub);
int rand_int(int);
void rdirichlet(double*, int, double*);
int rand_int_weighted(int, double*);
void logs_to_probs(double*, int);
void logs_to_probs_ind(double*, int, int*);
void weights_to_probs(double*, int);
