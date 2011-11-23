double gwish_nc_complete(int,int, Matrix);
//void make_sub_means(Matrix, int*, int, int, double*);
double j_g_decomposable(LPGraph, Matrix, Matrix, int, int);
Matrix rwish(int delta, Matrix D);
void gwish_blgibbs_iterate(LPGraph graph, int delta, Matrix D, Matrix K);
Matrix gwish_blgibbs_decomposable(LPGraph graph, int delta, Matrix D);

