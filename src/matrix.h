//A C++ matrix class that plays nice with R's blas and lapack libraries

//Alex Lenkoski <alex.lenkoski@gmail.com>

typedef class Mat* Matrix;
typedef class intMat* intMatrix;

class Mat{
 public:
  int p_r;
  int p_c;
  double *Data;
  
  Mat(int, int);
  Mat(int, int, double*);
  ~Mat();
  double get(int,int);
  void set(int,int,double);
  void set_iden();
  Matrix copy();
  Matrix invert();
  void print();
  double log_det();
  Matrix chol();
  Matrix transpose();
  Matrix sub_mat(int, int*);
  void set_sub_mat(int, int*, Matrix);
  void Add(Matrix);
  void mult_const(double);
};

class intMat{
 public:
  int p_r;
  int p_c;
  int* Data;

  intMat(int, int);
  intMat(int, int, int*);
  ~intMat();
  int get(int, int);
  void set(int, int, int);
  void print();

};
Matrix mult_square_mats(Matrix, Matrix);
Matrix mult_mats(Matrix, Matrix);
void get_complementary_set(int, int, int*,int*);
Matrix get_cond_matrix(int, int*, Matrix);
Matrix square_mat(Matrix);
double quad_form(Matrix A ,Matrix B); //returns ABA' when A 1 x p
