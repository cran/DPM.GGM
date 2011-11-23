typedef class iHMMState* State;

class iHMMState{

  //-------- Variables ---------------
 public:
  int n;
  int p;
  Matrix X;
  int *xi;
  int L;
  GGM* GGMlist;
  double alpha;
  double alpha0;
  double *gamma;
  
  //----------------------------------

  //--------- Functions --------------
  iHMMState(double*, int,int,int);
  ~iHMMState();
  void PrintX();
  void UpdateState();
  void UpdateAlpha();
  void UpdateGamma();
  void UpdateXi();
  void UpdateG();
  double PredictiveDistribution(int, int, GGM);
  void SampleNext(double*, double*, double*);
  void SampleAllK(double*);
  void Print();
 //----------------------------------

};

