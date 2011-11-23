typedef class DPMSliceState_CLASS* DPMSliceState;

class DPMSliceState_CLASS{

  //-------- Variables ---------------
 public:
  int n;
  int p;
  Matrix X;
  int *xi;
  int *r;
  int N;
  GGM* GGMlist;
  double *w;
  double *u;
  double alpha;
  //----------------------------------

  //--------- Functions --------------
  DPMSliceState_CLASS(double*, int,int,int);
  ~DPMSliceState_CLASS();
  void UpdateState();
  void UpdateAlpha();
  void UpdateW();
  void UpdateU();
  void UpdateXi();
  void UpdateGGM();
  void ExpandModel();
  double PredictiveDistribution(int, GGM);
  //----------------------------------

};
