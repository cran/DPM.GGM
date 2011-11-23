typedef class DPMState_CLASS* DPMState;

class DPMState_CLASS{

  //-------- Variables ---------------
 public:
  int n;
  int p;
  Matrix X;
  int *xi;
  int L;
  GGM* GGMlist;
  double alpha;
  //----------------------------------

  //--------- Functions --------------
  DPMState_CLASS(double*, int,int,int);
  ~DPMState_CLASS();
  void PrintX();
  void UpdateState();
  void UpdateAlpha();
  void UpdateXi();
  void UpdateG();
  double PredictiveDistribution(int, int, GGM);
  //----------------------------------

};
