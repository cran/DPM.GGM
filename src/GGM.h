//A Class for the Gaussian Graphical Model
typedef class GGM_CLASS* GGM;


class GGM_CLASS{
 public:

  
  //------ Data ------------------------------------------
  Matrix DplusU;//Our Posterior Parameter
  double* xbar;
  int n,p;
  int delta;
  //------------------------------------------------------

  //------ These are the actual parameters ---
  LPGraph graph;
  Matrix K;
  double* mu;
  //----------------------------------------------

  //-------- Functions -----------------
  GGM_CLASS(int);
  ~GGM_CLASS();
  void AddObs(Matrix, int);
  void DropObs(Matrix, int);
  void UpdateG();
  void SampleParams();
  void PrintInfo();
  //-------------------------------------

};


// Takes a graph and fills the vector "which" with an indication of whether a particular
// edge-away neighbor is decomposable.  Returns total number of decomposable neighbors
int FindDecomposableNeighbors(LPGraph graph, int *which);
void FlipEdge(LPGraph, int);

