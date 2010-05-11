#define nx 4096
#define ny 4096

#define oneau 1.49597870692e11

struct rotpar{
  int rotbad;
  int rotpairs;
  float rotcadence;
  int flatfield_version;
};

struct code_param{
  double convergence;
  int maxiter;
  double omega;
  double norm;
  double croprad;
  double rotcoef0;
  double rotcoef1;
  double rotcoef2;
};

struct list{
  int val;
  struct list *next;
};




