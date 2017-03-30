#include <math.h>

#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "tabix.h"

#include "Argument.h"
#include "Utils.h"
#include "VCFUtil.h"
#include "SimpleMatrix.h"
#include "IO.h"
#include "PlinkInputFile.h"
#include "Time.h"

/**
 *@return true when both @param marker1 and @param marker2 are the same.
 */
bool isSameVector(const std::vector< std::string >& marker1,
                  const std::vector< std::string >& marker2) {
  if (marker1.size() != marker2.size())
    return false;
  for (size_t i = 0; i < marker1.size(); ++i) {
    if (marker1[i] != marker2[i]) {
      return false;
    }
  }
  return true;
}

class Coef{
 public:
  float a1;
  float a2;
  float b;
};

class SPAModel{
 public:
  SPAModel(const std::string& fn) {
    LineReader lr(fn);
    std::vector<std::string> fd;
    std::string key;
    Coef value;
    std::string line;
    while(lr.readLine(&line)) {
      stringNaturalTokenize(line, "\t ", &fd);
      key = fd[0] + ":" + fd[3];
      minor[key] = fd[4];
      major[key] = fd[5];
      value.a1 = atof(fd[6]);
      value.a2 = atof(fd[7]);
      value.b  = atof(fd[8]);
      if (this->data.find(key) != this->data.end()) {
        fprintf(stderr, "Duplicated marker: %s\n", key.c_str() );
        continue;
      }
      this->data[key] = value;
      this->coef.push_back(value);
      this->markers.push_back(key);
    }
  }
  Coef& getData(const std::string& marker){
    if (data.find(marker) == data.end()) {
      fprintf(stderr, "%s:%d Non-exist marker: %s\nf", __FILE__, __LINE__, marker.c_str() );
    }
    return data[marker];
  }
  const Coef& operator [] (int idx) const {
    return coef[idx];
  }
  const std::string& getMinor(const std::string& markerName) {
    return minor[markerName];
  }
  const std::string& getMajor(const std::string& markerName) {
    return major[markerName];
  }
  size_t size() const {
    return this->coef.size();
  }
  const std::vector< std::string >& getMarkers() const{
    return this->markers;
  }
 private:
  std::map < std::string, Coef> data;
  std::vector< Coef > coef;
  std::vector< std::string > markers;
  std::map< std::string, std::string > minor;
  std::map< std::string, std::string > major;
};

/**
 * e.g. site file
 *  CHR     POS     ID      REF     ALT
 *  1       752566  rs3094315       G       A
 */
class SiteFile{
 public:
  SiteFile(const std::string& fn){
    LineReader lr(fn);
    std::vector<std::string> fd;
    std::string key;
    std::string line;
    while(lr.readLine(&line)) {
      stringNaturalTokenize(line, "\t ", &fd);
      if (fd[0] == "CHR") continue; // skip header
      key = fd[0] + ":" + fd[1];

      if (ref.count(key)) {
        fprintf(stderr, "Duplicated marker: %s\n", key.c_str() );
        continue;
      }
      ref[key] = fd[3];
      alt[key] = fd[4];
      markers.push_back(key);
    }
  }
  const std::string& getRef(const std::string& key) {
    if (ref.count(key) == 0) {
      fprintf(stderr, "%s:%d Non-exist marker: %s\n", __FILE__, __LINE__, key.c_str() );
    }
    return ref[key];
  }
  const std::string& getAlt(const std::string& key) {
    if (alt.count(key) == 0) {
      fprintf(stderr, "Non-exist marker: %s\n", key.c_str() );
    }
    return alt[key];
  }
  size_t size() const {
    return this->ref.size();
  }
  const std::vector< std::string >& getMarkers() const{
    return this->markers;
  }
 private:
  std::map<std::string, std::string> ref; // store marker->ref
  std::map<std::string, std::string> alt; // store marker->alt
  std::vector<std::string> markers;
};

class Param {
 public:
  SimpleMatrix* geno;
  SPAModel* model;
  double error;
  double extra[2]; // store x and f(x)
};

class RootFinder{
 public:
  RootFinder(gsl_function* func): F(*func) {
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    max_iter = 100;
  }
  virtual ~RootFinder() {
    gsl_root_fsolver_free (s);
  }
  /**
   * find a bracket @param out, such that (x, out) contains a rootj
   * @param x, start searching point
   * @param out, result
   * @param scale, scale > 0, search higher bracket; scale < 0, otherwise
   * @return 0 if succeed
   */
  int searchBracket(double x, double* out, double scale) {
    scale = scale > 0.0 ? 1.0 : -1.0;
    double l = GSL_FN_EVAL(&F,x);
    double step = 1.0;
    double& newX = *out;
    newX = x + step * scale;
    double newL = GSL_FN_EVAL(&F, newX);
    int time = 0;
    while ( (l>0 && newL>0) || (l<=0 && newL <= 0)) {
      step *= 2.0;
      newX = x + step * scale;
      newL = GSL_FN_EVAL(&F, newX);
      if (++time >= 10) {
        return -1;
      }
    }
    return  0;
  }
int calculate(double x_lo, double x_hi,
              double* out) {
  // F = func;
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

    int status;
    int iter = 0;
    double r;
    do {
      iter ++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 0.001);

      if (status == GSL_SUCCESS) {
        *out = 0.5 * (x_lo + x_hi);
#ifdef DEBUG
        fprintf(stderr, "solve at %g ( %g, %g )\n", *out, x_lo, x_hi);
#endif
        return 0;
      }
    } while (status == GSL_CONTINUE && iter < max_iter);
    return -1;
  }
 private:
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function& F;
  // double x_lo;
  // double x_hi;
  int max_iter;
};

/**
 * @return likelihood
 */
double computeLikelihood(const double x,
                         const double y,
                         void* params) {
  Param* p = (Param*) params;
  SimpleMatrix& geno = * (p->geno);
  SPAModel& model = * (p->model);
  double error = (p->error);
  const int numMarker = geno.ncol() / 2;

  double llk = 0;
  double l;
  for (int j = 0; j < numMarker; ++j) {
    if (geno[0][j] < 0) continue;  // skip missing genotypes
    const Coef& coef = model[j];
    const int& depth = geno[0][j * 2];
    if (depth == 0) continue;
    const int& refCount = geno[0][j * 2 + 1];
    assert(depth>= refCount);
    const int altCount = depth - refCount;
    l = 0.0;
    double f = 1.0 / (exp(-coef.a1 * x - coef.a2 * y - coef.b ) + 1.);
    l += pow(1.0 - error, refCount) * pow(error, altCount) * f * f;
    l += pow(0.5, depth) * 2.0 * f * (1.0 - f);
    l += pow(error, refCount) * pow(1.0 - error, altCount) * (1.0 - f) * (1.0 - f);
    l = log(l);
    if (std::isfinite(l)) {
      llk += l;
    } else {
      continue;
    }
  }
  return llk;
}

double computeLikelihood(const gsl_vector *v, void *params) {
  double x, y;
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
  
  return computeLikelihood(x, y, params);
}

double computeMinusLikelihood(const gsl_vector *v, void *params) {
  return -computeLikelihood(v, params);
}

double computeLikelihoodFixX(double y, void* params) {
  double x = ((Param*) params) -> extra[0];
  return computeLikelihood(x, y, params) - ((Param*) params) -> extra[1];
}

double computeLikelihoodFixY(double x, void* params) {
  double y = ((Param*) params) -> extra[0];
  return computeLikelihood(x, y, params) - ((Param*) params) -> extra[1];
}

/**
 * @return 0: if location are successfully inferred
 */
int resample(const SimpleMatrix& m,
             const SPAModel& model,
             const double error,
             const double x, const double y,
             const gsl_rng * r,
             SimpleMatrix* genoOut) {
  SimpleMatrix& out = *genoOut;
  out.resize(m.nrow(), m.ncol());
  const int numMarker = m.ncol() / 2;
  int minorAlleleCount;
  double p = 1.0 - error;
  for (int i = 0; i < numMarker; ++i ) {
    if (m[0][i*2] == 0) {
      out[0][i*2] = 0;
      out[0][i*2 + 1] = 0;
      continue;      
    }
    const int& depth = m[0][i*2];
    const Coef& coef = model[i];
    float p = 1.0 / (exp(-coef.a1 * x - coef.a2 * y - coef.b ) + 1.);
    minorAlleleCount = gsl_ran_binomial(r, p, depth);
    out[0][i*2] = depth;
    out[0][i*2 + 1] = minorAlleleCount;
  }
  
  return 0;
}
/**
 * @return 0: if location are successfully inferred
 */
int inferLocationSeq(SimpleMatrix& m, SPAModel& model,
                     double error,
                     double* x, double* y) {
  printf("m dim is %d x %d, and model has %zu markers\n", m.nrow(), m.ncol(), model.size());

  // initialize start points
  gsl_vector* init = gsl_vector_alloc (2);
  gsl_vector_set (init, 0, 0.0);
  gsl_vector_set (init, 1, 0.0);

  // step size
  gsl_vector* ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1.0);

  // init parameter
  Param param;
  param.geno = &m;
  param.model = &model;
  param.error = error;
  
  // initialize functions
  const gsl_multimin_fminimizer_type *T =
      gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;

  gsl_multimin_function minex_func;
  minex_func.n = 2;
  minex_func.f = computeMinusLikelihood;
  minex_func.params = &param;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, init, ss);

  // main loop
  size_t iter = 0;
  int status;
  double size;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status)
      break;
    
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1e-4);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }

    printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
            (int)iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            s->fval, size);
  }  while (status == GSL_CONTINUE && iter < 200);

  // assign optimization result
  if (status == GSL_SUCCESS) {
    *x = gsl_vector_get (s->x, 0);
    *y = gsl_vector_get (s->x, 1);
  } else {
    *x = *y = 0.0;
  }

  // free memory
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (init);

  return (status == GSL_SUCCESS ? 0 : 1);
}

/**
 * Calculate confidence interval (95% )
 */
int inferLocationSeqCI(SimpleMatrix& m, SPAModel& model,
                       double error,
                       const double x, const double y,  // mle estimator
                       double* xLow, double* xHigh,
                       double* yLow, double* yHigh) {
  // init parameter
  Param param;
  param.geno = &m;
  param.model = &model;
  param.error = error;

  *xLow = *xHigh = *yLow = *yHigh = 0;
  double llkFull = computeLikelihood(x, y, &param);
  double llkNull = llkFull - 0.5 * 5.991465; // 0.5 * qchisq(0.95, df = 2);

  // find xLow
  double llk = llkFull;
  double xTmp = x;
  double yTmp = y;
  double step = 10.0;
  int time = 0;
  
  gsl_function func;
  RootFinder f(&func);

  // find CI for x
  param.extra[0] = y;
  param.extra[1] = llkNull;
  func.function = computeLikelihoodFixY;
  func.params = &param;
  if (f.searchBracket(x, &xTmp, -1)) {
    fprintf(stderr, "cannot bracket!\n");
    return -1;
  }
  if (f.calculate(xTmp, x, xLow)) {
    fprintf(stderr, "cannot solve!\n");    
    return -1;
  }
  if (f.searchBracket(x, &xTmp, 1)) {
    fprintf(stderr, "cannot bracket!\n");
    return -1;    
  }
  if (f.calculate(x, xTmp, xHigh)) {
    fprintf(stderr, "cannot solve!\n");    
    return -1;
  }
  // find CI for y
  param.extra[0] = x;
  param.extra[1] = llkNull;
  func.function = computeLikelihoodFixX;
  func.params = &param;
  if (f.searchBracket(y, &yTmp, -1)) {
    fprintf(stderr, "cannot bracket!\n");
    return -1;
  }
  if (f.calculate(yTmp, y, yLow)) {
    fprintf(stderr, "cannot solve!\n");    
    return -1;
  }
  if (f.searchBracket(y, &yTmp, 1)) {
    fprintf(stderr, "cannot bracket!\n");
    return -1;    
  }
  if (f.calculate(y, yTmp, yHigh)) {
    fprintf(stderr, "cannot solve!\n");    
    return -1;
  }
  
  return 0;
}

int main(int argc, char** argv){
  time_t startTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&startTime));

  // set up GSL
  const gsl_rng_type * T;
  gsl_rng * r;

  int i, n = 10;
  double mu = 3.0;

  /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  
  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inSeq, "--inSeq", "input .geno file")
      ADD_STRING_PARAMETER(pl, inSite, "--inSite", "input .site file")
      ADD_STRING_PARAMETER(pl, inModel, "--inModel", "input SPA format output")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
      ADD_INT_PARAMETER(pl, bootstrap, "--bootstrap", "number of bootstrapped samples")
      ADD_BOOL_PARAMETER(pl, ci, "--ci", "calculate 95% confidence interval")            
      // ADD_PARAMETER_GROUP(pl, "Others")
      // ADD_DOUBLE_PARAMETER(pl, error, "--error", "error rate")
      // ADD_DOUBLE_PARAMETER(pl, depth, "--depth", "average depth")
      // ADD_INT_PARAMETER(pl, seed, "--seed", "average depth")
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);
  pl.Status();

  if (FLAG_REMAIN_ARG.size() > 0){
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    fprintf(stderr, "\n");
    abort();
  }

  REQUIRE_STRING_PARAMETER(FLAG_inSeq, "Please provide input file using: --inSeq");
  REQUIRE_STRING_PARAMETER(FLAG_inSite, "Please provide input file using: --inSite");
  REQUIRE_STRING_PARAMETER(FLAG_inModel, "Please provide input file using: --inModel");

  SPAModel spaModel(FLAG_inModel);
  LineReader genoReader(FLAG_inSeq);
  SiteFile siteFile(FLAG_inSite);

  FILE* fout = fopen((FLAG_outPrefix + ".loc").c_str(), "wt");

  // sanity check
  // 1. all sites in spa model outputs should be in siteFile
  // 2. the index of site file
  // record spa model order
  std::map<std::string, int> ord;
  for (int i = 0; i < spaModel.getMarkers().size(); ++i) {
    ord[spaModel.getMarkers()[i]] = i;
  }
  std::vector<int> siteIndex(siteFile.getMarkers().size(), -1);
  for (int i = 0; i < siteFile.getMarkers().size(); ++i) {
    const std::string& key = siteFile.getMarkers()[i];
    if (ord.count(key))  {
      if (siteIndex[i] >= 0) {
        fprintf(stderr, "duplicated marker [ %s ] ??\n", key.c_str());
      }else{
        siteIndex[i] = ord[key];
      }
    } else {
      continue;
    }
  }

  // if (!isSameVector(siteFile.getMarkers(), spaModel.getMarkers())) {
  //   fprintf(stderr, "Site file and marker does not match on rsid.\n");
  //   exit(1);
  // }
  int numMarker = spaModel.getMarkers().size();
  fprintf(stderr, "loaded %d markers\n", numMarker);
  
  // loop each marker in the plink file
  // const int numMarker = .getNumMarker();
  std::string key;
  fprintf(stderr, "PopId\tIndvId\tLoc1\tLoc2\n");
  fprintf(fout,   "PopId\tIndvId\tLoc1\tLoc2\n");

  SimpleMatrix m;
  std::vector<std::string> fd;
  int lineNo = 0;
  while(genoReader.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    const std::string& popId = fd[0];
    const std::string& indvId = fd[1];
    // const int numMarker = (fd.size() - 2) / 2;
    if (siteFile.size() * 2 + 2 != fd.size()) {
       fprintf(stderr, "Skip line %d due to incorrect column number!\n", lineNo);
       continue;
    }
    m.resize(1, numMarker * 2);
    for (int i = 0; i < siteFile.size(); ++i) {
      int idx = siteIndex[i];
      if (idx < 0) continue;
      assert( 0 <= idx && idx < numMarker);
      m[0][2 * idx]  = atof(fd[2 + 2 * i]);
      m[0][2 * idx + 1]  = atof(fd[2 + 2 * i + 1]);      
    }

    // flip ref/alt SNP
    for (int i = 0; i < numMarker; ++i) {
      const std::string& markerName = spaModel.getMarkers()[i]; 
      const std::string minor = spaModel.getMinor(markerName);
      const std::string ref = siteFile.getRef(markerName);

      if (m[0][i] < 0) continue; // skip NA allele
      if (minor == ref) {
        continue;
      } else {
        m[0][2 * i + 1] = m[0][2*i] - m[0][2 * i + 1];
      }
    }

    // infer location
    double x, y;
    double error = 0.01;
    if (!inferLocationSeq(m, spaModel, error, &x, &y)) {
      fprintf(stderr, "%s\t%s\t%g\t%g\n", popId.c_str(), indvId.c_str(), x, y);
      fprintf(fout,   "%s\t%s\t%g\t%g\n", popId.c_str(), indvId.c_str(), x, y);
      SimpleMatrix mBoot;
      double xBoot, yBoot;
      for (int i = 0; i < FLAG_bootstrap; ++i) {
        resample(m, spaModel, error, x, y, r, &mBoot);
        if (!inferLocationSeq(mBoot, spaModel, error, &xBoot, &yBoot)) {
          fprintf(stderr, "%s\t%s.Boot%d\t%g\t%g\n", popId.c_str(), indvId.c_str(), i+1, xBoot, yBoot);
          fprintf(fout,   "%s\t%s.Boot%d\t%g\t%g\n", popId.c_str(), indvId.c_str(), i+1, xBoot, yBoot);
        }  
      }
      if (FLAG_ci) {
        double xLow, xHigh, yLow, yHigh;
        if (!inferLocationSeqCI(m, spaModel, error, x, y, &xLow, &xHigh, &yLow, &yHigh)) {        
          fprintf(stderr, "%s\t%s.ConfInt95\t%g,%g\t%g,%g\n", popId.c_str(), indvId.c_str(), xLow, xHigh, yLow, yHigh);
          fprintf(fout,   "%s\t%s.ConfInt95\t%g,%g\t%g,%g\n", popId.c_str(), indvId.c_str(), xLow, xHigh, yLow, yHigh);
        } else {
          fprintf(stderr, "%s\t%s.ConfInt95\tNA\tNA\n", popId.c_str(), indvId.c_str());
          fprintf(fout,   "%s\t%s.ConfInt95\tNA\tNA\n", popId.c_str(), indvId.c_str());
        }
      }
    } else  {
      fprintf(stderr, "%s\t%s\tNA\tNA\n", popId.c_str(), indvId.c_str());
      fprintf(fout,   "%s\t%s\tNA\tNA\n", popId.c_str(), indvId.c_str());
    }
  }
  fclose(fout);
  time_t endTime = time(0);
  fprintf(stderr, "Analysis ends at: %s", ctime(&endTime));
  int elapsedSecond = (int) (endTime - startTime);
  fprintf(stderr, "Analysis took %d seconds", elapsedSecond);

  // free GSL
  gsl_rng_free (r);
  return 0;
};
