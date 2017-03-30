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

#include "tabix.h"

#include "Argument.h"
#include "Utils.h"
#include "VCFUtil.h"
#include "SimpleMatrix.h"
#include "IO.h"
#include "PlinkInputFile.h"
#include "Time.h"

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
      value.a1 = atof(fd[6]);
      value.a2 = atof(fd[7]);
      value.b  = atof(fd[8]);
      if (this->data.find(key) != this->data.end()) {
        fprintf(stderr, "Duplicated marker: %s\n", key.c_str() );
        continue;
      }
      this->data[key] = value;
      this->coef.push_back(value);
      this->markers.push_back(fd[1]);
    }
  }
  Coef& getData(const std::string& marker){
    if (data.find(marker) == data.end()) {
      fprintf(stderr, "Non-exist marker: %s\n", marker.c_str() );
    }
    return data[marker];
  }
  const Coef& operator [] (int idx) const {
    return coef[idx];
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
};

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

class MarkerFrequency{
 public:
  MarkerFrequency(const std::string& fn){
    LineReader lr(fn);
    std::vector<std::string> fd;
    std::string key;
    std::string line;
    while(lr.readLine(&line)) {
      stringNaturalTokenize(line, "\t ", &fd);
      if (fd[1] == "SNP") continue; // skip header
      key = fd[1];
      if (major.count(key)) {
        fprintf(stderr, "Duplicated marker: %s\n", key.c_str() );
        continue;
      }
      minor[key] = fd[2];
      major[key] = fd[3];
      markers.push_back(fd[1]);
    }
  }
  const std::string& getMajor(const std::string& key) {
    if (major.count(key) == 0) {
      fprintf(stderr, "Non-exist marker: %s\n", key.c_str() );
    }
    return major[key];
  }
  const std::string& getMinor(const std::string& key) {
    if (minor.count(key) == 0) {
      fprintf(stderr, "Non-exist marker: %s\n", key.c_str() );
    }
    return minor[key];
  }
  size_t size() const {
    return this->major.size();
  }
  const std::vector< std::string >& getMarkers() const{
    return this->markers;
  }
 private:
  std::map<std::string, std::string> major; // store major
  std::map<std::string, std::string> minor; // store minor
  std::vector<std::string> markers;
};

class Param {
 public:
  SimpleMatrix* geno;
  SPAModel* model;
};

double computeLikelihood(const gsl_vector *v, void *params) {
  double x, y;
  Param* p = (Param*) params;
  SimpleMatrix& geno = * (p->geno);
  SPAModel& model = * (p->model);
  const int numMarker = geno.ncol();

  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  double llk = 0;
  for (int j = 0; j < numMarker; ++j) {
    if (geno[0][j] < 0) continue;  // skip missing genotypes
    const Coef& coef = model[j];
    double f = 1.0 / (exp(-coef.a1 * x - coef.a2 * y - coef.b ) + 1);
    llk += geno[0][j] * log(f) + (2.0 - geno[0][j]) * log (1-f);
  }
  return -llk;
}

/**
 * @return 0: if location are successfully inferred
 */
int inferLocation(SimpleMatrix& m, SPAModel& model,
                  double* x, double* y) {
  printf("m dim is %d x %d, and model has %zu markers\n", m.nrow(), m.ncol(), model.size());

  // initialize start points
  gsl_vector* init = gsl_vector_alloc (2);
  gsl_vector_set (init, 0, 0.0);
  gsl_vector_set (init, 1, 0.0);

  gsl_vector* ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1.0);

  // init parameter
  Param param;
  param.geno = &m;
  param.model = &model;

  // initialize functions
  const gsl_multimin_fminimizer_type *T =
      gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;

  gsl_multimin_function minex_func;
  minex_func.n = 2;
  minex_func.f = computeLikelihood;
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
  }  while (status == GSL_CONTINUE && iter < 100);

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

int main(int argc, char** argv){
  time_t startTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&startTime));

  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inPlink, "--inPlink", "input PLINK file prefix")
      ADD_STRING_PARAMETER(pl, inModel, "--inModel", "input SPA format output")
      ADD_STRING_PARAMETER(pl, inFreq,  "--inFreq", "input PLINK frequency file")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
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

  REQUIRE_STRING_PARAMETER(FLAG_inPlink, "Please provide input file using: --inPlink");
  REQUIRE_STRING_PARAMETER(FLAG_inModel, "Please provide input file using: --inModel");
  REQUIRE_STRING_PARAMETER(FLAG_inFreq, "Please provide input file using: --inFreq");

  SPAModel spaModel(FLAG_inModel);
  MarkerFrequency markerFreq(FLAG_inFreq);
  PlinkInputFile pin(FLAG_inPlink);
  FILE* fout = fopen((FLAG_outPrefix + ".loc").c_str(), "wt");

  // sanity check
  if (!isSameVector(markerFreq.getMarkers(), spaModel.getMarkers())) {
    fprintf(stderr, "frequency file and marker does not match on rsid.\n");
    exit(1);
  }
  fprintf(stderr, "loaded %zu markers\n", markerFreq.size());

  // loop each marker in the plink file
  const int numMarker = pin.getNumMarker();
  std::string key;
  fprintf(stderr, "ID\tLoc1\tLoc2\n");
  fprintf(fout, "ID\tLoc1\tLoc2\n");

  SimpleMatrix m;
  std::vector<std::string> peopleToExtract;
  for (size_t i = 0; i < pin.indv.size(); ++i) {
    fprintf(stderr, "Process sample %zu of %zu\n", i,pin.indv.size());
    
    peopleToExtract.clear();
    const std::string& peopleName = pin.indv[i];
    peopleToExtract.push_back(peopleName);
    pin.readIntoMatrix(&m, &peopleToExtract, NULL);
    double x, y;

    // flip major/minor SNP
    for (int i = 0; i < numMarker; ++i) {
      key = pin.chrom[i] + ":" + toString(pin.pos[i]);
      const std::string& markerName = pin.snp[i];
      const std::string major = markerFreq.getMajor(markerName);
      const std::string minor = markerFreq.getMinor(markerName);

      if (m[0][i] < 0) continue; // skip NA allele
      if (pin.ref[i] == major[0]) {
        continue;
      } else {
        m[0][i] = 2 - m[0][i];
      }
    }

    // infer location
    if (!inferLocation(m, spaModel, &x, &y)) {
      fprintf(stderr, "%s\t%g\t%g\n", peopleName.c_str(), x, y);
      fprintf(fout, "%s\t%g\t%g\n", peopleName.c_str(), x, y);
    } else  {
      fprintf(stderr, "%s\tNA\tNA\n", peopleName.c_str());
      fprintf(fout, "%s\tNA\tNA\n", peopleName.c_str());
    }
  }
  fclose(fout);
  time_t endTime = time(0);
  fprintf(stderr, "Analysis ends at: %s", ctime(&endTime));
  int elapsedSecond = (int) (endTime - startTime);
  fprintf(stderr, "Analysis took %d seconds", elapsedSecond);

  return 0;
};
