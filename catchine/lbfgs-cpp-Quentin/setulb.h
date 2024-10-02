#ifndef SETULB_H
#define  SETULB_H

#include <vector>
#include <cassert>
#include <cstring>
#include <cmath>
#include <iostream>

typedef int integer;
typedef float real;
typedef unsigned int uinteger;
typedef int logical;
typedef int flag;
typedef int ftnlen;
typedef int ftnint;

typedef char *address;
typedef short int shortint;
typedef double doublereal;
typedef short int shortlogical;

int setulb_(integer *n, integer *m, doublereal *x, 
	    doublereal *l, doublereal *u, integer *nbd,
	    doublereal *f, doublereal 
	    *g, doublereal *factr, doublereal *pgtol,
	    doublereal *wa, integer * iwa,
	    char *task, integer *iprint, char *csave,
	    logical *lsave,  integer *isave, doublereal *dsave,
	    ftnlen task_len, ftnlen csave_len);

enum Lbfgs_result
  {
    LBFGS_FG,
    LBFGS_NEW_X,
    LBFGS_CONVERGENCE,
    LBFGS_ERROR
  };

typedef std::vector<double> uvector;

inline 
double norm_1(const uvector &u)
{
  double m = 0.0;
  for (size_t i = 0; i < u.size(); ++i)
    m += fabs(u[i]);
  return m;
}


inline 
double norm_2(const uvector &u)
{
  double m = 0.0;
  for (size_t i = 0; i < u.size(); ++i)
    m += u[i]*u[i];
  return sqrt(m);
}


class Lbfgs
{
  int _N, _m;
  char _csave[60], _task[60];
  double _dsave[29];
  int _isave[44],_lsave[4];
  double _factr, _pgtol;
  std::vector<int> _iwa, _nbd;
  std::vector<double> _wa;
  std::vector<double> _lc, _uc;
  
public:
  Lbfgs (size_t N, size_t m = 25) : _N(N), _m(m)
  {
    clear ();
  }

  void clear()
  {
    _nbd = std::vector<int> (_N, 0);
    _wa.resize((2*_m + 4) * _N + 12 * _m * _m + 12 * _m);
    _iwa.resize(3*_N);    
    _csave[0] = '\0';
    _lc.resize(_N, 0.0);
    _uc.resize(_N, 0.0);
    _factr = 1e7;
    _pgtol = 1e-5; 
    strcpy(_task, "START");
    for (size_t i = 5; i < 60; ++i)
      _task[i] = ' ';
  }

  Lbfgs_result iterate (uvector &x, double fx, uvector &g)
  {
    assert(x.size() == size_t(_N));
    assert(g.size() == size_t(_N));

    int iprint = -1;    
    setulb_(&_N, &_m, &x[0],
	    &_lc[0], &_uc[0], 
	    &_nbd[0], &fx, &g[0], &_factr, &_pgtol, &_wa[0],
	    &_iwa[0], _task, &iprint, _csave, _lsave, _isave,
	    _dsave, 60, 60);
    _task[59] = '\0';
    if(!strncmp(_task, "FG", 2))
      return LBFGS_FG;
    else if (!strncmp(_task, "NEW_X", 5))
      return LBFGS_NEW_X;
    else if (!strncmp(_task, "CONV", 4))
      {
	std::cerr << "setulb: " << _task << "\n";
	return LBFGS_CONVERGENCE;
      }
    else
      {
	std::cerr << "setulb: " << _task << "\n";
	return LBFGS_ERROR;
      }
  }
};


#endif
