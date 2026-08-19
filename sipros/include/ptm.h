#ifndef _PTM_H
#define _PTM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "proNovoConfig.h"

using namespace std;

class PTM_List {

  vector <char> _residue;
  vector <double> _mass_shift;
  vector <string> _symbol;

public:

  PTM_List() {}
  ~PTM_List() {}

  char residue(size_t n) const { return _residue[n]; }
  double mass_shift(size_t n) const { return _mass_shift[n]; }
  string symbol(size_t n) const { return _symbol[n]; }

  size_t size() const { return _residue.size(); }

  bool populateFromEnabledPtms();
  
  bool add_ptm(char, string, double);
  void dump();

};

#endif
