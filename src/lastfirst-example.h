#include <Rcpp.h>
#include "rankdistance.h"

using namespace Rcpp;

#include <functional>
#include <algorithm>

#include <iostream>
using namespace std;
using std::size_t;

void prettyprint_pts(map<int, vector<double>> s);
void prettyprint_q(map<int, vector<double>> Y_all, Type t);
void prettyprint_Nk(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all, Type t, Sign s);
void prettyprint_Q(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all, Type t, Sign s);
