#ifndef RANKDIST_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define RANKDIST_H

#include <Rcpp.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

enum Type {HAT, CHECK};
enum Sign {PLUS, MINUS};

int q_hat(std::vector<double> x, std::vector<double> y, std::map<int, std::vector<double>> Y_all);
int q_check(std::vector<double> x, std::vector<double> y, std::map<int, std::vector<double>> Y_all);

std::map<int, std::vector<double>> get_Nk(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all, Type t, Sign s);
std::map<int, std::vector<double>> Nk_hat_plus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all);
std::map<int, std::vector<double>> Nk_hat_minus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all);
std::map<int, std::vector<double>> Nk_check_plus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all);
std::map<int, std::vector<double>> Nk_check_minus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all);

std::vector<int> get_Q(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all, Type t, Sign s);
std::vector<int> Q_hat_plus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);
std::vector<int> Q_hat_minus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);
std::vector<int> Q_check_plus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);
std::vector<int> Q_check_minus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);

int colex(std::vector<int> an, std::vector<int> bn);
int revlex(std::vector<int> an, std::vector<int> bn);

std::map<int, std::vector<double>> firstlast(std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);
std::map<int, std::vector<double>> lastfirst(std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);

inline double dist_euc(std::vector<double> x, std::vector<double> y){
    double sum = 0;
    for(int i = 0; i < x.size(); i++){
        double diff = x[i] - y[i];
        sum = sum + pow(diff, 2.0);
    }
  return(sqrt(sum));
}

#endif
