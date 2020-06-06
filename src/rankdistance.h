#ifndef RANKDIST_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define RANKDIST_H

#include <Rcpp.h>
#include <vector>

int q_hat(std::vector<double> x, std::vector<double> y, std::map<int, std::vector<double>> Y_all);
int q_check(std::vector<double> x, std::vector<double> y, std::map<int, std::vector<double>> Y_all);

std::map<int, std::vector<double>> Nk_hat_plus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all);
std::map<int, std::vector<double>> Nk_check_plus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all);
std::map<int, std::vector<double>> Nk_check_minus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, int k, std::map<int, std::vector<double>> Y_all);

std::vector<int> Q_hat_plus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);
std::vector<int> Q_check_minus(std::pair<int, std::vector<double>> x, std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);

int colex(std::vector<int> an, std::vector<int> bn);
int revlex(std::vector<int> an, std::vector<int> bn);

std::map<int, std::vector<double>> firstlast(std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);
std::map<int, std::vector<double>> lastfirst(std::map<int, std::vector<double>> Yp, std::map<int, std::vector<double>> Y_all);

#endif
