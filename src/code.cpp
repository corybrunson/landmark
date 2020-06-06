#include <Rcpp.h>
#include "rankdistance.h"

using namespace Rcpp;

#include <functional>
#include <algorithm>

#include <iostream>
using namespace std;

using std::size_t;


////////////////////////////////////////////////////////////////////////////////

inline double sq_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = x-y;
  return(sum(pow(diff, 2.0)));
}

// inline double dist_euc(const NumericVector& x, const NumericVector& y){
//   return(sqrt(sq_dist(x,y)));
// }

inline double dist_euc(vector<double> x, vector<double> y){
    double sum = 0;
    for(int i = 0; i < x.size(); i++){
        double diff = x[i] - y[i];
        sum = sum + pow(diff, 2.0);
    }
  return(sqrt(sum));
}

// Uses the maxmin procedure to choose n landmarks for balls of radius eps.
// NOTE: Rcpp does now allow use of c++ constant (e.g. FLT_MAX) in parameters
//      -> must specify default eps=INF within function
//' @rdname landmarks_maxmin_cpp
//' @export
// [[Rcpp::export]]
IntegerVector landmarks_maxmin_cpp(const NumericMatrix& x, int n = 0, float eps = 0, const int seed_index = 0) {
    int num_pts = x.nrow();

    // error handling
    if(eps < 0){stop("Parameter 'eps' must be a positive number.");}
    if(n < 0){stop("Parameter 'n' must be >= 1.");}
    if(seed_index < 0 || seed_index >= num_pts){stop("Parameter 'seed_index' must be >=1 and <= number of data points.");}

    // additional parameter handling
    if(n > num_pts){
        warning("Warning: parameter 'n' was > max allowable value. Setting n = number of data points.");
        n = num_pts;
    }
    if(eps == 0){eps = FLT_MAX;}

    // store indices and values of X\L
    std::map<int, std::vector<double>> pts_left;
    for(int i = 0; i < num_pts; i++){
        NumericVector vec = x.row(i);
        std::vector<double> point(vec.begin(),vec.end());
        pts_left.emplace(i, point);
    }

    // store indices and values of landmark set L
    std::map<int, std::vector<double>> landmarks;
    landmarks.emplace(seed_index, pts_left.at(seed_index));
    pts_left.erase(seed_index); // remove seed landmark from X\L

    // compute remaining landmarks
    while(true){
        double d_max = 0;
        std::pair<int, std::vector<double>> max;
        // find max(d(x,L)) for x in X
        for(const auto& pt : pts_left){
            double d_min = DBL_MAX;
            // find min(d(x,l)) for l in L
            for(const auto& l : landmarks){
                double d = dist_euc(l.second, pt.second);
                if(d < d_min){d_min = d;}
            }
            if(d_min > d_max){
                d_max = d_min;
                max = pt;
            }
        }
        // done if farthest point is within epsilon and we have enough landmarks
        if(d_max <= eps && landmarks.size() >= n){break;} // TODO: this occurs one extra time because of c

        // otherwise add new max to L and remove from X\L
        landmarks.insert(max);
        pts_left.erase(max.first);
    }
    // only return the indices of landmarks (not the values)
    std::vector<int> landmark_idx;
    for (const auto& l : landmarks){landmark_idx.push_back(l.first);}

    IntegerVector ret = wrap(landmark_idx); // wrap into R data type
    return(ret+1); // switch to 1-based indexing for return
}


// Uses the euclidean lastfirst procedure to choose landmarks for nhds of size k.
//' @rdname landmarks_lastfirst_cpp
//' @export
// [[Rcpp::export]]
IntegerVector landmarks_lastfirst_cpp(const NumericMatrix& x, const int k, const int seed_index = 0) {
    int num_pts = x.nrow();

    // error handling
    if(k < 1 || k > num_pts){stop("Parameter 'k' must be >= 1 and <= number of data points.");}
    if(seed_index < 0 || seed_index >= num_pts){stop("Parameter 'seed_index' must be >=1 and <= number of data points.");}

    map<int, vector<double>> Y_all; // whole space Y
    map<int, vector<double>> landmarks; // landmark set L
    map<int, vector<double>> covered; // store union of Nk check plus over all landmarks

    // store indices and values of X\L
    for(int i = 0; i < num_pts; i++){
        NumericVector vec = x.row(i);
        std::vector<double> point(vec.begin(),vec.end());
        Y_all.emplace(i, point);
    }
    //map<int, vector<double>> pts_left(Y_all); // store indices and values of Y\L

    // choose l0 and add it to L
    pair<int, vector<double>> l_0(seed_index, Y_all.at(seed_index));
    landmarks.insert(l_0);
    //pts_left.erase(l_0.first); // remove l0 from list of points left

    pair<int, vector<double>> l_i = l_0;
    while(true){
        // update the list of covered points
        map<int, vector<double>> Nk = Nk_check_plus(l_i, Y_all, k, Y_all);
        for(const auto& x : Nk){ covered.insert(x); }
        if(covered.size() >= num_pts){break;}

        // compute lf(L)
        map<int, vector<double>> lf = lastfirst(landmarks, Y_all);

        // choose li from lf(L) and add it to L
        l_i = make_pair(lf.begin()->first, lf.begin()->second);
        landmarks.insert(l_i);
        //pts_left.erase(l_i.first);
    }

    // only return the indices of landmarks (not the values)
    std::vector<int> landmark_idx;
    for (const auto& l : landmarks){landmark_idx.push_back(l.first);}

    IntegerVector ret = wrap(landmark_idx); // wrap into R data type
    return(ret+1); // switch to 1-based indexing for return
}


// Uses the euclidean maxmin procedure to choose n landmarks.
//' @rdname landmarks_maxmin
//' @export
// [[Rcpp::export]]
IntegerVector landmark_maxmin(const NumericMatrix& x, const int n, const int seed_index = 0) {
  const size_t n_pts = x.nrow();
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());
  IntegerVector landmark_idx = no_init_vector(n);
  //landmark_idx[seed_index] = 0;
  landmark_idx[0] = seed_index;
  IntegerVector::iterator c_landmark = landmark_idx.begin();
  double new_min_dist;
  std::generate(landmark_idx.begin()+1, landmark_idx.end(), [&](){
    size_t i = 0;
    new_min_dist = std::numeric_limits<double>::infinity();

    // Replace nearest landmark distances
    std::replace_if(lm_dist.begin(), lm_dist.end(), [&i, &x, &c_landmark, &new_min_dist](const double c_dist){
      new_min_dist = sq_dist(x.row(*c_landmark), x.row(i++));
      return(new_min_dist < c_dist);
    }, new_min_dist);

    // Find the point that maximizes said distances, move to next landmark
    c_landmark++;
    auto max_landmark = std::max_element(lm_dist.begin(), lm_dist.end());
    return(std::distance(lm_dist.begin(), max_landmark));
  });
  return(landmark_idx+1); // 1-based return
}
