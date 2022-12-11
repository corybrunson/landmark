////////////////////////////////////////////////////////////////////////////////
// Maxmin-based Landmark Procedure
// Authors: Matt Piekenbrock, Jason Cory Brunson, Yara Skaf
// Description: Calculate a landmark set using the maxmin procedure.
////////////////////////////////////////////////////////////////////////////////

#include "rankdistance.h"
using namespace Rcpp;
using namespace std;
using std::size_t;

inline double sq_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = x-y;
  return(sum(pow(diff, 2.0)));
}

//' @title Maxmin in C++ (original)
//' @author Matt Piekenbrock
//' @description Use the Euclidean maxmin procedure to choose landmarks
//' @description Maxmin procedure to choose 'num' landmarks for balls of fixed
//' radius (supports euclidean distance only)
//' @param x a data matrix
//' @param num desired number of landmark points, or number of sets, in a ball
//'   cover (should be a positive integer)
//' @param seed_index index of the first landmark used to seed the algorithm
// [[Rcpp::export]]
IntegerVector landmark_maxmin(const NumericMatrix& x,
                              const int num,
                              const int seed_index = 1) {
  const size_t n_pts = x.nrow();
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());
  IntegerVector landmark_idx = no_init_vector(num);
  // landmark_idx[seed_index] = 0;
  landmark_idx[0] = seed_index-1;
  IntegerVector::iterator c_landmark = landmark_idx.begin();
  double new_min_dist;
  std::generate(landmark_idx.begin()+1, landmark_idx.end(), [&]() {
    size_t i = 0;
    new_min_dist = std::numeric_limits<double>::infinity();

    // Replace nearest landmark distances
    std::replace_if(lm_dist.begin(), lm_dist.end(),
                    [&i, &x, &c_landmark, &new_min_dist](const double c_dist) {
                      new_min_dist = sq_dist(x.row(*c_landmark), x.row(i++));
                      return(new_min_dist < c_dist);
                    }, new_min_dist);

    // Find the point that maximizes said distances, move to next landmark
    c_landmark++;
    auto max_landmark = std::max_element(lm_dist.begin(), lm_dist.end());
    return(std::distance(lm_dist.begin(), max_landmark));
  });
  // 1-based return
  return(landmark_idx + 1);
}

// define a tolerance for distances to avoid floating point rounding errors
#define EPS_TOL 0.001

//' @title Maxmin in C++ (extended)
//' @author Jason Cory Brunson
//' @author Yara Skaf
//' @description Maxmin procedure to choose landmarks for balls of fixed radius
//'   (supports Euclidean distance only)
//' @inheritParams landmark_maxmin
//' @param num desired number of landmark points, or number of sets, in a ball
//'   cover (should be a positive integer)
//' @param radius desired radius of a cover set (should be a positive real number)
//' @param cover boolean specifying whether to return cover sets in addition to
//'   the landmark points
// [[Rcpp::export]]
List landmarks_maxmin_cpp(const NumericMatrix& x,
                          int num = 0, float radius = -1,
                          const int seed_index = 1, const bool cover = false) {
  int num_pts = x.nrow();

  // error handling
  if (radius < 0 && radius != -1) {
    stop("Parameter `radius` must be positive.");
  }
  if (num < 0) {
    stop("Parameter `num` must be non-negative.");
  }
  if (seed_index < 1 || seed_index > num_pts) {
    stop("Parameter `seed_index` must be positive and at most `nrow(x)`.");
  }

  // additional parameter handling
  if (num > num_pts) {
    warning("Parameter `num` is too large; using `num = nrow(x)`.");
    num = num_pts;
  }
  // no parameters passed -> default behavior
  if (num == 0 && radius == -1) { num = std::min(num_pts, 24); }
  // Rcpp does now allow use of C++ constants (e.g. `FLT_MAX`) in parameters
  if (radius == -1) { radius = FLT_MAX; }

  // store indices and values of X\L
  map<int, vector<double>> pts_left;
  for (int i = 0; i < num_pts; i++) {
    NumericVector vec = x.row(i);
    vector<double> point(vec.begin(), vec.end());
    pts_left.emplace(i, point);
  }

  // remove seed landmark (and any duplicates) from X\L
  vector<double> seed_val = pts_left.at(seed_index - 1);
  for (auto pt = pts_left.cbegin(), next_pt = pt;
       pt != pts_left.cend();
       pt = next_pt) {
    ++next_pt;
    if (pt->second == seed_val) { pts_left.erase(pt); }
  }

  // store indices and values of landmark set L
  map<int, vector<double>> landmarks;
  // keep track of order landmarks are added
  vector<int> ordered_landmarks;
  landmarks.emplace(seed_index - 1, seed_val);
  ordered_landmarks.push_back(seed_index - 1);

  // compute remaining landmarks
  while (true) {
    map<int, vector<double>> maxmin;
    double d_max = 0;

    // find max(d(x,L)) for x in X
    for (const auto& pt : pts_left) {
      double d_min = DBL_MAX;

      // find min(d(x,l)) for l in L
      for (const auto& l : landmarks) {
        double d = dist_euc(l.second, pt.second);
        if (d < d_min){ d_min = d; }
      }

      // d_min is equal to the old max -> add this point to maxmin
      if (d_min == d_max){ maxmin.insert(pt); }

      // we have a new max -> clear out maxmin and add this point instead
      if (d_min > d_max) {
        d_max = d_min;
        maxmin.clear();
        maxmin.insert(pt);
      }
    }
    // done if farthest point is within radius and we have enough landmarks
    if (d_max <= radius && landmarks.size() >= num) {
      if (radius == FLT_MAX){ radius = d_max; }
      break;
    }

    // otherwise add new max to L and remove from X\L
    pair<int, vector<double>> l_i = make_pair(maxmin.begin()->first,
                                              maxmin.begin()->second);
    ordered_landmarks.push_back(l_i.first);
    landmarks.insert(l_i);

    // remove all points at this center & exit if all points are covered
    for (const auto& pt : maxmin) {
      if (pt.second == l_i.second) { pts_left.erase(pt.first); }
    }
    if (pts_left.size() <= 0) {
      // all points coincide with a point in L, so ball radius should be 0
      if (radius == FLT_MAX) { radius = 0; }
      break;
    }
  }

  List ret;
  // wrap into R data type
  IntegerVector landmarks_R = wrap(ordered_landmarks);
  ret.push_back(landmarks_R + 1);

  if (cover == true) {
    List cover_sets;
    for (const auto& l : ordered_landmarks) {
      vector<double> l_val = landmarks.at(l);
      IntegerVector ball;
      for (int i = 0; i < num_pts; i++) {
        // if pt is within the radius (+ tolerance), add it to the ball
        NumericVector vec = x.row(i);
        vector<double> pt(vec.begin(), vec.end());
        if (dist_euc(l_val, pt) <= radius+EPS_TOL) { ball.push_back(i); }
      }
      cover_sets.push_back(ball + 1);
    }
    ret.push_back(cover_sets);
  }
  return(ret);
}
