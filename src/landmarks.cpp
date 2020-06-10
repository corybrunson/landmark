////////////////////////////////////////////////////////////////////////////////
//' @title Maxmin + Neighborhood-based Landmark Sets
//' @author Matt Piekenbrock
//' @author Jason Cory Brunson
//' @author Yara Skaf
////////////////////////////////////////////////////////////////////////////////

#include "rankdistance.h"

using namespace Rcpp;
using namespace std;
using std::size_t;

inline double sq_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = x-y;
  return(sum(pow(diff, 2.0)));
}

inline double dist_euc(vector<double> x, vector<double> y){
    double sum = 0;
    for(int i = 0; i < x.size(); i++){
        double diff = x[i] - y[i];
        sum = sum + pow(diff, 2.0);
    }
  return(sqrt(sum));
}

// Uses the maxmin procedure to choose num_sets landmarks for balls of radius radius.
// NOTE: Rcpp does now allow use of c++ constant (e.g. FLT_MAX) in parameters
//      -> must specify default radius=INF within function
//' @rdname landmarks_maxmin
//' @description Compute landmark points using maxmin procedure.
//' @param x a data matrix.
//' @param num_sets a positive integer; the desired number of landmark points, or
//'   number of sets in a ball cover.
//' @param radius a positive real number; the desired radius of each cover set.
//' @param seed_index an integer (the first landmark to seed the algorithm)
//' @export
// [[Rcpp::export]]
IntegerVector landmarks_maxmin_cpp(const NumericMatrix& x, int num_sets = 0, float radius = -1, const int seed_index = 0) {
    int num_pts = x.nrow();

    // error handling
    if(radius < 0 && radius != -1){stop("Parameter 'radius' must be a positive number.");}
    if(num_sets < 0){stop("Parameter 'num_sets' must be >= 1.");}
    if(seed_index < 0 || seed_index >= num_pts){stop("Parameter 'seed_index' must be >=1 and <= number of data points.");}

    // additional parameter handling
    if(num_sets > num_pts){
        warning("Warning: parameter 'num_sets' was > max allowable value. Setting num_sets = number of data points.");
        num_sets = num_pts;
    }
    if(num_sets == 0 && radius == -1){num_sets = std::min(num_pts,24);} // no parameters passed -> default behavior
    if(radius == -1){radius = FLT_MAX;}

    // store indices and values of X\L
    map<int, vector<double>> pts_left;
    for(int i = 0; i < num_pts; i++){
        NumericVector vec = x.row(i);
        vector<double> point(vec.begin(),vec.end());
        pts_left.emplace(i, point);
    }

    // store indices and values of landmark set L
    map<int, vector<double>> landmarks;
    vector<double> seed_val = pts_left.at(seed_index);
    landmarks.emplace(seed_index, seed_val);

    // remove seed landmark (and any duplicates) from X\L
    for (const auto& x : pts_left){ if(x.second == seed_val){pts_left.erase(x.first);} }

    // keep track of order in which landmarks were added
    vector<int> ordered_landmarks;
    ordered_landmarks.push_back(seed_index);

    // compute remaining landmarks
    while(true){
        map<int, vector<double>> maxmin;
        double d_max = 0;

        // find max(d(x,L)) for x in X
        for(const auto& pt : pts_left){
            double d_min = DBL_MAX;
            // find min(d(x,l)) for l in L
            for(const auto& l : landmarks){
                double d = dist_euc(l.second, pt.second);
                if(d < d_min){d_min = d;}
            }

            // d_min is equal to the old max -> add this point to maxmin
            if(d_min == d_max){
                maxmin.insert(pt);
            }

            // we have a new max -> clear out maxmin and add this point instead
            if(d_min > d_max){
                d_max = d_min;
                maxmin.clear();
                maxmin.insert(pt);
            }
        }
        // done if farthest point is within radius and we have enough landmarks
        if(d_max <= radius && landmarks.size() >= num_sets){break;} // TODO: this occurs one extra time because of c

        // otherwise add new max to L and remove from X\L
        pair<int, vector<double>> l_i = make_pair(maxmin.begin()->first, maxmin.begin()->second);
        ordered_landmarks.push_back(l_i.first);
        landmarks.insert(l_i);

        // remove all points at this center
        for(const auto& pt : maxmin){
            if(pt.second == l_i.second){ pts_left.erase(pt.first); }
        }
        if(pts_left.size() <= 0){break;} // exit if all points are covered
    }
    // only return the indices of landmarks (not the values)
    IntegerVector ret = wrap(ordered_landmarks); // wrap into R data type
    return(ret+1); // switch to 1-based indexing for return
}


// Uses the euclidean lastfirst procedure to choose landmarks for nhds of size cardinality.
//' @rdname landmarks_lastfirst
//' @description Compute landmark points using maxmin procedure.
//' @param x a data matrix.
//' @param cardinality a positive integer; the desired cardinality of each
//'   landmark neighborhood, or of each set in a landmark cover.
//' @param num_sets a positive integer; the desired number of landmark points, or
//'   number of sets in a neighborhood cover.
//' @param seed_index an integer (the first landmark to seed the algorithm)
//' @export
// [[Rcpp::export]]
IntegerVector landmarks_lastfirst_cpp(const NumericMatrix& x, int num_sets = 0, int cardinality = 0, const int seed_index = 0) {
    int num_pts = x.nrow();

    // additional parameter handling
    if(num_sets > num_pts){
        warning("Warning: parameter 'num_sets' was > max allowable value. Setting num_sets = number of data points.");
        num_sets = num_pts;
    }
    if(num_sets == 0 && cardinality == 0){num_sets = std::min(num_pts,24);}
    if(cardinality == 0){cardinality = num_pts;}

    // error handling
    if(cardinality < 1 || cardinality > num_pts){stop("Parameter 'cardinality' must be >= 1 and <= number of data points.");}
    if(num_sets < 0){stop("Parameter 'num_sets' must be >= 1.");}
    if(seed_index < 0 || seed_index >= num_pts){stop("Parameter 'seed_index' must be >=1 and <= number of data points.");}

    map<int, vector<double>> Y_all; // whole space Y
    map<int, vector<double>> landmarks; // landmark set L
    map<int, vector<double>> covered; // store union of Nk check plus over all landmarks
    vector<int> ordered_landmarks; // keep track of order in which landmarks were added

    // store indices and values of X\L
    for(int i = 0; i < num_pts; i++){
        NumericVector vec = x.row(i);
        vector<double> point(vec.begin(),vec.end());
        Y_all.emplace(i, point);
    }

    // choose l0 and add it to L
    pair<int, vector<double>> l_0(seed_index, Y_all.at(seed_index));
    ordered_landmarks.push_back(seed_index);
    landmarks.insert(l_0);

    pair<int, vector<double>> l_i = l_0;
    while(true){
        // update the list of covered points
        map<int, vector<double>> Nk = Nk_check_plus(l_i, Y_all, cardinality, Y_all);
        for(const auto& x : Nk){ covered.insert(x); }

        // exit if all points in X are covered and enough landmarks have been chosen
        if(covered.size() >= num_pts && landmarks.size() >= num_sets){break;}

        // compute lf(L), then choose li from lf(L) and add it to L
        map<int, vector<double>> lf = lastfirst(landmarks, Y_all);
        l_i = make_pair(lf.begin()->first, lf.begin()->second);
        ordered_landmarks.push_back(l_i.first);
        landmarks.insert(l_i);
    }
    // only return the indices of landmarks (not the values)
    IntegerVector ret = wrap(ordered_landmarks); // wrap into R data type
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
