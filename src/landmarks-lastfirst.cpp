////////////////////////////////////////////////////////////////////////////////
// Neighborhood-based Landmark Procedure
// Authors: Yara Skaf, Jason Cory Brunson
// Description: Calculate a landmark set using the lastfirst procedure.
////////////////////////////////////////////////////////////////////////////////

#include "rankdistance.h"
using namespace Rcpp;
using namespace std;
using std::size_t;

//' @name landmarks_lastfirst_cpp
//' @title Lastfirst in C++
//' @author Yara Skaf
//' @author Jason Cory Brunson
//' @description Use the Euclidean lastfirst procedure to choose landmarks.
//' @details `landmark_lastfirst()` performs the lastfirst procedure to choose a
//'   given number landmarks for neighborhoods of at least a fixed cardinality.
//'   It supports Euclidean distances only and provides an option to collect
//'   covers.
//' @param x a data matrix
//' @param num desired number of landmark points, or number of sets, in a ball
//'   cover (should be a positive integer)
//' @param cardinality desired cardinality of a cover set (should be a positive
//'   integer)
//' @param seed_index index of the first landmark used to seed the algorithm
//' @param cover boolean specifying whether to return cover sets in addition to
//'   the landmark points

//' @rdname landmarks_lastfirst_cpp
// [[Rcpp::export]]
List landmarks_lastfirst_cpp(const NumericMatrix& x,
                             int num = 0, int cardinality = 0,
                             const int seed_index = 1,
                             const bool cover = false) {
    int num_pts = x.nrow();

    // default value to accomodate exit condition, without warning
    if (cardinality == 0 || cardinality > num_pts) { cardinality = num_pts; }
    // error handling
    if (cardinality < 0) { stop("Parameter `cardinality` must be non-negative."); }
    if (num < 0) { stop("Parameter `num` must be non-negative."); }
    if (seed_index < 1 || seed_index > num_pts) {
        stop("Parameter `seed_index` must be positive and at most `nrow(x)`.");
    }
    // additional parameter handling
    if (num > num_pts) {
        warning("Parameter `num` is too large; using `num = nrow(x)`.");
        num = num_pts;
    }
    if (num == 0 && cardinality == 0) { num = std::min(num_pts, 24); }

    // landmark set L
    map<int, vector<double>> landmarks;
    // keep track of order in which landmarks were added
    vector<int> ordered_landmarks;
    // store in-ranks from L
    map<int, map<int, int>> landmark_ranks;
    // whole space Y
    map<int, vector<double>> Y_all;
    for (int i = 0; i < num_pts; i++) {
        NumericVector vec = x.row(i);
        vector<double> point(vec.begin(),vec.end());
        Y_all.emplace(i, point);
    }
    // keep track of X\L
    map<int, vector<double>> pts_left(Y_all);

    // choose l0
    pair<int, vector<double>> l_i(seed_index - 1, Y_all.at(seed_index - 1));
    // minimum neighborhood size needed to cover X
    int min_k = num_pts;
    while (true) {
        // add landmark to landmark set, remove it (+ duplicates) from pts_left
        landmarks.insert(l_i);
        ordered_landmarks.push_back(l_i.first);
        for (auto pt = pts_left.cbegin(), next_pt = pt;
             pt != pts_left.cend(); pt = next_pt) {
            ++next_pt;
            if (pt->second == l_i.second) { pts_left.erase(pt); }
        }

        // maintain ranked k-nhd if we need to return cover sets
        if (cover == true) {
            // if cardinality was given, store that many neighbors
            // otherwise, must store min_k neighbors and remove extras later
            int num_neighbors = min(min_k, cardinality);
            map<int, int> ranked_nhd; // store index and q value of members
            for (const auto& y : Y_all) {
                int q = q_check(l_i.second, y.second, Y_all);
                if (q <= num_neighbors) { ranked_nhd.emplace(y.first, q); }
            }
            landmark_ranks.emplace(l_i.first, ranked_nhd);
        }

        // get the next landmark
        map<int, vector<double>> lf = lastfirst(landmarks, Y_all);
        if (lf.size() <= 0) { min_k = 1; break; }
        pair<int, vector<double>> next_l = make_pair(lf.begin()->first,
                                                     lf.begin()->second);

        // calculate the new minimum cardinality required to cover X
        for (const auto& l : landmarks) {
            int q = q_check(l.second, next_l.second, Y_all);
            if (q < min_k) { min_k = q; }
        }

        // exit if all points are covered and enough landmarks have been chosen
        if ( (min_k <= cardinality && landmarks.size() >= num) ||
             pts_left.size() <= 0 ) { break; }
        // otherwise continue with the next landmark
        l_i = next_l;
    }

    List ret;
    // wrap into R data type
    IntegerVector landmarks_R = wrap(ordered_landmarks);
    ret.push_back(landmarks_R + 1);

    if (cover == true) {
        List cover_sets;
        // if cardinality was given, nhds are already of appropriate size
        int cutoff = min_k;
        if (cardinality < num_pts) { cutoff = cardinality; }
        for (const auto& l : ordered_landmarks) {
            IntegerVector nhd;
            for (const auto& pt : landmark_ranks.at(l)) {
                if (pt.second <= cutoff) { nhd.push_back(pt.first); }
            }
            cover_sets.push_back(nhd + 1);
        }
        ret.push_back(cover_sets);
    }
    return(ret);
}

// original (more definitional, but slow) function to calculate lastfirst
// landmarks
List landmarks_lastfirst_cpp_deprecated(const NumericMatrix& x,
                                        int num = 0, int cardinality = 0,
                                        const int seed_index = 1,
                                        const bool cover = false) {
    int num_pts = x.nrow();

    // additional parameter handling
    if (num > num_pts) {
        warning("Parameter `num` is too large; using `num = nrow(x)`.");
        num = num_pts;
    }
    if (num == 0 && cardinality == 0) { num = std::min(num_pts, 24); }
    if (cardinality == 0) { cardinality = 1; }
    if (num == 0) { num = num_pts; }

    // error handling
    if (cardinality < 1 || cardinality > num_pts) {
        stop("Parameter `cardinality` must be >= 1 and <= number of data points.");
    }
    if (num < 0) { stop("Parameter `num` must be >= 1."); }
    if (seed_index < 1 || seed_index > num_pts) {
        stop("Parameter `seed_index` must be >=1 and <= number of data points.");
    }

    // whole space Y
    map<int, vector<double>> Y_all;
    // landmark set L
    map<int, vector<double>> landmarks;
    // store union of Nk check plus over all landmarks
    map<int, vector<double>> covered;
    // keep track of order in which landmarks were added
    vector<int> ordered_landmarks;

    // store k-nhds if cover == true
    map<int, vector<int>> nhds;


    // store indices and values of X\L
    for(int i = 0; i < num_pts; i++){
        NumericVector vec = x.row(i);
        vector<double> point(vec.begin(),vec.end());
        Y_all.emplace(i, point);
    }

    // choose l0 and add it to L
    pair<int, vector<double>> l_0(seed_index-1, Y_all.at(seed_index-1));
    ordered_landmarks.push_back(seed_index-1);
    landmarks.insert(l_0);

    // TODO: using this for now, but there is definitely a more efficient data
    // structure for this purpose
    // store in-ranks from L
    map<pair<int,int>, int> landmark_ranks;
    // minimum neighborhood size needed to cover X
    // int min_k = num_pts;

    pair<int, vector<double>> l_i = l_0;
    while (true) {
        vector<int> new_nhd;
        // update the list of covered points
        // expensive operation
        map<int, vector<double>> Nk = Nk_check_plus(l_i, Y_all,
                                                    cardinality, Y_all);
        for (const auto& x : Nk) {
            covered.insert(x);
            new_nhd.push_back(x.first);
        }
        // nhds.emplace(l_i.first, Nk);
        nhds.emplace(l_i.first, new_nhd);

        // exit if all points are covered and enough landmarks have been chosen
        if (covered.size() >= num_pts || landmarks.size() >= num) { break; }

        // compute lf(L), then choose li from lf(L) and add it to L
        map<int, vector<double>> lf = lastfirst(landmarks, Y_all);
        l_i = make_pair(lf.begin()->first, lf.begin()->second);
        ordered_landmarks.push_back(l_i.first);
        landmarks.insert(l_i);

        // cout << "\nlastfirst before insertion of landmark " << l_i.first;
        // for (const auto& x : lf) {
        //     cout << "\n" << x.first << ": ( ";
        //     for (const auto& v : x.second) { cout << v << " "; }
        //     cout << ")";
        // }
        // cout << "\n\n";
    }

    // cout << "landmark rank distances:";
    // for (const auto& l1 : landmarks) {
    //     for (const auto& l2 : landmarks) {
    //         int temp = q_hat(l1.second,l2.second,Y_all);
    //         cout << "\nq_hat(" << l1.first << "," << l2.first << ") = " << temp << "\n";
    //
    //         temp = q_check(l1.second,l2.second,Y_all);
    //         cout << "q_check(" << l1.first << "," << l2.first << ") = " << temp << "\n";
    //     }
    // }

    // wrap into R data type
    IntegerVector landmarks_R = wrap(ordered_landmarks);

    List cover_sets;
    if (cover == true) {
        for (const auto& x : ordered_landmarks) {
            IntegerVector nhd = wrap(nhds[x]);
            cover_sets.push_back(nhd + 1);
        }
    }
    // use null for cover_sets if cover==false
    else { cover_sets = R_NilValue; }

    List ret;
    ret.push_back(landmarks_R + 1);
    ret.push_back(cover_sets);
    return(ret);
}
