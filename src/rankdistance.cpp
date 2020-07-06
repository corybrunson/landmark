#include <Rcpp.h>
#include "rankdistance.h"

using namespace Rcpp;

#include <functional>
#include <algorithm>

#include <iostream>
using namespace std;

using std::size_t;


// q_hat(x,y) = |{z in Y | d(x,z) <= d(x,y)}|
// cardinality of smallest ball centered at 2 that contains y
int q_hat(vector<double> x, vector<double> y, map<int, vector<double>> Y_all){
    int q = 0;
    for(const auto& z : Y_all){
        if(dist_euc(x,z.second) <= dist_euc(x,y)){q++;}
    }
    return(q);
}

// q_check(x,y) = |{z in Y | d(x,z) < d(x,y)}| + 1
int q_check(vector<double> x, vector<double> y, map<int, vector<double>> Y_all){
    int q = 0;
    for(const auto& z : Y_all){
        if(dist_euc(x,z.second) < dist_euc(x,y)){q++;}
    }
    return(q+1);
}

// Nk_hat_plus(x,Y',k) = {y in Y' | q_hat(x,y) <= k}
// Nk_hat_minus(x,Y',k) = {y in Y' | q_hat(y,x) <= k}
// Nk_check_plus(x,Y',k) = {y in Y' | q_check(x,y) <= k}
// Nk_check_minus(x,Y',k) = {y in Y' | q_check(y,x) <= k}
map<int, vector<double>> get_Nk(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all, Type t, Sign s){
    map<int, vector<double>> Nk;
    for(const auto& y : Yp){
        int q;
        if(t == HAT && s == PLUS){ q = q_hat(x.second, y.second, Y_all); }
        if(t == HAT && s == MINUS){ q = q_hat(y.second, x.second, Y_all); }
        if(t == CHECK && s == PLUS){ q = q_check(x.second, y.second, Y_all); }
        if(t == CHECK && s == MINUS){ q = q_check(y.second, x.second, Y_all); }

        if(q <= k){Nk.insert(y);}
    }
    return(Nk);
}
map<int, vector<double>> Nk_hat_plus(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all){
    return(get_Nk(x, Yp, k, Y_all, HAT, PLUS));
}
map<int, vector<double>> Nk_hat_minus(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all){
    return(get_Nk(x, Yp, k, Y_all, HAT, MINUS));
}
map<int, vector<double>> Nk_check_plus(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all){
    return(get_Nk(x, Yp, k, Y_all, CHECK, PLUS));
}
map<int, vector<double>> Nk_check_minus(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all){
    return(get_Nk(x, Yp, k, Y_all, CHECK, MINUS));
}

// // Nk_hat_plus(x,Y',k) = {y in Y' | q_hat(x,y) <= k}
// map<int, vector<double>> Nk_hat_plus(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all){
//     map<int, vector<double>> Nk;
//     for(const auto& y : Yp){
//         int q = q_hat(x.second, y.second, Y_all);
//         if(q <= k){Nk.insert(y);}
//     }
//     return(get_Nk(x, Yp, k, Y_all, HAT, PLUS));
// }
//
// // Nk_check_plus(x,Y',k) = {y in Y' | q_check(x,y) <= k}
// map<int, vector<double>> Nk_check_plus(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all){
//     map<int, vector<double>> Nk;
//     for(const auto& y : Yp){
//         int q = q_check(x.second, y.second, Y_all);
//         if(q <= k){Nk.insert(y);}
//     }
//     return(Nk);
// }
//
// // Nk_check_minus(x,Y',k) = {y in Y' | q_check(y,x) <= k}
// map<int, vector<double>> Nk_check_minus(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k,  map<int, vector<double>> Y_all){
//     map<int, vector<double>> Nk;
//     for(const auto& y : Yp){
//         //int q = q_check(y.second, x.second, Yp);
//         int q = q_check(y.second, x.second, Y_all);
//         if(q <= k){Nk.insert(y);}
//     }
//     return(Nk);
// }

// Q_hat_plus(x,Y') = (|Nk_hat_plus(x,Y',1)|, ..., |Nk_hat_plus(x,Y',N)|)
// Q_hat_minus(x,Y') = (|Nk_hat_minus(x,Y',1)|, ..., |Nk_hat_minus(x,Y',N)|)
// Q_check_plus(x,Y') = (|Nk_check_plus(x,Y',1)|, ..., |Nk_check_plus(x,Y',N)|)
// Q_check_minus(x,Y') = (|Nk_check_minus(x,Y',1)|, ..., |Nk_check_minus(x,Y',N)|)
vector<int> get_Q(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all, Type t, Sign s){
    int bigN = Y_all.size();
    vector<int> q (bigN);
    // Compute Nk for k=0 up to N
    for(int i = 1; i <= bigN; i++){
        map<int, vector<double>> Nk;
        if(t == HAT && s == PLUS){ Nk = get_Nk(x, Yp, i, Y_all, HAT, PLUS); }
        if(t == HAT && s == MINUS){ Nk = get_Nk(x, Yp, i, Y_all, HAT, MINUS); }
        if(t == CHECK && s == PLUS){ Nk = get_Nk(x, Yp, i, Y_all, CHECK, PLUS); }
        if(t == CHECK && s == MINUS){ Nk = get_Nk(x, Yp, i, Y_all, CHECK, MINUS); }

        q[i-1] = Nk.size();
    }
    return(q);
}
vector<int> Q_hat_plus(pair<int, vector<double>> x, map<int, vector<double>> Yp,  map<int, vector<double>> Y_all){
    return(get_Q(x, Yp, Y_all, HAT, PLUS));
}
vector<int> Q_hat_minus(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all){
    return(get_Q(x, Yp, Y_all, HAT, MINUS));
}
vector<int> Q_check_plus(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all){
    return(get_Q(x, Yp, Y_all, CHECK, PLUS));
}
vector<int> Q_check_minus(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all){
    return(get_Q(x, Yp, Y_all, CHECK, MINUS));
}

// // Q_hat_plus(x,Y') = (|Nk_hat_plus(x,Y',1)|, ..., |Nk_hat_plus(x,Y',N)|)
// vector<int> Q_hat_plus(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all){
//     int bigN = Y_all.size();
//     vector<int> q (bigN);
//     // Compute Nk for k=0 up to N
//     for(int i = 1; i <= bigN; i++){
//         //map<int, vector<double>> Nk = Nk_hat_plus(x, Yp, i, Y_all);
//         map<int, vector<double>> Nk = get_Nk(x, Yp, i, Y_all, HAT, PLUS);
//         q[i-1] = Nk.size();
//     }
//     return(q);
// }
//
// // Q_check_minus(x,Y') = (|Nk_check_minus(x,Y',1)|, ..., |Nk_check_minus(x,Y',N)|)
// vector<int> Q_check_minus(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all){
//     int bigN = Y_all.size();
//     vector<int> q (bigN);
//     // Compute Nk for k=0 up to N
//     for(int i = 1; i <= bigN; i++){
//         //map<int, vector<double>> Nk = Nk_check_minus(x, Yp, i, Y_all);
//         map<int, vector<double>> Nk = get_Nk(x, Yp, i, Y_all, CHECK, MINUS);
//         q[i-1] = Nk.size();
//     }
//     return(q);
// }

// return -1 if an < bn
// return 0 if an == bn
// return +1 if an > bn
int colex(vector<int> an, vector<int> bn){
    int ret = 0;
    for(int i = an.size()-1; i >= 0; i--){
        if(an[i] != bn[i]){
            if(an[i] < bn[i]){
                ret = -1;
                break;
            }
            else{
                ret = 1;
                break;
            }
        }
    }
    return(ret);
}

// return -1 if an < bn
// return 0 if an == bn
// return +1 if an > bn
int revlex(vector<int> an, vector<int> bn){
    int ret = 0;
    for(int i = 0; i < an.size(); i++){
        if(an[i] != bn[i]){
            if(an[i] > bn[i]){
                ret = -1;
                break;
            }
            else{
                ret = 1;
                break;
            }
        }
    }
    return(ret);
}

// fl(Y') = {x in Y | Q_hat_plus(x,Y') = min Q_hat_plus(x',Y') for x' in Y}
map<int, vector<double>> firstlast(map<int, vector<double>> Yp, map<int, vector<double>> Y_all){
    map<int, vector<double>> fl;
    int length = Y_all.size();

    map<int, vector<double>> pts_left(Y_all);
    for(const auto& y : Yp){
        pts_left.erase(y.first);
    }

    // compute all Q(x,Y') and find the min
    vector<int> minQ(length, INT_MAX);
    for(const auto& y : pts_left){
        vector<int> newQ = get_Q(y, Yp, Y_all, HAT, PLUS);

        // this Q is equal to the old min -> add this point to fl
        if(colex(newQ, minQ) == 0){
            fl.emplace(y.first, Y_all.at(y.first));
        }

        // we have a new min -> clear out fl and add this point instead
        if(colex(newQ, minQ) < 0){
            minQ = newQ;
            fl.clear();
            fl.emplace(y.first, Y_all.at(y.first));
        }
    }
    return(fl);
}

// lf(Y') = {x in Y | Q_check_minus(x,Y') = max Q_check_minus(x',Y') for x' in Y}
map<int, vector<double>> lastfirst(map<int, vector<double>> Yp, map<int, vector<double>> Y_all){
    map<int, vector<double>> lf;
    int length = Y_all.size();

    map<int, vector<double>> pts_left(Y_all);
    for(const auto& y : Yp){
        pts_left.erase(y.first);
    }

    // compute all Q(x,Y') and find the min
    vector<int> maxQ(length, INT_MAX);
    for(const auto& y : pts_left){
        vector<int> newQ = get_Q(y, Yp, Y_all, CHECK, MINUS);

        // this Q is equal to the old max -> add this point to lf
        if(revlex(newQ, maxQ) == 0){
            lf.emplace(y.first, Y_all.at(y.first));
        }

        // we have a new max -> clear out lf and add this point instead
        if(revlex(newQ, maxQ) > 0){
            maxQ = newQ;
            lf.clear();
            lf.emplace(y.first, Y_all.at(y.first));
        }
    }
    return(lf);
}
