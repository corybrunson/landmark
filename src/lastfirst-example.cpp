#include "lastfirst-example.h"

//' @rdname run_lf_example
//' @export
// [[Rcpp::export]]
void run_lf_example(const NumericMatrix& y, int nhd_size){
    map<int, vector<double>> Y_all; // whole space Y
    map<int, vector<double>> landmarks; // landmark set L
    int length = y.nrow();

    // convert Y to cpp data structure
    for(int i = 0; i < length; i++){
        NumericVector vec = y.row(i);
        vector<double> point(vec.begin(),vec.end());
        Y_all.emplace(i, point);
    }
    map<int, vector<double>> pts_left(Y_all); // store indices and values of Y\L

    // print out all point in Y
    cout << "\nY";
    prettyprint_pts(Y_all);

    // first compute all q's for reference
    cout << "q check:\n--------\n";
    prettyprint_q(Y_all, CHECK);

    // choose l0 and add it to L
    pair<int, vector<double>> l_0(0, Y_all.at(0)); // for simplicity, choose the first pt
    landmarks.insert(l_0);
    pts_left.erase(l_0.first); // remove l0 from list of points left

    // keep track of covered points
    map<int, vector<double>> Nk = Nk_check_plus(l_0, Y_all, nhd_size, Y_all);
    map<int, vector<double>> covered(Nk); // store union of Nk check plus over all landmarks

    int i = 0;
    pair<int, vector<double>> l_i = l_0;
    while(true){
        cout << "------------------------------ i = " << i << " ------------------------------";
        // print landmark set
        cout << "\nL";
        prettyprint_pts(landmarks);

        // update the list of covered points
        map<int, vector<double>> Nk = Nk_check_plus(l_i, Y_all, nhd_size, Y_all);
        for(const auto& x : Nk){ covered.insert(x); }
        if(covered.size() >= length){break;}

        // compute Nk's + Q's for each point in Y\L
        for(const auto& y : pts_left){
            cout << "\ny = " << y.first << "\n";

            cout << "Nk check minus:\n k";
            for(int k = 1; k <= length; k++){
                prettyprint_Nk(y, landmarks, k, Y_all, CHECK, MINUS);
            }

            cout << "\nQ check minus: ";
            prettyprint_Q(y, landmarks, Y_all, CHECK, MINUS);
        }

        // compute lf(L)
        cout << "\nlf(L):";
        map<int, vector<double>> lf = lastfirst(landmarks, Y_all);
        prettyprint_pts(lf);

        // choose li from lf(L) and add it to L
        l_i = make_pair(lf.begin()->first, lf.begin()->second);
        landmarks.insert(l_i);
        pts_left.erase(l_i.first);

        i++;
    }

    cout << "\n\n\n-------------------------------------------------------------------\n";
    cout << "FINAL LANDMARK SET";
    cout << "\n-------------------------------------------------------------------\n";
    prettyprint_pts(landmarks);
}

void prettyprint_pts(map<int, vector<double>> s){
    for (const auto& x : s) {
        cout << "\n" << x.first << ": ( ";
        for(const auto& v : x.second){cout << v << " ";}
        cout << ")";
    }
    cout << "\n\n";
}

void prettyprint_q(map<int, vector<double>> Y_all, Type t){
    // table header
    for(const auto& x1 : Y_all){ cout << "\t" << x1.first; }
    cout << "\n";

    // table body
    for(const auto& x1 : Y_all){
        cout << "\n" << x1.first << ":";
        for(const auto& x2 : Y_all){
            int q;
            switch(t){
                case HAT    : q = q_hat(x1.second, x2.second, Y_all); break;
                case CHECK  : q = q_check(x1.second, x2.second, Y_all); break;
            }
            cout << "\t" << q;
        }
    }
    cout << "\n\n\n";
}

void prettyprint_Nk(pair<int, vector<double>> x, map<int, vector<double>> Yp, int k, map<int, vector<double>> Y_all, Type t, Sign s){
    map<int, vector<double>> Nk;
    if(t == HAT && s == PLUS){ Nk = Nk_hat_plus(x, Yp, k, Y_all); }
    if(t == CHECK && s == MINUS){ Nk = Nk_check_minus(x, Yp, k, Y_all); }

    cout << "\n " << k;
    for(const auto& pt : Nk){
        cout << "\t" << pt.first << " : ( ";
        for(const auto& v : pt.second){
            cout << v << " ";
        }
        cout << ")";
    }
}

void prettyprint_Q(pair<int, vector<double>> x, map<int, vector<double>> Yp, map<int, vector<double>> Y_all, Type t, Sign s){
    vector<int> q;
    if(t == HAT && s == PLUS){ q = Q_hat_plus(x, Yp, Y_all); }
    if(t == CHECK && s == MINUS){ q = Q_check_minus(x, Yp, Y_all); }

    cout << "( ";
    for(const auto& v : q){cout << v << " ";}
    cout << ")\n";
}
