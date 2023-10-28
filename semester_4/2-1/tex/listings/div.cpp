Vec div(const Vec& v1, const Vec& v2, double eps){
    if (v1.size() != v2.size()) throw invalid_argument("Vectors have to be of an equal length");
    Vec res;
    for (size_t i = 0; i < v1.size(); i++){
        if (fabs(v2[i]) > eps)
            res.push_back(v1[i]/v2[i]);
        else
            cerr << "Denomenator values that are less than set threshold, are ignored\n";
    }        
    return res;
}
