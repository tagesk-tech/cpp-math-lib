int main() {
    // Example 1: 2x + y – z =  8
    //           –3x – y + 2z = –11
    //           –2x + y + 2z = –3
    vector<vector<double>> A1 = {
        { 50,  89, -45},
        {-345, 45,  7},
        {-23433,  89,  9}
    };
    vector<double> b1 = {8, -11, -3};

    // Example 2: x + 2y + 3z = 14
    //            2x + 5y + 2z = 18
    //            3x + 2y + 4z = 20
    vector<vector<double>> A2 = {
        {1, 2, 3},
        {2, 5, 2},
        {3, 2, 4}
    };
    vector<double> b2 = {14, 18, 20};

    try {
        auto x1 = gauss_jordan(A1, b1);
        cout << "Solution for system 1:\n";
        for (size_t i = 0; i < x1.size(); ++i)
            cout << "  x[" << i << "] = " << x1[i] << "\n";

        auto x2 = gauss_jordan(A2, b2);
        cout << "\nSolution for system 2:\n";
        for (size_t i = 0; i < x2.size(); ++i)
            cout << "  x[" << i << "] = " << x2[i] << "\n";
    }
    catch (exception& e) {
        cerr << "Error: " << e.what() << "\n";
    }

    return 0;