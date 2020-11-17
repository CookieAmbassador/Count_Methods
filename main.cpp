#include <iostream>
#include <vector>
using namespace std;

vector<double> A_coeffs(vector<double> a, vector<double> b, vector<double> c, vector<double> d)
{
    vector<double> A_i;
    A_i.push_back(-c[0]/b[0]);
    for(int i=1;i<5;i++)
    {
        A_i.push_back(-c[i]/(a[i]*A_i[i-1] + b[i]));
    }
    return A_i;
};

vector<double> B_coeffs(vector<double> a, vector<double> b, vector<double> c, vector<double> d, vector<double> A_i)
{
    vector<double> B_i;
    B_i.push_back(d[0]/b[0]);
    for(int i=1;i<5;i++)
    {
        B_i.push_back((d[i]-a[i]*B_i[i-1])/(a[i]*A_i[i-1] + b[i]));
    }
    return B_i;
};


 int main()
{
    int n=5;
    double matrix[5][5] = {
        {2, -1, 0, 0, 0},
        {-3, 8, -1, 0, 0},
        {0, -5, 12, 2, 0},
        {0, 0, -6, 18, -4},
        {0, 0, 0, -5, 10}
    };
    vector<double> d = {-25, 72, -69, -156, 20};
    vector<double> a_i;
    vector<double> b_i;
    vector<double> c_i;
    for (int i=0;i<5;i++)
    {
        b_i.push_back(matrix[i][i]);
    }
    a_i.push_back(0);
    for (int i=1;i<n;i++)
    {
        a_i.push_back(matrix[i][i-1]);
    }
    for (int i=0;i<n-1;i++)
    {
        c_i.push_back(matrix[i][i+1]);
    }
    c_i.push_back(0);
    
    vector<double> A = A_coeffs(a_i, b_i, c_i, d);
    vector<double> B = B_coeffs(a_i, b_i, c_i, d, A);
    vector<double> x;
    x.push_back((d[n-1]-a_i[n-1]*B[n-2])/(a_i[n-1]*A[n-2] + b_i[n-1]));
    for (int i=n-2;i>-1;--i)
    {
        x.push_back(A[i]*x[n-i-2] + B[i]);
    }
    reverse(x.begin(), x.end());
    cout << "Ai: ";
    for(int i=0;i<A.size();i++)
    {
        cout <<"A["<<(i+1)<<"] "; printf("%.5f", A[i]); cout<<" ";
    }
    cout << endl;
    cout << "Bi: ";
    for(int i=0;i<B.size();i++)
    {
        cout <<"B["<<(i+1)<<"] "; printf("%.5f", B[i]); cout<<" ";
    }
    cout << endl;
    cout << "X: ";
    for(int i=0;i<x.size();i++)
       {
           cout << "x["<<(i+1)<<"] "; printf("%.5f", x[i]); cout<<" ";
       }
    cout << endl;
    vector<double> miss_value;
    double sum;
    for(int i=0;i<x.size();i++)
    {
        sum=0;
        for (int j=0;j<x.size();j++)
        {
            sum+=matrix[i][j]*x[j];
        }
        miss_value.push_back(sum);
    }
    for (int i=0;i<n;i++)
    {
        miss_value[i]-=d[i];
    }
    cout << "Компоненти вектору нев'язки: {";
    for (int i=0;i<n;i++)
    {
        printf("%.5f", miss_value[i]);
        if (i!!=miss_value.size()-1)
        {
            cout << ", ";
        }
    }
    cout <<"}" << endl<< "- без округлення - дуже мале число"<< endl;
    return 0;
}
