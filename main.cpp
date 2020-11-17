#include <iostream>
#include <cmath>
#include <vector>
#include "fstream"
using namespace std;

// stuff
const int n=6;
const double left_edge = 1.1;
const double right_edge=4;
const double step = (right_edge-left_edge)/n;

//Basic stuff
double interpolated_func(double x){
return sin(x) - 3*(1/tan(1/x));
}

vector<double> funct_results(vector<double> nodes){
    vector<double> f_i;
    for (int i=0;i<nodes.size();i++)
    {
        f_i.push_back(interpolated_func(nodes[i]));
    }
    return f_i;
}

//Gauss' method of solving equations
vector<double> Gauss(vector<vector<double>> matrix) {
    double diagonal_element, flag;
    int row_num = (int)matrix.size();
    int col_num = (int)matrix[0].size();
    for (int i = 0; i < row_num; i++) {
        for (int j = 0; j < col_num; j++) {
            if (i == j) {
                int large_pivot_row = i;
                for (int k = i; k < row_num; ++k) {
                    if (matrix[k][j] > matrix[large_pivot_row][j]) {
                        large_pivot_row = k;
                    }
                }
                if (i != large_pivot_row) {
                    for (int k = 0; k < col_num; k++) {
                        swap(matrix[large_pivot_row][k], matrix[i][k]);
                    }
                }
                diagonal_element = matrix[i][j];
                for (int l = 0; l < col_num; l++) {
                    matrix[i][l] /= diagonal_element;
                }
                for (int k = 0; k < row_num; k++) {
                    flag = matrix[k][j];
                    if (k != i) {
                        for (int l = 0; l < col_num; l++) {
                            matrix[k][l] = (matrix[k][l]) - flag * (matrix[i][l]);
                        }
                    }

                }
            }
        }
    }
    vector<double> result;
    for (int i = 0; i < row_num; i++) {
        result.push_back(matrix[i][col_num - 1]);
    }
    return result;
}

//Functions for Lagrange polynomial
vector<double> Lagrange(int n, vector<double>nodes){
    vector<double> Lagrange_results;
    double multiplication_result=1;
    vector<double> Lagrange_pol;
    double result=0;
    for (int k=0;k<=nodes.size();k++)
    {
       result=0;
        for (int i=0;i<=n;i++)
        {
            
                   
            for (int j=0;j<=n;j++)
            {
                
                if (j!=i)
                {
                    multiplication_result*=(nodes[k]-nodes[j])/(nodes[i]-nodes[j]);
                }
            }
            Lagrange_pol.push_back(interpolated_func(nodes[i])*multiplication_result);
            multiplication_result=1;
        }
        for (int i=0;i<=n;i++)
        {
            result+=Lagrange_pol[i];
        }
        Lagrange_pol.clear();
        Lagrange_results.push_back(result);
    }
    return Lagrange_results;
}

vector<double> Lagrange_Rest(vector<double>f_i, vector<double>Li){
    vector<double>rest;
    for(int i=0;i<=f_i.size();i++)
    {
        rest.push_back(f_i[i]-Li[i]);
    }
    return rest;
}

//Functions for Newton polynomial
vector<double> Newton_polynomial(vector<double> nodes){
    vector<double>Pn;
    double sum=0;
    double difference=0;
    double divisor=0;
    for(int k=0;k<=nodes.size();k++)
    {
        sum = interpolated_func(nodes[0]);
        for(int i=1;i<=n;i++)
        {
            difference=0;
            for (int j=0;j<=i;j++)
            {
                divisor=1;
                for(int m=0;m<=i;m++)
                {
                    if(m!=j)
                    {
                        divisor*=(nodes[j]-nodes[m]);
                    }
                }
                difference+=interpolated_func(nodes[j])/divisor;
            }
            for (int j=0;j<i;j++)
            {
                difference*=(nodes[k]-nodes[j]);
            }
            sum+=difference;
        }
        Pn.push_back(sum);
    }
    return Pn;
}

vector<double> Newton_rest(vector<double> f_i, vector<double> Pn){
    vector<double> rest;
    for (int i=0;i<f_i.size();i++)
    {
        rest.push_back(f_i[i]-Pn[i]);
    }
    return rest;
}

//Coefficients for Smallest Squares Method
vector<vector<double>> Cj_coeffitients(vector<double> nodes){
    vector<vector<double>>Cj;
    double temp=0;
    for(int i=0;i<=n;i++)
    {
        Cj.push_back(vector<double>());
        for (int j=0;j<=n;j++)
        {
            temp=0;
            for (int k=0;k<=n;k++)
            {
                temp+=pow(nodes[k], j+i);
            }
            Cj[i].push_back(temp);
        }
    }
    return Cj;
}

vector<double> dk_value(vector<double> nodes, vector<double>func_resutls){
    vector<double> dk_values;
    double dk=0;
    for (int k=0;k<=n;k++)
    {
        dk=0;
        for(int i=0;i<=n;i++)
        {
            dk+=func_resutls[i]*pow(nodes[i], k);
        }
        dk_values.push_back(dk);
    }
    return dk_values;
}

vector<double> smallest_squares_polynomial(vector<double> nodes, vector<double>coefficients){
    vector<double> F;
    double sum=0;
    for(int i=0;i<=nodes.size();i++)
    {
        sum=0;
        for(int j=0;j<=n;j++)
        {
            sum+=coefficients[j]*pow(nodes[i], j);
        }
        F.push_back(sum);
        
    }
    return F;
}

vector<double> smallest_squares_rest(vector<double>f_i, vector<double>F_i){
    vector<double> rest;
    for(int i=0;i<=f_i.size();i++)
    {
        rest.push_back(f_i[i]-F_i[i]);
    }

    return rest;
}

//Nodes to extrtapolate and interpolate between nodes
vector<double>& widen_nodes(vector<double>& nodes){
    for(int j=0;j<=n+1;j++)
    {
        nodes.push_back(left_edge+(2*j-1)*step/2);
    }
    return nodes;
}

//Functions for sweep method
vector<double> A_coeffs(vector<double> a, vector<double> b, vector<double> c, vector<double> d){
    vector<double> A_i;
    A_i.push_back(-c[0]/b[0]);
    double p=0;
    for(int i=1;i<=n-2;i++)
    {
        p=(a[i]*A_i[i-1] + b[i]);
        A_i.push_back(-c[i]/p);
        p=0;
    }
    return A_i;
};

vector<double> B_coeffs(vector<double> a, vector<double> b, vector<double> c, vector<double> d, vector<double> A_i){
    vector<double> B_i;
    B_i.push_back(d[0]/b[0]);
    double p=0;
    for(int i=1;i<n-1;i++)
    {
        p=(a[i]*A_i[i-1] + b[i]);
        B_i.push_back((d[i]-a[i]*B_i[i-1])/p);
        p=0;
    }
    return B_i;
};

vector<double> C_coeffs(vector<double> a, vector<double> b, vector<double> d, vector<double> A, vector<double> B){
    vector<double> x(n+1,0);
    x[x.size()-2] = (d[d.size()-1]-a[a.size()-1]*B[B.size()-2])/(a[a.size()-1]*A[A.size()-2] + b[b.size()-1]);
    for (int i= n-2;i>=1;i--)
    {
        x[i]=A[i]*x[i+1]+B[i];
    }
    return x;
}

//Spline
double Spline_values(vector<double> a, vector<double> b, vector<double> c, vector<double> d, vector<double> nodes, double x){
    double spline_i = 0;
    int k=0;
    for (int i=1;i<=n;i++)
    {
        if (x>=nodes[i-1] && x<=nodes[i])
        {
            k=i;
            break;
        }
        else if(x<left_edge)
        {
            k=1;
            break;
        }
        else if(x>right_edge)
        {
            k=n;
            break;
        }
        else
        {
            continue;
        }
    }
    spline_i = a[k] + b[k]*(x-nodes[k]) + c[k]*(x-nodes[k])*(x-nodes[k])/2 + d[k]*(x-nodes[k])*(x-nodes[k])*(x-nodes[k])/6;
    return spline_i;
}

int main()
{
//Creating vectors;
    vector<double> interpolation_nodes;
    vector<double> function_values;
    
    //Newton,Lagrange,SQM
    vector<double> d_ks;
    vector<vector<double>> cj_s;
    vector<double> ai_s;
    vector<double> lagrange_polynomial;
    vector<double> lagrange_rest;
    vector<double> divided_differences_vector;
    vector<double> newton_polynomial;
    vector<double> newton_rest;
    vector<double> smallest_squares;
    vector<double> smallest_squares_rest_polynomial;
    
    //Splines
    vector<double> spline_values;
    vector<double> spline_a(n+1,0);
    vector<double> spline_b(n+1,0);
    vector<double> spline_c;
    vector<double> spline_d(n+1,0);
    vector<double> a_i;
    vector<double> b_i;
    vector<double> c_i;
    vector<double> d_i;
    vector<double> A;
    vector<double> B;
//
    //Getting data
    for (int i=0;i<=n;i++){
        interpolation_nodes.push_back(left_edge+i*step);
    }
    widen_nodes(interpolation_nodes);
    function_values=funct_results(interpolation_nodes);
    
    //Sweeping input data
    for (int i=1;i<=n-1;i++){
        b_i.push_back(4*step);
    }
    a_i.push_back(0);
    for(int i=2;i<=n-1;i++){
        a_i.push_back(step);
    }
    for(int i=0;i<n-2;i++){
        c_i.push_back(step);
    }
    c_i.push_back(0);
    for (int i=1;i<=n-1;i++){
        d_i.push_back(6*(function_values[i+1]-2*function_values[i]+function_values[i-1])/step);
    }
    
    //Sweeping coefficients
    A = A_coeffs(a_i, b_i, c_i, d_i);
    B = B_coeffs(a_i, b_i, c_i, d_i, A);
    
    //Spline coefficients
    spline_c = C_coeffs(a_i, b_i, d_i, A, B);
    for (int i=1;i<=n;i++){
        spline_d[i]=(spline_c[i]-spline_c[i-1])/step;
    }
    for (int i=1;i<=n;i++){
        spline_b[i]=(step*spline_c[i]/2)-(step*step)*spline_d[i]/6 + (function_values[i]-function_values[i-1])/step;
    }
    for (int i=0;i<=n;i++){
        spline_a[i]=function_values[i];
    }
    for (int i=0;i<=interpolation_nodes.size();i++){
        spline_values.push_back(Spline_values(spline_a, spline_b, spline_c, spline_d, interpolation_nodes, interpolation_nodes[i]));
    }
    
    //Getting results on screen
    cout  << "Table of spline coefficients: "<< endl;
    for (int i=1;i<=n;i++){
        printf("i = %d| ai =%10.6f| bi =%10.6f| ci = %10.6f| di = %10.6f| \n", i, spline_a[i], spline_b[i], spline_c[i], spline_d[i]);
    }
    cout << "__ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __" <<endl;
    cout << "Spline - function table" << endl;
    for (int i=0;i<=interpolation_nodes.size()-1;i++){
        if (i==n+1)
        {
            cout << "                            x~" << endl;
        }
        printf("x%2i = %11.6f|f(x) = %10.6f| S(x) = %10.6f| Rn(x) = %10.6f| \n",i, interpolation_nodes[i],function_values[i], spline_values[i], function_values[i]-spline_values[i]);
    }
    cout << "__ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __" <<endl;
    
    //Calling functions
    lagrange_polynomial = Lagrange(n, interpolation_nodes);
    lagrange_rest=Lagrange_Rest(function_values, lagrange_polynomial);
    newton_polynomial = Newton_polynomial(interpolation_nodes);
    newton_rest=Newton_rest(funct_results(interpolation_nodes), newton_polynomial);
    d_ks=dk_value(interpolation_nodes, function_values);
    cj_s=Cj_coeffitients(interpolation_nodes);
    
        //Getting a coefficients via Gauss' method
    for (int i=0;i<cj_s.size();i++)
    {
        cj_s[i].push_back(d_ks[i]);
    }
    ai_s=Gauss(cj_s);
        //
    smallest_squares=smallest_squares_polynomial(interpolation_nodes, ai_s);
    smallest_squares_rest_polynomial=smallest_squares_rest(function_values, smallest_squares);
    
    //Printing values
    cout << "Polynomial table: " << endl;
    for (int i=0;i<=interpolation_nodes.size()-1;i++)
    {
        if (i==n+1)
        {
              cout << "                            x~" << endl;
        }
        printf("x%2i =%10.4f| f(x) = %10.6f| S(x) = %10.6f| Ln(x) = %10.6f| Pn(x) = %10.6f| F(x) = %10.6f| \n",i, interpolation_nodes[i],function_values[i], spline_values[i], lagrange_polynomial[i], newton_polynomial[i], smallest_squares[i]);
    }
    cout <<endl << "__ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __" << endl << "Rest table: " << endl;
    cout.precision(9);
    for (int i=0;i<=interpolation_nodes.size()-1;i++)
    {
        if (i==n+1)
        {
              cout << "                            x~" << endl;
        }
        printf("x%2i = %10.6f| S(x)_rest = %10.6f| Ln(x)_rest = %10.6f| Pn(x)_rest = %10.6f| F(x)_rest =%10.6f| \n",i, interpolation_nodes[i], function_values[i]-spline_values[i], lagrange_rest[i], newton_rest[i], smallest_squares_rest_polynomial[i]);
    }
    //And the curtains fall
    return 0;
}
