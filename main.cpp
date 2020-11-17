#include "vector"
#include <iostream>
#include <cmath>
using namespace std;

const int n = 3;
const double epsilon=0.0001;
//determinant value = 12.796

//functions for square root method
double determinant(double matrix[n][n], int m) {
        double result = 0;
        double temp;
        int* transposition = new int[n];
        for (int i = 0; i < n; ++i) {
            transposition[i] = i;
        }
    
        while (true) {
            temp = 1;
            for (int i = 0; i < n; ++i) {
                temp *= matrix[i][transposition[i]];
            }
            for (int i = 0; i < n - 1; ++i) {  // sign
                for (int x = i + 1; x < n; ++x) {
                    if (transposition[i] > transposition[x]) {
                        temp *= -1;
                    }
                }
            }
            result += temp;
    
            // next transposition
            int j = n - 2;
            while (j != -1 && transposition[j] > transposition[j + 1]) {
                --j;
            }
            if (j == -1) {
                break;
            }
            int k = n - 1;
            while (transposition[j] > transposition[k]) {
                --k;
            }
            swap(transposition[j], transposition[k]);
            int l = j + 1, r = n - 1;
            while (l < r) {
                swap(transposition[l++], transposition[r--]);
            }
        }
    
        delete[] transposition;
        return result;
}

void transponedMatrix(double matrix[n][n], double result[n][n]){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[j][i] = matrix[i][j];
        }
    }
}

void matrix_multiplication(double matrix[n][n], double another_one[n][n], double symmetric_matrix[n][n]){
  
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                symmetric_matrix[i][j] += matrix[i][k] * another_one[k][j];
            }
        }
    }
}

void make_matrix_s(double matrix[n][n], double s[n][n]){
     for (int i = 0; i < n; ++i) {
          {
              double sum = 0;
              for (int k = 0; k < i; ++k) {
                  sum += s[k][i] * s[k][i];
              }
              s[i][i] = sqrt(matrix[i][i] - sum);
          }
          for (int j = i + 1; j < n; ++j) {
              double sum = 0;
              for (int k = 0; k < i; ++k) {
                  sum += s[k][i] * s[k][j];
              }
              s[i][j] = (matrix[i][j] - sum) / s[i][i];
          }
      }
    
}

void calculating_y_and_x_vector(double s[n][n], double f[n], double y[n], double x[n]){
    y[0]=f[0]/s[0][0];
    for (int i=1;i<n;++i)
    {
        double sum=0;
        for (int k=0;k<i;++k)
        {
            sum+=s[k][i]*y[k];
        }
        y[i]=f[i]-sum;
        y[i]/=s[i][i];
    }
    x[n-1]=y[n-1]/s[n-1][n-1];
    for (int i = n-2;i>=0;--i)
    {
        double sum=0;
        for (int k=i+1; k<n;++k)
        {
            sum+=s[i][k]*x[k];
        }
        x[i]=y[i]-sum;
        x[i]/=s[i][i];
    }
}

void calculating_rest(double matrix[n][n], double x[n], double f[n], double rest[n]){
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            rest[i]+=matrix[i][j]*x[j];
        }
        rest[i]-=f[i];
    }
}
//

//functions for simple iterration methods
double matrixs_norm(double matrix[n][n]){
    double row_norm=0;
    for (int i=0;i<n;i++)
    {
        row_norm+=abs(matrix[0][i]);
    }
    double max=0;
    for (int i=1;i<n;i++)
    {
        max=0;
        for (int j=0;j<n;j++)
        {
            max+=abs(matrix[i][j]);
        }
        if (max>row_norm)
        {
            row_norm=max;
        }
    }
    return row_norm;
}

void matrix_B(double matrix[n][n], double alpha, double E[n][n], double B[n][n]){
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            B[i][j]=E[i][j]-alpha*matrix[i][j];
        }
    }
}

void vector_g(double f[n], double alpha, double g[n]){
    for (int i=0;i<n;i++)
    {
        g[i]=f[i]*alpha;
    }
}

void calculating_x(double B[n][n], double x[n], double Bx[n], double g[n]){
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            Bx[i]+=B[i][j]*x[j];
        }
    }
    for (int i=0;i<n;i++)
    {
        x[i]=Bx[i];
    }
    for (int i=0;i<n;i++)
    {
        x[i]+=g[i];
        
    }
}

void transposed_f(double transposed_matrix[n][n], double f[n], double transposed_f[n]){
    for (int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            transposed_f[i]+=transposed_matrix[i][j]*f[j];
        }
    }
    for (int i=0;i<n;i++)
    {
        f[i] = transposed_f[i];
    }
}
//

int main()
{
    double matrix[n][n] = {{1.7, -2.2, 3.0}, {2.1, 1.9, -2.3}, {4.2, 3.9, -3.1}};
    double f[n] = {1.8, 2.8, 5.1};
    double transponed_matrix[n][n];
    double symmetric_matrix[n][n];
    double y[n];
    double x[n];
    double rest_srm[n]={0,0,0};
    double s[n][n];
    double new_f[n];
    transponedMatrix(matrix, transponed_matrix);
    transposed_f(transponed_matrix, f, new_f);
    matrix_multiplication(transponed_matrix, matrix, symmetric_matrix);
    for (int i =0;i<n;i++){
    for (int j=0;j<n;j++)
    {
        s[i][j]=0;
    }
}
    double iteration_x[n]={0,0,0};
    double E[n][n] = {{1,0,0}, {0, 1, 0}, {0, 0, 1}};
    double g[n]={0,0,0};
    double d = determinant(matrix, n);
    double matrix_norm = matrixs_norm(symmetric_matrix);
    double alpha = (1.0)/matrix_norm;
    double iterations_amount=0;
    double B[n][n];
    for (int i =0;i<n;i++){
    for (int j=0;j<n;j++)
    {
        B[i][j]=0;
    }
}
    matrix_B(symmetric_matrix, alpha, E, B);
    
    if (d!=0){
    cout << "Determinant is not equal to 0" << endl;
}

    make_matrix_s(symmetric_matrix, s);
    calculating_y_and_x_vector(s, f, y, x);
    calculating_rest(symmetric_matrix, x, f, rest_srm);
    vector_g(f, alpha, g);
    cout << "symmetric matrix: " << endl;
    for (int i =0;i<n;i++){
    for (int j=0;j<n;j++)
    {
        cout << symmetric_matrix[i][j] << " " ;
    }
    cout << endl;
}
    cout << endl << "s:" << endl;
    for (int i =0;i<n;i++){
    for (int j=0;j<n;j++){
        cout << s[i][j] << " " ;
    }
    cout << endl;
}
    cout << endl << "Square root method" << endl;
    cout << endl << "x:" << endl;
    for (int i=0;i<n;i++){
    cout << "x"<<i+1<<"="<<x[i]<<" ";
}
    cout << endl;
    cout << endl << "rest vector:" << endl << "{ ";
    for (int i=0;i<n;i++){
        cout <<rest_srm[i]<<" ";
    }
    cout <<"}" << endl;
    
    bool is_accurate=true;;
    double iteration_rest[n]={0,0,0};
    while(true)
    {
        iterations_amount++;
        is_accurate=true;
        double temp[n]={0,0,0};
        calculating_x(B, iteration_x, temp, g);
        for (int i=0;i<n;i++)
        {
            iteration_rest[i]=0;
        }
        calculating_rest(symmetric_matrix, iteration_x, f, iteration_rest);
        for (int i=0;i<n;i++)
        {
            if (abs(iteration_rest[i])>epsilon)
            {
                is_accurate=false;
            }
        }
        if (is_accurate==true)
        {
            break;
        }
    }
    cout << "simple iteration method" << endl;
     cout << endl << "x:" << endl;
        for (int i=0;i<n;i++){
        cout << "x"<<i+1<<"="<<iteration_x[i]<<" ";
    }
        cout << endl;
        cout << endl << "rest vector:" << endl << "{ ";
        for (int i=0;i<n;i++){
            cout <<iteration_rest[i]<<" ";
        }
        cout <<"}" << endl << endl;
        cout <<"Iterations needed "<< iterations_amount << endl;
    return 0;
}
