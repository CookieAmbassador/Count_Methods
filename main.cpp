#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

const int n=6;
const double left_edge = 1.1;
const double right_edge=4;
const double step = (right_edge-left_edge)/n;
double original_value  = -19.76016;

double interpolated_func(double x)
{
    return sin(x) - 3*(1/tan(1/x));
}
double checking_function_0(double x)
{
    return pow(x,0);
}
double checking_function_1(double x)
{
    return x;
}
double checking_function_2(double x)
{
    return pow(x,2);
}
double checking_function_3(double x)
{
    return pow(x,3);
}
double checking_function_4(double x)
{
    return pow(x,4);
}
double checking_function_2np1(double x)
{
    return pow(x,2*n+1);
}
double checking_function_2np2(double x)
{
    return pow(x,2*n+2);
}
vector<double> funct_results(vector<double> nodes)
{
    vector<double> f_i;
    for (int i=0;i<nodes.size();i++)
    {
            f_i.push_back(interpolated_func(nodes[i]));
    }
    return f_i;
}

double left_rectangles(vector<double> nodes, double (*function)(double))
{
    double integral_of_function = 0;
    for (int k=0;k<n;k++)
    {
        integral_of_function+=function(nodes[k]);
    }
    integral_of_function*=step;
    return integral_of_function;
}

double right_rectangles(vector<double> nodes,double (*function)(double))
{
    double integral_of_function = 0;
    for (int k=1;k<n+1;k++)
    {
        integral_of_function+=function(nodes[k]);
    }
    integral_of_function*=step;
    return integral_of_function;
}

double middle_rectangles(vector<double> nodes, double (*function)(double))
{
    double integral_of_function = 0;
    for (int k=0;k<n;k++)
    {
        integral_of_function+=function((nodes[k]+nodes[k+1])/2);
    }
    integral_of_function*=step;
    return integral_of_function;
}

double trapezoids(vector<double> nodes, double (*function)(double))
{
    double first_addition=step*(function(nodes[0]) + function(nodes[n]))/2;
    double second_addition=0;
    for(int k=1;k<n;k++)
    {
        second_addition+=function(nodes[k]);
    }
    second_addition*=step;
    return first_addition+second_addition;
}

double simpson(vector<double> nodes, double (*function)(double))
{
    double integral_of_function=0;
    for(int k=0;k<n+1;k++)
    {
        if(k==0)
        {
            integral_of_function+=function(nodes[k]);
        }
        else if(k==n)
        {
            integral_of_function+=function(nodes[k]);
        }
        else if(k%2==1)
        {
            integral_of_function+=4*function(nodes[k]);
        }
        else
        {
            integral_of_function+=2*function(nodes[k]);
        }
    }
    integral_of_function*=step/3;
    return integral_of_function;
}

double Gauss(vector<double>A_k, vector<double>t_k, double (*function)(double))
{
    double integral_of_function=0;
    double x=0;
    for(int k=0;k<n+1;k++)
    {
        x=(left_edge+right_edge)/2 + t_k[k]*(right_edge-left_edge)/2;
        integral_of_function+=A_k[k]*function(x);
    }
    integral_of_function*=(right_edge-left_edge)/2;
    return integral_of_function;
}

int main()
{
    vector<double> interpolation_nodes;
    vector<double> function_values;
    double left_rectangles_formulae=0;
    double right_rectangles_formulae=0;
    double middle_rectangles_formulae=0;
    double trapezoids_formulae=0;
    double simpsons_formulae=0;
    double gauss_formulae=0;
    vector<double> t_k={0,-0.40584515137739716690,0.40584515137739716690,-0.74153118559939443986,0.74153118559939443986,-0.94910791234275852452, 0.94910791234275852452};
    vector<double> A_k={0.41795918367346938775, 0.38183005050511894495,0.38183005050511894495,0.27970539148927666790,0.27970539148927666790,0.12948496616886969327, 0.12948496616886969327};
    for (int i=0;i<=n;i++)
    {
         interpolation_nodes.push_back(left_edge+i*step);
    }
    cout << "interpolation nodes are: {";
    for (int i=0;i<=n;i++)
    {
        cout << interpolation_nodes[i];
        if (i!=n)
        {
            cout << ", ";
        }
    }
    cout << "}" << endl << endl;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,interpolated_func);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,interpolated_func);
    middle_rectangles_formulae= middle_rectangles( interpolation_nodes, interpolated_func);
    trapezoids_formulae=trapezoids(interpolation_nodes,interpolated_func);
    simpsons_formulae=simpson(interpolation_nodes,interpolated_func);
    gauss_formulae=Gauss(A_k, t_k,interpolated_func);
    cout << endl << endl<<"Integral of f(x)=sin(x)-3ctg(1/x) values:"<< endl;
    cout << "Original_value: " << original_value << endl;
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-trapezoids_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl<<"Integral of x^N with following N values:" << endl;
    cout << "__ __N=0__ __"<<endl;
    original_value=2.9;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,checking_function_0);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,checking_function_0);
    middle_rectangles_formulae= middle_rectangles(interpolation_nodes,checking_function_0);
    trapezoids_formulae=trapezoids(interpolation_nodes,checking_function_0);
    simpsons_formulae=simpson(interpolation_nodes,checking_function_0);
    gauss_formulae=Gauss(A_k, t_k,checking_function_0);
    cout << "Original_value: " << original_value << endl;
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-simpsons_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl;
    cout << "__ __N=1__ __"<<endl;
    original_value=7.395;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,checking_function_1);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,checking_function_1);
    middle_rectangles_formulae= middle_rectangles(interpolation_nodes,checking_function_1);
    trapezoids_formulae=trapezoids(interpolation_nodes,checking_function_1);
    simpsons_formulae=simpson(interpolation_nodes,checking_function_1);
    gauss_formulae=Gauss(A_k, t_k,checking_function_1);
    cout << "Original_value: " << original_value << endl;
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-simpsons_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl;
    cout << "__ __N=2__ __"<<endl;
    original_value=20.88966666666667;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,checking_function_2);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,checking_function_2);
    middle_rectangles_formulae= middle_rectangles(interpolation_nodes,checking_function_2);
    trapezoids_formulae=trapezoids(interpolation_nodes,checking_function_2);
    simpsons_formulae=simpson(interpolation_nodes,checking_function_2);
    gauss_formulae=Gauss(A_k, t_k,checking_function_2);
    cout << "Original_value: " << original_value << endl;
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-simpsons_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl;
    cout << "__ __N=3__ __"<<endl;
    original_value=63.633975;
    cout << "Original_value: " << original_value << endl;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,checking_function_3);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,checking_function_3);
    middle_rectangles_formulae= middle_rectangles(interpolation_nodes,checking_function_3);
    trapezoids_formulae=trapezoids(interpolation_nodes,checking_function_3);
    simpsons_formulae=simpson(interpolation_nodes,checking_function_3);
    gauss_formulae=Gauss(A_k, t_k,checking_function_3);
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-simpsons_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl;
    cout << "__ __N=4__ __"<<endl;
    original_value=204.477898;
    cout << "Original_value: " << original_value << endl;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,checking_function_4);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,checking_function_4);
    middle_rectangles_formulae= middle_rectangles(interpolation_nodes,checking_function_4);
    trapezoids_formulae=trapezoids(interpolation_nodes,checking_function_4);
    simpsons_formulae=simpson(interpolation_nodes,checking_function_4);
    gauss_formulae=Gauss(A_k, t_k,checking_function_4);
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-simpsons_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl;
    cout << "__ __N=2N+1__ __"<<endl;
    original_value=1.917396087160726*pow(10,7);
    cout << "Original_value: " << original_value << endl;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,checking_function_2np1);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,checking_function_2np1);
    middle_rectangles_formulae= middle_rectangles(interpolation_nodes,checking_function_2np1);
    trapezoids_formulae=trapezoids(interpolation_nodes,checking_function_2np1);
    simpsons_formulae=simpson(interpolation_nodes,checking_function_2np1);
    gauss_formulae=Gauss(A_k, t_k,checking_function_2np1);
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-simpsons_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl;
    cout << "__ __N=2N+2__ __"<<endl;
    original_value=7.158278798818346*pow(10,7);
    cout << "Original_value: " << original_value << endl;
    left_rectangles_formulae= left_rectangles(interpolation_nodes,checking_function_2np2);
    right_rectangles_formulae= right_rectangles(interpolation_nodes,checking_function_2np2);
    middle_rectangles_formulae= middle_rectangles(interpolation_nodes,checking_function_2np2);
    trapezoids_formulae=trapezoids(interpolation_nodes,checking_function_2np2);
    simpsons_formulae=simpson(interpolation_nodes,checking_function_2np2);
    gauss_formulae=Gauss(A_k, t_k,checking_function_2np2);
    printf("Left rectangles formulae: %.7g, and it's error is: %.6g \n", left_rectangles_formulae, -(original_value-left_rectangles_formulae));
    printf("Right rectangles formulae: %.7g, and it's error is: %.6g\n", right_rectangles_formulae,-(original_value-right_rectangles_formulae));
    printf("Middle rectangles formulae: %.7g, and it's error is: %.6g\n", middle_rectangles_formulae,-(original_value-middle_rectangles_formulae));
    printf("Trapezoids formulae: %.7g, and it's error is: %.6g\n", trapezoids_formulae,-(original_value-trapezoids_formulae));
    printf("Simpson's formulae: %.7g, and it's error is: %.6g\n", simpsons_formulae, -(original_value-simpsons_formulae));
    printf("Gauss' formulae: %.7g, and it's error is: %.6g\n", gauss_formulae,-(original_value-gauss_formulae));
    cout << endl;
    return 0;
}
