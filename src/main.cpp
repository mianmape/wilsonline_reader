#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include "ipglasma.hpp"
#include <gsl/gsl_integration.h>
#include <array>
#include <cmath>

#define PI 3.14159265

using namespace std;

int LeviCivita(int i, int j, int k) {
    if ((i == j) || (i == k) || (j == k)) return 0;
    if ((i == 1 && j == 2 && k == 3) || (i == 2 && j == 3 && k == 1) || (i == 3 && j == 1 && k == 2)) return +1;
    if ((i == 3 && j == 2 && k == 1) || (i == 1 && j == 3 && k == 2) || (i == 2 && j == 1 && k == 3)) return -1;
    return 0;
}


std::complex<double> MultiplyElements(const std::complex<double>& a, const std::complex<double>& b) {
    return a * b;
}

/*
std::array<double, 2> RotateVector(double x, double y, double angle, std::array<double, 2>& cp) {
   std::array<double, 2> rotVec;
   double rad = angle * PI / 180.0;
   double x_rot = cos(rad) * x - sin(rad) * y + cp[0];
   double y_rot = sin(rad) * x + cos(rad) * y + cp[1];
   rotVec[0] = x_rot;
   rotVec[1] = y_rot;
   return rotVec;
}
*/

std::array<double, 6> CalcTriangleCoord(double d, std::array<double, 2>& cp) {
    double height = (sqrt(3) / 2) * d;
    std::array<double, 6> coord;
    coord[0] = cp[0];
    coord[1] = cp[1] + height / 2;
    coord[2] = cp[0] - d / 2;
    coord[3] = cp[1] - height / 2;
    coord[4] = cp[0] + d / 2;
    coord[5] = cp[1] - height / 2;

   return coord;
}


std::complex<double> ComputeBaryonOperator(const WilsonLine& w0, const WilsonLine& w1, const WilsonLine& w2) {
    std::complex<double> sum = {0.0, 0.0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                int eps_ijk = LeviCivita(i + 1, j + 1, k + 1);

                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        for (int n = 0; n < 3; n++) {
                            int eps_lmn = LeviCivita(l + 1, m + 1, n + 1);

                            auto v0 = w0.Element(i, m);
                            auto v1 = w1.Element(j, l);
                            auto v2 = w2.Element(k, n);

                            auto product1 = MultiplyElements(v0, v1);
                            auto product2 = MultiplyElements(product1, v2);

                            sum += static_cast<double>(eps_ijk * eps_lmn) * product2;
                       }
                    }
                }
            }
        }
    }
    return sum;
}


struct integrand_params {
    double x, y;
    IPGlasma* ipglasma;
};


// This function calculates the baryon operator
double baryonOperator(double x, double y, IPGlasma* ipglasma) {
    std::array<double, 2> cp = {x, y};

    // Tringle side length [1/GeV]
    double d = 0.4;

    std::array<double, 6> coord = CalcTriangleCoord(d, cp);

    WilsonLine w0 = ipglasma->GetWilsonLine(coord[0], coord[1]);
    WilsonLine w1 = ipglasma->GetWilsonLine(coord[2], coord[3]);
    WilsonLine w2 = ipglasma->GetWilsonLine(coord[4], coord[5]);

    std::complex<double> B123 = ComputeBaryonOperator(w0, w1, w2);

    double result = (1.0 + std::real(B123) / 6.0);
    return  2 * result;
}


// This computes the dipole operator
double dipoleOperator(double x, double y, IPGlasma* ipglasma) {
    std::array<double, 2> cp = {x, y};

    // Distance between the quarks [1/GeV]
    double d = 0.4;
    
    WilsonLine Vx = ipglasma->GetWilsonLine(x, y + d / 2.0);
    WilsonLine Vy = ipglasma->GetWilsonLine(x, y - d / 2.0);
    
    WilsonLine product = Vx.MultiplyByHermitianConjugate(Vy);
    std::complex<double> trace = product.Trace();
    
    return 2 * (1.0 - std::real(trace) / 3.0); 
}


double integrand(double y, void* p) {
    integrand_params* params = (integrand_params*) p;
    return dipoleOperator(params->x, y, params->ipglasma);
}


double integrand2(double x, void* p) {
    integrand_params* params = (integrand_params*) p;
    params->x = x;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &integrand;
    F.params = p;

    double result, error;
    gsl_integration_qags(&F, -24, 24, 0, 1e-2, 1000, w, &result, &error);

    return result;
}


double integrand_y_numerator(double y, void* p) {
    integrand_params* params = (integrand_params*) p;
    double x = params->x;
    double value = baryonOperator(x, y, params->ipglasma);
    return (x * x + y * y) / 2.0 * value;
}


double integrand_y_denominator(double y, void* p) {
    integrand_params* params = (integrand_params*) p;
    double x = params->x;
    return baryonOperator(x, y, params->ipglasma);
}


double integrate_y_numerator(double x, IPGlasma* ipglasma) {
    integrand_params params;
    params.x = x;
    params.ipglasma = ipglasma;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &integrand_y_numerator;
    F.params = &params;
    double result, error;
    gsl_integration_qags(&F, -24, 24, 0, 1e-2, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}


double integrate_y_denominator(double x, IPGlasma* ipglasma) {
    integrand_params params;
    params.x = x;
    params.ipglasma = ipglasma;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &integrand_y_denominator;
    F.params = &params;
    double result, error;
    gsl_integration_qags(&F, -24, 24, 0, 1e-2, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}


double integrand_x_numerator(double x, void* p) {
    integrand_params* params = (integrand_params*) p;
    IPGlasma* ipglasma = params->ipglasma;
    double inner_result = integrate_y_numerator(x, ipglasma);
    return inner_result;
}


double integrand_x_denominator(double x, void* p) {
    integrand_params* params = (integrand_params*) p;
    IPGlasma* ipglasma = params->ipglasma;
    double inner_result = integrate_y_denominator(x, ipglasma);
    return inner_result;
}


double integrate_outer(gsl_function* F) {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_integration_qags(F, -24, 24, 0, 1e-2, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}


// This function calculates the B-slope
double calculate_Bt0(IPGlasma* ipglasma) {
    integrand_params params;
    params.ipglasma = ipglasma;


    gsl_function F_numerator;
    F_numerator.function = &integrand_x_numerator;
    F_numerator.params = &params;
    double numerator = integrate_outer(&F_numerator);

    gsl_function F_denominator;
    F_denominator.function = &integrand_x_denominator;
    F_denominator.params = &params;
    double denominator = integrate_outer(&F_denominator);
    
    double Bt0 = numerator / denominator;
    return Bt0;
}
    

int main(int argc, char* argv[])
{
    // Arguments: filaname
    string fname = argv[1];
    IPGlasma ipglasma(fname, 0, BINARY);
   
    
    
    std::complex<double> N;
    /*
    std::array<double, 6> coord;
    std::array<double, 2> cp = {0.0, 0.0};
    */

    
    double min_d=0.01; double max_d=10; double factor=std::pow(10, 0.1);
    for (double d=min_d; d <= max_d; d *= factor)
    {

        std::array<double, 2> cp = {0.0, 0.0};
    
        WilsonLine Vx = ipglasma.GetWilsonLine(cp[0], cp[1] + d / 2.0);
        WilsonLine Vy = ipglasma.GetWilsonLine(cp[0], cp[1] - d / 2.0);
    
        WilsonLine product = Vx.MultiplyByHermitianConjugate(Vy);
        std::complex<double> trace = product.Trace();
    
        N = 1.0 - trace / 3.0;

        /*    
        coord = CalcTriangleCoord(d, cp, coord);

        WilsonLine w0 = ipglasma.GetWilsonLine(coord[0], coord[1]);
        WilsonLine w1 = ipglasma.GetWilsonLine(coord[2], coord[3]);
        WilsonLine w2 = ipglasma.GetWilsonLine(coord[4], coord[5]);

        B123 = ComputeBaryonOperator(w0, w1, w2);
        */

        std::cout << "d = " << fixed << setprecision(6) << d << setprecision(6) << " " << real(N) << " " << imag(N) << std::endl;
    }
    
    
    /* 
    // Calculate the cross section 
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    double result, error;

    integrand_params params; 
    params.ipglasma = &ipglasma;

    gsl_function F;
    F.function = &integrand2;
    F.params = &params;
    gsl_integration_qags(&F, -24, 24, 0, 1e-2, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    cout << "Result: " << result << " Error: " << error << endl;
    
    
    // Calculate the B-slope 
    double Bt0 = calculate_Bt0(&ipglasma);
    cout << "Result: " << Bt0 << endl;
    */
   
    return 0;
}   

