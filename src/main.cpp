#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <complex>
#include <cmath>
#include "ipglasma.hpp"

#define PI 3.141592653589793
using namespace std;

int LeviCivita(int i, int j, int k) {
    if ((i == j) || (i == k) || (j == k)) return 0;
    if ((i == 1 && j == 2 && k == 3) || (i == 2 && j == 3 && k == 1) || (i == 3 && j == 1 && k == 2)) return +1;
    if ((i == 3 && j == 2 && k == 1) || (i == 1 && j == 3 && k == 2) || (i == 2 && j == 1 && k == 3)) return -1;
    return 0;
}


// Calculate the equilateral tirangle coordinates
array<double, 6> CalcTriangleCoord(double d, const array<double, 2>& cp) {
    double h = sqrt(3.0) / 2.0 * d;

    // Vertices relative to centroid
    array<double, 6> coord;
    coord[0] = cp[0];           coord[1] = cp[1] + 2.0 * h / 3.0; 
    coord[2] = cp[0] - d / 2.0; coord[3] = cp[1] - h / 3.0;       
    coord[4] = cp[0] + d / 2.0; coord[5] = cp[1] - h / 3.0;       
    return coord;
}

// Multiply two complex numbers
inline complex<double> MultiplyElements(const complex<double>& a, const complex<double>& b) {
    return a * b;
}

// Compute baryon Wilson line operator
complex<double> ComputeBaryonOperator(const WilsonLine& w0, const WilsonLine& w1, const WilsonLine& w2) {
    complex<double> sum = {0.0, 0.0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                int eps_ijk = LeviCivita(i + 1, j + 1, k + 1);
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        for (int n = 0; n < 3; n++) {
                            int eps_lmn = LeviCivita(l + 1, m + 1, n + 1);
                            complex<double> prod =
                                MultiplyElements(w0.Element(i, l),
                                MultiplyElements(w1.Element(j, m), w2.Element(k, n)));
                            sum += static_cast<double>(eps_ijk * eps_lmn) * prod;
                        }
                    }
                }
            }
        }
    }
    return sum;
}

// Compute baryon the amplitude 
double BaryonAmplitude(double x, double y, IPGlasma& ipglasma) {
    const double l = 6.5; // triangle side length
    array<double, 2> cp = {x, y};
    array<double, 6> coord = CalcTriangleCoord(l, cp);

    WilsonLine w0 = ipglasma.GetWilsonLine(coord[0], coord[1]);
    WilsonLine w1 = ipglasma.GetWilsonLine(coord[2], coord[3]);
    WilsonLine w2 = ipglasma.GetWilsonLine(coord[4], coord[5]);

    complex<double> B123 = ComputeBaryonOperator(w0, w1, w2);
    complex<double> amplitude = 1.0 - B123 / 6.0; 
    return std::real(amplitude); 
}


// Dipole amplitude
double dipoleOperator(double x, double y, IPGlasma* ipglasma) {
    std::array<double, 2> cp = {x, y};
    double d = 0.4; // distance between the quarks
    WilsonLine Vx = ipglasma->GetWilsonLine(x, y + d / 2.0);
    WilsonLine Vy = ipglasma->GetWilsonLine(x, y - d / 2.0);
    WilsonLine product = Vx.MultiplyByHermitianConjugate(Vy);
    std::complex<double> trace = product.Trace();
    return 2.0 * (1.0 - std::real(trace) / 3.0);
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <ipglasma-file>" << endl;
        return 1;
    }

    string fname = argv[1];
    IPGlasma ipglasma(fname, 0, BINARY);
    
    /*
    double d = 0.4;
    double min_d=0.01; double max_d=10; double factor=std::pow(10, 0.1);
    for (double d=min_d; d <= max_d; d *= factor)
    {

        std::array<double, 2> cp = {0.0, 0.0};
        double d = 0.4;
        WilsonLine Vx = ipglasma.GetWilsonLine(cp[0], cp[1] + d / 2.0);
        WilsonLine Vy = ipglasma.GetWilsonLine(cp[0], cp[1] - d / 2.0);

        WilsonLine product = Vx.MultiplyByHermitianConjugate(Vy);
        std::complex<double> trace = product.Trace();

        N = 1.0 - trace / 3.0;

            
        double l = 0.4;
        coord = CalcTriangleCoord(l, cp, coord);

        WilsonLine w0 = ipglasma.GetWilsonLine(coord[0], coord[1]);
        WilsonLine w1 = ipglasma.GetWilsonLine(coord[2], coord[3]);
        WilsonLine w2 = ipglasma.GetWilsonLine(coord[4], coord[5]);

        B123 = ComputeBaryonOperator(w0, w1, w2);
        

        std::cout << "d = " << fixed << setprecision(6) << d << setprecision(6) << " " << real(N) << " " << imag(N) << std::endl;
    }
    */
    
    // Number of lattice points
    const int Nx = 512;
    const int Ny = 512;

    // Size of the lattice
    const double Lx = 50.67731; 
    const double Ly = 50.67731;
 
    const double dx = Lx / Nx;
    const double dy = Ly / Ny;
    
    /*
    cout << fixed << setprecision(10);
    for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
            double x = -Lx / 2.0 + (ix + 0.5) * dx;
            double y = -Ly / 2.0 + (iy + 0.5) * dy;
            double N_val = BaryonAmplitude(x, y, ipglasma);
            cout << N_val << "\n";
        }
    }
    */


    double numerator = 0.0;
    double denominator = 0.0;

    for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
            double x = -Lx / 2.0 + (ix + 0.5) * dx;
            double y = -Ly / 2.0 + (iy + 0.5) * dy;
            double N_val = BaryonAmplitude(x, y, ipglasma);

            if (std::isnan(N_val) || std::isinf(N_val)) continue;          

            double b2 = x * x + y * y;
            numerator   += 0.5 * b2 * N_val * dx * dy;
            denominator += N_val * dx * dy;
        }
    }

    
    cout << fixed << setprecision(10) << numerator << " " << denominator << endl;

    return 0;
}

