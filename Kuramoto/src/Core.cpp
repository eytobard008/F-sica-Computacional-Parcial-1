#include "Kuramoto.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <limits>

using namespace std;

// ========== FUNCIÓN SAFE_SYSTEM ==========
int safe_system(const std::string& command) {
    int result = system(command.c_str());
    (void)result;
    return result;
}

// ========== CLASE KURAMOTO2 ==========

Kuramoto2::Kuramoto2(double dt, double omega1, double omega2) 
    : dt_(dt), omega1_(omega1), omega2_(omega2), theta1_(0), theta2_(0), K_(0) {}

void Kuramoto2::Inicie(double theta1, double theta2) {
    theta1_ = theta1;
    theta2_ = theta2;
}

void Kuramoto2::set_coupling(double K) {
    K_ = K;
}

void Kuramoto2::Paso() {
    // RK4 para 2 osciladores - ecuaciones explícitas
    double k1_1, k1_2, k2_1, k2_2, k3_1, k3_2, k4_1, k4_2;
    
    // k1
    k1_1 = omega1_ + (K_/2.0) * sin(theta2_ - theta1_);
    k1_2 = omega2_ + (K_/2.0) * sin(theta1_ - theta2_);
    
    // k2
    double theta1_temp = theta1_ + dt_ * k1_1 / 2.0;
    double theta2_temp = theta2_ + dt_ * k1_2 / 2.0;
    k2_1 = omega1_ + (K_/2.0) * sin(theta2_temp - theta1_temp);
    k2_2 = omega2_ + (K_/2.0) * sin(theta1_temp - theta2_temp);
    
    // k3
    theta1_temp = theta1_ + dt_ * k2_1 / 2.0;
    theta2_temp = theta2_ + dt_ * k2_2 / 2.0;
    k3_1 = omega1_ + (K_/2.0) * sin(theta2_temp - theta1_temp);
    k3_2 = omega2_ + (K_/2.0) * sin(theta1_temp - theta2_temp);
    
    // k4
    theta1_temp = theta1_ + dt_ * k3_1;
    theta2_temp = theta2_ + dt_ * k3_2;
    k4_1 = omega1_ + (K_/2.0) * sin(theta2_temp - theta1_temp);
    k4_2 = omega2_ + (K_/2.0) * sin(theta1_temp - theta2_temp);
    
    // Combinar
    theta1_ += dt_ * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1) / 6.0;
    theta2_ += dt_ * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2) / 6.0;
}

double Kuramoto2::calculate_order_parameter() const {
    complex<double> sum = polar(1.0, theta1_) + polar(1.0, theta2_);
    return abs(sum) / 2.0;
}

double Kuramoto2::GetTheta1() const { return theta1_; }
double Kuramoto2::GetTheta2() const { return theta2_; }

// ========== CLASE KURAMOTO3 ==========

Kuramoto3::Kuramoto3(double dt, double omega1, double omega2, double omega3) 
    : dt_(dt), omega1_(omega1), omega2_(omega2), omega3_(omega3), 
      theta1_(0), theta2_(0), theta3_(0), K_(0) {}

void Kuramoto3::Inicie(double theta1, double theta2, double theta3) {
    theta1_ = theta1;
    theta2_ = theta2;
    theta3_ = theta3;
}

void Kuramoto3::set_coupling(double K) {
    K_ = K;
}

void Kuramoto3::Paso() {
    // RK4 para 3 osciladores - ecuaciones explícitas
    double k1_1, k1_2, k1_3, k2_1, k2_2, k2_3, k3_1, k3_2, k3_3, k4_1, k4_2, k4_3;
    
    // k1
    k1_1 = omega1_ + (K_/3.0) * (sin(theta2_ - theta1_) + sin(theta3_ - theta1_));
    k1_2 = omega2_ + (K_/3.0) * (sin(theta1_ - theta2_) + sin(theta3_ - theta2_));
    k1_3 = omega3_ + (K_/3.0) * (sin(theta1_ - theta3_) + sin(theta2_ - theta3_));
    
    // k2
    double theta1_temp = theta1_ + dt_ * k1_1 / 2.0;
    double theta2_temp = theta2_ + dt_ * k1_2 / 2.0;
    double theta3_temp = theta3_ + dt_ * k1_3 / 2.0;
    
    k2_1 = omega1_ + (K_/3.0) * (sin(theta2_temp - theta1_temp) + sin(theta3_temp - theta1_temp));
    k2_2 = omega2_ + (K_/3.0) * (sin(theta1_temp - theta2_temp) + sin(theta3_temp - theta2_temp));
    k2_3 = omega3_ + (K_/3.0) * (sin(theta1_temp - theta3_temp) + sin(theta2_temp - theta3_temp));
    
    // k3
    theta1_temp = theta1_ + dt_ * k2_1 / 2.0;
    theta2_temp = theta2_ + dt_ * k2_2 / 2.0;
    theta3_temp = theta3_ + dt_ * k2_3 / 2.0;
    
    k3_1 = omega1_ + (K_/3.0) * (sin(theta2_temp - theta1_temp) + sin(theta3_temp - theta1_temp));
    k3_2 = omega2_ + (K_/3.0) * (sin(theta1_temp - theta2_temp) + sin(theta3_temp - theta2_temp));
    k3_3 = omega3_ + (K_/3.0) * (sin(theta1_temp - theta3_temp) + sin(theta2_temp - theta3_temp));
    
    // k4
    theta1_temp = theta1_ + dt_ * k3_1;
    theta2_temp = theta2_ + dt_ * k3_2;
    theta3_temp = theta3_ + dt_ * k3_3;
    
    k4_1 = omega1_ + (K_/3.0) * (sin(theta2_temp - theta1_temp) + sin(theta3_temp - theta1_temp));
    k4_2 = omega2_ + (K_/3.0) * (sin(theta1_temp - theta2_temp) + sin(theta3_temp - theta2_temp));
    k4_3 = omega3_ + (K_/3.0) * (sin(theta1_temp - theta3_temp) + sin(theta2_temp - theta3_temp));
    
    // Combinar
    theta1_ += dt_ * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1) / 6.0;
    theta2_ += dt_ * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2) / 6.0;
    theta3_ += dt_ * (k1_3 + 2*k2_3 + 2*k3_3 + k4_3) / 6.0;
}

double Kuramoto3::calculate_order_parameter() const {
    complex<double> sum = polar(1.0, theta1_) + polar(1.0, theta2_) + polar(1.0, theta3_);
    return abs(sum) / 3.0;
}

double Kuramoto3::GetTheta1() const { return theta1_; }
double Kuramoto3::GetTheta2() const { return theta2_; }
double Kuramoto3::GetTheta3() const { return theta3_; }

// ========== FUNCIONES BÁSICAS DE CÁLCULO ==========

string format_number(double value) {
    ostringstream oss;
    oss << fixed << setprecision(3) << value;
    string result = oss.str();
    size_t pos = result.find(".000");
    if (pos != string::npos && pos == result.length() - 4) {
        result = result.substr(0, result.length() - 4);
    }
    return result;
}

double calculate_critical_K_2osc(const vector<double>& omega) {
    return abs(omega[1] - omega[0]);
}

double calcular_K_umbral_3osc(const vector<double>& omega) {
    return *max_element(omega.begin(), omega.end()) - *min_element(omega.begin(), omega.end());
}

void ensure_directories() {
    safe_system("mkdir -p results/datos/2osc results/datos/3osc");
    safe_system("mkdir -p results/visual/2osc results/visual/3osc");
    safe_system("mkdir -p scripts/gnuplot scripts/python");
}