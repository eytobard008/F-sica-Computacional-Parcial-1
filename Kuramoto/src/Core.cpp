#include "Kuramoto.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <limits>

using namespace std;

/**
 * @brief Ejecuta un comando de sistema de forma segura
 * @param comando Comando a ejecutar
 * @return Resultado de la ejecución
 */
int ejecutarComando(const std::string& comando) {
    int resultado = system(comando.c_str());
    return resultado;
}

// ========== IMPLEMENTACIÓN KURAMOTO2 ==========

Kuramoto2::Kuramoto2(double dt, double frecuencia1, double frecuencia2) {
    dt_ = dt;
    frecuencia1_ = frecuencia1;
    frecuencia2_ = frecuencia2;
    theta1_ = 0;
    theta2_ = 0;
    K_ = 0;
}

void Kuramoto2::iniciar(double theta1, double theta2) {
    theta1_ = theta1;
    theta2_ = theta2;
}

void Kuramoto2::configurarAcoplamiento(double K) {
    K_ = K;
}

void Kuramoto2::pasoIntegracion() {
    Oscilador o1(frecuencia1_, theta1_);
    Oscilador o2(frecuencia2_, theta2_);
    
    double t1 = o1.getTheta();
    double t2 = o2.getTheta();
    
    double k1_1 = o1.getFrecuencia() + o1.calcularInfluencia(o2, K_, 2);
    double k1_2 = o2.getFrecuencia() + o2.calcularInfluencia(o1, K_, 2);
    
    Oscilador o1_temp(o1.getFrecuencia(), t1 + dt_ * k1_1 / 2.0);
    Oscilador o2_temp(o2.getFrecuencia(), t2 + dt_ * k1_2 / 2.0);
    
    double k2_1 = o1.getFrecuencia() + o1_temp.calcularInfluencia(o2_temp, K_, 2);
    double k2_2 = o2.getFrecuencia() + o2_temp.calcularInfluencia(o1_temp, K_, 2);
    
    o1_temp.setTheta(t1 + dt_ * k2_1 / 2.0);
    o2_temp.setTheta(t2 + dt_ * k2_2 / 2.0);
    
    double k3_1 = o1.getFrecuencia() + o1_temp.calcularInfluencia(o2_temp, K_, 2);
    double k3_2 = o2.getFrecuencia() + o2_temp.calcularInfluencia(o1_temp, K_, 2);
    
    o1_temp.setTheta(t1 + dt_ * k3_1);
    o2_temp.setTheta(t2 + dt_ * k3_2);
    
    double k4_1 = o1.getFrecuencia() + o1_temp.calcularInfluencia(o2_temp, K_, 2);
    double k4_2 = o2.getFrecuencia() + o2_temp.calcularInfluencia(o1_temp, K_, 2);
    
    theta1_ = t1 + dt_ * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1) / 6.0;
    theta2_ = t2 + dt_ * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2) / 6.0;
}

double Kuramoto2::calcularParametroOrden() const {
    complex<double> suma = polar(1.0, theta1_) + polar(1.0, theta2_);
    return abs(suma) / 2.0;
}

double Kuramoto2::getTheta1() const { return theta1_; }
double Kuramoto2::getTheta2() const { return theta2_; }

// ========== IMPLEMENTACIÓN KURAMOTO3 ==========

Kuramoto3::Kuramoto3(double dt, double frecuencia1, double frecuencia2, double frecuencia3) {
    dt_ = dt;
    frecuencia1_ = frecuencia1;
    frecuencia2_ = frecuencia2;
    frecuencia3_ = frecuencia3;
    theta1_ = 0;
    theta2_ = 0;
    theta3_ = 0;
    K_ = 0;
}

void Kuramoto3::iniciar(double theta1, double theta2, double theta3) {
    theta1_ = theta1;
    theta2_ = theta2;
    theta3_ = theta3;
}

void Kuramoto3::configurarAcoplamiento(double K) {
    K_ = K;
}

void Kuramoto3::pasoIntegracion() {
    Oscilador o1(frecuencia1_, theta1_);
    Oscilador o2(frecuencia2_, theta2_);
    Oscilador o3(frecuencia3_, theta3_);
    
    double t1 = o1.getTheta();
    double t2 = o2.getTheta();
    double t3 = o3.getTheta();
    
    double k1_1 = o1.getFrecuencia() + o1.calcularInfluencia(o2, K_, 3) + o1.calcularInfluencia(o3, K_, 3);
    double k1_2 = o2.getFrecuencia() + o2.calcularInfluencia(o1, K_, 3) + o2.calcularInfluencia(o3, K_, 3);
    double k1_3 = o3.getFrecuencia() + o3.calcularInfluencia(o1, K_, 3) + o3.calcularInfluencia(o2, K_, 3);
    
    Oscilador o1_temp(o1.getFrecuencia(), t1 + dt_ * k1_1 / 2.0);
    Oscilador o2_temp(o2.getFrecuencia(), t2 + dt_ * k1_2 / 2.0);
    Oscilador o3_temp(o3.getFrecuencia(), t3 + dt_ * k1_3 / 2.0);
    
    double k2_1 = o1.getFrecuencia() + o1_temp.calcularInfluencia(o2_temp, K_, 3) + o1_temp.calcularInfluencia(o3_temp, K_, 3);
    double k2_2 = o2.getFrecuencia() + o2_temp.calcularInfluencia(o1_temp, K_, 3) + o2_temp.calcularInfluencia(o3_temp, K_, 3);
    double k2_3 = o3.getFrecuencia() + o3_temp.calcularInfluencia(o1_temp, K_, 3) + o3_temp.calcularInfluencia(o2_temp, K_, 3);
    
    o1_temp.setTheta(t1 + dt_ * k2_1 / 2.0);
    o2_temp.setTheta(t2 + dt_ * k2_2 / 2.0);
    o3_temp.setTheta(t3 + dt_ * k2_3 / 2.0);
    
    double k3_1 = o1.getFrecuencia() + o1_temp.calcularInfluencia(o2_temp, K_, 3) + o1_temp.calcularInfluencia(o3_temp, K_, 3);
    double k3_2 = o2.getFrecuencia() + o2_temp.calcularInfluencia(o1_temp, K_, 3) + o2_temp.calcularInfluencia(o3_temp, K_, 3);
    double k3_3 = o3.getFrecuencia() + o3_temp.calcularInfluencia(o1_temp, K_, 3) + o3_temp.calcularInfluencia(o2_temp, K_, 3);
    
    o1_temp.setTheta(t1 + dt_ * k3_1);
    o2_temp.setTheta(t2 + dt_ * k3_2);
    o3_temp.setTheta(t3 + dt_ * k3_3);
    
    double k4_1 = o1.getFrecuencia() + o1_temp.calcularInfluencia(o2_temp, K_, 3) + o1_temp.calcularInfluencia(o3_temp, K_, 3);
    double k4_2 = o2.getFrecuencia() + o2_temp.calcularInfluencia(o1_temp, K_, 3) + o2_temp.calcularInfluencia(o3_temp, K_, 3);
    double k4_3 = o3.getFrecuencia() + o3_temp.calcularInfluencia(o1_temp, K_, 3) + o3_temp.calcularInfluencia(o2_temp, K_, 3);
    
    theta1_ = t1 + dt_ * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1) / 6.0;
    theta2_ = t2 + dt_ * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2) / 6.0;
    theta3_ = t3 + dt_ * (k1_3 + 2*k2_3 + 2*k3_3 + k4_3) / 6.0;
}

double Kuramoto3::calcularParametroOrden() const {
    complex<double> suma = polar(1.0, theta1_) + 
                         polar(1.0, theta2_) + 
                         polar(1.0, theta3_);
    return abs(suma) / 3.0;
}

double Kuramoto3::getTheta1() const { return theta1_; }
double Kuramoto3::getTheta2() const { return theta2_; }
double Kuramoto3::getTheta3() const { return theta3_; }

// ========== FUNCIONES BÁSICAS ==========

/**
 * @brief Formatea un número para su uso en nombres de archivo
 * @param valor Número a formatear
 * @return String formateado
 */
string formatearNumero(double valor) {
    ostringstream oss;
    oss << fixed << setprecision(3) << valor;
    string resultado = oss.str();
    size_t pos = resultado.find(".000");
    if (pos != string::npos && pos == resultado.length() - 4) {
        resultado = resultado.substr(0, resultado.length() - 4);
    }
    return resultado;
}

/**
 * @brief Calcula el umbral de acoplamiento para 2 osciladores
 * @param frecuencia1 Frecuencia del primer oscilador
 * @param frecuencia2 Frecuencia del segundo oscilador
 * @return Valor del umbral K
 */
double calcularUmbralK2osc(double frecuencia1, double frecuencia2) {
    return abs(frecuencia2 - frecuencia1);
}

/**
 * @brief Calcula el umbral de acoplamiento para 3 osciladores
 * @param frecuencia1 Frecuencia del primer oscilador
 * @param frecuencia2 Frecuencia del segundo oscilador
 * @param frecuencia3 Frecuencia del tercer oscilador
 * @return Valor del umbral K
 */
double calcularUmbralK3osc(double frecuencia1, double frecuencia2, double frecuencia3) {
    double max_frecuencia = max(frecuencia1, max(frecuencia2, frecuencia3));
    double min_frecuencia = min(frecuencia1, min(frecuencia2, frecuencia3));
    return max_frecuencia - min_frecuencia;
}

/**
 * @brief Crea la estructura de directorios necesaria para los resultados
 */
void crearDirectorios() {
    ejecutarComando("mkdir -p results/datos/2osc results/datos/3osc");
    ejecutarComando("mkdir -p results/visual/2osc results/visual/3osc");
    ejecutarComando("mkdir -p scripts/gnuplot scripts/python");
}
