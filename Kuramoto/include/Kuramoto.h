#ifndef KURAMOTO_H
#define KURAMOTO_H

#include <string>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @file Kuramoto.h
 * @brief Sistema de osciladores acoplados modelo Kuramoto
 * @author Esteban Tobar
 * @date 29 de octubre de 2025
 */

// ========== CLASE OSCILADOR ==========

/**
 * @brief Representa un oscilador individual en el sistema Kuramoto
 * 
 * Cada oscilador tiene una frecuencia natural y una fase que evoluciona en el tiempo
 * según las ecuaciones del modelo Kuramoto con acoplamiento sinusoidal.
 */
class Oscilador {
private:
    double frecuencia_;
    double theta_;

public:
    /**
     * @brief Constructor del oscilador
     * @param frecuencia Frecuencia natural del oscilador
     * @param theta0 Fase inicial del oscilador
     */
    Oscilador(double frecuencia, double theta0) {
        frecuencia_ = frecuencia;
        theta_ = theta0;
    }
    
    double getTheta() const { return theta_; }
    double getFrecuencia() const { return frecuencia_; }
    void setTheta(double theta) { theta_ = theta; }
    
    /**
     * @brief Calcula la influencia de otro oscilador sobre este
     * @param otro Oscilador que ejerce la influencia
     * @param K Constante de acoplamiento
     * @param n Número total de osciladores en el sistema
     * @return Valor de la influencia calculada
     */
    double calcularInfluencia(const Oscilador& otro, double K, int n) const {
        return K * sin(otro.theta_ - theta_) / n;
    }
};

// ========== CLASE SISTEMA 2 OSCILADORES ==========

/**
 * @brief Implementa el sistema de Kuramoto para 2 osciladores
 * 
 * Resuelve las ecuaciones diferenciales usando el método Runge-Kutta de 4to orden
 * para un sistema de dos osciladores acoplados.
 */
class Kuramoto2 {
private:
    double dt_;
    double frecuencia1_, frecuencia2_;
    double theta1_, theta2_;
    double K_;

public:
    Kuramoto2(double dt, double frecuencia1, double frecuencia2);
    void iniciar(double theta1, double theta2);
    void configurarAcoplamiento(double K);
    void pasoIntegracion();
    double calcularParametroOrden() const;
    double getTheta1() const;
    double getTheta2() const;
};

// ========== CLASE SISTEMA 3 OSCILADORES ==========

/**
 * @brief Implementa el sistema de Kuramoto para 3 osciladores
 * 
 * Extensión del modelo a tres osciladores con acoplamiento global.
 */
class Kuramoto3 {
private:
    double dt_;
    double frecuencia1_, frecuencia2_, frecuencia3_;
    double theta1_, theta2_, theta3_;
    double K_;

public:
    Kuramoto3(double dt, double frecuencia1, double frecuencia2, double frecuencia3);
    void iniciar(double theta1, double theta2, double theta3);
    void configurarAcoplamiento(double K);
    void pasoIntegracion();
    double calcularParametroOrden() const;
    double getTheta1() const;
    double getTheta2() const;
    double getTheta3() const;
};

// ========== ESTRUCTURA DE CONFIGURACIÓN ==========

/**
 * @brief Almacena todos los parámetros de configuración para las simulaciones
 * 
 * Contiene parámetros para sistemas de 2 y 3 osciladores, así como configuraciones
 * para análisis de bifurcación y mapas de sensibilidad.
 */
struct ConfiguracionSimulacion {
    double dt;
    int pasos;
    
    // Parámetros para 2 osciladores
    double frecuencia1_2osc, frecuencia2_2osc;
    double theta1_2osc, theta2_2osc;
    double K1_2osc, K2_2osc, K3_2osc, K4_2osc;
    
    // Parámetros para 3 osciladores  
    double frecuencia1_3osc, frecuencia2_3osc, frecuencia3_3osc;
    double theta1_3osc, theta2_3osc, theta3_3osc;
    double K1_3osc, K2_3osc, K3_3osc, K4_3osc;
    
    // Parámetros para análisis de bifurcación
    int puntos_K_3osc;
    double K_max_3osc;
    double K_min_3osc;
    double K_umbral_2osc;
    double K_umbral_3osc;
};

// ========== DECLARACIONES DE FUNCIONES ==========

// Core.cpp
int ejecutarComando(const std::string& comando);
double calcularUmbralK2osc(double frecuencia1, double frecuencia2);
double calcularUmbralK3osc(double frecuencia1, double frecuencia2, double frecuencia3);
void crearDirectorios();
std::string formatearNumero(double valor);

// Analisis.cpp  
void analizarSistema2osc(const ConfiguracionSimulacion& config);
void analizarSistema3osc(const ConfiguracionSimulacion& config);

// Análisis de bifurcaciones para 2 osciladores
void analizarBifurcacion2oscFase(const ConfiguracionSimulacion& config);
void analizarBifurcacion2oscDerivada(const ConfiguracionSimulacion& config);

// Análisis de bifurcaciones para 3 osciladores  
void analizarBifurcacion3oscFase(const ConfiguracionSimulacion& config);
void analizarBifurcacion3oscDerivada(const ConfiguracionSimulacion& config);

// Generación de mapas de sensibilidad
void generarMapasSensibilidadFase(const ConfiguracionSimulacion& config);
void generarMapasSensibilidadAmplitud(const ConfiguracionSimulacion& config);
void generarMapasFrecuenciaFase(const ConfiguracionSimulacion& config);
void generarMapasFrecuenciaVelocidad(const ConfiguracionSimulacion& config);

// Análisis de puntos fijos
void generarPuntosFijos3D(const ConfiguracionSimulacion& config);

// Visualizacion.cpp
void generarAnimaciones2osc(const ConfiguracionSimulacion& config);
void generarAnimaciones3osc(const ConfiguracionSimulacion& config);
void graficarSeriesTemporales(const ConfiguracionSimulacion& config);
void graficarFasesTiempo(const ConfiguracionSimulacion& config);
void graficarDiagramasEstabilidad(const ConfiguracionSimulacion& config);
void graficarEspaciosFase(const ConfiguracionSimulacion& config);
void graficarParametrosOrden(const ConfiguracionSimulacion& config);
void generarGraficas3D();
void generarEspacioFase3DAlpha(const ConfiguracionSimulacion& config);

#endif
