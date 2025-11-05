#ifndef KURAMOTO_H
#define KURAMOTO_H

#include <vector>
#include <string>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ========== CLASE KURAMOTO2 ==========

class Kuramoto2 {
private:
    double dt_;
    double omega1_, omega2_;
    double theta1_, theta2_;
    double K_;

public:
    Kuramoto2(double dt, double omega1, double omega2);
    void Inicie(double theta1, double theta2);
    void set_coupling(double K);
    void Paso();
    double calculate_order_parameter() const;
    double GetTheta1() const;
    double GetTheta2() const;
};

// ========== CLASE KURAMOTO3 ==========

class Kuramoto3 {
private:
    double dt_;
    double omega1_, omega2_, omega3_;
    double theta1_, theta2_, theta3_;
    double K_;

public:
    Kuramoto3(double dt, double omega1, double omega2, double omega3);
    void Inicie(double theta1, double theta2, double theta3);
    void set_coupling(double K);
    void Paso();
    double calculate_order_parameter() const;
    double GetTheta1() const;
    double GetTheta2() const;
    double GetTheta3() const;
};

// Estructura de configuración
struct SimulationConfig {
    // Parámetros generales
    double dt;
    int steps;
    
    // 2 Osciladores
    std::vector<double> omega_2osc;
    std::vector<double> initial_phases_2osc;
    std::vector<double> K_values_2osc;
    
    // 3 Osciladores  
    std::vector<double> omega_3osc;
    std::vector<double> initial_phases_3osc;
    std::vector<double> K_values_3osc;
    
    // Parámetros configurables para bifurcación 3OSC
    int n_K_points_3osc;
    double K_max_3osc;
    double K_min_3osc;
    double K_umbral_3osc;
};

// ========== DECLARACIONES DE FUNCIONES ==========

// Funciones en Core.cpp
int safe_system(const std::string& command);
double calculate_critical_K_2osc(const std::vector<double>& omega);
double calcular_K_umbral_3osc(const std::vector<double>& omega);
void ensure_directories();
std::string format_number(double value);

// Funciones en Analisis.cpp  
void analyze_2osc_system(const SimulationConfig& config);
void analyze_3osc_system(const SimulationConfig& config);
void analyze_3osc_bifurcation(const SimulationConfig& config);
void generate_sensitivity_maps(const SimulationConfig& config);
void generate_frequency_maps(const SimulationConfig& config);
void analyze_2osc_bifurcation(const SimulationConfig& config);
void generate_3d_fixed_points(const SimulationConfig& config);

// Funciones en Visualizacion.cpp
void generate_2osc_animations(const SimulationConfig& config);
void generate_3osc_animations(const SimulationConfig& config);
void plot_time_series(const SimulationConfig& config);
void plot_stability_diagrams(const SimulationConfig& config);
void plot_phase_spaces(const SimulationConfig& config);
void plot_order_parameters(const SimulationConfig& config);
void generate_3d_plots();
void generate_3d_alpha_phase_space(const SimulationConfig& config);

#endif