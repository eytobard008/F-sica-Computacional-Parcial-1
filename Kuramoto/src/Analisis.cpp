#include "Kuramoto.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

// ========== ANÁLISIS DE SISTEMAS ==========

/**
 * @brief Analiza el sistema de 2 osciladores para diferentes valores de K
 * @param config Configuración de la simulación
 */
void analizarSistema2osc(const ConfiguracionSimulacion& config) {
    cout << "ANALIZANDO SISTEMA DE 2 OSCILADORES" << endl;
    
    double K_umbral_2osc = calcularUmbralK2osc(config.frecuencia1_2osc, config.frecuencia2_2osc);
    cout << "  K umbral: " << K_umbral_2osc << endl;
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        cout << "  Simulando K = " << K << endl;
        
        Kuramoto2 modelo(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
        modelo.iniciar(config.theta1_2osc, config.theta2_2osc);
        modelo.configurarAcoplamiento(K);
        
        string K_str = formatearNumero(K);
        ofstream archivo_sim("results/datos/2osc/simulation_2osc_K" + K_str + ".dat");
        archivo_sim << "# time theta1 theta2 Delta cos1 cos2 sin1 sin2 r\n";
        
        for (int paso = 0; paso < config.pasos; paso++) {
            double tiempo = paso * config.dt;
            modelo.pasoIntegracion();
            
            if (paso % 2 == 0) {
                double theta1 = modelo.getTheta1();
                double theta2 = modelo.getTheta2();
                double Delta = theta2 - theta1;
                double cos1 = cos(theta1);
                double cos2 = cos(theta2);
                double sin1 = sin(theta1);
                double sin2 = sin(theta2);
                double r = modelo.calcularParametroOrden();
                
                archivo_sim << tiempo << " " << theta1 << " " << theta2 << " " << Delta << " "
                        << cos1 << " " << cos2 << " " << sin1 << " " << sin2 << " " << r << "\n";
            }
        }
        archivo_sim.close();
    }
}

/**
 * @brief Analiza el sistema de 3 osciladores para diferentes valores de K
 * @param config Configuración de la simulación
 */
void analizarSistema3osc(const ConfiguracionSimulacion& config) {
    cout << "ANALIZANDO SISTEMA DE 3 OSCILADORES" << endl;
    
    double valores_K[4] = {config.K1_3osc, config.K2_3osc, config.K3_3osc, config.K4_3osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        cout << "  Simulando K = " << K << endl;
        
        Kuramoto3 modelo(config.dt, config.frecuencia1_3osc, config.frecuencia2_3osc, config.frecuencia3_3osc);
        modelo.iniciar(config.theta1_3osc, config.theta2_3osc, config.theta3_3osc);
        modelo.configurarAcoplamiento(K);
        
        string K_str = formatearNumero(K);
        ofstream archivo_sim("results/datos/3osc/simulation_3osc_K" + K_str + ".dat");
        archivo_sim << "# time theta1 theta2 theta3 cos1 cos2 cos3 sin1 sin2 sin3 r\n";
        
        for (int paso = 0; paso < config.pasos; paso++) {
            double tiempo = paso * config.dt;
            modelo.pasoIntegracion();
            
            if (paso % 2 == 0) {
                double theta1 = modelo.getTheta1();
                double theta2 = modelo.getTheta2();
                double theta3 = modelo.getTheta3();
                double cos1 = cos(theta1);
                double cos2 = cos(theta2);
                double cos3 = cos(theta3);
                double sin1 = sin(theta1);
                double sin2 = sin(theta2);
                double sin3 = sin(theta3);
                double r = modelo.calcularParametroOrden();
                
                archivo_sim << tiempo << " " << theta1 << " " << theta2 << " " << theta3 << " "
                        << cos1 << " " << cos2 << " " << cos3 << " "
                        << sin1 << " " << sin2 << " " << sin3 << " " << r << "\n";
            }
        }
        archivo_sim.close();
    }
}

// ========== BIFURCACIONES 2 OSCILADORES ==========

/**
 * @brief Analiza la bifurcación en sistema de 2 osciladores usando diferencias de fase
 * @param config Configuración de la simulación
 */
void analizarBifurcacion2oscFase(const ConfiguracionSimulacion& config) {
    cout << "BIFURCACION 2 OSCILADORES - DIFERENCIAS DE FASES" << endl;
    
    double K_umbral_2osc = calcularUmbralK2osc(config.frecuencia1_2osc, config.frecuencia2_2osc);
    
    ofstream archivo_datos("results/datos/2osc/bifurcation_2osc_phase.dat");
    archivo_datos << "# K phase_diff\n";
    
    int n_puntos = 200;
    for (int i = 0; i <= n_puntos; ++i) {
        double K = 3.0 * K_umbral_2osc * i / n_puntos;
        
        Kuramoto2 modelo(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
        modelo.iniciar(config.theta1_2osc, config.theta2_2osc);
        modelo.configurarAcoplamiento(K);
        
        for (int paso = 0; paso < config.pasos; paso++) {
            modelo.pasoIntegracion();
        }
        
        double diff_fase = abs(modelo.getTheta2() - modelo.getTheta1());
        archivo_datos << K << " " << diff_fase << "\n";
    }
    archivo_datos.close();
    
    ofstream script("scripts/gnuplot/plot_bifurcation_2osc_phase.gnu");
    script << "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/2osc/bifurcation_2osc_phase.png'\n";
    script << "set xlabel 'Acoplamiento K'\n";
    script << "set ylabel '|θ₂ - θ₁| final (rad)'\n";
    script << "set title 'Bifurcacion 2 Osciladores - Diferencias de Fase'\n";
    script << "set grid\n";
    script << "plot 'results/datos/2osc/bifurcation_2osc_phase.dat' u 1:2 w l lw 2 lc 'red' notitle\n";
    script.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_bifurcation_2osc_phase.gnu");
}

/**
 * @brief Analiza la bifurcación en sistema de 2 osciladores usando derivadas de fase
 * @param config Configuración de la simulación
 */
void analizarBifurcacion2oscDerivada(const ConfiguracionSimulacion& config) {
    cout << "BIFURCACION 2 OSCILADORES - DERIVADAS DE FASES" << endl;
    
    double K_umbral_2osc = calcularUmbralK2osc(config.frecuencia1_2osc, config.frecuencia2_2osc);
    
    ofstream archivo_datos("results/datos/2osc/bifurcation_2osc_deriv.dat");
    archivo_datos << "# K deriv_sum\n";
    
    int n_puntos = 200;
    for (int i = 0; i <= n_puntos; ++i) {
        double K = 3.0 * K_umbral_2osc * i / n_puntos;
        
        Kuramoto2 modelo(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
        modelo.iniciar(config.theta1_2osc, config.theta2_2osc);
        modelo.configurarAcoplamiento(K);
        
        for (int paso = 0; paso < config.pasos; paso++) {
            modelo.pasoIntegracion();
        }
        
        double theta1 = modelo.getTheta1();
        double theta2 = modelo.getTheta2();
        double Delta = theta2 - theta1;
        
        double dDelta_dt = (config.frecuencia2_2osc - config.frecuencia1_2osc) - K * sin(Delta);
        double suma_deriv = abs(dDelta_dt);
        
        archivo_datos << K << " " << suma_deriv << "\n";
    }
    archivo_datos.close();
    
    ofstream script("scripts/gnuplot/plot_bifurcation_2osc_deriv.gnu");
    script << "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/2osc/bifurcation_2osc_deriv.png'\n";
    script << "set xlabel 'Acoplamiento K'\n";
    script << "set ylabel '|dΔ/dt|'\n";
    script << "set title 'Bifurcacion 2 Osciladores - Derivadas de Fase'\n";
    script << "set grid\n";
    script << "plot 'results/datos/2osc/bifurcation_2osc_deriv.dat' u 1:2 w l lw 2 lc 'blue' notitle\n";
    script.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_bifurcation_2osc_deriv.gnu");
}

// ========== BIFURCACIONES 3 OSCILADORES ==========

/**
 * @brief Analiza la bifurcación en sistema de 3 osciladores usando diferencias de fase
 * @param config Configuración de la simulación
 */
void analizarBifurcacion3oscFase(const ConfiguracionSimulacion& config) {
    cout << "BIFURCACION 3 OSCILADORES - DIFERENCIAS DE FASES" << endl;
    
    ofstream archivo_datos("results/datos/3osc/bifurcation_3osc_phase.dat");
    archivo_datos << "# K phase_diffs_sum\n";
    
    int n_puntos = 200;
    for (int i = 0; i <= n_puntos; ++i) {
        double K = config.K_max_3osc * i / n_puntos;
        
        Kuramoto3 modelo(config.dt, config.frecuencia1_3osc, config.frecuencia2_3osc, config.frecuencia3_3osc);
        modelo.iniciar(config.theta1_3osc, config.theta2_3osc, config.theta3_3osc);
        modelo.configurarAcoplamiento(K);
        
        for (int paso = 0; paso < config.pasos; paso++) {
            modelo.pasoIntegracion();
        }
        
        double theta1 = modelo.getTheta1();
        double theta2 = modelo.getTheta2();
        double theta3 = modelo.getTheta3();
        
        double delta21 = theta2 - theta1;
        double delta31 = theta3 - theta1;
        double delta32 = theta3 - theta2;
        
        double suma_diferencias_fase = abs(delta21) + abs(delta31) + abs(delta32);
        archivo_datos << K << " " << suma_diferencias_fase << "\n";
    }
    archivo_datos.close();
    
    ofstream script("scripts/gnuplot/plot_bifurcation_3osc_phase.gnu");
    script << "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/3osc/bifurcation_3osc_phase.png'\n";
    script << "set xlabel 'Acoplamiento K'\n";
    script << "set ylabel 'Σ|Δθ| final'\n";
    script << "set title 'Bifurcacion 3 Osciladores - Diferencias de Fase'\n";
    script << "set grid\n";
    script << "plot 'results/datos/3osc/bifurcation_3osc_phase.dat' u 1:2 w l lw 2 lc 'green' notitle\n";
    script.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_bifurcation_3osc_phase.gnu");
}

/**
 * @brief Analiza la bifurcación en sistema de 3 osciladores usando derivadas de fase
 * @param config Configuración de la simulación
 */
void analizarBifurcacion3oscDerivada(const ConfiguracionSimulacion& config) {
    cout << "BIFURCACION 3 OSCILADORES - DERIVADAS DE FASES" << endl;
    
    ofstream archivo_datos("results/datos/3osc/bifurcation_3osc_deriv.dat");
    archivo_datos << "# K deriv_sum\n";
    
    double min_suma = numeric_limits<double>::max();
    double K_min = 0.0;
    
    int n_puntos = config.puntos_K_3osc;
    for (int i = 0; i <= n_puntos; ++i) {
        double K = config.K_max_3osc * i / n_puntos;
        
        Kuramoto3 modelo(config.dt, config.frecuencia1_3osc, config.frecuencia2_3osc, config.frecuencia3_3osc);
        modelo.iniciar(config.theta1_3osc, config.theta2_3osc, config.theta3_3osc);
        modelo.configurarAcoplamiento(K);
        
        for (int paso = 0; paso < config.pasos; paso++) {
            modelo.pasoIntegracion();
        }
        
        double theta1 = modelo.getTheta1();
        double theta2 = modelo.getTheta2();
        double theta3 = modelo.getTheta3();
        
        double Delta21 = theta2 - theta1;
        double Delta31 = theta3 - theta1;
        double Delta32 = theta3 - theta2;
        
        double d21_dt = (config.frecuencia2_3osc - config.frecuencia1_3osc) - (K/3.0) * (sin(Delta21) + sin(Delta31 - Delta21));
        double d31_dt = (config.frecuencia3_3osc - config.frecuencia1_3osc) - (K/3.0) * (sin(Delta31) + sin(Delta21 - Delta31));
        double d32_dt = (config.frecuencia3_3osc - config.frecuencia2_3osc) - (K/3.0) * (sin(Delta32) + sin(Delta31 - Delta32));
        
        double suma_deriv = abs(d21_dt) + abs(d31_dt) + abs(d32_dt);
        
        archivo_datos << K << " " << suma_deriv << "\n";
        
        if (suma_deriv < min_suma) {
            min_suma = suma_deriv;
            K_min = K;
        }
    }
    archivo_datos.close();
    
    const_cast<ConfiguracionSimulacion&>(config).K_min_3osc = K_min;
    
    cout << "  K_umbral teorico: " << config.K_umbral_3osc << endl;
    cout << "  K_minimo encontrado: " << K_min << endl;
    
    ofstream script("scripts/gnuplot/plot_bifurcation_3osc_deriv.gnu");
    script << "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/3osc/bifurcation_3osc_deriv.png'\n";
    script << "set xlabel 'Acoplamiento K'\n";
    script << "set ylabel 'Σ|dΔ/dt|'\n";
    script << "set title 'Bifurcacion 3 Osciladores - Derivadas de Fase'\n";
    script << "set grid\n";
    script << "plot 'results/datos/3osc/bifurcation_3osc_deriv.dat' u 1:2 w l lw 2 lc 'blue' notitle\n";
    script.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_bifurcation_3osc_deriv.gnu");
}

// ========== MAPAS DE SENSIBILIDAD ==========

/**
 * @brief Genera mapas de sensibilidad basados en diferencias de fase
 * @param config Configuración de la simulación
 */
void generarMapasSensibilidadFase(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO MAPAS DE SENSIBILIDAD - FASES" << endl;
    
    const int tamano_grid = 50;
    const double epsilon = 0.01;
    
    double minimo_global = numeric_limits<double>::max();
    double maximo_global = numeric_limits<double>::min();
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        
        string K_str = formatearNumero(K);
        ofstream archivo_sens("results/datos/2osc/sensitivity_phase_K" + K_str + ".dat");
        archivo_sens << "# theta1_0 theta2_0 sensitivity\n";
        
        for (int i = 0; i < tamano_grid; ++i) {
            double theta1_0 = 2 * M_PI * i / tamano_grid;
            
            for (int j = 0; j < tamano_grid; ++j) {
                double theta2_0 = 2 * M_PI * j / tamano_grid;
                
                Kuramoto2 modelo_centro(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo_centro.iniciar(theta1_0, theta2_0);
                modelo_centro.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo_centro.pasoIntegracion();
                double theta1_centro = modelo_centro.getTheta1();
                double theta2_centro = modelo_centro.getTheta2();
                
                Kuramoto2 modelo_dtheta1(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo_dtheta1.iniciar(theta1_0 + epsilon, theta2_0);
                modelo_dtheta1.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo_dtheta1.pasoIntegracion();
                double theta1_d1 = modelo_dtheta1.getTheta1();
                double theta2_d1 = modelo_dtheta1.getTheta2();
                
                Kuramoto2 modelo_dtheta2(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo_dtheta2.iniciar(theta1_0, theta2_0 + epsilon);
                modelo_dtheta2.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo_dtheta2.pasoIntegracion();
                double theta1_d2 = modelo_dtheta2.getTheta1();
                double theta2_d2 = modelo_dtheta2.getTheta2();
                
                double sensibilidad_fase = (abs(theta1_d1 - theta1_centro) + abs(theta2_d1 - theta2_centro) +
                                   abs(theta1_d2 - theta1_centro) + abs(theta2_d2 - theta2_centro)) / (2.0 * epsilon);
                
                archivo_sens << theta1_0 << " " << theta2_0 << " " << sensibilidad_fase << "\n";
                
                minimo_global = min(minimo_global, sensibilidad_fase);
                maximo_global = max(maximo_global, sensibilidad_fase);
            }
            archivo_sens << "\n";
        }
        archivo_sens.close();
    }
}

/**
 * @brief Genera mapas de sensibilidad basados en amplitudes
 * @param config Configuración de la simulación
 */
void generarMapasSensibilidadAmplitud(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO MAPAS DE SENSIBILIDAD - AMPLITUDES" << endl;
    
    const int tamano_grid = 50;
    const double epsilon = 0.01;
    
    double minimo_global = numeric_limits<double>::max();
    double maximo_global = numeric_limits<double>::min();
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        
        string K_str = formatearNumero(K);
        ofstream archivo_sens("results/datos/2osc/sensitivity_amplitude_K" + K_str + ".dat");
        archivo_sens << "# theta1_0 theta2_0 sensitivity\n";
        
        for (int i = 0; i < tamano_grid; ++i) {
            double theta1_0 = 2 * M_PI * i / tamano_grid;
            
            for (int j = 0; j < tamano_grid; ++j) {
                double theta2_0 = 2 * M_PI * j / tamano_grid;
                
                Kuramoto2 modelo_centro(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo_centro.iniciar(theta1_0, theta2_0);
                modelo_centro.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo_centro.pasoIntegracion();
                double theta1_centro = modelo_centro.getTheta1();
                double theta2_centro = modelo_centro.getTheta2();
                
                Kuramoto2 modelo_dtheta1(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo_dtheta1.iniciar(theta1_0 + epsilon, theta2_0);
                modelo_dtheta1.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo_dtheta1.pasoIntegracion();
                double theta1_d1 = modelo_dtheta1.getTheta1();
                double theta2_d1 = modelo_dtheta1.getTheta2();
                
                Kuramoto2 modelo_dtheta2(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo_dtheta2.iniciar(theta1_0, theta2_0 + epsilon);
                modelo_dtheta2.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo_dtheta2.pasoIntegracion();
                double theta1_d2 = modelo_dtheta2.getTheta1();
                double theta2_d2 = modelo_dtheta2.getTheta2();
                
                double sensibilidad_amplitud = (abs(cos(theta1_d1) - cos(theta1_centro)) + abs(cos(theta2_d1) - cos(theta2_centro)) +
                                      abs(cos(theta1_d2) - cos(theta1_centro)) + abs(cos(theta2_d2) - cos(theta2_centro))) / (2.0 * epsilon);
                
                archivo_sens << theta1_0 << " " << theta2_0 << " " << sensibilidad_amplitud << "\n";
                
                minimo_global = min(minimo_global, sensibilidad_amplitud);
                maximo_global = max(maximo_global, sensibilidad_amplitud);
            }
            archivo_sens << "\n";
        }
        archivo_sens.close();
    }
}

// ========== MAPAS DE FRECUENCIA ==========

/**
 * @brief Genera mapas de frecuencia basados en diferencias de fase
 * @param config Configuración de la simulación
 */
void generarMapasFrecuenciaFase(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO MAPAS DE FRECUENCIA - FASES" << endl;
    
    const int tamano_grid = 60;
    
    double minimo_global = numeric_limits<double>::max();
    double maximo_global = numeric_limits<double>::min();
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        
        string K_str = formatearNumero(K);
        ofstream archivo_frec("results/datos/2osc/frequency_phase_K" + K_str + ".dat");
        archivo_frec << "# theta1_0 theta2_0 freq_diff\n";
        
        for (int i = 0; i < tamano_grid; ++i) {
            double theta1_0 = 2 * M_PI * i / tamano_grid;
            
            for (int j = 0; j < tamano_grid; ++j) {
                double theta2_0 = 2 * M_PI * j / tamano_grid;
                
                Kuramoto2 modelo(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo.iniciar(theta1_0, theta2_0);
                modelo.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo.pasoIntegracion();
                
                double theta1 = modelo.getTheta1();
                double theta2 = modelo.getTheta2();
                
                double freq1 = config.frecuencia1_2osc + (K/2.0) * sin(theta2 - theta1);
                double freq2 = config.frecuencia2_2osc + (K/2.0) * sin(theta1 - theta2);
                double diff_frec = abs(freq1 - freq2);
                
                archivo_frec << theta1_0 << " " << theta2_0 << " " << diff_frec << "\n";
                
                minimo_global = min(minimo_global, diff_frec);
                maximo_global = max(maximo_global, diff_frec);
            }
            archivo_frec << "\n";
        }
        archivo_frec.close();
    }
}

/**
 * @brief Genera mapas de frecuencia basados en velocidades
 * @param config Configuración de la simulación
 */
void generarMapasFrecuenciaVelocidad(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO MAPAS DE FRECUENCIA - VELOCIDADES" << endl;
    
    const int tamano_grid = 60;
    
    double minimo_global = numeric_limits<double>::max();
    double maximo_global = numeric_limits<double>::min();
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        
        string K_str = formatearNumero(K);
        ofstream archivo_frec("results/datos/2osc/frequency_velocity_K" + K_str + ".dat");
        archivo_frec << "# theta1_0 theta2_0 velocity_diff\n";
        
        for (int i = 0; i < tamano_grid; ++i) {
            double theta1_0 = 2 * M_PI * i / tamano_grid;
            
            for (int j = 0; j < tamano_grid; ++j) {
                double theta2_0 = 2 * M_PI * j / tamano_grid;
                
                Kuramoto2 modelo(config.dt, config.frecuencia1_2osc, config.frecuencia2_2osc);
                modelo.iniciar(theta1_0, theta2_0);
                modelo.configurarAcoplamiento(K);
                for (int paso = 0; paso < 2000; paso++) modelo.pasoIntegracion();
                
                double theta1 = modelo.getTheta1();
                double theta2 = modelo.getTheta2();
                
                double omega1_inst = config.frecuencia1_2osc + (K/2.0) * sin(theta2 - theta1);
                double omega2_inst = config.frecuencia2_2osc + (K/2.0) * sin(theta1 - theta2);
                double diff_velocidad = abs(omega1_inst - omega2_inst);
                
                archivo_frec << theta1_0 << " " << theta2_0 << " " << diff_velocidad << "\n";
                
                minimo_global = min(minimo_global, diff_velocidad);
                maximo_global = max(maximo_global, diff_velocidad);
            }
            archivo_frec << "\n";
        }
        archivo_frec.close();
    }
}

// ========== PUNTOS FIJOS 3D ==========

/**
 * @brief Genera análisis de puntos fijos en 3D para sistema de 3 osciladores
 * @param config Configuración de la simulación
 */
void generarPuntosFijos3D(const ConfiguracionSimulacion& config) {
    cout << "ANALIZANDO PUNTOS FIJOS 3D" << endl;
    
    const int tamano_grid = 50;
    
    double valores_K[4] = {config.K1_3osc, config.K2_3osc, config.K3_3osc, config.K4_3osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        cout << "  K = " << K << endl;
        
        string K_str = formatearNumero(K);
        ofstream archivo_fijos("results/datos/3osc/fixed_points_3d_K" + K_str + ".dat");
        archivo_fijos << "# Delta21 Delta31 derivative_magnitude\n";
        
        for (int i = 0; i < tamano_grid; ++i) {
            double Delta21 = -M_PI + 2*M_PI * i / tamano_grid;
            
            for (int j = 0; j < tamano_grid; ++j) {
                double Delta31 = -M_PI + 2*M_PI * j / tamano_grid;
                
                double d21_dt = - K * (sin(Delta21) + sin(Delta31 - Delta21));
                double d31_dt = - K * (sin(Delta31) + sin(Delta21 - Delta31));
                
                double magnitud_derivada = sqrt(d21_dt*d21_dt + d31_dt*d31_dt);
                
                archivo_fijos << Delta21 << " " << Delta31 << " " << magnitud_derivada << "\n";
            }
            archivo_fijos << "\n";
        }
        archivo_fijos.close();
    }
}
