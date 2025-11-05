#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Kuramoto.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

int main() {
    cout << "=========================================\n";
    cout << "     SISTEMA DE KURAMOTO - ANALISIS COMPLETO\n";
    cout << "=========================================\n\n";
    
    SimulationConfig config;
    
    // ========== CONFIGURACIÓN DE PARÁMETROS GENERALES ==========
    
    // Parámetros generales
    config.dt = 0.01;
    config.steps = 5000;
    
    // 2 Osciladores
    config.omega_2osc = {1.0, 1.7};
    config.initial_phases_2osc = {0.0, M_PI/4};
    
    // 3 Osciladores  
    config.omega_3osc = {0.3, 1.0, 1.7};
    config.initial_phases_3osc = {0.0, M_PI/3, 2*M_PI/3};
    
    // ========== CONFIGURACIÓN DE PARÁMETROS K ==========
    
    // 2 OSCILADORES - K crítico y valores
    double K_c_2osc = calculate_critical_K_2osc(config.omega_2osc);
    config.K_values_2osc = {0.0, 0.93*K_c_2osc, K_c_2osc, 3.0*K_c_2osc};
    
    // 3 OSCILADORES - Parámetros configurables
    config.n_K_points_3osc = 2000;           // Número de puntos K
    config.K_max_3osc = 4.0;                // K máximo para barrido
    config.K_umbral_3osc = calcular_K_umbral_3osc(config.omega_3osc);  // K umbral calculado una vez
    
    cout << "PARAMETROS CONFIGURADOS:" << endl;
    cout << "  2 Osciladores - K umbral: " << K_c_2osc << endl;
    cout << "  3 Osciladores - K umbral: " << config.K_umbral_3osc << endl;
    cout << "  3 Osciladores - Puntos K: " << config.n_K_points_3osc << endl;
    cout << "  3 Osciladores - K maximo: " << config.K_max_3osc << endl;
    
    // ========== EJECUCIÓN PRINCIPAL ==========
    
    ensure_directories();
    
    cout << "\n1. ANALISIS 2 OSCILADORES..." << endl;
    analyze_2osc_system(config);
    
    cout << "\n2. BIFURCACION 3 OSCILADORES..." << endl;
    analyze_3osc_bifurcation(config);
    
    // Construir vector de K para 3 osciladores después del análisis de bifurcación
    config.K_values_3osc = {0.0, config.K_min_3osc, config.K_umbral_3osc, 3.0 * config.K_umbral_3osc};
    
    cout << "   Valores K 3osc configurados: ";
    for (double K : config.K_values_3osc) cout << K << " ";
    cout << endl;
    
    cout << "\n3. ANALISIS 3 OSCILADORES..." << endl;
    analyze_3osc_system(config);
    
    cout << "\n4. MAPAS DE SENSIBILIDAD..." << endl;
    generate_sensitivity_maps(config);
    
    cout << "\n5. MAPAS DE FRECUENCIAS..." << endl;
    generate_frequency_maps(config);
    
    cout << "\n6. ANIMACIONES 2OSC..." << endl;
    generate_2osc_animations(config);
    
    cout << "\n7. ANIMACIONES 3OSC..." << endl;
    generate_3osc_animations(config);
    
    cout << "\n8. PUNTOS FIJOS 3D..." << endl;
    generate_3d_fixed_points(config);
    
    cout << "\n9. GRAFICAS VARIAS..." << endl;
    plot_time_series(config);
    plot_phase_spaces(config);
    plot_order_parameters(config);
    plot_stability_diagrams(config);
    generate_3d_alpha_phase_space(config);
    
    cout << "\n10. BIFURCACION 2OSC..." << endl;
    analyze_2osc_bifurcation(config);
    
    cout << "\n11. GRAFICAS 3D INTERACTIVAS..." << endl;
    generate_3d_plots();
    
    cout << "\n=========================================\n";
    cout << "     ✅ ANALISIS COMPLETADO\n";
    cout << "=========================================\n";
    
    cout << "\nRESUMEN DE PARAMETROS CALCULADOS:" << endl;
    cout << "  2 Osciladores - K umbral: " << K_c_2osc << endl;
    cout << "  3 Osciladores - K umbral: " << config.K_umbral_3osc << endl;
    cout << "  3 Osciladores - K minimo: " << config.K_min_3osc << endl;
    
    return 0;
}