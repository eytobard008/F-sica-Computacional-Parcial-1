#include "Kuramoto.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

// ========== ANÁLISIS DE 2 OSCILADORES ==========

void analyze_2osc_system(const SimulationConfig& config) {
    cout << "ANALIZANDO SISTEMA DE 2 OSCILADORES" << endl;
    cout << "  Pasos: " << config.steps << " | dt: " << config.dt << " | Tiempo total: " << config.steps*config.dt << "s" << endl;
    
    double K_c = calculate_critical_K_2osc(config.omega_2osc);
    cout << "  K umbral: " << K_c << endl;
    cout << "  Valores de K: ";
    for (double K : config.K_values_2osc) cout << K << " ";
    cout << endl;
    
    for (double K : config.K_values_2osc) {
        cout << "  Simulando K = " << K;
        if (abs(K - K_c) < 0.01) cout << " [UMBRAL]";
        cout << endl;
        
        Kuramoto2 model(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
        model.Inicie(config.initial_phases_2osc[0], config.initial_phases_2osc[1]);
        model.set_coupling(K);
        
        string K_str = format_number(K);
        
        // Abrir archivos para este K
        ofstream time_file("results/datos/2osc/time_series_K" + K_str + ".dat");
        ofstream circle_file("results/datos/2osc/circle_coords_K" + K_str + ".dat");
        ofstream phase_file("results/datos/2osc/phase_space_K" + K_str + ".dat");
        ofstream order_file("results/datos/2osc/order_parameter_K" + K_str + ".dat");
        
        time_file << "time theta1 theta2 phase_diff\n";
        circle_file << "time cos1 cos2 sin1 sin2 r\n";
        phase_file << "theta1 theta2\n";
        order_file << "time r\n";
        
        // Simular y escribir datos
        for (int step = 0; step < config.steps; step++) {
            double t = step * config.dt;
            model.Paso();
            
            // Series temporales cada 5 pasos
            if (step % 5 == 0) {
                double phase_diff = model.GetTheta2() - model.GetTheta1();
                time_file << t << " " << model.GetTheta1() << " " << model.GetTheta2() << " " << phase_diff << "\n";
                phase_file << model.GetTheta1() << " " << model.GetTheta2() << "\n";
            }
            
            // Coordenadas para animaciones cada 2 pasos
            if (step % 2 == 0) {
                double cos1 = cos(model.GetTheta1());
                double cos2 = cos(model.GetTheta2());
                double sin1 = sin(model.GetTheta1());
                double sin2 = sin(model.GetTheta2());
                double r = model.calculate_order_parameter();
                circle_file << t << " " << cos1 << " " << cos2 << " " << sin1 << " " << sin2 << " " << r << "\n";
            }
            
            // Parámetro de orden cada 5 pasos
            if (step % 5 == 0) {
                double r = model.calculate_order_parameter();
                order_file << t << " " << r << "\n";
            }
        }
        
        time_file.close();
        circle_file.close();
        phase_file.close();
        order_file.close();
    }
}

// ========== ANÁLISIS DE 3 OSCILADORES ==========

void analyze_3osc_system(const SimulationConfig& config) {
    cout << "ANALIZANDO SISTEMA DE 3 OSCILADORES" << endl;
    
    cout << "  Valores de K: ";
    for (double K : config.K_values_3osc) cout << K << " ";
    cout << endl;
    
    for (double K : config.K_values_3osc) {
        cout << "  Simulando K = " << K;
        if (abs(K - config.K_umbral_3osc) < 0.01) cout << " [UMBRAL]";
        cout << endl;
        
        Kuramoto3 model(config.dt, config.omega_3osc[0], config.omega_3osc[1], config.omega_3osc[2]);
        model.Inicie(config.initial_phases_3osc[0], config.initial_phases_3osc[1], config.initial_phases_3osc[2]);
        model.set_coupling(K);
        
        string K_str = format_number(K);
        
        // Abrir archivos para este K
        ofstream circle_file("results/datos/3osc/circle_coords_K" + K_str + ".dat");
        ofstream phase_3d_file("results/datos/3osc/phase_space_3d_K" + K_str + ".dat");
        ofstream order_file("results/datos/3osc/order_parameter_K" + K_str + ".dat");
        
        circle_file << "time cos1 cos2 cos3 sin1 sin2 sin3 r\n";
        phase_3d_file << "theta1 theta2 theta3\n";
        order_file << "time r\n";
        
        // Simular y escribir datos
        for (int step = 0; step < config.steps; step++) {
            double t = step * config.dt;
            model.Paso();
            
            // Coordenadas para animaciones cada 2 pasos
            if (step % 2 == 0) {
                double cos1 = cos(model.GetTheta1());
                double cos2 = cos(model.GetTheta2());
                double cos3 = cos(model.GetTheta3());
                double sin1 = sin(model.GetTheta1());
                double sin2 = sin(model.GetTheta2());
                double sin3 = sin(model.GetTheta3());
                double r = model.calculate_order_parameter();
                circle_file << t << " " << cos1 << " " << cos2 << " " << cos3 << " "
                          << sin1 << " " << sin2 << " " << sin3 << " " << r << "\n";
            }
            
            // Espacio de fases 3D
            phase_3d_file << model.GetTheta1() << " " << model.GetTheta2() << " " << model.GetTheta3() << "\n";
            
            // Parámetro de orden cada 5 pasos
            if (step % 5 == 0) {
                double r = model.calculate_order_parameter();
                order_file << t << " " << r << "\n";
            }
        }
        
        circle_file.close();
        phase_3d_file.close();
        order_file.close();
    }
}

// ========== ANÁLISIS DE BIFURCACIÓN 3 OSCILADORES ==========

void analyze_3osc_bifurcation(const SimulationConfig& config) {
    cout << "ANALISIS DE BIFURCACION 3 OSCILADORES" << endl;
    
    ofstream data_file("results/datos/3osc/bifurcation_data.dat");
    data_file << "K sum_derivatives\n";
    
    double min_sum = numeric_limits<double>::max();
    double K_min = 0.0;
    vector<pair<double, double>> data_points;
    
    // Barrido de K
    for (int i = 0; i <= config.n_K_points_3osc; ++i) {
        double K = (config.K_max_3osc * i) / config.n_K_points_3osc;
        
        // Simular hasta estado estacionario
        Kuramoto3 model(config.dt, config.omega_3osc[0], config.omega_3osc[1], config.omega_3osc[2]);
        model.Inicie(config.initial_phases_3osc[0], config.initial_phases_3osc[1], config.initial_phases_3osc[2]);
        model.set_coupling(K);
        
        for (int step = 0; step < config.steps; step++) {
            model.Paso();
        }
        
        // Obtener fases finales
        double theta1 = model.GetTheta1();
        double theta2 = model.GetTheta2();
        double theta3 = model.GetTheta3();
        
        // Calcular diferencias de fase
        double delta21 = theta2 - theta1;
        double delta31 = theta3 - theta1;
        double delta32 = theta3 - theta2;
        
        // Calcular derivadas
        double d21_dt = config.omega_3osc[1] - config.omega_3osc[0] - (K/3.0) * (sin(delta21) + sin(delta31 - delta21));
        double d31_dt = config.omega_3osc[2] - config.omega_3osc[0] - (K/3.0) * (sin(delta31) + sin(delta21 - delta31));
        double d32_dt = config.omega_3osc[2] - config.omega_3osc[1] - (K/3.0) * (sin(delta32) + sin(delta31 - delta32));
        
        // Suma de valores absolutos
        double sum_abs = abs(d21_dt) + abs(d31_dt) + abs(d32_dt);
        
        data_file << K << " " << sum_abs << "\n";
        data_points.push_back({K, sum_abs});
        
        // Encontrar mínimo
        if (sum_abs < min_sum) {
            min_sum = sum_abs;
            K_min = K;
        }
    }
    data_file.close();
    
    // Actualizar K_min en la configuración
    const_cast<SimulationConfig&>(config).K_min_3osc = K_min;
    
    // Encontrar punto más cercano a K_umbral
    double closest_K_to_umbral = data_points[0].first;
    double closest_sum_to_umbral = data_points[0].second;
    double min_distance = abs(config.K_umbral_3osc - closest_K_to_umbral);
    
    for (const auto& point : data_points) {
        double distance = abs(config.K_umbral_3osc - point.first);
        if (distance < min_distance) {
            min_distance = distance;
            closest_K_to_umbral = point.first;
            closest_sum_to_umbral = point.second;
        }
    }
    
    cout << "  K_umbral teorico: " << config.K_umbral_3osc << endl;
    cout << "  K_umbral en grafica: " << closest_K_to_umbral << endl;
    cout << "  K_minimo: " << K_min << endl;
    
    // Gráfica de bifurcación
    ofstream script("scripts/gnuplot/plot_bifurcation_3osc.gnu");
    script << "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/3osc/bifurcation_3osc.png'\n";
    script << "set xlabel 'Acoplamiento K'\n";
    script << "set ylabel 'Σ|dΔ/dt|'\n";
    script << "set title 'Bifurcacion 3 Osciladores'\n";
    script << "set grid\n";
    script << "set key top right\n";
    
    script << "plot 'results/datos/3osc/bifurcation_data.dat' u 1:2 w l lw 2 lc 'blue' t 'Σ|dΔ/dt|', \\\n";
    script << "     '-' w p pt 7 ps 2 lc 'red' t 'K_{umbral}', \\\n";
    script << "     '-' w p pt 7 ps 2 lc 'green' t 'K_{min}'\n";
    script << closest_K_to_umbral << " " << closest_sum_to_umbral << "\ne\n";
    script << K_min << " " << min_sum << "\ne\n";
    
    script.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_bifurcation_3osc.gnu");
    cout << "  Grafica de bifurcacion generada" << endl;
}

// ========== MAPAS DE SENSIBILIDAD ==========

void generate_sensitivity_maps(const SimulationConfig& config) {
    cout << "GENERANDO MAPAS DE SENSIBILIDAD" << endl;
    
    const int grid_size = 50;
    const double epsilon = 0.01;
    
    // Primera pasada: calcular min/max global para escala unificada
    double global_min_phase = numeric_limits<double>::max();
    double global_max_phase = numeric_limits<double>::min();
    double global_min_osc = numeric_limits<double>::max();
    double global_max_osc = numeric_limits<double>::min();
    
    for (double K : config.K_values_2osc) {
        for (int i = 0; i < grid_size; ++i) {
            double theta1_0 = 2 * M_PI * i / grid_size;
            for (int j = 0; j < grid_size; ++j) {
                double theta2_0 = 2 * M_PI * j / grid_size;
                
                // Simulación central
                Kuramoto2 model_center(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model_center.Inicie(theta1_0, theta2_0);
                model_center.set_coupling(K);
                for (int step = 0; step < 2000; step++) model_center.Paso();
                double theta1_center = model_center.GetTheta1();
                double theta2_center = model_center.GetTheta2();
                
                // Simulación con perturbación en theta1
                Kuramoto2 model_dtheta1(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model_dtheta1.Inicie(theta1_0 + epsilon, theta2_0);
                model_dtheta1.set_coupling(K);
                for (int step = 0; step < 2000; step++) model_dtheta1.Paso();
                double theta1_d1 = model_dtheta1.GetTheta1();
                double theta2_d1 = model_dtheta1.GetTheta2();
                
                // Simulación con perturbación en theta2
                Kuramoto2 model_dtheta2(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model_dtheta2.Inicie(theta1_0, theta2_0 + epsilon);
                model_dtheta2.set_coupling(K);
                for (int step = 0; step < 2000; step++) model_dtheta2.Paso();
                double theta1_d2 = model_dtheta2.GetTheta1();
                double theta2_d2 = model_dtheta2.GetTheta2();
                
                double sens_phase = (abs(theta1_d1 - theta1_center) + abs(theta2_d1 - theta2_center) +
                                   abs(theta1_d2 - theta1_center) + abs(theta2_d2 - theta2_center)) / (2.0 * epsilon);
                
                double sens_osc = (abs(cos(theta1_d1) - cos(theta1_center)) + abs(cos(theta2_d1) - cos(theta2_center)) +
                                 abs(cos(theta1_d2) - cos(theta1_center)) + abs(cos(theta2_d2) - cos(theta2_center))) / (2.0 * epsilon);
                
                global_min_phase = min(global_min_phase, sens_phase);
                global_max_phase = max(global_max_phase, sens_phase);
                global_min_osc = min(global_min_osc, sens_osc);
                global_max_osc = max(global_max_osc, sens_osc);
            }
        }
    }
    
    // Segunda pasada: generar datos con escala unificada
    for (double K : config.K_values_2osc) {
        cout << "  K = " << K << endl;
        
        string K_str = format_number(K);
        ofstream sens_phase_file("results/datos/2osc/sensitivity_phase_K" + K_str + ".dat");
        ofstream sens_osc_file("results/datos/2osc/sensitivity_osc_K" + K_str + ".dat");
        
        for (int i = 0; i < grid_size; ++i) {
            double theta1_0 = 2 * M_PI * i / grid_size;
            
            for (int j = 0; j < grid_size; ++j) {
                double theta2_0 = 2 * M_PI * j / grid_size;
                
                // Simulación central
                Kuramoto2 model_center(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model_center.Inicie(theta1_0, theta2_0);
                model_center.set_coupling(K);
                for (int step = 0; step < 2000; step++) model_center.Paso();
                double theta1_center = model_center.GetTheta1();
                double theta2_center = model_center.GetTheta2();
                
                // Simulación con perturbación en theta1
                Kuramoto2 model_dtheta1(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model_dtheta1.Inicie(theta1_0 + epsilon, theta2_0);
                model_dtheta1.set_coupling(K);
                for (int step = 0; step < 2000; step++) model_dtheta1.Paso();
                double theta1_d1 = model_dtheta1.GetTheta1();
                double theta2_d1 = model_dtheta1.GetTheta2();
                
                // Simulación con perturbación en theta2
                Kuramoto2 model_dtheta2(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model_dtheta2.Inicie(theta1_0, theta2_0 + epsilon);
                model_dtheta2.set_coupling(K);
                for (int step = 0; step < 2000; step++) model_dtheta2.Paso();
                double theta1_d2 = model_dtheta2.GetTheta1();
                double theta2_d2 = model_dtheta2.GetTheta2();
                
                double sens_phase = (abs(theta1_d1 - theta1_center) + abs(theta2_d1 - theta2_center) +
                                   abs(theta1_d2 - theta1_center) + abs(theta2_d2 - theta2_center)) / (2.0 * epsilon);
                
                double sens_osc = (abs(cos(theta1_d1) - cos(theta1_center)) + abs(cos(theta2_d1) - cos(theta2_center)) +
                                 abs(cos(theta1_d2) - cos(theta1_center)) + abs(cos(theta2_d2) - cos(theta2_center))) / (2.0 * epsilon);
                
                sens_phase_file << theta1_0 << " " << theta2_0 << " " << sens_phase << "\n";
                sens_osc_file << theta1_0 << " " << theta2_0 << " " << sens_osc << "\n";
            }
            sens_phase_file << "\n";
            sens_osc_file << "\n";
        }
        
        sens_phase_file.close();
        sens_osc_file.close();
    }
    
    // Gráfica de mapas de sensibilidad
    ofstream script("scripts/gnuplot/plot_sensitivity_maps.gnu");
    script << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/2osc/sensitivity_maps.png'\n";
    script << "set multiplot layout 2," << config.K_values_2osc.size() << "\n";
    script << "set pm3d map\n";
    script << "set xlabel 'θ₁⁰'\n";
    script << "set ylabel 'θ₂⁰'\n";
    script << "set cbrange [" << global_min_phase << ":" << global_max_phase << "]\n";
    
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script << "set title 'Sensibilidad Fases K=" << K << "'\n";
        script << "splot 'results/datos/2osc/sensitivity_phase_K" << K_str << ".dat' u 1:2:3 with pm3d\n";
    }
    
    script << "set cbrange [" << global_min_osc << ":" << global_max_osc << "]\n";
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script << "set title 'Sensibilidad Oscilaciones K=" << K << "'\n";
        script << "splot 'results/datos/2osc/sensitivity_osc_K" << K_str << ".dat' u 1:2:3 with pm3d\n";
    }
    script << "unset multiplot\n";
    script.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_sensitivity_maps.gnu");
}

// ========== MAPAS DE FRECUENCIAS ==========

void generate_frequency_maps(const SimulationConfig& config) {
    cout << "GENERANDO MAPAS DE FRECUENCIAS" << endl;
    
    const int grid_size = 60;
    
    // Primera pasada: calcular min/max global
    double global_min_phase = numeric_limits<double>::max();
    double global_max_phase = numeric_limits<double>::min();
    double global_min_osc = numeric_limits<double>::max();
    double global_max_osc = numeric_limits<double>::min();
    
    // Calcular todos los mínimos/máximos globales
    for (double K : config.K_values_2osc) {
        for (int i = 0; i < grid_size; ++i) {
            double theta1_0 = 2 * M_PI * i / grid_size;
            for (int j = 0; j < grid_size; ++j) {
                double theta2_0 = 2 * M_PI * j / grid_size;
                
                Kuramoto2 model(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model.Inicie(theta1_0, theta2_0);
                model.set_coupling(K);
                for (int step = 0; step < 2000; step++) model.Paso();
                
                double theta1 = model.GetTheta1();
                double theta2 = model.GetTheta2();
                
                // Calcular frecuencias instantáneas
                double freq1 = config.omega_2osc[0] + (K/2.0) * sin(theta2 - theta1);
                double freq2 = config.omega_2osc[1] + (K/2.0) * sin(theta1 - theta2);
                double freq_diff = abs(freq1 - freq2);
                
                double omega1_inst = config.omega_2osc[0] + (K/2.0) * sin(theta2 - theta1);
                double omega2_inst = config.omega_2osc[1] + (K/2.0) * sin(theta1 - theta2);
                double osc_vel_diff = abs(omega1_inst - omega2_inst);
                
                global_min_phase = min(global_min_phase, freq_diff);
                global_max_phase = max(global_max_phase, freq_diff);
                global_min_osc = min(global_min_osc, osc_vel_diff);
                global_max_osc = max(global_max_osc, osc_vel_diff);
            }
        }
    }
    
    // Segunda pasada: generar datos con escala unificada
    for (double K : config.K_values_2osc) {
        cout << "  K = " << K << endl;
        
        string K_str = format_number(K);
        ofstream freq_phase_file("results/datos/2osc/frequency_diff_K" + K_str + ".dat");
        ofstream freq_osc_file("results/datos/2osc/frequency_vel_K" + K_str + ".dat");
        
        for (int i = 0; i < grid_size; ++i) {
            double theta1_0 = 2 * M_PI * i / grid_size;
            
            for (int j = 0; j < grid_size; ++j) {
                double theta2_0 = 2 * M_PI * j / grid_size;
                
                Kuramoto2 model(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
                model.Inicie(theta1_0, theta2_0);
                model.set_coupling(K);
                for (int step = 0; step < 2000; step++) model.Paso();
                
                double theta1 = model.GetTheta1();
                double theta2 = model.GetTheta2();
                
                // Calcular frecuencias instantáneas
                double freq1 = config.omega_2osc[0] + (K/2.0) * sin(theta2 - theta1);
                double freq2 = config.omega_2osc[1] + (K/2.0) * sin(theta1 - theta2);
                double freq_diff = abs(freq1 - freq2);
                
                double omega1_inst = config.omega_2osc[0] + (K/2.0) * sin(theta2 - theta1);
                double omega2_inst = config.omega_2osc[1] + (K/2.0) * sin(theta1 - theta2);
                double osc_vel_diff = abs(omega1_inst - omega2_inst);
                
                freq_phase_file << theta1_0 << " " << theta2_0 << " " << freq_diff << "\n";
                freq_osc_file << theta1_0 << " " << theta2_0 << " " << osc_vel_diff << "\n";
            }
            freq_phase_file << "\n";
            freq_osc_file << "\n";
        }
        
        freq_phase_file.close();
        freq_osc_file.close();
    }
    
    // Gráfica con escala unificada
    ofstream script("scripts/gnuplot/plot_frequency_maps.gnu");
    script << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/2osc/frequency_maps.png'\n";
    script << "set multiplot layout 2," << config.K_values_2osc.size() << "\n";
    script << "set pm3d map\n";
    script << "set xlabel 'θ₁⁰'\n";
    script << "set ylabel 'θ₂⁰'\n";
    
    // Escala unificada para primera fila
    script << "set cbrange [" << global_min_phase << ":" << global_max_phase << "]\n";
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script << "set title 'Dif. Frecuencias K=" << K << "'\n";
        script << "splot 'results/datos/2osc/frequency_diff_K" << K_str << ".dat' u 1:2:3 with pm3d\n";
    }
    
    // Escala unificada para segunda fila  
    script << "set cbrange [" << global_min_osc << ":" << global_max_osc << "]\n";
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script << "set title 'Dif. Velocidades K=" << K << "'\n";
        script << "splot 'results/datos/2osc/frequency_vel_K" << K_str << ".dat' u 1:2:3 with pm3d\n";
    }
    script << "unset multiplot\n";
    script.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_frequency_maps.gnu");
}

// ========== ANÁLISIS DE BIFURCACIÓN 2 OSCILADORES ==========

void analyze_2osc_bifurcation(const SimulationConfig& config) {
    cout << "ANALIZANDO BIFURCACION 2 OSCILADORES" << endl;
    
    double K_c = calculate_critical_K_2osc(config.omega_2osc);
    
    // Curva continua con muchos puntos de K
    ofstream bifurcation_file("results/datos/2osc/bifurcation_data.dat");
    bifurcation_file << "K phase_diff_final\n";
    
    int n_points = 100;
    for (int i = 0; i <= n_points; ++i) {
        double K = 3.0 * K_c * i / n_points;
        
        Kuramoto2 model(config.dt, config.omega_2osc[0], config.omega_2osc[1]);
        model.Inicie(config.initial_phases_2osc[0], config.initial_phases_2osc[1]);
        model.set_coupling(K);
        
        // Simular hasta estado estacionario
        for (int step = 0; step < config.steps; step++) {
            model.Paso();
        }
        
        double phase_diff_final = abs(model.GetTheta2() - model.GetTheta1());
        bifurcation_file << K << " " << phase_diff_final << "\n";
        
        if (i % 20 == 0) {
            cout << "  K = " << K << ", |Δθ| = " << phase_diff_final << endl;
        }
    }
    bifurcation_file.close();
    
    // Gráfica de bifurcación
    ofstream script("scripts/gnuplot/plot_bifurcation_2osc.gnu");
    script << "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n";
    script << "set output 'results/visual/2osc/bifurcation_2osc.png'\n";
    script << "set xlabel 'Acoplamiento K'\n";
    script << "set ylabel '|θ₂ - θ₁| final (rad)'\n";
    script << "set title 'Diagrama de Bifurcacion 2 Osciladores'\n";
    script << "set grid\n";
    script << "plot 'results/datos/2osc/bifurcation_data.dat' u 1:2 w l lw 2 lc 'red' notitle\n";
    script.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_bifurcation_2osc.gnu");
}

// ========== ANÁLISIS DE PUNTOS FIJOS 3D ==========

void generate_3d_fixed_points(const SimulationConfig& config) {
    cout << "ANALIZANDO PUNTOS FIJOS 3D" << endl;
    
    const int grid_size = 50;
    
    for (double K : config.K_values_3osc) {
        cout << "  K = " << K << endl;
        
        string K_str = format_number(K);
        ofstream fixed_file("results/datos/3osc/fixed_points_3d_K" + K_str + ".dat");
        fixed_file << "delta21 delta31 derivative_magnitude\n";
        
        for (int i = 0; i < grid_size; ++i) {
            double delta21 = -M_PI + 2*M_PI * i / grid_size;
            
            for (int j = 0; j < grid_size; ++j) {
                double delta31 = -M_PI + 2*M_PI * j / grid_size;
                
                // Calcular derivadas usando ecuaciones de Kuramoto
                double d21_dt = - K * (sin(delta21) + sin(delta31 - delta21));
                double d31_dt = - K * (sin(delta31) + sin(delta21 - delta31));
                
                double derivative_magnitude = sqrt(d21_dt*d21_dt + d31_dt*d31_dt);
                
                fixed_file << delta21 << " " << delta31 << " " << derivative_magnitude << "\n";
            }
            fixed_file << "\n";
        }
        fixed_file.close();
    }
}