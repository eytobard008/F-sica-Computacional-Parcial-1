#include "Kuramoto.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;

// ========== ANIMACIONES CIRCULARES 2 OSCILADORES ==========

void generate_2osc_animations(const SimulationConfig& config) {
    cout << "GENERANDO ANIMACIONES 2 OSCILADORES" << endl;
    
    for (double K : config.K_values_2osc) {
        cout << "  K = " << K << endl;
        
        string K_str = format_number(K);
        
        ofstream script("scripts/gnuplot/animate_2osc_K" + K_str + ".gnu");
        script << "set terminal gif animate delay 10 size 800,800 enhanced font 'Arial,12'\n";
        script << "set output 'results/visual/2osc/animation_2osc_K" << K_str << ".gif'\n";
        script << "set size square\n";
        script << "set xrange [-1.2:1.2]\n";
        script << "set yrange [-1.2:1.2]\n";
        script << "set xlabel 'cos(θ)'\n";
        script << "set ylabel 'sin(θ)'\n";
        script << "set title 'Movimiento en Circulo Unitario (K=" << K << ")'\n";
        script << "set grid\n";
        script << "set key outside\n";
        
        // Leer número de líneas
        ifstream data_file("results/datos/2osc/circle_coords_K" + K_str + ".dat");
        string line;
        int total_lines = 0;
        getline(data_file, line); // skip header
        while (getline(data_file, line)) total_lines++;
        data_file.close();
        
        int step = max(1, total_lines / 80);
        
        script << "do for [i=1:" << total_lines-1 << ":" << step << "] {\n";
        script << "  # Circulo unitario como fondo\n";
        script << "  set object 1 circle at 0,0 size 1 fc rgb 'white' fs transparent solid 0.05 border rgb 'gray'\n";
        script << "  \n";
        script << "  plot 'results/datos/2osc/circle_coords_K" << K_str << ".dat' every ::1::i u 2:4 w l lc 'blue' lw 2 t 'Oscilador 1', \\\n";
        script << "       '' every ::1::i u 3:5 w l lc 'red' lw 2 t 'Oscilador 2', \\\n";
        script << "       '' every ::i::i u 2:4 w p pt 7 ps 3 lc 'blue' t '', \\\n";
        script << "       '' every ::i::i u 3:5 w p pt 7 ps 3 lc 'red' t '', \\\n";
        script << "       '' every ::i::i u (0):(0):(($2+$3)/2):(($4+$5)/2) w vectors head filled lc 'green' lw 3 t 'Vector r', \\\n";
        script << "  \n";
        script << "  unset object 1\n";
        script << "}\n";
        script.close();
        
        string gnuplot_cmd = "gnuplot scripts/gnuplot/animate_2osc_K" + K_str + ".gnu";
        int result = safe_system(gnuplot_cmd.c_str());
        
        if (result == 0) {
            cout << "  Animacion creada para K=" << K << endl;
        } else {
            cout << "  Error creando animacion para K=" << K << endl;
        }
    }
}

// ========== ANIMACIONES CIRCULARES 3 OSCILADORES ==========

void generate_3osc_animations(const SimulationConfig& config) {
    cout << "GENERANDO ANIMACIONES 3 OSCILADORES" << endl;
    
    for (double K : config.K_values_3osc) {
        cout << "  K = " << K << endl;
        
        string K_str = format_number(K);
        
        ofstream script("scripts/gnuplot/animate_3osc_K" + K_str + ".gnu");
        script << "set terminal gif animate delay 10 size 800,800 enhanced font 'Arial,12'\n";
        script << "set output 'results/visual/3osc/animation_3osc_K" << K_str << ".gif'\n";
        script << "set size square\n";
        script << "set xrange [-1.2:1.2]\n";
        script << "set yrange [-1.2:1.2]\n";
        script << "set xlabel 'cos(θ)'\n";
        script << "set ylabel 'sin(θ)'\n";
        script << "set title 'Movimiento en Circulo Unitario 3OSC (K=" << K << ")'\n";
        script << "set grid\n";
        script << "set key outside\n";
        
        // Leer número de líneas
        ifstream data_file("results/datos/3osc/circle_coords_K" + K_str + ".dat");
        string line;
        int total_lines = 0;
        getline(data_file, line); // skip header
        while (getline(data_file, line)) total_lines++;
        data_file.close();
        
        int step = max(1, total_lines / 80);
        
        script << "do for [i=1:" << total_lines-1 << ":" << step << "] {\n";
        script << "  # Circulo unitario como fondo\n";
        script << "  set object 1 circle at 0,0 size 1 fc rgb 'white' fs transparent solid 0.05 border rgb 'gray'\n";
        script << "  \n";
        // USAR COLUMNAS CORRECTAS (2:5, 3:6, 4:7 para cos:sin)
        script << "  plot 'results/datos/3osc/circle_coords_K" << K_str << ".dat' every ::1::i u 2:5 w l lc 'blue' lw 2 t 'Oscilador 1', \\\n";
        script << "       '' every ::1::i u 3:6 w l lc 'red' lw 2 t 'Oscilador 2', \\\n";
        script << "       '' every ::1::i u 4:7 w l lc 'green' lw 2 t 'Oscilador 3', \\\n";
        script << "       '' every ::i::i u 2:5 w p pt 7 ps 3 lc 'blue' t '', \\\n";
        script << "       '' every ::i::i u 3:6 w p pt 7 ps 3 lc 'red' t '', \\\n";
        script << "       '' every ::i::i u 4:7 w p pt 7 ps 3 lc 'green' t '', \\\n";
        script << "       '' every ::i::i u (0):(0):(($2+$3+$4)/3):(($5+$6+$7)/3) w vectors head filled lc 'purple' lw 3 t 'Vector r'\n";
        script << "  \n";
        script << "  unset object 1\n";
        script << "}\n";
        script.close();
        
        safe_system(("gnuplot scripts/gnuplot/animate_3osc_K" + K_str + ".gnu").c_str());
        cout << "  Animacion 3OSC creada para K=" << K << endl;
    }
}

// ========== GRÁFICAS DE SERIES TEMPORALES ==========

void plot_time_series(const SimulationConfig& config) {
    cout << "GENERANDO GRAFICAS DE SERIES TEMPORALES" << endl;
    
    // Para 2 osciladores
    ofstream script_2osc("scripts/gnuplot/plot_time_series_2osc.gnu");
    script_2osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_2osc << "set output 'results/visual/2osc/time_series_2osc.png'\n";
    script_2osc << "set multiplot layout 2,2\n";
    script_2osc << "set xlabel 'Tiempo (s)'\n";
    script_2osc << "set ylabel 'α = cos(θ)'\n";
    script_2osc << "set yrange [-1.2:1.2]\n";
    script_2osc << "set grid\n";
    
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script_2osc << "set title 'Oscilaciones K=" << K << "'\n";
        script_2osc << "plot 'results/datos/2osc/circle_coords_K" << K_str << ".dat' u 1:2 w l lw 2 lc 'blue' t 'α₁', ";
        script_2osc << "'' u 1:3 w l lw 2 lc 'red' t 'α₂'\n";
    }
    script_2osc << "unset multiplot\n";
    script_2osc.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_time_series_2osc.gnu");
    
    // Para 3 osciladores
    ofstream script_3osc("scripts/gnuplot/plot_time_series_3osc.gnu");
    script_3osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_3osc << "set output 'results/visual/3osc/time_series_3osc.png'\n";
    script_3osc << "set multiplot layout 2,2\n";
    script_3osc << "set xlabel 'Tiempo (s)'\n";
    script_3osc << "set ylabel 'α = cos(θ)'\n";
    script_3osc << "set yrange [-1.2:1.2]\n";
    script_3osc << "set grid\n";
    
    for (double K : config.K_values_3osc) {
        string K_str = format_number(K);
        script_3osc << "set title 'Oscilaciones K=" << K << "'\n";
        script_3osc << "plot 'results/datos/3osc/circle_coords_K" << K_str << ".dat' u 1:2 w l lw 2 lc 'blue' t 'α₁', ";
        script_3osc << "'' u 1:3 w l lw 2 lc 'green' t 'α₂', ";
        script_3osc << "'' u 1:4 w l lw 2 lc 'red' t 'α₃'\n";
    }
    script_3osc << "unset multiplot\n";
    script_3osc.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_time_series_3osc.gnu");
}

// ========== DIAGRAMAS DE ESTABILIDAD ==========

void plot_stability_diagrams(const SimulationConfig& config) {
    cout << "GENERANDO DIAGRAMAS DE ESTABILIDAD" << endl;
    
    double delta_omega = config.omega_2osc[1] - config.omega_2osc[0];
    double K_c = calculate_critical_K_2osc(config.omega_2osc);
    (void)K_c; // Silenciar warning variable no usada
    
    for (double K : config.K_values_2osc) {
        cout << "  K = " << K << endl;
        
        string K_str = format_number(K);
        
        ofstream stability_file("results/datos/2osc/stability_K" + K_str + ".dat");
        stability_file << "phi dphi_dt V\n";
        
        // Calcular ambas funciones: dφ/dt y V(φ)
        int n_points = 200;
        for (int i = 0; i < n_points; ++i) {
            double phi = 2.0 * M_PI * i / n_points;
            double dphi_dt = delta_omega - K * sin(phi);
            double V = -delta_omega * phi - K * cos(phi); // Potencial efectivo
            
            stability_file << phi << " " << dphi_dt << " " << V << "\n";
        }
        stability_file.close();
        
        // Guardar puntos fijos
        ofstream fixed_file("results/datos/2osc/fixed_points_K" + K_str + ".dat");
        fixed_file << "phi type\n";
        
        if (K > fabs(delta_omega)) {
            // Dos puntos fijos
            double phi_stable = asin(delta_omega / K);
            double phi_unstable = M_PI - phi_stable;
            
            fixed_file << phi_stable << " stable\n";
            fixed_file << phi_unstable << " unstable\n";
        } else if (fabs(K - fabs(delta_omega)) < 1e-10) {
            // Un punto marginal (semiestable)
            double phi_marginal = (delta_omega > 0) ? M_PI/2 : 3*M_PI/2;
            fixed_file << phi_marginal << " marginal\n";
        }
        fixed_file.close();
    }
    
    // Script MEJORADO - con líneas verticales
    ofstream script_stab("scripts/gnuplot/plot_stability_diagrams.gnu");
    script_stab << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_stab << "set output 'results/visual/2osc/stability_diagrams.png'\n";
    script_stab << "set multiplot layout 2,2\n";
    script_stab << "set xlabel 'Diferencia de Fase φ (rad)'\n";
    script_stab << "set ylabel 'dφ/dt'\n";
    script_stab << "set xrange [0:6.283]\n";
    script_stab << "set xtics ('0' 0, 'π/2' 1.57, 'π' 3.14, '3π/2' 4.71, '2π' 6.28)\n";
    script_stab << "set grid\n";
    
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script_stab << "set title 'Analisis de Estabilidad K=" << K << "'\n";
        script_stab << "set y2label 'Potencial V(φ)'\n";
        script_stab << "set y2tics\n";
        
        // Leer puntos fijos para agregar líneas verticales
        ifstream fixed_file("results/datos/2osc/fixed_points_K" + K_str + ".dat");
        string line;
        vector<double> stable_points, unstable_points;
        
        getline(fixed_file, line); // skip header
        while (getline(fixed_file, line)) {
            istringstream iss(line);
            double phi;
            string type;
            if (iss >> phi >> type) {
                if (type == "stable") stable_points.push_back(phi);
                else if (type == "unstable") unstable_points.push_back(phi);
            }
        }
        fixed_file.close();
        
        script_stab << "plot 'results/datos/2osc/stability_K" << K_str << ".dat' u 1:2 w l lw 3 lc 'blue' axis x1y1 t 'dφ/dt', ";
        script_stab << "'' u 1:3 w l lw 2 lc 'red' axis x1y2 t 'V(φ)', ";
        script_stab << "'results/datos/2osc/fixed_points_K" << K_str << ".dat' u 1:2 w p pt 7 ps 2 lc 'green' axis x1y1 t 'Puntos fijos', ";
        
        // Agregar líneas verticales para puntos fijos
        for (double phi : stable_points) {
            script_stab << "'' u (" << phi << "):(0) w l lw 2 lc 'green' dt 2 notitle, ";
        }
        for (double phi : unstable_points) {
            script_stab << "'' u (" << phi << "):(0) w l lw 2 lc 'orange' dt 2 notitle, ";
        }
        
        script_stab << "0 w l lc 'black' lw 1 axis x1y1 t ''\n";
    }
    script_stab << "unset multiplot\n";
    script_stab.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_stability_diagrams.gnu");
}

// ========== DIAGRAMAS DE ESPACIO DE FASES ==========

void plot_phase_spaces(const SimulationConfig& config) {
    cout << "GENERANDO DIAGRAMAS DE ESPACIO DE FASES" << endl;
    
    // Diagrama con todas las K en misma gráfica (theta1 vs theta2)
    ofstream script_combined("scripts/gnuplot/plot_phase_spaces_combined.gnu");
    script_combined << "set terminal pngcairo size 1000,800 enhanced font 'Arial,12'\n";
    script_combined << "set output 'results/visual/2osc/phase_spaces_combined.png'\n";
    script_combined << "set xlabel 'θ₁ (rad)'\n";
    script_combined << "set ylabel 'θ₂ (rad)'\n";
    script_combined << "set title 'Espacio de Fases - Todas las K'\n";
    script_combined << "set grid\n";
    script_combined << "set key outside right\n";
    
    script_combined << "plot ";
    bool first = true;
    vector<string> colors = {"red", "blue", "green", "orange", "purple", "brown"};
    
    for (size_t i = 0; i < config.K_values_2osc.size(); ++i) {
        double K = config.K_values_2osc[i];
        string K_str = format_number(K);
        string color = colors[i % colors.size()];
        
        if (!first) script_combined << ", ";
        script_combined << "'results/datos/2osc/phase_space_K" << K_str << ".dat' u 1:2 w l lw 1.5 lc '" 
                       << color << "' t 'K=" << K << "'";
        first = false;
    }
    script_combined << "\n";
    script_combined.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_phase_spaces_combined.gnu");
    
    // Diagramas individuales (alpha1 vs alpha2) - uno por cada K
    ofstream script_individual("scripts/gnuplot/plot_phase_spaces_individual.gnu");
    script_individual << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_individual << "set output 'results/visual/2osc/phase_spaces_individual.png'\n";
    script_individual << "set multiplot layout 2,2\n";
    script_individual << "set xlabel 'α₁ = cos(θ₁)'\n";
    script_individual << "set ylabel 'α₂ = cos(θ₂)'\n";
    script_individual << "set xrange [-1.1:1.1]\n";
    script_individual << "set yrange [-1.1:1.1]\n";
    script_individual << "set grid\n";
    
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script_individual << "set title 'K=" << K << "'\n";
        script_individual << "plot 'results/datos/2osc/circle_coords_K" << K_str << ".dat' u 2:3 w l lw 1.5 lc 'blue' notitle\n";
    }
    script_individual << "unset multiplot\n";
    script_individual.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_phase_spaces_individual.gnu");
}

// ========== GRÁFICAS DE PARÁMETRO DE ORDEN ==========

void plot_order_parameters(const SimulationConfig& config) {
    cout << "GENERANDO GRAFICAS DE PARAMETRO DE ORDEN" << endl;
    
    // Para 2 osciladores
    ofstream script_2osc("scripts/gnuplot/plot_order_parameters_2osc.gnu");
    script_2osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_2osc << "set output 'results/visual/2osc/order_parameters_2osc.png'\n";
    script_2osc << "set multiplot layout 2,2\n";
    script_2osc << "set xlabel 'Tiempo (s)'\n";
    script_2osc << "set ylabel 'Parametro de Orden r'\n";
    script_2osc << "set yrange [0:1.1]\n";
    script_2osc << "set grid\n";
    
    for (double K : config.K_values_2osc) {
        string K_str = format_number(K);
        script_2osc << "set title 'K=" << K << "'\n";
        script_2osc << "plot 'results/datos/2osc/order_parameter_K" << K_str << ".dat' u 1:2 w l lw 2 lc 'red' notitle\n";
    }
    script_2osc << "unset multiplot\n";
    script_2osc.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_order_parameters_2osc.gnu");
    
    // Para 3 osciladores
    ofstream script_3osc("scripts/gnuplot/plot_order_parameters_3osc.gnu");
    script_3osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_3osc << "set output 'results/visual/3osc/order_parameters_3osc.png'\n";
    script_3osc << "set multiplot layout 2,2\n";
    script_3osc << "set xlabel 'Tiempo (s)'\n";
    script_3osc << "set ylabel 'Parametro de Orden r'\n";
    script_3osc << "set yrange [0:1.1]\n";
    script_3osc << "set grid\n";
    
    for (double K : config.K_values_3osc) {
        string K_str = format_number(K);
        script_3osc << "set title 'K=" << K << "'\n";
        script_3osc << "plot 'results/datos/3osc/order_parameter_K" << K_str << ".dat' u 1:2 w l lw 2 lc 'red' notitle\n";
    }
    script_3osc << "unset multiplot\n";
    script_3osc.close();
    
    safe_system("gnuplot scripts/gnuplot/plot_order_parameters_3osc.gnu");
}

// ========== GRÁFICAS 3D INTERACTIVAS ==========

void generate_3d_plots() {
    cout << "GENERANDO GRAFICAS 3D INTERACTIVAS HTML" << endl;
    
    ofstream py_script("scripts/python/interactive_3d.py");
    py_script << "import plotly.graph_objects as go\n";
    py_script << "import numpy as np\n";
    py_script << "import os\n\n";
    
    py_script << "# 1. Espacio de fases 3D superpuesto\n";
    py_script << "fig = go.Figure()\n";
    py_script << "data_files = [f for f in os.listdir('results/datos/3osc/') if f.startswith('phase_space_3d_K')]\n";
    py_script << "colors = ['red', 'blue', 'green', 'orange']\n\n";
    
    py_script << "for idx, data_file in enumerate(data_files):\n";
    py_script << "    K_value = data_file.split('_K')[1].replace('.dat', '')\n";
    py_script << "    try:\n";
    py_script << "        data = np.loadtxt('results/datos/3osc/' + data_file, skiprows=1)\n";
    py_script << "        if len(data) > 0:\n";
    py_script << "            fig.add_trace(go.Scatter3d(x=data[:,0], y=data[:,1], z=data[:,2], mode='lines', name='K='+K_value, line=dict(color=colors[idx % len(colors)], width=4)))\n";
    py_script << "    except: pass\n\n";
    
    py_script << "fig.update_layout(title='Espacio de Fases 3D', scene=dict(xaxis_title='θ₁', yaxis_title='θ₂', zaxis_title='θ₃'))\n";
    py_script << "fig.write_html('results/visual/3osc/phase_space_3d.html')\n\n";
    
    py_script << "# 2. Puntos fijos 3D con marcadores en z=0\n";
    py_script << "fixed_files = [f for f in os.listdir('results/datos/3osc/') if f.startswith('fixed_points_3d_K')]\n";
    py_script << "for fixed_file in fixed_files:\n";
    py_script << "    K_value = fixed_file.split('_K')[1].replace('.dat', '')\n";
    py_script << "    try:\n";
    py_script << "        data = np.loadtxt('results/datos/3osc/' + fixed_file, skiprows=1)\n";
    py_script << "        n = int(len(data)**0.5)\n";
    py_script << "        X = data[:,0].reshape(n, n)\n";
    py_script << "        Y = data[:,1].reshape(n, n)\n";
    py_script << "        Z = data[:,2].reshape(n, n)\n";
    py_script << "        \n";
    py_script << "        # Marcar puntos donde z ≈ 0\n";
    py_script << "        zero_mask = np.abs(Z) < 0.1\n";
    py_script << "        zero_indices = np.where(zero_mask)\n";
    py_script << "        \n";
    py_script << "        fig2 = go.Figure()\n";
    py_script << "        fig2.add_trace(go.Surface(z=Z, x=X, y=Y, colorscale='Viridis', name='Superficie'))\n";
    py_script << "        \n";
    py_script << "        if len(zero_indices[0]) > 0:\n";
    py_script << "            x_zeros = X[zero_indices]\n";
    py_script << "            y_zeros = Y[zero_indices]\n";
    py_script << "            z_zeros = Z[zero_indices]\n";
    py_script << "            fig2.add_trace(go.Scatter3d(x=x_zeros, y=y_zeros, z=z_zeros, mode='markers', marker=dict(size=3, color='red'), name='Puntos fijos'))\n";
    py_script << "        \n";
    py_script << "        fig2.update_layout(title='Puntos Fijos 3D - K='+K_value, scene=dict(xaxis_title='Δ₂₁', yaxis_title='Δ₃₁', zaxis_title='|dΔ/dt|'))\n";
    py_script << "        fig2.write_html('results/visual/3osc/fixed_points_3d_K'+K_value+'.html')\n";
    py_script << "    except: pass\n\n";
    
    py_script << "print('Gráficas 3D HTML guardadas en results/visual/3osc/')\n";
    py_script.close();
    
    safe_system("python scripts/python/interactive_3d.py");
    cout << "  Gráficas 3D HTML creadas en results/visual/3osc/" << endl;

    py_script << "# 3. Espacio de fases con coordenadas alpha\n";
py_script << "alpha_files = [f for f in os.listdir('results/datos/3osc/') if f.startswith('alpha_phase_space_3d_K')]\n";
py_script << "for alpha_file in alpha_files:\n";
py_script << "    K_value = alpha_file.split('_K')[1].replace('.dat', '')\n";
py_script << "    try:\n";
py_script << "        data = np.loadtxt('results/datos/3osc/' + alpha_file, skiprows=1)\n";
py_script << "        if len(data) > 0:\n";
py_script << "            fig3 = go.Figure(data=[go.Scatter3d(x=data[:,0], y=data[:,1], z=data[:,2], mode='lines', line=dict(color='purple', width=4))])\n";
py_script << "            fig3.update_layout(title='Espacio de Fases Alpha - K='+K_value, scene=dict(xaxis_title='α₁', yaxis_title='α₂', zaxis_title='α₃', xaxis_range=[-1.1,1.1], yaxis_range=[-1.1,1.1], zaxis_range=[-1.1,1.1]))\n";
py_script << "            fig3.write_html('results/visual/3osc/alpha_phase_space_K'+K_value+'.html')\n";
py_script << "    except: pass\n\n";

py_script << "# 4. Espacio de fases de fases (theta) por separado\n";
py_script << "theta_files = [f for f in os.listdir('results/datos/3osc/') if f.startswith('phase_space_3d_K')]\n";
py_script << "for theta_file in theta_files:\n";
py_script << "    K_value = theta_file.split('_K')[1].replace('.dat', '')\n";
py_script << "    try:\n";
py_script << "        data = np.loadtxt('results/datos/3osc/' + theta_file, skiprows=1)\n";
py_script << "        if len(data) > 0:\n";
py_script << "            fig4 = go.Figure(data=[go.Scatter3d(x=data[:,0], y=data[:,1], z=data[:,2], mode='lines', line=dict(color='blue', width=4))])\n";
py_script << "            fig4.update_layout(title='Espacio de Fases Theta - K='+K_value, scene=dict(xaxis_title='θ₁', yaxis_title='θ₂', zaxis_title='θ₃'))\n";
py_script << "            fig4.write_html('results/visual/3osc/theta_phase_space_K'+K_value+'.html')\n";
py_script << "    except: pass\n\n";

py_script << "print('Gráficas 3D HTML guardadas en results/visual/3osc/')\n";
py_script.close();  // ESTO VA AL FINAL

safe_system("python scripts/python/interactive_3d.py");
cout << "  Gráficas 3D HTML creadas en results/visual/3osc/" << endl;
}

// ========== GRÁFICAS 3D DE ESPACIO DE FASES CON COORDENADAS ALPHA ==========

void generate_3d_alpha_phase_space(const SimulationConfig& config) {
    cout << "GENERANDO DATOS PARA GRAFICAS 3D ALPHA" << endl;
    
    // Generar datos para gráficas 3D de alpha
    for (double K : config.K_values_3osc) {
        string K_str = format_number(K);
        
        ifstream circle_file("results/datos/3osc/circle_coords_K" + K_str + ".dat");
        ofstream alpha_3d_file("results/datos/3osc/alpha_phase_space_3d_K" + K_str + ".dat");
        
        alpha_3d_file << "alpha1 alpha2 alpha3\n";
        
        string line;
        getline(circle_file, line); // skip header
        
        while (getline(circle_file, line)) {
            istringstream iss(line);
            double time, alpha1, alpha2, alpha3, sin1, sin2, sin3, r;
            if (iss >> time >> alpha1 >> alpha2 >> alpha3 >> sin1 >> sin2 >> sin3 >> r) {
                alpha_3d_file << alpha1 << " " << alpha2 << " " << alpha3 << "\n";
            }
        }
        
        circle_file.close();
        alpha_3d_file.close();
        cout << "  Datos 3D alpha generados para K=" << K << endl;
    }
}