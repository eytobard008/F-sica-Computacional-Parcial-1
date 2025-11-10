#include "Kuramoto.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;

// ========== ANIMACIONES CIRCULARES ==========

/**
 * @brief Genera animaciones GIF para el sistema de 2 osciladores
 * @param config Configuración de la simulación
 */
void generarAnimaciones2osc(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO ANIMACIONES 2 OSCILADORES" << endl;
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        cout << "  K = " << K << endl;
        
        string K_str = formatearNumero(K);
        
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
        
        ifstream archivo_datos("results/datos/2osc/simulation_2osc_K" + K_str + ".dat");
        string linea;
        int total_lineas = 0;
        getline(archivo_datos, linea);
        while (getline(archivo_datos, linea)) total_lineas++;
        archivo_datos.close();
        
        int paso_anim = max(1, total_lineas / 80);
        
        script << "do for [i=1:" << total_lineas-1 << ":" << paso_anim << "] {\n";
        script << "  set object 1 circle at 0,0 size 1 fc rgb 'white' fs transparent solid 0.05 border rgb 'gray'\n";
        script << "  plot 'results/datos/2osc/simulation_2osc_K" << K_str << ".dat' every ::1::i u 5:7 w l lc 'blue' lw 2 t 'Oscilador 1', \\\n";
        script << "       '' every ::1::i u 6:8 w l lc 'red' lw 2 t 'Oscilador 2', \\\n";
        script << "       '' every ::i::i u 5:7 w p pt 7 ps 3 lc 'blue' t '', \\\n";
        script << "       '' every ::i::i u 6:8 w p pt 7 ps 3 lc 'red' t '', \\\n";
        script << "       '' every ::i::i u (0):(0):(($5+$6)/2):(($7+$8)/2) w vectors head filled lc 'green' lw 3 t 'Vector r'\n";
        script << "  unset object 1\n";
        script << "}\n";
        script.close();
        
        ejecutarComando(("gnuplot scripts/gnuplot/animate_2osc_K" + K_str + ".gnu").c_str());
        cout << "  Animacion creada para K=" << K << endl;
    }
}

/**
 * @brief Genera animaciones GIF para el sistema de 3 osciladores
 * @param config Configuración de la simulación
 */
void generarAnimaciones3osc(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO ANIMACIONES 3 OSCILADORES" << endl;
    
    double valores_K[4] = {config.K1_3osc, config.K2_3osc, config.K3_3osc, config.K4_3osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        cout << "  K = " << K << endl;
        
        string K_str = formatearNumero(K);
        
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
        
        ifstream archivo_datos("results/datos/3osc/simulation_3osc_K" + K_str + ".dat");
        string linea;
        int total_lineas = 0;
        getline(archivo_datos, linea);
        while (getline(archivo_datos, linea)) total_lineas++;
        archivo_datos.close();
        
        int paso_anim = max(1, total_lineas / 80);
        
        script << "do for [i=1:" << total_lineas-1 << ":" << paso_anim << "] {\n";
        script << "  set object 1 circle at 0,0 size 1 fc rgb 'white' fs transparent solid 0.05 border rgb 'gray'\n";
        script << "  plot 'results/datos/3osc/simulation_3osc_K" << K_str << ".dat' every ::1::i u 5:8 w l lc 'blue' lw 2 t 'Oscilador 1', \\\n";
        script << "       '' every ::1::i u 6:9 w l lc 'red' lw 2 t 'Oscilador 2', \\\n";
        script << "       '' every ::1::i u 7:10 w l lc 'green' lw 2 t 'Oscilador 3', \\\n";
        script << "       '' every ::i::i u 5:8 w p pt 7 ps 3 lc 'blue' t '', \\\n";
        script << "       '' every ::i::i u 6:9 w p pt 7 ps 3 lc 'red' t '', \\\n";
        script << "       '' every ::i::i u 7:10 w p pt 7 ps 3 lc 'green' t '', \\\n";
        script << "       '' every ::i::i u (0):(0):(($5+$6+$7)/3):(($8+$9+$10)/3) w vectors head filled lc 'purple' lw 3 t 'Vector r'\n";
        script << "  unset object 1\n";
        script << "}\n";
        script.close();
        
        ejecutarComando(("gnuplot scripts/gnuplot/animate_3osc_K" + K_str + ".gnu").c_str());
        cout << "  Animacion 3OSC creada para K=" << K << endl;
    }
}

// ========== GRÁFICAS DE SERIES TEMPORALES ==========

/**
 * @brief Genera gráficas de series temporales para las coordenadas alpha
 * @param config Configuración de la simulación
 */
void graficarSeriesTemporales(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO GRAFICAS DE SERIES TEMPORALES (ALPHA)" << endl;
    
    ofstream script_2osc("scripts/gnuplot/plot_time_series_2osc.gnu");
    script_2osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_2osc << "set output 'results/visual/2osc/time_series_2osc.png'\n";
    script_2osc << "set multiplot layout 2,2\n";
    script_2osc << "set xlabel 'Tiempo (s)'\n";
    script_2osc << "set ylabel 'α = cos(θ)'\n";
    script_2osc << "set yrange [-1.2:1.2]\n";
    script_2osc << "set grid\n";
    
    double valores_K_2osc[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K_2osc[idx_K];
        string K_str = formatearNumero(K);
        script_2osc << "set title 'Oscilaciones K=" << K << "'\n";
        script_2osc << "plot 'results/datos/2osc/simulation_2osc_K" << K_str << ".dat' u 1:5 w l lw 2 lc 'blue' t 'α₁', ";
        script_2osc << "'' u 1:6 w l lw 2 lc 'red' t 'α₂'\n";
    }
    script_2osc << "unset multiplot\n";
    script_2osc.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_time_series_2osc.gnu");
    
    ofstream script_3osc("scripts/gnuplot/plot_time_series_3osc.gnu");
    script_3osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_3osc << "set output 'results/visual/3osc/time_series_3osc.png'\n";
    script_3osc << "set multiplot layout 2,2\n";
    script_3osc << "set xlabel 'Tiempo (s)'\n";
    script_3osc << "set ylabel 'α = cos(θ)'\n";
    script_3osc << "set yrange [-1.2:1.2]\n";
    script_3osc << "set grid\n";
    
    double valores_K_3osc[4] = {config.K1_3osc, config.K2_3osc, config.K3_3osc, config.K4_3osc};
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K_3osc[idx_K];
        string K_str = formatearNumero(K);
        script_3osc << "set title 'Oscilaciones K=" << K << "'\n";
        script_3osc << "plot 'results/datos/3osc/simulation_3osc_K" << K_str << ".dat' u 1:5 w l lw 2 lc 'blue' t 'α₁', ";
        script_3osc << "'' u 1:6 w l lw 2 lc 'green' t 'α₂', ";
        script_3osc << "'' u 1:7 w l lw 2 lc 'red' t 'α₃'\n";
    }
    script_3osc << "unset multiplot\n";
    script_3osc.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_time_series_3osc.gnu");
}

// ========== GRÁFICAS DE FASES vs TIEMPO ==========

/**
 * @brief Genera gráficas de fases en función del tiempo
 * @param config Configuración de la simulación
 */
void graficarFasesTiempo(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO GRAFICAS DE FASES VS TIEMPO" << endl;
    
    ofstream script_2osc("scripts/gnuplot/plot_phases_time_2osc.gnu");
    script_2osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_2osc << "set output 'results/visual/2osc/phases_time_2osc.png'\n";
    script_2osc << "set multiplot layout 2,2\n";
    script_2osc << "set xlabel 'Tiempo (s)'\n";
    script_2osc << "set ylabel 'Fase θ (rad)'\n";
    script_2osc << "set grid\n";
    
    double valores_K_2osc[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K_2osc[idx_K];
        string K_str = formatearNumero(K);
        script_2osc << "set title 'Fases vs Tiempo K=" << K << "'\n";
        script_2osc << "plot 'results/datos/2osc/simulation_2osc_K" << K_str << ".dat' u 1:2 w l lw 2 lc 'blue' t 'θ₁', ";
        script_2osc << "'' u 1:3 w l lw 2 lc 'red' t 'θ₂'\n";
    }
    script_2osc << "unset multiplot\n";
    script_2osc.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_phases_time_2osc.gnu");
    
    ofstream script_3osc("scripts/gnuplot/plot_phases_time_3osc.gnu");
    script_3osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_3osc << "set output 'results/visual/3osc/phases_time_3osc.png'\n";
    script_3osc << "set multiplot layout 2,2\n";
    script_3osc << "set xlabel 'Tiempo (s)'\n";
    script_3osc << "set ylabel 'Fase θ (rad)'\n";
    script_3osc << "set grid\n";
    
    double valores_K_3osc[4] = {config.K1_3osc, config.K2_3osc, config.K3_3osc, config.K4_3osc};
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K_3osc[idx_K];
        string K_str = formatearNumero(K);
        script_3osc << "set title 'Fases vs Tiempo K=" << K << "'\n";
        script_3osc << "plot 'results/datos/3osc/simulation_3osc_K" << K_str << ".dat' u 1:2 w l lw 2 lc 'blue' t 'θ₁', ";
        script_3osc << "'' u 1:3 w l lw 2 lc 'green' t 'θ₂', ";
        script_3osc << "'' u 1:4 w l lw 2 lc 'red' t 'θ₃'\n";
    }
    script_3osc << "unset multiplot\n";
    script_3osc.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_phases_time_3osc.gnu");
}

// ========== DIAGRAMAS DE ESTABILIDAD ==========

/**
 * @brief Genera diagramas de estabilidad para el sistema de 2 osciladores
 * @param config Configuración de la simulación
 */
void graficarDiagramasEstabilidad(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO DIAGRAMAS DE ESTABILIDAD" << endl;
    
    double delta_frecuencia = config.frecuencia2_2osc - config.frecuencia1_2osc;
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        
        string K_str = formatearNumero(K);
        
        ofstream archivo_estabilidad("results/datos/2osc/stability_K" + K_str + ".dat");
        archivo_estabilidad << "# phi dphi_dt V\n";
        
        int n_puntos = 200;
        for (int i = 0; i < n_puntos; ++i) {
            double phi = 2.0 * M_PI * i / n_puntos;
            double dphi_dt = delta_frecuencia - K * sin(phi);
            double V = -delta_frecuencia * phi - K * cos(phi);
            
            archivo_estabilidad << phi << " " << dphi_dt << " " << V << "\n";
        }
        archivo_estabilidad.close();
        
        ofstream archivo_fijos("results/datos/2osc/fixed_points_K" + K_str + ".dat");
        archivo_fijos << "# phi type\n";
        
        if (K > fabs(delta_frecuencia)) {
            double phi_estable = asin(delta_frecuencia / K);
            double phi_inestable = M_PI - phi_estable;
            
            archivo_fijos << phi_estable << " stable\n";
            archivo_fijos << phi_inestable << " unstable\n";
        } else if (fabs(K - fabs(delta_frecuencia)) < 1e-10) {
            double phi_marginal = (delta_frecuencia > 0) ? M_PI/2 : 3*M_PI/2;
            archivo_fijos << phi_marginal << " marginal\n";
        } else {
            archivo_fijos << "0 no_fixed_points\n";
        }
        archivo_fijos.close();
    }
    
    ofstream script_estab("scripts/gnuplot/plot_stability_diagrams.gnu");
    script_estab << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_estab << "set output 'results/visual/2osc/stability_diagrams.png'\n";
    script_estab << "set multiplot layout 2,2\n";
    script_estab << "set xlabel 'Diferencia de Fase φ (rad)'\n";
    script_estab << "set ylabel 'dφ/dt'\n";
    script_estab << "set xrange [0:6.283]\n";
    script_estab << "set xtics ('0' 0, 'π/2' 1.57, 'π' 3.14, '3π/2' 4.71, '2π' 6.28)\n";
    script_estab << "set grid\n";
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        string K_str = formatearNumero(K);
        script_estab << "set title 'Analisis de Estabilidad K=" << K << "'\n";
        script_estab << "set y2label 'Potencial V(φ)'\n";
        script_estab << "set y2tics\n";
        
        ifstream archivo_fijos("results/datos/2osc/fixed_points_K" + K_str + ".dat");
        string linea;
        vector<double> puntos_estables, puntos_inestables;
        
        getline(archivo_fijos, linea);
        while (getline(archivo_fijos, linea)) {
            istringstream iss(linea);
            double phi;
            string tipo;
            if (iss >> phi >> tipo) {
                if (tipo == "stable") puntos_estables.push_back(phi);
                else if (tipo == "unstable") puntos_inestables.push_back(phi);
            }
        }
        archivo_fijos.close();
        
        script_estab << "plot 'results/datos/2osc/stability_K" << K_str << ".dat' u 1:2 w l lw 3 lc 'blue' axis x1y1 t 'dφ/dt', ";
        script_estab << "'' u 1:3 w l lw 2 lc 'red' axis x1y2 t 'V(φ)', ";
        script_estab << "'results/datos/2osc/fixed_points_K" << K_str << ".dat' u 1:2 w p pt 7 ps 2 lc 'green' axis x1y1 t 'Puntos fijos', ";
        
        for (double phi : puntos_estables) {
            script_estab << "'' u (" << phi << "):(0) w l lw 2 lc 'green' dt 2 notitle, ";
        }
        for (double phi : puntos_inestables) {
            script_estab << "'' u (" << phi << "):(0) w l lw 2 lc 'orange' dt 2 notitle, ";
        }
        
        script_estab << "0 w l lc 'black' lw 1 axis x1y1 t ''\n";
    }
    script_estab << "unset multiplot\n";
    script_estab.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_stability_diagrams.gnu");
}

// ========== DIAGRAMAS DE ESPACIO DE FASES ==========

/**
 * @brief Genera diagramas de espacio de fases para diferentes configuraciones
 * @param config Configuración de la simulación
 */
void graficarEspaciosFase(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO DIAGRAMAS DE ESPACIO DE FASES" << endl;
    
    ofstream script_combinado("scripts/gnuplot/plot_phase_spaces_combined.gnu");
    script_combinado << "set terminal pngcairo size 1000,800 enhanced font 'Arial,12'\n";
    script_combinado << "set output 'results/visual/2osc/phase_spaces_combined.png'\n";
    script_combinado << "set xlabel 'θ₁ (rad)'\n";
    script_combinado << "set ylabel 'θ₂ (rad)'\n";
    script_combinado << "set title 'Espacio de Fases - Todas las K'\n";
    script_combinado << "set grid\n";
    script_combinado << "set key outside right\n";
    
    script_combinado << "plot ";
    bool primero = true;
    vector<string> colores = {"red", "blue", "green", "orange"};
    
    double valores_K[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        string K_str = formatearNumero(K);
        string color = colores[idx_K % colores.size()];
        
        if (!primero) script_combinado << ", ";
        script_combinado << "'results/datos/2osc/simulation_2osc_K" << K_str << ".dat' u 2:3 w l lw 1.5 lc '" 
                       << color << "' t 'K=" << K << "'";
        primero = false;
    }
    script_combinado << "\n";
    script_combinado.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_phase_spaces_combined.gnu");
    
    ofstream script_individual("scripts/gnuplot/plot_phase_spaces_individual.gnu");
    script_individual << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_individual << "set output 'results/visual/2osc/phase_spaces_individual.png'\n";
    script_individual << "set multiplot layout 2,2\n";
    script_individual << "set xlabel 'α₁ = cos(θ₁)'\n";
    script_individual << "set ylabel 'α₂ = cos(θ₂)'\n";
    script_individual << "set xrange [-1.1:1.1]\n";
    script_individual << "set yrange [-1.1:1.1]\n";
    script_individual << "set grid\n";
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        string K_str = formatearNumero(K);
        script_individual << "set title 'K=" << K << "'\n";
        script_individual << "plot 'results/datos/2osc/simulation_2osc_K" << K_str << ".dat' u 5:6 w l lw 1.5 lc 'blue' notitle\n";
    }
    script_individual << "unset multiplot\n";
    script_individual.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_phase_spaces_individual.gnu");
}

// ========== GRÁFICAS DE PARÁMETRO DE ORDEN ==========

/**
 * @brief Genera gráficas del parámetro de orden en función del tiempo
 * @param config Configuración de la simulación
 */
void graficarParametrosOrden(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO GRAFICAS DE PARAMETRO DE ORDEN" << endl;
    
    ofstream script_2osc("scripts/gnuplot/plot_order_parameters_2osc.gnu");
    script_2osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_2osc << "set output 'results/visual/2osc/order_parameters_2osc.png'\n";
    script_2osc << "set multiplot layout 2,2\n";
    script_2osc << "set xlabel 'Tiempo (s)'\n";
    script_2osc << "set ylabel 'Parametro de Orden r'\n";
    script_2osc << "set yrange [0:1.1]\n";
    script_2osc << "set grid\n";
    
    double valores_K_2osc[4] = {config.K1_2osc, config.K2_2osc, config.K3_2osc, config.K4_2osc};
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K_2osc[idx_K];
        string K_str = formatearNumero(K);
        script_2osc << "set title 'K=" << K << "'\n";
        script_2osc << "plot 'results/datos/2osc/simulation_2osc_K" << K_str << ".dat' u 1:9 w l lw 2 lc 'red' notitle\n";
    }
    script_2osc << "unset multiplot\n";
    script_2osc.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_order_parameters_2osc.gnu");
    
    ofstream script_3osc("scripts/gnuplot/plot_order_parameters_3osc.gnu");
    script_3osc << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n";
    script_3osc << "set output 'results/visual/3osc/order_parameters_3osc.png'\n";
    script_3osc << "set multiplot layout 2,2\n";
    script_3osc << "set xlabel 'Tiempo (s)'\n";
    script_3osc << "set ylabel 'Parametro de Orden r'\n";
    script_3osc << "set yrange [0:1.1]\n";
    script_3osc << "set grid\n";
    
    double valores_K_3osc[4] = {config.K1_3osc, config.K2_3osc, config.K3_3osc, config.K4_3osc};
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K_3osc[idx_K];
        string K_str = formatearNumero(K);
        script_3osc << "set title 'K=" << K << "'\n";
        script_3osc << "plot 'results/datos/3osc/simulation_3osc_K" << K_str << ".dat' u 1:11 w l lw 2 lc 'red' notitle\n";
    }
    script_3osc << "unset multiplot\n";
    script_3osc.close();
    
    ejecutarComando("gnuplot scripts/gnuplot/plot_order_parameters_3osc.gnu");
}

// ========== GRÁFICAS 3D INTERACTIVAS ==========

/**
 * @brief Genera gráficas 3D interactivas en formato HTML
 */
void generarGraficas3D() {
    cout << "GENERANDO GRAFICAS 3D INTERACTIVAS HTML" << endl;
    
    ofstream script_python("scripts/python/interactive_3d.py");
    script_python << "import plotly.graph_objects as go\n";
    script_python << "import numpy as np\n";
    script_python << "import os\n\n";
    
    script_python << "fig = go.Figure()\n";
    script_python << "data_files = [f for f in os.listdir('results/datos/3osc/') if f.startswith('simulation_3osc_K')]\n";
    script_python << "colors = ['red', 'blue', 'green', 'orange']\n\n";
    
    script_python << "for idx, data_file in enumerate(data_files):\n";
    script_python << "    K_value = data_file.split('_K')[1].replace('.dat', '')\n";
    script_python << "    try:\n";
    script_python << "        data = np.loadtxt('results/datos/3osc/' + data_file, skiprows=1)\n";
    script_python << "        if len(data) > 0:\n";
    script_python << "            fig.add_trace(go.Scatter3d(x=data[:,1], y=data[:,2], z=data[:,3], mode='lines', name='K='+K_value, line=dict(color=colors[idx % len(colors)], width=4)))\n";
    script_python << "    except: pass\n\n";
    
    script_python << "fig.update_layout(title='Espacio de Fases 3D (Theta)', scene=dict(xaxis_title='θ₁', yaxis_title='θ₂', zaxis_title='θ₃'))\n";
    script_python << "fig.write_html('results/visual/3osc/phase_space_3d_theta.html')\n\n";
    
    script_python << "alpha_files = [f for f in os.listdir('results/datos/3osc/') if f.startswith('alpha_phase_space_3d_K')]\n";
    script_python << "for alpha_file in alpha_files:\n";
    script_python << "    K_value = alpha_file.split('_K')[1].replace('.dat', '')\n";
    script_python << "    try:\n";
    script_python << "        data = np.loadtxt('results/datos/3osc/' + alpha_file, skiprows=1)\n";
    script_python << "        if len(data) > 0:\n";
    script_python << "            fig_alpha = go.Figure(data=[go.Scatter3d(x=data[:,0], y=data[:,1], z=data[:,2], mode='lines', line=dict(color='purple', width=4))])\n";
    script_python << "            fig_alpha.update_layout(title='Espacio de Fases 3D (Alpha) - K='+K_value, scene=dict(xaxis_title='α₁', yaxis_title='α₂', zaxis_title='α₃', xaxis_range=[-1.1,1.1], yaxis_range=[-1.1,1.1], zaxis_range=[-1.1,1.1]))\n";
    script_python << "            fig_alpha.write_html('results/visual/3osc/phase_space_3d_alpha_K'+K_value+'.html')\n";
    script_python << "    except: pass\n\n";
    
    script_python << "fixed_files = [f for f in os.listdir('results/datos/3osc/') if f.startswith('fixed_points_3d_K')]\n";
    script_python << "for fixed_file in fixed_files:\n";
    script_python << "    K_value = fixed_file.split('_K')[1].replace('.dat', '')\n";
    script_python << "    try:\n";
    script_python << "        data = np.loadtxt('results/datos/3osc/' + fixed_file, skiprows=1)\n";
    script_python << "        if len(data) > 0:\n";
    script_python << "            n = int(len(data)**0.5)\n";
    script_python << "            X = data[:,0].reshape(n, n)\n";
    script_python << "            Y = data[:,1].reshape(n, n)\n";
    script_python << "            Z = data[:,2].reshape(n, n)\n";
    script_python << "            \n";
    script_python << "            zero_mask = np.abs(Z) < 0.1\n";
    script_python << "            zero_indices = np.where(zero_mask)\n";
    script_python << "            \n";
    script_python << "            fig2 = go.Figure()\n";
    script_python << "            fig2.add_trace(go.Surface(z=Z, x=X, y=Y, colorscale='Viridis', name='Superficie'))\n";
    script_python << "            \n";
    script_python << "            if len(zero_indices[0]) > 0:\n";
    script_python << "                x_zeros = X[zero_indices]\n";
    script_python << "                y_zeros = Y[zero_indices]\n";
    script_python << "                z_zeros = Z[zero_indices]\n";
    script_python << "                fig2.add_trace(go.Scatter3d(x=x_zeros, y=y_zeros, z=z_zeros, mode='markers', marker=dict(size=3, color='red'), name='Puntos fijos'))\n";
    script_python << "            \n";
    script_python << "            fig2.update_layout(title='Puntos Fijos 3D - K='+K_value, \n";
    script_python << "                              scene=dict(xaxis_title='Δ₂₁', yaxis_title='Δ₃₁', zaxis_title='|dΔ/dt|',\n";
    script_python << "                                        xaxis=dict(tickvals=[-3.14, -1.57, 0, 1.57, 3.14], ticktext=['-π', '-π/2', '0', 'π/2', 'π']),\n";
    script_python << "                                        yaxis=dict(tickvals=[-3.14, -1.57, 0, 1.57, 3.14], ticktext=['-π', '-π/2', '0', 'π/2', 'π'])))\n";
    script_python << "            fig2.write_html('results/visual/3osc/fixed_points_3d_K'+K_value+'.html')\n";
    script_python << "    except: pass\n\n";
    
    script_python << "print('Gráficas 3D HTML guardadas en results/visual/3osc/')\n";
    script_python.close();
    
    if (system("python3 --version > /dev/null 2>&1") == 0) {
        ejecutarComando("python3 scripts/python/interactive_3d.py");
    } else {
        cout << "  Python3 no disponible - omitiendo gráficas 3D interactivas" << endl;
    }
    cout << "  Gráficas 3D HTML creadas en results/visual/3osc/" << endl;
}

// ========== GRÁFICAS 3D DE ESPACIO DE FASES CON COORDENADAS ALPHA ==========

/**
 * @brief Genera datos para gráficas 3D del espacio de fases usando coordenadas alpha
 * @param config Configuración de la simulación
 */
void generarEspacioFase3DAlpha(const ConfiguracionSimulacion& config) {
    cout << "GENERANDO DATOS PARA GRAFICAS 3D ALPHA" << endl;
    
    double valores_K[4] = {config.K1_3osc, config.K2_3osc, config.K3_3osc, config.K4_3osc};
    
    for (int idx_K = 0; idx_K < 4; idx_K++) {
        double K = valores_K[idx_K];
        string K_str = formatearNumero(K);
        
        ifstream archivo_sim("results/datos/3osc/simulation_3osc_K" + K_str + ".dat");
        ofstream archivo_alpha_3d("results/datos/3osc/alpha_phase_space_3d_K" + K_str + ".dat");
        
        archivo_alpha_3d << "# alpha1 alpha2 alpha3\n";
        
        string linea;
        getline(archivo_sim, linea);
        
        while (getline(archivo_sim, linea)) {
            istringstream iss(linea);
            double tiempo, theta1, theta2, theta3, alpha1, alpha2, alpha3, sin1, sin2, sin3, r;
            if (iss >> tiempo >> theta1 >> theta2 >> theta3 >> alpha1 >> alpha2 >> alpha3 >> sin1 >> sin2 >> sin3 >> r) {
                archivo_alpha_3d << alpha1 << " " << alpha2 << " " << alpha3 << "\n";
            }
        }
        
        archivo_sim.close();
        archivo_alpha_3d.close();
        cout << "  Datos 3D alpha generados para K=" << K << endl;
    }
}
