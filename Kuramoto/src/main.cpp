#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Kuramoto.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

/**
 * @brief Función principal del sistema de análisis Kuramoto
 * @return Código de salida del programa
 */
int main() {
    cout << "=========================================\n";
    cout << "     SISTEMA DE KURAMOTO - ANALISIS COMPLETO\n";
    cout << "=========================================\n\n";
    
    ConfiguracionSimulacion config;
    
    // ========== CONFIGURACIÓN DE PARÁMETROS GENERALES ==========
    
    config.dt = 0.01;
    config.pasos = 5000;
    
    // 2 Osciladores
    config.frecuencia1_2osc = 1.0;
    config.frecuencia2_2osc = 1.7;
    config.theta1_2osc = 0.0;
    config.theta2_2osc = M_PI/4;
    
    // 3 Osciladores  
    config.frecuencia1_3osc = 0.3;
    config.frecuencia2_3osc = 1.0;
    config.frecuencia3_3osc = 1.7;
    config.theta1_3osc = 0.0;
    config.theta2_3osc = M_PI/3;
    config.theta3_3osc = 2*M_PI/3;
    
    // ========== CONFIGURACIÓN DE PARÁMETROS K ==========
    
    config.K_umbral_2osc = calcularUmbralK2osc(config.frecuencia1_2osc, config.frecuencia2_2osc);
    config.K1_2osc = 0.0;
    config.K2_2osc = 0.93 * config.K_umbral_2osc;
    config.K3_2osc = config.K_umbral_2osc;
    config.K4_2osc = 3.0 * config.K_umbral_2osc;
    
    config.puntos_K_3osc = 2000;
    config.K_max_3osc = 4.0;
    config.K_umbral_3osc = calcularUmbralK3osc(config.frecuencia1_3osc, config.frecuencia2_3osc, config.frecuencia3_3osc);
    
    cout << "PARAMETROS CONFIGURADOS:" << endl;
    cout << "  2 Osciladores - K umbral: " << config.K_umbral_2osc << endl;
    cout << "  3 Osciladores - K umbral: " << config.K_umbral_3osc << endl;
    cout << "  3 Osciladores - Puntos K: " << config.puntos_K_3osc << endl;
    cout << "  3 Osciladores - K maximo: " << config.K_max_3osc << endl;
    
    // ========== EJECUCIÓN PRINCIPAL ==========
    
    crearDirectorios();
    
    cout << "\n1. ANALISIS 2 OSCILADORES..." << endl;
    analizarSistema2osc(config);
    
    cout << "\n2. BIFURCACION 3 OSCILADORES - DERIVADAS..." << endl;
    analizarBifurcacion3oscDerivada(config);
    
    config.K1_3osc = 0.0;
    config.K2_3osc = config.K_min_3osc;
    config.K3_3osc = config.K_umbral_3osc;
    config.K4_3osc = 3.0 * config.K_umbral_3osc;
    
    cout << "   Valores K 3osc configurados: " << config.K1_3osc << " " << config.K2_3osc << " " << config.K3_3osc << " " << config.K4_3osc << endl;
    
    cout << "\n3. ANALISIS 3 OSCILADORES..." << endl;
    analizarSistema3osc(config);
    
    cout << "\n4. MAPAS DE SENSIBILIDAD - FASES..." << endl;
    generarMapasSensibilidadFase(config);
    
    cout << "\n5. MAPAS DE SENSIBILIDAD - AMPLITUDES..." << endl;
    generarMapasSensibilidadAmplitud(config);
    
    cout << "\n6. MAPAS DE FRECUENCIA - FASES..." << endl;
    generarMapasFrecuenciaFase(config);
    
    cout << "\n7. MAPAS DE FRECUENCIA - VELOCIDADES..." << endl;
    generarMapasFrecuenciaVelocidad(config);
    
    cout << "\n8. ANIMACIONES 2OSC..." << endl;
    generarAnimaciones2osc(config);
    
    cout << "\n9. ANIMACIONES 3OSC..." << endl;
    generarAnimaciones3osc(config);
    
    cout << "\n10. PUNTOS FIJOS 3D..." << endl;
    generarPuntosFijos3D(config);
    
    cout << "\n11. BIFURCACION 2OSC - FASES..." << endl;
    analizarBifurcacion2oscFase(config);
    
    cout << "\n12. BIFURCACION 2OSC - DERIVADAS..." << endl;
    analizarBifurcacion2oscDerivada(config);
    
    cout << "\n13. BIFURCACION 3OSC - FASES..." << endl;
    analizarBifurcacion3oscFase(config);
    
    cout << "\n14. GRAFICAS DE SERIES TEMPORALES..." << endl;
    graficarSeriesTemporales(config);
    
    cout << "\n15. GRAFICAS DE FASES VS TIEMPO..." << endl;
    graficarFasesTiempo(config);
    
    cout << "\n16. DIAGRAMAS DE ESTABILIDAD..." << endl;
    graficarDiagramasEstabilidad(config);
    
    cout << "\n17. DIAGRAMAS DE ESPACIO DE FASES..." << endl;
    graficarEspaciosFase(config);
    
    cout << "\n18. PARAMETROS DE ORDEN..." << endl;
    graficarParametrosOrden(config);
    
    cout << "\n19. PREPARANDO DATOS 3D ALPHA..." << endl;
    generarEspacioFase3DAlpha(config);
    
    cout << "\n20. GRAFICAS 3D INTERACTIVAS..." << endl;
    generarGraficas3D();
    
    cout << "\n=========================================\n";
    cout << "     ✅ ANALISIS COMPLETADO\n";
    cout << "=========================================\n";
    
    cout << "\nRESUMEN DE PARAMETROS CALCULADOS:" << endl;
    cout << "  2 Osciladores - K umbral: " << config.K_umbral_2osc << endl;
    cout << "  3 Osciladores - K umbral: " << config.K_umbral_3osc << endl;
    cout << "  3 Osciladores - K minimo: " << config.K_min_3osc << endl;
    
    return 0;
}
