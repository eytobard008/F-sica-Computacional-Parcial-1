/**
 * @file main.cpp
 * @brief Programa principal para la simulación de N partículas en una caja 2D
 * @author Esteban Tobar
 * @date 29 de octubre de 2025
 */

#include "../include/Npepas.h"
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <vector>
#include <ctime>

using namespace std;

int main() {
    // Parámetros fijos de la simulación
    double masa = 1.0;
    double dt = 0.01;
    double W = 20.0, H = 10.0;
    double R, tmax, velocidadMedia, dispersionVelocidad;
    int N;
    bool usarMallaOpt;
    bool contarCalculos = true;

    // Solicitar parámetros al usuario
    SolicitarDatos(R, N, tmax, velocidadMedia, dispersionVelocidad, usarMallaOpt);
    
    // Validar los datos ingresados
    if (!VerificarDatos(R, N, velocidadMedia, dispersionVelocidad)) {
        return 1;
    }

    // Inicializar generador de números aleatorios
    mt19937 gen(time(0));
    vector<Bola> bolas(N);

    // Inicializar sistema (usar disposición aleatoria por defecto)
    bool usarMalla = false;
    
    if (!InicializarSistema(bolas.data(), N, R, W, H, gen, velocidadMedia, dispersionVelocidad, usarMalla)) {
        return 1;
    }

    // Ejecutar simulación principal
    SimularSistema(bolas.data(), N, W, H, dt, tmax, "results/Npepas.dat", "results/energia.dat", "results/presion.dat", masa, velocidadMedia, dispersionVelocidad, usarMallaOpt, contarCalculos);

    // Generar visualizaciones y análisis
    int pasos = static_cast<int>(tmax / dt);
    GenerarAnimacion("results/Npepas.dat", "results/presion.dat", W, H, N, pasos, dt, R, usarMallaOpt, velocidadMedia, dispersionVelocidad, "Npepas.gif");
    GenerarHistogramaAnimado("results/Npepas.dat", N, pasos, dt, velocidadMedia, "histograma.gif", "histograma_zoom.gif");
    GenerarGraficaEnergia("results/energia.dat", "energia.png");
    GenerarGraficaPresion("results/presion.dat", "presion.png");

    // Mostrar estadísticas finales
    MostrarEstadisticas(bolas.data(), N, dt, tmax);

    cout << "\nSimulación completada.\n";
    cout << "Usa 'make view' para ver las animaciones y gráficas.\n";

    return 0;
}
