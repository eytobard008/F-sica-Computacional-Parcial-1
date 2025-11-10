/**
 * @file Core.cpp
 * @brief Implementación de la física fundamental y métodos de la clase Bola
 * @author Esteban Tobar
 * @date 29 de octubre de 2025
 */

#include "../include/Npepas.h"
#include <cmath>
#include <vector>
#include <ctime>

using namespace std;

// Variables globales para conteo de cálculos
long long calculosDistancias = 0;
long long calculosMalla = 0;
long long choquesEntreParticulas = 0;
long long calculosVecindad = 0;

// ==========================
// Implementación de la clase Bola
// ==========================

void Bola::Inicie(double x0, double y0, double vx0, double vy0, double R0) {
    x = x0;
    y = y0;
    vx = vx0;
    vy = vy0;
    R = R0;
    choquesPared = 0;
}

void Bola::Mover(double dt) {
    x += vx * dt;
    y += vy * dt;
}

void Bola::RebotarEnCaja(double W, double H, double& suma_vx_colisiones, double& suma_vy_colisiones) {
    // Rebote en pared izquierda
    if (x < R) {
        x = R;
        suma_vx_colisiones += 2.0 * abs(vx);
        vx = abs(vx);
        choquesPared++;
    }
    // Rebote en pared derecha
    if (x > W - R) {
        x = W - R;
        suma_vx_colisiones += 2.0 * abs(vx);
        vx = -abs(vx);
        choquesPared++;
    }
    // Rebote en pared inferior
    if (y < R) {
        y = R;
        suma_vy_colisiones += 2.0 * abs(vy);
        vy = abs(vy);
        choquesPared++;
    }
    // Rebote en pared superior
    if (y > H - R) {
        y = H - R;
        suma_vy_colisiones += 2.0 * abs(vy);
        vy = -abs(vy);
        choquesPared++;
    }
}

bool Bola::SolapaCon(const Bola &otra) const {
    double dx = x - otra.x;
    double dy = y - otra.y;
    double distancia2 = dx * dx + dy * dy;
    double sumaR = R + otra.R;
    return distancia2 < sumaR * sumaR;
}

void Bola::Colisionar(Bola &otra) {
    double dx = otra.x - x;
    double dy = otra.y - y;
    double dist2 = dx * dx + dy * dy;
    double sumaR = R + otra.R;
    double sumaR2 = sumaR * sumaR;

    if (dist2 <= sumaR2 && dist2 > 1e-12) {
        double dist = sqrt(dist2);

        // Vector normal unitario (apunta de esta partícula a la otra)
        double nx = dx / dist;
        double ny = dy / dist;

        // Vector tangente unitario (perpendicular al normal)
        double tx = -ny;
        double ty = nx;

        // Proyectar velocidades en componentes normal y tangencial
        double v1n = vx * nx + vy * ny;
        double v1t = vx * tx + vy * ty;
        double v2n = otra.vx * nx + otra.vy * ny;
        double v2t = otra.vx * tx + otra.vy * ty;

        // En colisión elástica, las componentes normales se intercambian
        // y las tangenciales se mantienen
        double v1n_nueva = v2n;
        double v2n_nueva = v1n;

        // Convertir de vuelta a coordenadas cartesianas
        vx = v1n_nueva * nx + v1t * tx;
        vy = v1n_nueva * ny + v1t * ty;
        otra.vx = v2n_nueva * nx + v2t * tx;
        otra.vy = v2n_nueva * ny + v2t * ty;

        // Separar partículas para evitar que se queden pegadas
        double separacion = 0.5 * (sumaR - dist + 0.001);
        x -= separacion * nx;
        y -= separacion * ny;
        otra.x += separacion * nx;
        otra.y += separacion * ny;

        choquesEntreParticulas++;
    }
}

// ==========================
// Implementación de malla espacial optimizada
// ==========================

void InicializarMalla(vector<Celda>& malla, int& celdasX, int& celdasY, vector<int>& offsetsVecinos, double W, double H, double tamanoCelda) {
    // Calcular número de celdas en cada dirección
    celdasX = static_cast<int>(ceil(W / tamanoCelda));
    celdasY = static_cast<int>(ceil(H / tamanoCelda));
    int totalCeldas = celdasX * celdasY;
    
    malla.resize(totalCeldas);
    
    // Precalcular offsets para los 9 vecinos (incluyendo celda central)
    // Esto evita calcular repetidamente los índices de celdas vecinas
    offsetsVecinos.clear();
    for (int dj = -1; dj <= 1; dj++) {
        for (int di = -1; di <= 1; di++) {
            offsetsVecinos.push_back(dj * celdasX + di);
        }
    }
}

void ActualizarMalla(Bola* bolas, int N, vector<Celda>& malla, int celdasX, int celdasY, vector<int>& celdasUsadas, double tamanoCelda, double W, double H, bool contarCalculos) {
    // Limpiar solo las celdas que fueron usadas en el frame anterior
    // Esto es más eficiente que limpiar todas las celdas
    for (int indice : celdasUsadas) {
        malla[indice].particulas.clear();
    }
    celdasUsadas.clear();

    // Asignar cada partícula a su celda correspondiente
    for (int i = 0; i < N; i++) {
        double x = bolas[i].GetX();
        double y = bolas[i].GetY();

        // Calcular índices de celda
        int celdaX = static_cast<int>(x / tamanoCelda);
        int celdaY = static_cast<int>(y / tamanoCelda);

        if (contarCalculos) {
            calculosMalla += 2; // Por los cálculos de celdaX y celdaY
        }

        // Verificar que la partícula esté dentro de los límites de la malla
        if (celdaX >= 0 && celdaX < celdasX && celdaY >= 0 && celdaY < celdasY) {
            int indice = celdaY * celdasX + celdaX;
            malla[indice].particulas.push_back(i);
            celdasUsadas.push_back(indice);
            
            if (contarCalculos) {
                calculosMalla += 2; // Por el cálculo del índice y seguimiento de celdas usadas
            }
        }
    }
}

void DetectarColisionesConMalla(Bola* bolas, int N, vector<Celda>& malla, int celdasX, int celdasY, const vector<int>& offsetsVecinos, const vector<int>& celdasUsadas, bool contarCalculos) {
    int totalCeldas = celdasX * celdasY;
    
    // Procesar solo las celdas que contienen partículas
    for (int indiceCelda : celdasUsadas) {
        Celda& celdaActual = malla[indiceCelda];
        int numParticulas = celdaActual.particulas.size();
        
        if (numParticulas == 0) continue;
        
        // Revisar colisiones entre partículas dentro de la misma celda
        // Se usa i+1 para evitar verificar el mismo par dos veces
        for (int i = 0; i < numParticulas; i++) {
            for (int j = i + 1; j < numParticulas; j++) {
                int idx1 = celdaActual.particulas[i];
                int idx2 = celdaActual.particulas[j];
                
                // Verificación rápida de colisión usando distancia al cuadrado
                double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                double dist2 = dx * dx + dy * dy;
                double sumaR = bolas[idx1].GetR() + bolas[idx2].GetR();

                if (contarCalculos) {
                    calculosDistancias++;
                }

                if (dist2 <= sumaR * sumaR) {
                    bolas[idx1].Colisionar(bolas[idx2]);
                }
            }
        }
        
        // Revisar colisiones con partículas en celdas vecinas
        // Esto es necesario porque partículas en celdas adyacentes pueden colisionar
        for (int offset : offsetsVecinos) {
            int indiceVecino = indiceCelda + offset;
            
            // Verificar que la celda vecina esté dentro de los límites
            if (indiceVecino < 0 || indiceVecino >= totalCeldas) {
                if (contarCalculos) {
                    calculosVecindad++; // 1 cálculo por verificación de límites
                }
                continue;
            }
            
            // Verificación de offset cero
            if (contarCalculos) {
                calculosVecindad++; // Por la verificación offset == 0
            }
            
            // No comparar con los mismos
            if (offset == 0) continue;
            
            Celda& celdaVecina = malla[indiceVecino];
            int numParticulasVecina = celdaVecina.particulas.size();
            
            if (numParticulasVecina == 0) continue;
            
            // Comparar todas las partículas de la celda actual con todas las de la vecina
            for (int i = 0; i < numParticulas; i++) {
                for (int j = 0; j < numParticulasVecina; j++) {
                    int idx1 = celdaActual.particulas[i];
                    int idx2 = celdaVecina.particulas[j];
                    
                    // Asegurar que idx1 < idx2 para evitar verificar el mismo par dos veces
                    // cuando se procese la celda vecina
                    if (idx1 >= idx2) continue;
                    
                    // Verificación rápida de colisión
                    double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                    double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                    double dist2 = dx * dx + dy * dy;
                    double sumaR = bolas[idx1].GetR() + bolas[idx2].GetR();

                    if (contarCalculos) {
                        calculosDistancias++;
                    }

                    if (dist2 <= sumaR * sumaR) {
                        bolas[idx1].Colisionar(bolas[idx2]);
                    }
                }
            }
        }
        
        // Contar operaciones de gestión de vecindad por celda activa
        if (contarCalculos) {
            calculosVecindad += 9; // Por los 9 offsets procesados por cada celda activa
            calculosVecindad += numParticulas; // Por el procesamiento interno de la celda
        }
    }
}
