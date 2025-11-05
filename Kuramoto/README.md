# Sistema de Kuramoto - Análisis Completo

Sistema de simulación y análisis del modelo de Kuramoto para osciladores acoplados.

## Características

- **2 Osciladores**: Simulación completa, mapas de sensibilidad, frecuencias, animaciones
- **3 Osciladores**: Espacios de fase 3D, bifurcaciones, puntos fijos, animaciones
- **Método Numérico**: Runge-Kutta de 4to orden (RK4)
- **Visualizaciones**: Gráficas 2D/3D, animaciones GIF, scripts interactivos

## Instalación y Uso

```bash
# Compilar
make

# Ejecutar análisis completo
make run

# Limpiar archivos generados
make clean

# Abrir visualizaciones en 2D
make view

# Abrir Gráficas 3D interactivas
make 3d

```

## Estructura del Proyecto

```bash

Kuramoto/
├── include/
│   └── Kuramoto.h
├── src/
│   ├── Core.cpp
│   ├── Analisis.cpp
│   ├── Visualizacion.cpp
│   └── main.cpp
├── documents/
│   ├── analisis.tex
│   └── analisis.pdf
├── results/
│   ├── datos/
│   └── visual/
├── scripts/
│   ├── gnuplot/
│   └── python/
├── Makefile
├── documents/
│ ├── Kuramoto.tex
│ └── Kuramoto.pdf
└── README.md

```

Configuración
Modificar parámetros en src/main.cpp:

```bash
cpp
config.dt = 0.01;
config.steps = 5000;
config.omega_2osc = {1.0, 1.7};
config.omega_3osc = {0.3, 1.0, 1.7};
```

## Programación Orientada a Objetos

- **Kuramoto2**: 2 osciladores con RK4
- **Kuramoto3**: 3 osciladores con RK4  
- **SimulationConfig**: Configuración centralizada

## Resultados

En `results/visual/` se generan:

- Series temporales
- Mapas de sensibilidad de fases, frecuencias, coordenadas angulares y velocidades angulares
- Animaciones circulares
- Diagramas de bifurcación
- Espacios de fase 2D/3D
- Análisis de estabilidad

## Requisitos

- g++ (C++11)
- Gnuplot
- Python 3 con Matplotlib (para gráficas 3D)
