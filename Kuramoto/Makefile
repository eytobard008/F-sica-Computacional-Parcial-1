CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm

SRCDIR = src
INCDIR = include
BINDIR = bin

SOURCES = $(SRCDIR)/Core.cpp $(SRCDIR)/Analisis.cpp $(SRCDIR)/Visualizacion.cpp $(SRCDIR)/main.cpp
TARGET = $(BINDIR)/kuramoto_system

all: $(TARGET)

$(TARGET): $(SOURCES)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET) $(LDFLAGS)
	@echo "✅ Sistema compilado: $(TARGET)"

run: $(TARGET)
	./$(TARGET)

view:
	@echo "📊 Abriendo carpeta de visualizaciones..."
	@if command -v xdg-open > /dev/null; then \
		xdg-open results/visual/2osc/; \
	elif command -v open > /dev/null; then \
		open results/visual/2osc/; \
	else \
		echo "Abre manualmente: results/visual/2osc/"; \
	fi

3d:
	@echo "🔄 Ejecutando visualizaciones 3D interactivas..."
	@python scripts/python/interactive_3d.py

clean:
	@echo "🧹 Limpiando archivos generados..."
	rm -rf $(BINDIR)
	find results/ -name "*.dat" -delete
	find results/ -name "*.png" -delete
	find results/ -name "*.gif" -delete
	find scripts/ -name "*.gnu" -delete
	find scripts/ -name "*.py" -delete
	find results/ -name "*.html" -delete
	@echo "✅ Limpieza completada (carpetas preservadas)"

.PHONY: all run clean view 3d