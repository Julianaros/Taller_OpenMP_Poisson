# Makefile para Taller OpenMP Poisson
# Compilador y flags
CXX = g++
CXXFLAGS = -fopenmp -O2 -Wall -std=c++17
LDFLAGS = -fopenmp

# Directorios
SRC_DIR = src
BIN_DIR = bin

# Archivos fuente
SOURCES = $(wildcard $(SRC_DIR)/actividad_*.cpp)
# Generar nombres de ejecutables: actividad_0, actividad_1, actividad_2, etc.
EXECUTABLES = $(patsubst $(SRC_DIR)/%.cpp,%,$(SOURCES))

# Carpetas que se generan al ejecutar (incluyendo actividad0)
OUTPUT_DIRS = actividad0 actividad1 actividad2 actividad3 actividad4 actividad5 actividad6 actividad7

# Regla por defecto
all: $(BIN_DIR) $(EXECUTABLES)

# Crear directorio bin si no existe
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Regla genérica para compilar cada actividad
actividad_%: $(SRC_DIR)/actividad_%.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)
	@echo "Compilado: $@ -> $(BIN_DIR)/$@"

# Reglas individuales (alternativa explícita)
actividad_0: $(SRC_DIR)/actividad_0.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

actividad_1: $(SRC_DIR)/actividad_1.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

actividad_2: $(SRC_DIR)/actividad_2.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

actividad_3: $(SRC_DIR)/actividad_3.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

actividad_4: $(SRC_DIR)/actividad_4.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

actividad_5: $(SRC_DIR)/actividad_5.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

actividad_6: $(SRC_DIR)/actividad_6.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

actividad_7: $(SRC_DIR)/actividad_7.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

# Ejecutar actividades
run-0: actividad_0
	cd $(BIN_DIR) && ./actividad_0

run-1: actividad_1
	cd $(BIN_DIR) && ./actividad_1

run-2: actividad_2
	cd $(BIN_DIR) && ./actividad_2

run-3: actividad_3
	cd $(BIN_DIR) && ./actividad_3

run-4: actividad_4
	cd $(BIN_DIR) && ./actividad_4

run-5: actividad_5
	cd $(BIN_DIR) && ./actividad_5

run-6: actividad_6
	cd $(BIN_DIR) && ./actividad_6

run-7: actividad_7
	cd $(BIN_DIR) && ./actividad_7

# Ejecutar todas las actividades en secuencia
run-all: all
	@echo "Ejecutando todas las actividades..."
	@$(MAKE) run-0
	@$(MAKE) run-1
	@$(MAKE) run-2
	@$(MAKE) run-3
	@$(MAKE) run-4
	@$(MAKE) run-5
	@$(MAKE) run-6
	@$(MAKE) run-7

# Limpiar archivos compilados y carpetas generadas
clean:
	@echo "Limpiando archivos compilados..."
	rm -rf $(BIN_DIR)
	@echo "Limpiando carpetas de salida..."
	rm -rf $(OUTPUT_DIRS)
	@echo "Limpieza completada."

# Limpiar solo ejecutables
clean-bin:
	rm -rf $(BIN_DIR)

# Limpiar solo carpetas de salida
clean-output:
	rm -rf $(OUTPUT_DIRS)

# Mostrar información del proyecto
info:
	@echo "=== INFORMACIÓN DEL PROYECTO ==="
	@echo "Directorio fuente: $(SRC_DIR)"
	@echo "Directorio binarios: $(BIN_DIR)"
	@echo "Compilador: $(CXX)"
	@echo "Flags de compilación: $(CXXFLAGS)"
	@echo "Flags de enlazado: $(LDFLAGS)"
	@echo ""
	@echo "Archivos fuente encontrados:"
	@for file in $(SOURCES); do echo "  $$file"; done
	@echo ""
	@echo "Ejecutables que se generarán:"
	@for exe in $(EXECUTABLES); do echo "  $(BIN_DIR)/$$exe"; done
	@echo ""
	@echo "Carpetas de salida:"
	@for dir in $(OUTPUT_DIRS); do echo "  $$dir"; done

# Verificar dependencias del sistema
check-deps:
	@echo "=== VERIFICANDO DEPENDENCIAS ==="
	@echo -n "Verificando g++... "
	@which g++ > /dev/null && echo "✓ Encontrado" || echo "✗ No encontrado"
	@echo -n "Verificando soporte OpenMP... "
	@echo '#include <omp.h>' | g++ -fopenmp -x c++ - -o /tmp/test_omp 2>/dev/null && echo "✓ Disponible" || echo "✗ No disponible"
	@rm -f /tmp/test_omp
	@echo -n "Verificando soporte C++17... "
	@echo 'int main(){return 0;}' | g++ -std=c++17 -x c++ - -o /tmp/test_cpp17 2>/dev/null && echo "✓ Disponible" || echo "✗ No disponible"
	@rm -f /tmp/test_cpp17

# Compilar en modo debug
debug: CXXFLAGS += -g -DDEBUG
debug: all

# Compilar en modo release (optimización máxima)
release: CXXFLAGS += -O3 -DNDEBUG
release: all

# Crear estructura de directorios del proyecto
setup:
	@echo "=== CREANDO ESTRUCTURA DEL PROYECTO ==="
	@mkdir -p $(SRC_DIR)
	@mkdir -p $(BIN_DIR)
	@echo "Directorios creados:"
	@echo "  $(SRC_DIR)/ - Para archivos fuente (.cpp)"
	@echo "  $(BIN_DIR)/ - Para ejecutables compilados"
	@echo ""
	@echo "Las carpetas de salida se crearán automáticamente al ejecutar:"
	@for dir in $(OUTPUT_DIRS); do echo "  $dir/ - Resultados de $dir"; done

# Declarar targets que no son archivos
.PHONY: all clean clean-bin clean-output run-0 run-1 run-2 run-3 run-4 run-5 run-6 run-7 run-all info check-deps debug release setup
