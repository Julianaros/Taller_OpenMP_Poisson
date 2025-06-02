// Simulador Unificado de la Ecuación de Poisson en 2D
// Combina 5 ejemplos diferentes con sus respectivas condiciones de frontera y términos fuente
// Utiliza diferencias finitas y método iterativo de Jacobi/Gauss-Seidel

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <string>
#include <filesystem> // Para crear directorios

const double TOL = 1e-6;
const double e = 8.85e-12; // Permitividad eléctrica (para ejemplo original)

enum class Ejemplo {
    ORIGINAL = 0,    // Término fuente gaussiano
    EJEMPLO_1 = 1,   // ∇²V = (x² + y²)e^xy
    EJEMPLO_2 = 2,   // ∇²V = 0 (Laplace)
    EJEMPLO_3 = 3,   // ∇²V = 4
    EJEMPLO_4 = 4    // ∇²V = x/y + y/x
};

struct DominioConfig {
    double x_min, x_max;
    double y_min, y_max;
    std::string descripcion;
    std::string solucion_analitica;
};

// Función para crear el directorio data si no existe
void create_data_directory() {
    try {
        if (!std::filesystem::exists("data")) {
            std::filesystem::create_directory("data");
            std::cout << "Directorio 'data' creado exitosamente." << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error al crear el directorio 'data': " << ex.what() << std::endl;
    }
}

// Configuraciones de dominio para cada ejemplo
DominioConfig getDominioConfig(Ejemplo ejemplo) {
    switch (ejemplo) {
        case Ejemplo::ORIGINAL:
            return {0.0, 1.0, 0.0, 1.0, "Término fuente gaussiano", "No disponible"};
        case Ejemplo::EJEMPLO_1:
            return {0.0, 2.0, 0.0, 1.0, "∇²V = (x² + y²)e^xy", "V(x,y) = e^xy"};
        case Ejemplo::EJEMPLO_2:
            return {1.0, 2.0, 0.0, 1.0, "∇²V = 0 (Laplace)", "V(x,y) = ln(y² + x²)"};
        case Ejemplo::EJEMPLO_3:
            return {1.0, 2.0, 0.0, 2.0, "∇²V = 4", "V(x,y) = (x-y)²"};
        case Ejemplo::EJEMPLO_4:
            return {1.0, 2.0, 1.0, 2.0, "∇²V = x/y + y/x", "V(x,y) = xy ln(yx)"};
        default:
            return {0.0, 1.0, 0.0, 1.0, "Desconocido", "No disponible"};
    }
}

// Inicializa las condiciones de frontera según el ejemplo seleccionado
void initialize_boundary_conditions(Ejemplo ejemplo, int M, int N, std::vector<std::vector<double>>& V, 
                                   double h, double k, const DominioConfig& config) {
    V.resize(M + 1, std::vector<double>(N + 1, 0.0));
    
    switch (ejemplo) {
        case Ejemplo::ORIGINAL: {
            // Condiciones originales
            for (int j = 0; j <= N; ++j) {
                V[M][j] = 0.0;  // x = xf
                V[0][j] = 0.0;  // x = x0
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = 0.0;                    // y = y0
                V[i][N] = std::pow(x, 1.0);       // y = yf
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_1: {
            // Ejemplo 1: ∇²V = (x² + y²)e^xy
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = 1.0;                    // V(0,y) = 1
                V[M][j] = std::exp(2.0 * y);      // V(2,y) = e^(2y)
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = 1.0;                    // V(x,0) = 1
                V[i][N] = std::exp(x);            // V(x,1) = e^x
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_2: {
            // Ejemplo 2: Ecuación de Laplace
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = std::log(y*y + 1.0);         // V(1,y) = ln(y² + 1)
                V[M][j] = std::log(y*y + 4.0);         // V(2,y) = ln(y² + 4)
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = 2.0 * std::log(x);           // V(x,0) = 2ln(x)
                V[i][N] = std::log(x*x + 1.0);         // V(x,1) = ln(x² + 1)
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_3: {
            // Ejemplo 3: ∇²V = 4
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = (1.0 - y) * (1.0 - y);      // V(1,y) = (1-y)²
                V[M][j] = (2.0 - y) * (2.0 - y);      // V(2,y) = (2-y)²
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = x * x;                       // V(x,0) = x²
                V[i][N] = (x - 2.0) * (x - 2.0);       // V(x,2) = (x-2)²
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_4: {
            // Ejemplo 4: ∇²V = x/y + y/x
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = y * std::log(y);                    // V(1,y) = y ln(y)
                V[M][j] = 2.0 * y * std::log(2.0 * y);       // V(2,y) = 2y ln(2y)
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = x * std::log(x);                    // V(x,1) = x ln(x)
                V[i][N] = x * std::log(4.0 * x * x);          // V(x,2) = x ln(4x²)
            }
            break;
        }
    }
}

// Calcula el término fuente según el ejemplo seleccionado
void calculate_source_term(Ejemplo ejemplo, int M, int N, std::vector<std::vector<double>>& source, 
                          double h, double k, const DominioConfig& config) {
    source.resize(M + 1, std::vector<double>(N + 1, 0.0));
    
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            double x = config.x_min + i * h;
            double y = config.y_min + j * k;
            
            switch (ejemplo) {
                case Ejemplo::ORIGINAL: {
                    // Término fuente gaussiano
                    double mu = 0.5, sigma = 0.1;
                    source[i][j] = std::exp(-((x - mu)*(x - mu) + (y - mu)*(y - mu)) / 
                                          (2 * sigma * sigma)) / (std::sqrt(2 * M_PI * sigma));
                    break;
                }
                
                case Ejemplo::EJEMPLO_1: {
                    // f(x,y) = (x² + y²)e^xy
                    source[i][j] = (x*x + y*y) * std::exp(x * y);
                    break;
                }
                
                case Ejemplo::EJEMPLO_2: {
                    // f(x,y) = 0 (Ecuación de Laplace)
                    source[i][j] = 0.0;
                    break;
                }
                
                case Ejemplo::EJEMPLO_3: {
                    // f(x,y) = 4 (constante)
                    source[i][j] = 4.0;
                    break;
                }
                
                case Ejemplo::EJEMPLO_4: {
                    // f(x,y) = x/y + y/x
                    if (x > 0 && y > 0) {
                        source[i][j] = x/y + y/x;
                    } else {
                        source[i][j] = 0.0;
                    }
                    break;
                }
            }
            
            // Condiciones de contorno para la fuente (siempre 0 en los bordes)
            if (i == 0 || i == M || j == 0 || j == N) {
                source[i][j] = 0.0;
            }
        }
    }
}

// Resuelve la ecuación de Poisson iterativamente
int solve_poisson(Ejemplo ejemplo, std::vector<std::vector<double>>& V, 
                 const std::vector<std::vector<double>>& source, int M, int N, double h, double k) {
    double delta = 1.0;
    int iterations = 0;
    
    // Factor para el ejemplo original (con permitividad eléctrica)
    double factor = (ejemplo == Ejemplo::ORIGINAL) ? (1.0 / e) : 1.0;
    
    while (delta > TOL) {
        delta = 0.0;
        iterations++;
        
        for (int i = 1; i < M; ++i) {
            for (int j = 1; j < N; ++j) {
                double V_new;
                
                if (ejemplo == Ejemplo::ORIGINAL) {
                    // Fórmula original con signo negativo
                    V_new = (
                        ((V[i + 1][j] + V[i - 1][j]) * k * k) +
                        ((V[i][j + 1] + V[i][j - 1]) * h * h) -
                        (factor * source[i][j] * h * h * k * k)) /
                        (2.0 * (h * h + k * k));
                } else if (ejemplo == Ejemplo::EJEMPLO_3) {
                    // Ejemplo 3 con signo negativo
                    V_new = (
                        ((V[i + 1][j] + V[i - 1][j]) * k * k) +
                        ((V[i][j + 1] + V[i][j - 1]) * h * h) -
                        (source[i][j] * h * h * k * k)) /
                        (2.0 * (h * h + k * k));
                } else {
                    // Otros ejemplos con signo positivo
                    V_new = (
                        ((V[i + 1][j] + V[i - 1][j]) * k * k) +
                        ((V[i][j + 1] + V[i][j - 1]) * h * h) +
                        (source[i][j] * h * h * k * k)) /
                        (2.0 * (h * h + k * k));
                }
                
                delta = std::max(delta, std::abs(V_new - V[i][j]));
                V[i][j] = V_new;
            }
        }
    }
    
    return iterations;
}

// Calcula la solución analítica (donde esté disponible)
void calculate_analytical_solution(Ejemplo ejemplo, int M, int N, std::vector<std::vector<double>>& V_analytical, 
                                  double h, double k, const DominioConfig& config) {
    V_analytical.resize(M + 1, std::vector<double>(N + 1, 0.0));
    
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            double x = config.x_min + i * h;
            double y = config.y_min + j * k;
            
            switch (ejemplo) {
                case Ejemplo::EJEMPLO_1:
                    V_analytical[i][j] = std::exp(x * y);
                    break;
                case Ejemplo::EJEMPLO_2:
                    V_analytical[i][j] = std::log(y*y + x*x);
                    break;
                case Ejemplo::EJEMPLO_3:
                    V_analytical[i][j] = (x - y) * (x - y);
                    break;
                case Ejemplo::EJEMPLO_4:
                    if (x > 0 && y > 0) {
                        V_analytical[i][j] = x * y * std::log(y * x);
                    }
                    break;
                default:
                    V_analytical[i][j] = 0.0; // No hay solución analítica disponible
                    break;
            }
        }
    }
}

// Calcula el error entre la solución numérica y analítica
double calculate_error(const std::vector<std::vector<double>>& V_numerical, 
                      const std::vector<std::vector<double>>& V_analytical, 
                      int M, int N) {
    double max_error = 0.0;
    double sum_squared_error = 0.0;
    int count = 0;
    
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            double error = std::abs(V_numerical[i][j] - V_analytical[i][j]);
            max_error = std::max(max_error, error);
            sum_squared_error += error * error;
            count++;
        }
    }
    
    double rms_error = std::sqrt(sum_squared_error / count);
    
    std::cout << "Error máximo: " << max_error << std::endl;
    std::cout << "Error RMS: " << rms_error << std::endl;
    
    return max_error;
}

// Exporta los resultados a un archivo .dat en la carpeta data
void export_to_file(const std::vector<std::vector<double>>& V, double h, double k, int M, int N, 
                   const std::string& filename, const DominioConfig& config) {
    // Crear la ruta completa con la carpeta data
    std::string full_path = "data/" + filename;
    
    std::ofstream file(full_path);
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo " << full_path << " para escritura." << std::endl;
        return;
    }
    
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            double x = config.x_min + i * h;
            double y = config.y_min + j * k;
            file << x << "\t" << y << "\t" << V[i][j] << "\n";
        }
    }
    file.close();
    std::cout << "Resultados exportados a " << full_path << std::endl;
}

// Función principal que ejecuta la simulación para el ejemplo seleccionado
void run_simulation(Ejemplo ejemplo, int M = 50, int N = 50) {
    DominioConfig config = getDominioConfig(ejemplo);
    
    double h = (config.x_max - config.x_min) / M;
    double k = (config.y_max - config.y_min) / N;
    
    std::vector<std::vector<double>> V, source, V_analytical;
    
    // Mostrar información del ejemplo
    std::cout << "\n=== EJEMPLO " << static_cast<int>(ejemplo) << ": " << config.descripcion << " ===" << std::endl;
    std::cout << "Dominio: x ∈ [" << config.x_min << ", " << config.x_max << "], y ∈ [" 
              << config.y_min << ", " << config.y_max << "]" << std::endl;
    std::cout << "Grilla: " << (M+1) << " x " << (N+1) << " puntos" << std::endl;
    std::cout << "Tolerancia: " << TOL << std::endl;
    if (config.solucion_analitica != "No disponible") {
        std::cout << "Solución analítica: " << config.solucion_analitica << std::endl;
    }
    std::cout << std::endl;
    
    // Inicialización y resolución
    auto start = std::chrono::high_resolution_clock::now();
    initialize_boundary_conditions(ejemplo, M, N, V, h, k, config);
    calculate_source_term(ejemplo, M, N, source, h, k, config);
    
    std::cout << "Resolviendo la ecuación..." << std::endl;
    int iterations = solve_poisson(ejemplo, V, source, M, N, h, k);
    auto end = std::chrono::high_resolution_clock::now();
    
    // Cálculo del tiempo de ejecución
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double tiempo_segundos = duration.count() / 1000.0;
    
    // Mostrar resultados
    std::cout << "\n=== RESULTADOS ===" << std::endl;
    std::cout << "Tiempo de ejecución: " << tiempo_segundos << " segundos" << std::endl;
    std::cout << "Número de iteraciones: " << iterations << std::endl;
    std::cout << "Número de threads: 1 (secuencial)" << std::endl;
    
    // Análisis de error (si hay solución analítica disponible)
    if (config.solucion_analitica != "No disponible") {
        calculate_analytical_solution(ejemplo, M, N, V_analytical, h, k, config);
        std::cout << "\n=== ANÁLISIS DE ERROR ===" << std::endl;
        calculate_error(V, V_analytical, M, N);
        
        // Exportar solución analítica también
        std::string analytical_filename = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + "_analytical.dat";
        export_to_file(V_analytical, h, k, M, N, analytical_filename, config);
    }
    
    // Exportar solución numérica
    std::string numerical_filename = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + "_numerical.dat";
    export_to_file(V, h, k, M, N, numerical_filename, config);
    
    std::cout << "\nSimulación completada." << std::endl;
}

int main() {
    std::cout << "=== SIMULADOR DE ECUACIÓN DE POISSON 2D (SERIAL) ===" << std::endl;
    
    // Crear el directorio data al inicio del programa
    create_data_directory();
    
    std::cout << "\nEjemplos disponibles:" << std::endl;
    std::cout << "0 - Ejemplo Original: Término fuente gaussiano" << std::endl;
    std::cout << "1 - Ejemplo 1: ∇²V = (x² + y²)e^xy" << std::endl;
    std::cout << "2 - Ejemplo 2: ∇²V = 0 (Ecuación de Laplace)" << std::endl;
    std::cout << "3 - Ejemplo 3: ∇²V = 4" << std::endl;
    std::cout << "4 - Ejemplo 4: ∇²V = x/y + y/x" << std::endl;
    std::cout << "5 - Ejecutar todos los ejemplos" << std::endl;
    
    int opcion;
    std::cout << "\nSeleccione una opción (0-5): ";
    std::cin >> opcion;
    
    if (opcion >= 0 && opcion <= 4) {
        // Ejecutar ejemplo específico
        run_simulation(static_cast<Ejemplo>(opcion));
    } else if (opcion == 5) {
        // Ejecutar todos los ejemplos
        std::cout << "\n=== EJECUTANDO TODOS LOS EJEMPLOS ===\n" << std::endl;
        for (int i = 0; i <= 4; ++i) {
            run_simulation(static_cast<Ejemplo>(i));
            std::cout << "\n" << std::string(60, '=') << std::endl;
        }
    } else {
        std::cout << "Opción inválida. Ejecutando ejemplo 1 por defecto." << std::endl;
        run_simulation(Ejemplo::EJEMPLO_1);
    }
    
    return 0;
}