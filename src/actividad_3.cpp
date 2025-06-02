// Actividad 3: sections
// Ejecutar funciones independientes en paralelo con #pragma omp sections

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <string>
#include <filesystem>
#include <omp.h>

const double TOL = 1e-4;
const double e = 8.85e-12;

enum class Ejemplo {
    ORIGINAL = 0,
    EJEMPLO_1 = 1,
    EJEMPLO_2 = 2,
    EJEMPLO_3 = 3,
    EJEMPLO_4 = 4
};

struct DominioConfig {
    double x_min, x_max;
    double y_min, y_max;
    std::string descripcion;
    std::string solucion_analitica;
};

void create_data_directory() {
    try {
        if (!std::filesystem::exists("actividad3")) {
            std::filesystem::create_directory("actividad3");
            std::cout << "Directorio 'actividad3' creado exitosamente." << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error al crear el directorio 'actividad3': " << ex.what() << std::endl;
    }
}

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

void initialize_boundary_conditions(Ejemplo ejemplo, int M, int N, std::vector<std::vector<double>>& V, 
                                   double h, double k, const DominioConfig& config) {
    V.resize(M + 1, std::vector<double>(N + 1, 0.0));
    
    switch (ejemplo) {
        case Ejemplo::ORIGINAL: {
            for (int j = 0; j <= N; ++j) {
                V[M][j] = 0.0;
                V[0][j] = 0.0;
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = 0.0;
                V[i][N] = std::pow(x, 1.0);
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_1: {
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = 1.0;
                V[M][j] = std::exp(2.0 * y);
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = 1.0;
                V[i][N] = std::exp(x);
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_2: {
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = std::log(y*y + 1.0);
                V[M][j] = std::log(y*y + 4.0);
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = 2.0 * std::log(x);
                V[i][N] = std::log(x*x + 1.0);
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_3: {
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = (1.0 - y) * (1.0 - y);
                V[M][j] = (2.0 - y) * (2.0 - y);
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = x * x;
                V[i][N] = (x - 2.0) * (x - 2.0);
            }
            break;
        }
        
        case Ejemplo::EJEMPLO_4: {
            for (int j = 0; j <= N; ++j) {
                double y = config.y_min + j * k;
                V[0][j] = y * std::log(y);
                V[M][j] = 2.0 * y * std::log(2.0 * y);
            }
            for (int i = 0; i <= M; ++i) {
                double x = config.x_min + i * h;
                V[i][0] = x * std::log(x);
                V[i][N] = x * std::log(4.0 * x * x);
            }
            break;
        }
    }
}

void calculate_source_term(Ejemplo ejemplo, int M, int N, std::vector<std::vector<double>>& source, 
                          double h, double k, const DominioConfig& config) {
    source.resize(M + 1, std::vector<double>(N + 1, 0.0));
    
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            double x = config.x_min + i * h;
            double y = config.y_min + j * k;
            
            switch (ejemplo) {
                case Ejemplo::ORIGINAL: {
                    double mu = 0.5, sigma = 0.1;
                    source[i][j] = std::exp(-((x - mu)*(x - mu) + (y - mu)*(y - mu)) / 
                                          (2 * sigma * sigma)) / (std::sqrt(2 * M_PI * sigma));
                    break;
                }
                
                case Ejemplo::EJEMPLO_1: {
                    source[i][j] = (x*x + y*y) * std::exp(x * y);
                    break;
                }
                
                case Ejemplo::EJEMPLO_2: {
                    source[i][j] = 0.0;
                    break;
                }
                
                case Ejemplo::EJEMPLO_3: {
                    source[i][j] = 4.0;
                    break;
                }
                
                case Ejemplo::EJEMPLO_4: {
                    if (x > 0 && y > 0) {
                        source[i][j] = x/y + y/x;
                    } else {
                        source[i][j] = 0.0;
                    }
                    break;
                }
            }
            
            if (i == 0 || i == M || j == 0 || j == N) {
                source[i][j] = 0.0;
            }
        }
    }
}

// ACTIVIDAD 3: Función principal que utiliza sections para paralelizar la inicialización
void initialize_grid_and_source_parallel(Ejemplo ejemplo, int M, int N, 
                                        std::vector<std::vector<double>>& V,
                                        std::vector<std::vector<double>>& source,
                                        double h, double k, const DominioConfig& config) {
    
    auto start_init = std::chrono::high_resolution_clock::now();
    
    // PARALELIZACIÓN CON SECTIONS: Ejecutar funciones independientes en paralelo
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            std::cout << "Thread " << omp_get_thread_num() << " ejecutando initialize_boundary_conditions" << std::endl;
            initialize_boundary_conditions(ejemplo, M, N, V, h, k, config);
        }
        
        #pragma omp section
        {
            std::cout << "Thread " << omp_get_thread_num() << " ejecutando calculate_source_term" << std::endl;
            calculate_source_term(ejemplo, M, N, source, h, k, config);
        }
    }
    
    auto end_init = std::chrono::high_resolution_clock::now();
    auto duration_init = std::chrono::duration_cast<std::chrono::milliseconds>(end_init - start_init);
    double tiempo_init = duration_init.count() / 1000.0;
    
    std::cout << "Tiempo de inicialización paralela (sections): " << tiempo_init << " segundos" << std::endl;
}

int solve_poisson(Ejemplo ejemplo, std::vector<std::vector<double>>& V, 
                 const std::vector<std::vector<double>>& source, int M, int N, double h, double k) {
    double delta = 1.0;
    int iterations = 0;
    
    double factor = (ejemplo == Ejemplo::ORIGINAL) ? (1.0 / e) : 1.0;
    
    while (delta > TOL) {
        delta = 0.0;
        iterations++;
        
        // Usar parallel for simple para la resolución principal
        #pragma omp parallel for reduction(max:delta)
        for (int i = 1; i < M; ++i) {
            for (int j = 1; j < N; ++j) {
                double V_new;
                
                if (ejemplo == Ejemplo::ORIGINAL) {
                    V_new = (
                        ((V[i + 1][j] + V[i - 1][j]) * k * k) +
                        ((V[i][j + 1] + V[i][j - 1]) * h * h) -
                        (factor * source[i][j] * h * h * k * k)) /
                        (2.0 * (h * h + k * k));
                } else if (ejemplo == Ejemplo::EJEMPLO_3) {
                    V_new = (
                        ((V[i + 1][j] + V[i - 1][j]) * k * k) +
                        ((V[i][j + 1] + V[i][j - 1]) * h * h) -
                        (source[i][j] * h * h * k * k)) /
                        (2.0 * (h * h + k * k));
                } else {
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
                    V_analytical[i][j] = 0.0;
                    break;
            }
        }
    }
}

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

void export_to_file(const std::vector<std::vector<double>>& V, double h, double k, int M, int N, 
                   const std::string& filename, const DominioConfig& config) {
    std::string full_path = "actividad3/" + filename;
    
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

void run_simulation(Ejemplo ejemplo, int M = 500, int N = 500) {
    DominioConfig config = getDominioConfig(ejemplo);
    
    double h = (config.x_max - config.x_min) / M;
    double k = (config.y_max - config.y_min) / N;
    
    std::vector<std::vector<double>> V, source, V_analytical;
    
    std::cout << "\n=== EJEMPLO " << static_cast<int>(ejemplo) << ": " << config.descripcion << " ===" << std::endl;
    std::cout << "Dominio: x ∈ [" << config.x_min << ", " << config.x_max << "], y ∈ [" 
              << config.y_min << ", " << config.y_max << "]" << std::endl;
    std::cout << "Grilla: " << (M+1) << " x " << (N+1) << " puntos" << std::endl;
    std::cout << "Tolerancia: " << TOL << std::endl;
    if (config.solucion_analitica != "No disponible") {
        std::cout << "Solución analítica: " << config.solucion_analitica << std::endl;
    }
    std::cout << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // ACTIVIDAD 3: Inicialización paralela con sections
    std::cout << "=== INICIALIZACIÓN PARALELA CON SECTIONS ===" << std::endl;
    initialize_grid_and_source_parallel(ejemplo, M, N, V, source, h, k, config);
    
    std::cout << "\nResolviendo la ecuación..." << std::endl;
    int iterations = solve_poisson(ejemplo, V, source, M, N, h, k);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double tiempo_segundos = duration.count() / 1000.0;
    
    std::cout << "\n=== RESULTADOS ACTIVIDAD 3: sections ===" << std::endl;
    std::cout << "Tiempo de ejecución total: " << tiempo_segundos << " segundos" << std::endl;
    std::cout << "Número de iteraciones: " << iterations << std::endl;
    std::cout << "Número de threads: " << omp_get_max_threads() << std::endl;
    
    if (config.solucion_analitica != "No disponible") {
        calculate_analytical_solution(ejemplo, M, N, V_analytical, h, k, config);
        std::cout << "\n=== ANÁLISIS DE ERROR ===" << std::endl;
        calculate_error(V, V_analytical, M, N);
        
        std::string analytical_filename = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + "_sections_analytical.dat";
        export_to_file(V_analytical, h, k, M, N, analytical_filename, config);
    }
    
    std::string numerical_filename = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + "_sections_numerical.dat";
    export_to_file(V, h, k, M, N, numerical_filename, config);
    
    std::cout << "\nSimulación completada." << std::endl;
}

int main() {
    std::cout << "=== ACTIVIDAD 3: SIMULADOR DE ECUACIÓN DE POISSON 2D (sections) ===" << std::endl;
    std::cout << "Esta implementación utiliza #pragma omp parallel sections" << std::endl;
    std::cout << "para ejecutar funciones independientes en paralelo.\n" << std::endl;
    
    create_data_directory();
    
    std::cout << "Ejemplos disponibles:" << std::endl;
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
        run_simulation(static_cast<Ejemplo>(opcion));
    } else if (opcion == 5) {
        std::cout << "\n=== EJECUTANDO TODOS LOS EJEMPLOS CON SECTIONS ===\n" << std::endl;
        
        std::cout << "REFLEXIÓN SOBRE SECTIONS:" << std::endl;
        std::cout << "• #pragma omp sections permite ejecutar diferentes tareas en paralelo" << std::endl;
        std::cout << "• Útil cuando tenemos funciones independientes que pueden ejecutarse simultáneamente" << std::endl;
        std::cout << "• En este caso, paraleliza initialize_boundary_conditions y calculate_source_term" << std::endl;
        std::cout << "• Para funciones muy rápidas, el overhead puede superar el beneficio" << std::endl;
        std::cout << "• Cada section se ejecuta en un thread diferente\n" << std::endl;
        
        for (int i = 0; i <= 4; ++i) {
            run_simulation(static_cast<Ejemplo>(i));
            std::cout << "\n" << std::string(60, '=') << std::endl;
        }
        
        std::cout << "\n=== CONCLUSIONES ACTIVIDAD 3 ===" << std::endl;
        std::cout << "¿Vale la pena paralelizar funciones tan rápidas?" << std::endl;
        std::cout << "- Para funciones de inicialización rápidas, el overhead de OpenMP" << std::endl;
        std::cout << "  puede ser mayor que el beneficio obtenido" << std::endl;
        std::cout << "- Es más útil cuando las secciones tienen trabajo computacional significativo" << std::endl;
        std::cout << "- El reordenamiento de funciones fue necesario para evitar dependencias" << std::endl;

    } else {
        std::cout << "Opción inválida. Ejecutando ejemplo 1 por defecto." << std::endl;
        run_simulation(Ejemplo::EJEMPLO_1);
    }
    
    return 0;
}
