// Actividad 6: Contar el número de iteraciones con acceso compartido
// Simulador de la Ecuación de Poisson en 2D con OpenMP
// Comparación entre #pragma omp critical y #pragma omp atomic

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

// Variable global para contar iteraciones
int iterations_global = 0;

void create_directories() {
    try {
        if (!std::filesystem::exists("actividad6")) {
            std::filesystem::create_directory("actividad6");
            std::cout << "Directorio 'actividad6' creado exitosamente." << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error al crear el directorio: " << ex.what() << std::endl;
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

// Versión con #pragma omp critical
int solve_poisson_critical(Ejemplo ejemplo, std::vector<std::vector<double>>& V, 
                          const std::vector<std::vector<double>>& source, int M, int N, double h, double k) {
    double delta = 1.0;
    int iterations = 0;
    double factor = (ejemplo == Ejemplo::ORIGINAL) ? (1.0 / e) : 1.0;
    
    while (delta > TOL) {
        delta = 0.0;
        
        // Incrementar contador de iteraciones usando critical
        #pragma omp critical
        {
            iterations++;
        }
        
        #pragma omp parallel for collapse(2) reduction(max:delta)
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

// Versión con #pragma omp atomic
int solve_poisson_atomic(Ejemplo ejemplo, std::vector<std::vector<double>>& V, 
                        const std::vector<std::vector<double>>& source, int M, int N, double h, double k) {
    double delta = 1.0;
    iterations_global = 0; // Resetear variable global
    double factor = (ejemplo == Ejemplo::ORIGINAL) ? (1.0 / e) : 1.0;
    
    while (delta > TOL) {
        delta = 0.0;
        
        // Incrementar contador de iteraciones usando atomic
        #pragma omp atomic
        iterations_global++;
        
        #pragma omp parallel for collapse(2) reduction(max:delta)
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
    
    return iterations_global;
}

void export_to_file(const std::vector<std::vector<double>>& V, double h, double k, int M, int N, 
                   const std::string& filename, const DominioConfig& config) {
    std::string full_path = "actividad6/" + filename;
    
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

void run_comparison(Ejemplo ejemplo, int M = 50, int N = 50) {
    DominioConfig config = getDominioConfig(ejemplo);
    
    double h = (config.x_max - config.x_min) / M;
    double k = (config.y_max - config.y_min) / N;
    
    std::vector<std::vector<double>> V_critical, V_atomic, source;
    
    std::cout << "\n=== ACTIVIDAD 6: EJEMPLO " << static_cast<int>(ejemplo) << " ===" << std::endl;
    std::cout << "Descripción: " << config.descripcion << std::endl;
    std::cout << "Número de threads disponibles: " << omp_get_max_threads() << std::endl;
    std::cout << std::endl;
    
    // Inicialización común
    initialize_boundary_conditions(ejemplo, M, N, V_critical, h, k, config);
    calculate_source_term(ejemplo, M, N, source, h, k, config);
    
    // Copiar condiciones iniciales para la segunda prueba
    V_atomic = V_critical;
    
    // ===== PRUEBA CON CRITICAL =====
    std::cout << "=== PRUEBA CON #pragma omp critical ===" << std::endl;
    auto start_critical = std::chrono::high_resolution_clock::now();
    
    int iterations_critical = solve_poisson_critical(ejemplo, V_critical, source, M, N, h, k);
    
    auto end_critical = std::chrono::high_resolution_clock::now();
    auto duration_critical = std::chrono::duration_cast<std::chrono::milliseconds>(end_critical - start_critical);
    double tiempo_critical = duration_critical.count() / 1000.0;
    
    std::cout << "Tiempo de ejecución: " << tiempo_critical << " segundos" << std::endl;
    std::cout << "Número de iteraciones: " << iterations_critical << std::endl;
    
    // ===== PRUEBA CON ATOMIC =====
    std::cout << "\n=== PRUEBA CON #pragma omp atomic ===" << std::endl;
    auto start_atomic = std::chrono::high_resolution_clock::now();
    
    int iterations_atomic = solve_poisson_atomic(ejemplo, V_atomic, source, M, N, h, k);
    
    auto end_atomic = std::chrono::high_resolution_clock::now();
    auto duration_atomic = std::chrono::duration_cast<std::chrono::milliseconds>(end_atomic - start_atomic);
    double tiempo_atomic = duration_atomic.count() / 1000.0;
    
    std::cout << "Tiempo de ejecución: " << tiempo_atomic << " segundos" << std::endl;
    std::cout << "Número de iteraciones: " << iterations_atomic << std::endl;
    
    // ===== COMPARACIÓN =====
    std::cout << "\n=== COMPARACIÓN DE RENDIMIENTO ===" << std::endl;
    std::cout << "Diferencia de tiempo: " << std::abs(tiempo_critical - tiempo_atomic) << " segundos" << std::endl;
    
    if (tiempo_critical < tiempo_atomic) {
        double mejora = ((tiempo_atomic - tiempo_critical) / tiempo_atomic) * 100.0;
        std::cout << "CRITICAL es " << mejora << "% más rápido que ATOMIC" << std::endl;
    } else if (tiempo_atomic < tiempo_critical) {
        double mejora = ((tiempo_critical - tiempo_atomic) / tiempo_critical) * 100.0;
        std::cout << "ATOMIC es " << mejora << "% más rápido que CRITICAL" << std::endl;
    } else {
        std::cout << "Ambos métodos tienen el mismo rendimiento" << std::endl;
    }
    
    // Verificar que ambos métodos llegaron al mismo resultado
    std::cout << "Verificación: Ambos métodos convergen en " << iterations_critical << " iteraciones: " 
              << (iterations_critical == iterations_atomic ? "SÍ" : "NO") << std::endl;
    
    // Exportar resultados
    std::string filename_critical = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + "_critical.dat";
    std::string filename_atomic = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + "_atomic.dat";
    
    export_to_file(V_critical, h, k, M, N, filename_critical, config);
    export_to_file(V_atomic, h, k, M, N, filename_atomic, config);
    
    // Exportar resumen de comparación
    std::string summary_filename = "actividad6/comparacion_ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + ".txt";
    std::ofstream summary_file(summary_filename);
    if (summary_file.is_open()) {
        summary_file << "=== ACTIVIDAD 6: COMPARACIÓN CRITICAL vs ATOMIC ===" << std::endl;
        summary_file << "Ejemplo: " << static_cast<int>(ejemplo) << " - " << config.descripcion << std::endl;
        summary_file << "Threads utilizados: " << omp_get_max_threads() << std::endl;
        summary_file << std::endl;
        summary_file << "RESULTADOS CON #pragma omp critical:" << std::endl;
        summary_file << "  Tiempo: " << tiempo_critical << " segundos" << std::endl;
        summary_file << "  Iteraciones: " << iterations_critical << std::endl;
        summary_file << std::endl;
        summary_file << "RESULTADOS CON #pragma omp atomic:" << std::endl;
        summary_file << "  Tiempo: " << tiempo_atomic << " segundos" << std::endl;
        summary_file << "  Iteraciones: " << iterations_atomic << std::endl;
        summary_file << std::endl;
        summary_file << "CONCLUSIONES:" << std::endl;
        if (tiempo_critical < tiempo_atomic) {
            double mejora = ((tiempo_atomic - tiempo_critical) / tiempo_atomic) * 100.0;
            summary_file << "  - CRITICAL es " << mejora << "% más rápido que ATOMIC" << std::endl;
        } else if (tiempo_atomic < tiempo_critical) {
            double mejora = ((tiempo_critical - tiempo_atomic) / tiempo_critical) * 100.0;
            summary_file << "  - ATOMIC es " << mejora << "% más rápido que CRITICAL" << std::endl;
        } else {
            summary_file << "  - Ambos métodos tienen rendimiento similar" << std::endl;
        }
        summary_file.close();
        std::cout << "Resumen exportado a " << summary_filename << std::endl;
    }
}

int main() {
    std::cout << "=== ACTIVIDAD 6: CONTADOR DE ITERACIONES CON CRITICAL Y ATOMIC ===" << std::endl;
    
    create_directories();
    
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
        run_comparison(static_cast<Ejemplo>(opcion));
    } else if (opcion == 5) {
        std::cout << "\n=== EJECUTANDO COMPARACIÓN PARA TODOS LOS EJEMPLOS ===\n" << std::endl;
        for (int i = 0; i <= 4; ++i) {
            run_comparison(static_cast<Ejemplo>(i));
            std::cout << "\n" << std::string(60, '=') << std::endl;
        }
    } else {
        std::cout << "Opción inválida. Ejecutando ejemplo 1 por defecto." << std::endl;
        run_comparison(Ejemplo::EJEMPLO_1);
    }
    
    return 0;
}
