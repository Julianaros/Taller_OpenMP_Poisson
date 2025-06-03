// Actividad 7: Ecuación de Poisson 2D con Tasks de OpenMP
// Divide el dominio en bloques y asigna tareas a distintos hilos

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

// Función para crear el directorio necesario
void create_directories() {
    try {
        if (!std::filesystem::exists("actividad7")) {
            std::filesystem::create_directory("actividad7");
            std::cout << "Directorio 'actividad7' creado exitosamente." << std::endl;
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

// Función que resuelve usando tasks
int solve_poisson_with_tasks(Ejemplo ejemplo, std::vector<std::vector<double>>& V, 
                           const std::vector<std::vector<double>>& source, int M, int N, 
                           double h, double k, int block_size = 50) {
    double delta = 1.0;
    int iterations = 0;
    
    double factor = (ejemplo == Ejemplo::ORIGINAL) ? (1.0 / e) : 1.0;
    
    while (delta > TOL) {
        delta = 0.0;
        iterations++;
        
        // Array para almacenar los deltas locales de cada tarea
        std::vector<double> local_deltas;
        
        #pragma omp parallel
        {
            #pragma omp single
            {
                // Crear tareas para cada bloque
                for (int bi = 1; bi < M; bi += block_size) {
                    for (int bj = 1; bj < N; bj += block_size) {
                        #pragma omp task
                        {
                            double local_delta = 0.0;
                            int end_i = std::min(bi + block_size, M);
                            int end_j = std::min(bj + block_size, N);
                            
                            for (int i = bi; i < end_i; ++i) {
                                for (int j = bj; j < end_j; ++j) {
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
                                    
                                    local_delta = std::max(local_delta, std::abs(V_new - V[i][j]));
                                    V[i][j] = V_new;
                                }
                            }
                            
                            // Actualizar delta global de forma crítica
                            #pragma omp critical
                            {
                                delta = std::max(delta, local_delta);
                            }
                        }
                    }
                }
                
                // Esperar a que todas las tareas terminen
                #pragma omp taskwait
            }
        }
    }
    
    return iterations;
}

void export_to_file(const std::vector<std::vector<double>>& V, double h, double k, int M, int N, 
                   const std::string& filename, const DominioConfig& config) {
    std::string full_path = "actividad7/" + filename;
    
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

void run_simulation_with_tasks(Ejemplo ejemplo, int M, int N, int block_size) {
    DominioConfig config = getDominioConfig(ejemplo);
    
    double h = (config.x_max - config.x_min) / M;
    double k = (config.y_max - config.y_min) / N;
    
    std::vector<std::vector<double>> V, source;
    
    std::cout << "\n=== ACTIVIDAD 7 - EJEMPLO " << static_cast<int>(ejemplo) << " CON TASKS ===" << std::endl;
    std::cout << "Descripción: " << config.descripcion << std::endl;
    std::cout << "Dominio: x ∈ [" << config.x_min << ", " << config.x_max << "], y ∈ [" 
              << config.y_min << ", " << config.y_max << "]" << std::endl;
    std::cout << "Grilla: " << (M+1) << " x " << (N+1) << " puntos" << std::endl;
    std::cout << "Tamaño de bloque: " << block_size << " x " << block_size << std::endl;
    std::cout << "Número de threads: " << omp_get_max_threads() << std::endl;
    std::cout << "Tolerancia: " << TOL << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    initialize_boundary_conditions(ejemplo, M, N, V, h, k, config);
    calculate_source_term(ejemplo, M, N, source, h, k, config);
    
    std::cout << "\nResolviendo con tasks..." << std::endl;
    int iterations = solve_poisson_with_tasks(ejemplo, V, source, M, N, h, k, block_size);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double tiempo_segundos = duration.count() / 1000.0;
    
    std::cout << "\n=== RESULTADOS CON TASKS ===" << std::endl;
    std::cout << "Tiempo de ejecución: " << tiempo_segundos << " segundos" << std::endl;
    std::cout << "Número de iteraciones: " << iterations << std::endl;
    std::cout << "Número de threads utilizados: " << omp_get_max_threads() << std::endl;
    
    // Exportar resultados
    std::string filename = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + "_tasks.dat";
    export_to_file(V, h, k, M, N, filename, config);
}

int main() {
    std::cout << "=== ACTIVIDAD 7: ECUACIÓN DE POISSON 2D CON TASKS ===\n" << std::endl;
    
    create_directories();
    
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
    
    // Configurar el número de threads
    int num_threads;
    std::cout << "Ingrese el número de threads a usar: ";
    std::cin >> num_threads;
    omp_set_num_threads(num_threads);
    
    // Configurar tamaño de grilla
    int M, N;
    std::cout << "Ingrese el número de divisiones en X (M): ";
    std::cin >> M;
    std::cout << "Ingrese el número de divisiones en Y (N): ";
    std::cin >> N;
    
    // Configurar tamaño de bloque
    int block_size;
    std::cout << "Ingrese el tamaño de bloque (recomendado: 25-100): ";
    std::cin >> block_size;
    
    if (opcion >= 0 && opcion <= 4) {
        run_simulation_with_tasks(static_cast<Ejemplo>(opcion), M, N, block_size);
    } else if (opcion == 5) {
        std::cout << "\n=== EJECUTANDO TODOS LOS EJEMPLOS CON TASKS ===\n" << std::endl;
        for (int i = 0; i <= 4; ++i) {
            run_simulation_with_tasks(static_cast<Ejemplo>(i), M, N, block_size);
            std::cout << "\n" << std::string(60, '=') << std::endl;
        }
    } else {
        std::cout << "Opción inválida. Ejecutando ejemplo 1 por defecto." << std::endl;
        run_simulation_with_tasks(Ejemplo::EJEMPLO_1, M, N, block_size);
    }
    
    std::cout << "\n=== ANÁLISIS DE RENDIMIENTO ===" << std::endl;
    std::cout << "Para comparar el rendimiento:" << std::endl;
    std::cout << "1. Compare el tiempo obtenido con el tiempo de parallel for" << std::endl;
    std::cout << "2. Observe cómo diferentes tamaños de bloque afectan el rendimiento" << std::endl;
    std::cout << "3. Note que tasks permite mayor flexibilidad en la distribución de trabajo" << std::endl;
    std::cout << "4. Tasks puede ser útil para balanceo de carga dinámico" << std::endl;
    
    return 0;
}
