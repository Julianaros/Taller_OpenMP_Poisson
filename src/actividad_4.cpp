// Actividad 4: parallel y for separados
// Controlar la región paralela con diferentes estrategias de scheduling

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <string>
#include <filesystem>
#include <omp.h>

const double TOL = 1e-6;
const double e = 8.85e-12;

enum class Ejemplo {
    ORIGINAL = 0,
    EJEMPLO_1 = 1,
    EJEMPLO_2 = 2,
    EJEMPLO_3 = 3,
    EJEMPLO_4 = 4
};

enum class ScheduleStrategy {
    DEFAULT = 0,
    STATIC = 1,
    DYNAMIC = 2
};

struct DominioConfig {
    double x_min, x_max;
    double y_min, y_max;
    std::string descripcion;
    std::string solucion_analitica;
};

void create_data_directory() {
    try {
        if (!std::filesystem::exists("actividad4")) {
            std::filesystem::create_directory("actividad4");
            std::cout << "Directorio 'actividad4' creado exitosamente." << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error al crear el directorio 'actividad4': " << ex.what() << std::endl;
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

// ACTIVIDAD 4: Función que resuelve con parallel y for separados
int solve_poisson_parallel_explicit(Ejemplo ejemplo, std::vector<std::vector<double>>& V, 
                                   const std::vector<std::vector<double>>& source, 
                                   int M, int N, double h, double k, ScheduleStrategy strategy) {
    double delta = 1.0;
    int iterations = 0;
    
    double factor = (ejemplo == Ejemplo::ORIGINAL) ? (1.0 / e) : 1.0;
    
    std::string strategy_name;
    switch (strategy) {
        case ScheduleStrategy::DEFAULT:
            strategy_name = "DEFAULT";
            break;
        case ScheduleStrategy::STATIC:
            strategy_name = "STATIC";
            break;
        case ScheduleStrategy::DYNAMIC:
            strategy_name = "DYNAMIC";
            break;
    }
    
    std::cout << "Ejecutando con estrategia: " << strategy_name << std::endl;
    
    while (delta > TOL) {
        delta = 0.0;
        iterations++;
        
        // ACTIVIDAD 4: Usar parallel y for separados con control explícito
        #pragma omp parallel
        {
            double local_delta = 0.0;
            
            // Aplicar diferentes estrategias de scheduling
            switch (strategy) {
                case ScheduleStrategy::DEFAULT:
                    #pragma omp for reduction(max:delta)
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
                    break;
                    
                case ScheduleStrategy::STATIC:
                    #pragma omp for schedule(static) reduction(max:delta)
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
                    break;
                    
                case ScheduleStrategy::DYNAMIC:
                    #pragma omp for schedule(dynamic) reduction(max:delta)
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
                    break;
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
    std::string full_path = "actividad4/" + filename;
    
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

void run_simulation_with_strategy(Ejemplo ejemplo, ScheduleStrategy strategy, int M = 50, int N = 50) {
    DominioConfig config = getDominioConfig(ejemplo);
    
    double h = (config.x_max - config.x_min) / M;
    double k = (config.y_max - config.y_min) / N;
    
    std::vector<std::vector<double>> V, source, V_analytical;
    
    std::string strategy_name;
    std::string strategy_suffix;
    switch (strategy) {
        case ScheduleStrategy::DEFAULT:
            strategy_name = "DEFAULT";
            strategy_suffix = "default";
            break;
        case ScheduleStrategy::STATIC:
            strategy_name = "STATIC";
            strategy_suffix = "static";
            break;
        case ScheduleStrategy::DYNAMIC:
            strategy_name = "DYNAMIC";
            strategy_suffix = "dynamic";
            break;
    }
    
    std::cout << "\n=== EJEMPLO " << static_cast<int>(ejemplo) << " - ESTRATEGIA " << strategy_name << " ===" << std::endl;
    std::cout << config.descripcion << std::endl;
    std::cout << "Dominio: x ∈ [" << config.x_min << ", " << config.x_max << "], y ∈ [" 
              << config.y_min << ", " << config.y_max << "]" << std::endl;
    std::cout << "Grilla: " << (M+1) << " x " << (N+1) << " puntos" << std::endl;
    std::cout << "Tolerancia: " << TOL << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Inicialización
    initialize_boundary_conditions(ejemplo, M, N, V, h, k, config);
    calculate_source_term(ejemplo, M, N, source, h, k, config);
    
    std::cout << "\nResolviendo la ecuación con control explícito..." << std::endl;
    int iterations = solve_poisson_parallel_explicit(ejemplo, V, source, M, N, h, k, strategy);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double tiempo_segundos = duration.count() / 1000.0;
    
    std::cout << "\n=== RESULTADOS ACTIVIDAD 4: " << strategy_name << " ===" << std::endl;
    std::cout << "Tiempo de ejecución: " << tiempo_segundos << " segundos" << std::endl;
    std::cout << "Número de iteraciones: " << iterations << std::endl;
    std::cout << "Número de threads: " << omp_get_max_threads() << std::endl;
    
    if (config.solucion_analitica != "No disponible") {
        calculate_analytical_solution(ejemplo, M, N, V_analytical, h, k, config);
        std::cout << "\n=== ANÁLISIS DE ERROR ===" << std::endl;
        calculate_error(V, V_analytical, M, N);
        
        std::string analytical_filename = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + 
                                        "_" + strategy_suffix + "_analytical.dat";
        export_to_file(V_analytical, h, k, M, N, analytical_filename, config);
    }
    
    std::string numerical_filename = "ejemplo_" + std::to_string(static_cast<int>(ejemplo)) + 
                                   "_" + strategy_suffix + "_numerical.dat";
    export_to_file(V, h, k, M, N, numerical_filename, config);
    
    std::cout << "Simulación completada." << std::endl;
}

void run_all_strategies_for_example(Ejemplo ejemplo) {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "EJECUTANDO TODAS LAS ESTRATEGIAS PARA EJEMPLO " << static_cast<int>(ejemplo) << std::endl;
    
    // Ejecutar con todas las estrategias
    run_simulation_with_strategy(ejemplo, ScheduleStrategy::DEFAULT);
    std::cout << "\n" << std::string(40, '-') << std::endl;
    
    run_simulation_with_strategy(ejemplo, ScheduleStrategy::STATIC);
    std::cout << "\n" << std::string(40, '-') << std::endl;
    
    run_simulation_with_strategy(ejemplo, ScheduleStrategy::DYNAMIC);
}

int main() {
    std::cout << "=== ACTIVIDAD 4: SIMULADOR DE ECUACIÓN DE POISSON 2D ===" << std::endl;
    std::cout << "Control explícito de región paralela con diferentes estrategias de scheduling" << std::endl;
    std::cout << "\nEsta implementación utiliza:" << std::endl;
    std::cout << "• #pragma omp parallel" << std::endl;
    std::cout << "• #pragma omp for reduction(max:delta)" << std::endl;
    std::cout << "• schedule(static) y schedule(dynamic)\n" << std::endl;
    
    create_data_directory();
    
    std::cout << "Ejemplos disponibles:" << std::endl;
    std::cout << "0 - Ejemplo Original: Término fuente gaussiano" << std::endl;
    std::cout << "1 - Ejemplo 1: ∇²V = (x² + y²)e^xy" << std::endl;
    std::cout << "2 - Ejemplo 2: ∇²V = 0 (Ecuación de Laplace)" << std::endl;
    std::cout << "3 - Ejemplo 3: ∇²V = 4" << std::endl;
    std::cout << "4 - Ejemplo 4: ∇²V = x/y + y/x" << std::endl;
    std::cout << "5 - Ejecutar todos los ejemplos con todas las estrategias" << std::endl;
    
    int opcion;
    std::cout << "\nSeleccione una opción (0-5): ";
    std::cin >> opcion;
    
    if (opcion >= 0 && opcion <= 4) {
        run_all_strategies_for_example(static_cast<Ejemplo>(opcion));
    } else if (opcion == 5) {
        std::cout << "\n=== EJECUTANDO TODOS LOS EJEMPLOS CON TODAS LAS ESTRATEGIAS ===\n" << std::endl;
        
        std::cout << "INFORMACIÓN SOBRE ESTRATEGIAS DE SCHEDULING:" << std::endl;
        std::cout << "• DEFAULT: OpenMP decide la distribución de trabajo" << std::endl;
        std::cout << "• STATIC: Distribuye iteraciones en bloques fijos entre threads" << std::endl;
        std::cout << "• DYNAMIC: Distribuye iteraciones dinámicamente durante la ejecución" << std::endl;
        std::cout << "• STATIC es mejor para cargas de trabajo uniformes" << std::endl;
        std::cout << "• DYNAMIC es mejor para cargas de trabajo irregulares" << std::endl;
        
        for (int i = 0; i <= 4; ++i) {
            run_all_strategies_for_example(static_cast<Ejemplo>(i));
        }
        
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "=== CONCLUSIONES ACTIVIDAD 4 ===" << std::endl;
        std::cout << "\nVentajas del control explícito con #pragma omp parallel:" << std::endl;
        std::cout << "• Mayor control sobre la región paralela" << std::endl;
        std::cout << "• Posibilidad de combinar múltiples directivas for dentro de la misma región" << std::endl;
        std::cout << "• Mejor control de variables privadas y compartidas" << std::endl;
        
        std::cout << "\nComparación de estrategias de scheduling:" << std::endl;
        std::cout << "• STATIC: Menor overhead, mejor para trabajo uniforme" << std::endl;
        std::cout << "• DYNAMIC: Mayor overhead, mejor balance de carga para trabajo irregular" << std::endl;
        std::cout << "• DEFAULT: OpenMP elige automáticamente (usualmente static)" << std::endl;
        
        std::cout << "\nPara este problema (ecuación de Poisson):" << std::endl;
        std::cout << "• El trabajo es bastante uniforme por iteración" << std::endl;
        std::cout << "• STATIC debería ser más eficiente que DYNAMIC" << std::endl;
        std::cout << "• La diferencia puede ser pequeña para grillas regulares" << std::endl;

    } else {
        std::cout << "Opción inválida. Ejecutando ejemplo 1 con todas las estrategias." << std::endl;
        run_all_strategies_for_example(Ejemplo::EJEMPLO_1);
    }
    
    return 0;
}
