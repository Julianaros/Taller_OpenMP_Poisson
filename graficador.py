import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def load_data(filename):
    """Carga los datos de un archivo .dat"""
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    
    # Determinar las dimensiones de la malla
    unique_x = np.unique(x)
    unique_y = np.unique(y)
    nx = len(unique_x)
    ny = len(unique_y)
    
    # Reorganizar los datos en una malla 2D
    X = x.reshape(nx, ny)
    Y = y.reshape(nx, ny)
    Z = z.reshape(nx, ny)
    
    return X, Y, Z

def plot_2d_comparison(X, Y, Z_num, Z_anal, example_num):
    """Crea gráficos 2D comparando soluciones numérica y analítica"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Gráfico de la solución numérica
    c1 = ax1.contourf(X, Y, Z_num, levels=20, cmap='viridis')
    ax1.set_title(f'Ejemplo {example_num} - Solución Numérica')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    fig.colorbar(c1, ax=ax1)
    
    # Gráfico de la solución analítica
    c2 = ax2.contourf(X, Y, Z_anal, levels=20, cmap='viridis')
    ax2.set_title(f'Ejemplo {example_num} - Solución Analítica')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    fig.colorbar(c2, ax=ax2)
    
    plt.tight_layout()
    plt.savefig(f'figuras/ejemplo_{example_num}_2d_comparison.png')
    plt.close()

def plot_3d_comparison(X, Y, Z_num, Z_anal, example_num):
    """Crea gráficos 3D comparando soluciones numérica y analítica"""
    fig = plt.figure(figsize=(14, 6))
    
    # Gráfico 3D de la solución numérica
    ax1 = fig.add_subplot(121, projection='3d')
    surf1 = ax1.plot_surface(X, Y, Z_num, cmap='viridis', alpha=0.8)
    ax1.set_title(f'Ejemplo {example_num} - Solución Numérica (3D)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('V(x,y)')
    fig.colorbar(surf1, ax=ax1, shrink=0.5)
    
    # Gráfico 3D de la solución analítica (si está disponible)
    if Z_anal is not None:
        ax2 = fig.add_subplot(122, projection='3d')
        surf2 = ax2.plot_surface(X, Y, Z_anal, cmap='viridis', alpha=0.8)
        ax2.set_title(f'Ejemplo {example_num} - Solución Analítica (3D)')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel('V(x,y)')
        fig.colorbar(surf2, ax=ax2, shrink=0.5)
    
    plt.tight_layout()
    plt.savefig(f'figuras/ejemplo_{example_num}_3d_comparison.png')
    plt.close()

def plot_3d_single(X, Y, Z, example_num, title):
    """Crea un gráfico 3D individual"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('V(x,y)')
    fig.colorbar(surf, ax=ax, shrink=0.5)
    plt.tight_layout()
    plt.savefig(f'figuras/ejemplo_{example_num}_3d_numerical.png')
    plt.close()

def plot_error(X, Y, error, example_num):
    """Crea gráficos del error entre soluciones numérica y analítica"""
    fig = plt.figure(figsize=(12, 6))
    
    # Gráfico de contorno del error
    ax1 = fig.add_subplot(121)
    c1 = ax1.contourf(X, Y, error, levels=20, cmap='hot')
    ax1.set_title(f'Ejemplo {example_num} - Error (Contorno)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    fig.colorbar(c1, ax=ax1, label='Error absoluto')
    
    # Gráfico 3D del error
    ax2 = fig.add_subplot(122, projection='3d')
    surf = ax2.plot_surface(X, Y, error, cmap='hot', alpha=0.8)
    ax2.set_title(f'Ejemplo {example_num} - Error (3D)')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('Error absoluto')
    fig.colorbar(surf, ax=ax2, shrink=0.5, label='Error absoluto')
    
    plt.tight_layout()
    plt.savefig(f'figuras/ejemplo_{example_num}_error.png')
    plt.close()

def process_example(example_num):
    """Procesa un ejemplo específico"""
    num_file = f'data/ejemplo_{example_num}_numerical.dat'
    anal_file = f'data/ejemplo_{example_num}_analytical.dat'
    
    if not os.path.exists(num_file):
        print(f"Archivo numérico para ejemplo {example_num} no encontrado")
        return
    
    # Cargar datos numéricos
    X_num, Y_num, Z_num = load_data(num_file)
    
    # Graficar siempre la solución numérica en 3D (especialmente para el ejemplo 0)
    plot_3d_single(X_num, Y_num, Z_num, example_num, f'Ejemplo {example_num} - Solución Numérica (3D)')
    
    # Verificar si existe solución analítica
    if os.path.exists(anal_file):
        # Cargar datos analíticos
        X_anal, Y_anal, Z_anal = load_data(anal_file)
        
        # Verificar que las mallas coincidan
        assert np.allclose(X_num, X_anal), "Las mallas X no coinciden"
        assert np.allclose(Y_num, Y_anal), "Las mallas Y no coinciden"
        
        # Calcular error absoluto
        error = np.abs(Z_num - Z_anal)
        
        # Crear gráficos comparativos
        plot_2d_comparison(X_num, Y_num, Z_num, Z_anal, example_num)
        plot_3d_comparison(X_num, Y_num, Z_num, Z_anal, example_num)
        plot_error(X_num, Y_num, error, example_num)
    else:
        print(f"No hay solución analítica para el ejemplo {example_num}")
        
        # Graficar solo solución numérica en 2D si no hay analítica
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        c = ax.contourf(X_num, Y_num, Z_num, levels=20, cmap='viridis')
        ax.set_title(f'Ejemplo {example_num} - Solución Numérica')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        fig.colorbar(c, ax=ax)
        plt.tight_layout()
        plt.savefig(f'figuras/ejemplo_{example_num}_numerical_only.png')
        plt.close()
    
    print(f"Gráficos para ejemplo {example_num} generados exitosamente")

def main():
    # Crear directorio de figuras si no existe
    if not os.path.exists('figuras'):
        os.makedirs('figuras')
    
    # Preguntar qué ejemplos procesar
    print("Ejemplos disponibles:")
    print("0 - Ejemplo Original")
    print("1 - Ejemplo 1")
    print("2 - Ejemplo 2")
    print("3 - Ejemplo 3")
    print("4 - Ejemplo 4")
    print("5 - Todos los ejemplos")
    
    try:
        option = int(input("Seleccione una opción (0-5): "))
    except ValueError:
        print("Opción inválida. Usando opción por defecto (1).")
        option = 1
    
    if option == 5:
        # Procesar todos los ejemplos
        for i in range(5):
            process_example(i)
    elif 0 <= option <= 4:
        # Procesar un ejemplo específico
        process_example(option)
    else:
        print("Opción inválida. Procesando ejemplo 1 por defecto.")
        process_example(1)

if __name__ == "__main__":
    main()