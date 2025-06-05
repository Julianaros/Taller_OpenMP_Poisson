# Proyecto de Paralelización con OpenMP: Ecuación de Poisson 2D

Este repositorio contiene un proyecto desarrollado en la asignatura de **Sistemas Distribuidos**, cuyo objetivo es resolver numéricamente distintos casos de la ecuación de Poisson en dos dimensiones usando el método de diferencias finitas, y paralelizar su solución con **OpenMP**.

---

## Información General

- **Docente:** Carlos Andrés Gómez Vasco  
- **Estudiantes:** Julián Aros, Andrés Gómez, Laura Oliveros  

---

##  Objetivo del Proyecto

Aplicar distintas directivas de OpenMP sobre un código secuencial que resuelve la ecuación de Poisson 2D mediante diferencias finitas, analizando el impacto de cada técnica sobre:

- Tiempo de ejecución  
- Convergencia  
- Escalabilidad y rendimiento  

---

##  Ejemplos Analizados

A continuación se listan los 5 ejemplos tratados, incluyendo su dominio, condiciones de frontera, término fuente y solución analítica (si está disponible):

### Ejemplo 0: Fuente Gaussiana

**Ecuación diferencial:**  
∇²V = f(x, y) = (1 / 2πσ²) * exp(−((x−μ)² + (y−μ)²) / 2σ²)  
con μ = 0.5, σ = 0.1

**Dominio:** [0, 1] × [0, 1]  
**Condiciones de frontera:**  
- V(0, y) = V(1, y) = V(x, 0) = 0  
- V(x, 1) = x

**Interpretación:**  
Campo eléctrico generado por una carga puntual con distribución gaussiana.

**Solución analítica:** No disponible.

---

### Ejemplo 1: Fuente Exponencial

**Ecuación diferencial:**  
∇²V = (x² + y²)·e^(xy)

**Dominio:** [0, 2] × [0, 1]  
**Condiciones de frontera:**
- V(0, y) = 1  
- V(2, y) = e^(2y)  
- V(x, 0) = 1  
- V(x, 1) = e^x

**Solución exacta:**  
V(x, y) = e^(xy)

**Interpretación:**  
Caso de verificación con solución conocida.

---

### Ejemplo 2: Ecuación de Laplace

**Ecuación diferencial:**  
∇²V = 0

**Dominio:** [1, 2] × [0, 1]  
**Condiciones de frontera:**
- V(1, y) = ln(y² + 1)  
- V(2, y) = ln(y² + 4)  
- V(x, 0) = 2ln(x)  
- V(x, 1) = ln(x² + 1)

**Solución exacta:**  
V(x, y) = ln(x² + y²)

**Interpretación:**  
Problema clásico de equilibrio, útil en electrostática y flujo potencial.

---

### Ejemplo 3: Fuente Constante

**Ecuación diferencial:**  
∇²V = 4

**Dominio:** [1, 2] × [0, 2]  
**Condiciones de frontera:**
- V(1, y) = (1−y)²  
- V(2, y) = (2−y)²  
- V(x, 0) = x²  
- V(x, 2) = (x−2)²

**Solución exacta:**  
V(x, y) = (x − y)²

**Interpretación:**  
Generación uniforme de calor, carga constante o presión constante sobre una membrana.

---

### Ejemplo 4: Fuente Racional

**Ecuación diferencial:**  
∇²V = x/y + y/x

**Dominio:** [1, 2] × [1, 2]  
**Condiciones de frontera:**
- V(1, y) = y·ln(y)  
- V(2, y) = 2y·ln(2y)  
- V(x, 1) = x·ln(x)  
- V(x, 2) = x·ln(4x²)

**Solución exacta:**  
V(x, y) = xy·ln(xy)

**Interpretación:**  
Dependencias racionales en el término fuente, comportamiento logarítmico.

---

## Métodos Numéricos

- **Discretización:** Diferencias finitas centradas.
- **Iteración:** Jacobi / Gauss-Seidel.
- **Criterio de convergencia:**  
  \[
  \max_{i,j} |V^{(n+1)}_{i,j} - V^{(n)}_{i,j}| < 10^{-6}
  \]

---

## Paralelización con OpenMP

Se probaron diferentes estrategias:

| Actividad | Estrategia                             |
|----------|----------------------------------------|
| 1        | `#pragma omp parallel for`             |
| 2        | `collapse(2)`                          |
| 3        | `sections` para inicialización         |
| 4        | `parallel` + `for` + `schedule`        |
| 5        | `nowait`, `barrier`, `single`, `critical` |
| 6        | Contador de iteraciones con `atomic` y `critical` |
| 7        | División del dominio con `task` y `taskwait` |

**Observaciones clave:**
- La directiva `parallel for` con `reduction(max:delta)` fue la más eficiente y sencilla.
- El uso de `schedule(dynamic)` aumentó el overhead en la mayoría de los casos.
- `critical` es útil para contadores compartidos, pero puede ralentizar si se usa mal.

---

## Comparación de Rendimiento

Se midió el rendimiento en todos los ejemplos. En general, se observó:

- Reducción drástica del tiempo respecto al código secuencial.
- Cada estrategia tiene ventajas en diferentes contextos.
- Para cargas balanceadas, `schedule(static)` fue más eficiente.

---
