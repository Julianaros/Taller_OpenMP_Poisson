# Proyecto de Paralelización con OpenMP: Ecuación de Poisson 2D

Este repositorio contiene un proyecto desarrollado en la asignatura de **Sistemas Distribuidos**, cuyo objetivo es resolver numéricamente distintos casos de la ecuación de Poisson en dos dimensiones usando el método de diferencias finitas, y paralelizar su solución con **OpenMP**.

---

## Información General

- **Docente:** Carlos Andrés Gómez Vasco  
- **Estudiantes:** Julián Aros, Andrés Gómez, Laura Oliveros  

---

## Objetivo del Proyecto

Aplicar distintas directivas de OpenMP sobre un código secuencial que resuelve la ecuación de Poisson 2D mediante diferencias finitas, analizando el impacto de cada técnica sobre:

- Tiempo de ejecución  
- Convergencia  
- Escalabilidad y rendimiento  

---

## Ejemplos Analizados

### Ejemplo 0: Fuente Gaussiana

**Ecuación diferencial:**

$$
\nabla^2 V = f(x, y) = \frac{1}{2\pi\sigma^2} \exp\left( -\frac{(x - \mu)^2 + (y - \mu)^2}{2\sigma^2} \right)
$$

donde $\mu = 0.5$, $\sigma = 0.1$

**Dominio:** $[0, 1] \times [0, 1]$  
**Condiciones de frontera:**  
- $V(0, y) = V(1, y) = V(x, 0) = 0$  
- $V(x, 1) = x$

**Solución analítica:** No disponible.

---

### Ejemplo 1: Fuente Exponencial

**Ecuación diferencial:**

$$
\nabla^2 V = (x^2 + y^2)e^{xy}
$$

**Solución exacta:**

$$
V(x, y) = e^{xy}
$$

**Dominio:** $[0, 2] \times [0, 1]$  
**Condiciones de frontera:**
- $V(0, y) = 1$  
- $V(2, y) = e^{2y}$  
- $V(x, 0) = 1$  
- $V(x, 1) = e^x$

---

### Ejemplo 2: Ecuación de Laplace

**Ecuación diferencial:**

$$
\nabla^2 V = 0
$$

**Solución exacta:**

$$
V(x, y) = \ln(x^2 + y^2)
$$

**Dominio:** $[1, 2] \times [0, 1]$  
**Condiciones de frontera:**
- $V(1, y) = \ln(y^2 + 1)$  
- $V(2, y) = \ln(y^2 + 4)$  
- $V(x, 0) = 2\ln(x)$  
- $V(x, 1) = \ln(x^2 + 1)$

---

### Ejemplo 3: Fuente Constante

**Ecuación diferencial:**

$$
\nabla^2 V = 4
$$

**Solución exacta:**

$$
V(x, y) = (x - y)^2
$$

**Dominio:** $[1, 2] \times [0, 2]$  
**Condiciones de frontera:**
- $V(1, y) = (1 - y)^2$  
- $V(2, y) = (2 - y)^2$  
- $V(x, 0) = x^2$  
- $V(x, 2) = (x - 2)^2$

---

### Ejemplo 4: Fuente Racional

**Ecuación diferencial:**

$$
\nabla^2 V = \frac{x}{y} + \frac{y}{x}
$$

**Solución exacta:**

$$
V(x, y) = xy \ln(xy)
$$

**Dominio:** $[1, 2] \times [1, 2]$  
**Condiciones de frontera:**
- $V(1, y) = y\ln(y)$  
- $V(2, y) = 2y\ln(2y)$  
- $V(x, 1) = x\ln(x)$  
- $V(x, 2) = x\ln(4x^2)$

---

## Métodos Numéricos

- **Discretización:** Diferencias finitas centradas.

\[
\nabla^2 V \approx \frac{V_{i+1,j} - 2V_{i,j} + V_{i-1,j}}{h^2} + \frac{V_{i,j+1} - 2V_{i,j} + V_{i,j-1}}{k^2}
\]

- **Iteración:** Método de Jacobi / Gauss-Seidel.

\[
V^{(n+1)}_{i,j} = \frac{
(V^{(n)}_{i+1,j} + V^{(n)}_{i-1,j})k^2 + (V^{(n)}_{i,j+1} + V^{(n)}_{i,j-1})h^2 \pm f_{i,j} h^2 k^2
}{
2(h^2 + k^2)
}
\]

- **Criterio de convergencia:**

\[
\max_{i,j} |V^{(n+1)}_{i,j} - V^{(n)}_{i,j}| < 10^{-6}
\]

---

## ⚙️ Paralelización con OpenMP

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

- `reduction(max:delta)` es fundamental para evitar condiciones de carrera.
- `schedule(static)` ofrece mejor rendimiento en dominios regulares.
- `task` es más flexible pero introduce mayor overhead.

---

## 📊 Comparación de Rendimiento

Se observaron mejoras notables en el rendimiento, siendo el uso de `parallel for` con reducción la mejor opción en la mayoría de los casos.

---

