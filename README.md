## Código proyecto semestral IIQ2003 - G8

### Modelamiento de un sistema TES sometido a radiación solar y por reflejo de espejos

### 1. Introducción y planteamiento del modelo
En el siguiente código se modela el gradiente de temperaturas dentro de un cilindro vertical de largo L y radio R, que en su interior transporta sal fundida de densidad $\rho$ con un perfil de velocidad axial $v_x$ constante mediante una bomba ($\frac{\delta v_z}{\delta z}=0$). El líquido de trabajo además entra al cilindro a una temperatura uniforme $T_0$. En el contexto de TES aplicada a energía solar, existe una fuente de radiación $\dot{q_{rad} }$ proveniente de una torre de concentración solar. Para este modelo se toman en consideración los siguientes supuestos:
* Estado estacionario -> $\frac{\delta}{\delta t} = 0$
* Simetría axial -> $\frac{\delta}{\delta \theta}=0$
* Fluido incompresible -> $∇ \cdot \vec{v} =0$
* Fluido newtoniano
* Disipación viscosa despreciable
* Flujo unidireccional
* Propiedades fisicoquímicas constantes
* Grosor de la pared ínfimo

Entonces, de las ecuaciones de Navier-Stokes se plantea la siguiente EDO:

$$ \frac{\delta p}{\delta z} = \mu (\frac{1}{r}\frac{\delta v_x}{\delta r} + \frac{\delta^2 v_x}{\delta r^2}) +\rho g $$

Sujeto a la siguiente condición de borde:
* CB1: Condición de no resbalar (Dirichlet)
    $$v_x \bigg|_{r=R}=0$$

Lo que resulta en la siguiente ecuación para el perfil de velocidad (perfil parabólico de Poiseuille):
$$v_z(r) = \frac{\nabla P-\rho g}{4\mu}(r^2-R^2)$$

Con la definición del gradiente de presión:
$$\nabla P=\frac{\delta p}{\delta z}=\frac{p|_L-p|_0}{L}$$

Para simplificar las ecuaciones posteriormente, se define el perfil de velocidad como $v_z(r) = \Phi (r^2-R^2)$, con $\Phi = \frac{\nabla P-\rho g}{4\mu}$.

Por otro lado, se tiene el modelo de la temperatura dentro del cilindro, que se obtiene de las ecuaciones de conservación de la energía (en coordenadas cilíndricas):

$$\rho c_p v_z \frac{\delta T}{\delta z} = k\left[\frac{1}{r}\frac{\delta}{\delta r}\left(r \frac{\delta T}{\delta r}\right) + \frac{\delta^2 T}{\delta z^2}\right]$$

Lo que resulta en la siguiente Ecuación Diferencial Ordinaria (EDO) si expandimos el término radial:
$$\rho c_p v_z \frac{\delta T}{\delta z} = k\left[\frac{1}{r}\frac{\delta T}{\delta r} + \frac{\delta^2 T}{\delta r^2}+\frac{\delta^2 T}{\delta z^2}\right]$$

La EDO se somete a las siguientes condiciones de borde:

* CB2: Temperatura de entrada (Dirichlet)
    $$T|_{z=0,r} = T_0$$

* CB3: Cambio en la frontera (Neumann)
    $$\frac{\delta T}{\delta z}|_{z=0,r} = 0$$

* CB4: Simetría del perfil radialmente (Neumann)
    $$\frac{\delta T}{\delta r}|_{z,r=0} = 0$$

* CB5: Calor por radiación (Neumann)
    $$\dot{q}_{\text{rad}} = k\frac{\delta T}{\delta r}|_{z,r=R}$$


## 1.2 Discretización del modelo y condiciones de borde

Reemplazando el perfil de velocidades en el balance de energía:

$$\rho c_p \Phi (r^2-R^2) \frac{\delta T}{\delta z} = k\left[\frac{1}{r}\frac{\delta T}{\delta r} + \frac{\delta^2 T}{\delta r^2}+\frac{\delta^2 T}{\delta z^2}\right]$$

Y aplicando diferencias finitas centrales de segundo orden (donde $i$ es radial y $j$ es axial):

$$
\rho c_p \Phi ((i\Delta r)^2-R^2)\left(\frac{T_{i,j+1}-T_{i,j-1}}{2\Delta z}\right) = k\left[\frac{1}{i\Delta r}\left(\frac{T_{i+1,j}-T_{i-1,j}}{2\Delta r}\right)+\left(\frac{T_{i+1,j}-2T_{i,j}+T_{i-1,j}}{\Delta r^2}\right)+\left(\frac{T_{i,j+1}-2T_{i,j}+T_{i,j-1}}{\Delta z^2}\right)\right]
$$

Reordenando los términos:

$$
\frac{\rho c_p \Phi}{2\Delta z}[(i\Delta r)^2-R^2]-\frac{k}{\Delta z^2}T_{i,j+1} - (\frac{\rho c_p \Phi}{2\Delta z}[(i\Delta r)^2-R^2]+\frac{k}{\Delta z^2})T_{i,j-1} + k(\frac{2}{\Delta z^2}+\frac{2}{\Delta r^2})T_{i,j} - \frac{k}{\Delta r^{2}}(1+\frac{1}{2i})T_{i+1,j} - \frac{k}{\Delta r^{2}}(1 -\frac{1}{2i})T_{i-1,j} = 0
$$

Y discretizando las condiciones de borde:

* CB2: Temperatura de entrada (Dirichlet)
    $$T_{i,0} = T_0$$

* CB3: Cambio en la frontera (Neumann, diferencia finita adelantada de segundo orden para $j=0$)
    $$-3T_{i,0} + 4T_{i,1} - T_{i,2} = 0$$
    $$4T_{i,1} - T_{i,2} = 3T_0$$

* CB4: Simetría del perfil radialmente (Neumann, usando el punto fantasma $T_{-1,j} = T_{1,j}$ o diferencia finita adelantada para $i=0$)
    $$T_{1,j} - T_{-1,j}=0$$
    $$-3T_{0,j} + 4T_{1,j} - T_{2,j} = 0$$

* CB5: Calor por radiación (Neumann)

$$q_{\text{rad}} = k\frac{T_{N_r,j}-T_{N_r-1,j}}{\Delta r}$$

$$ \frac{q_{\text{rad}}\Delta r}{k} = T_{N_r,j}-T_{N_r-1,j} $$

## 2. Resolución por SOR

### 2.1 Descripción del método
El método SOR es un método iterativo para resolver sistemas del tipo 

$$ a_{ij}T_{i+1,j}+b_{ij}T_{i-1,j}+c_{ij}T_{i,j+1}+d_{ij}T_{i,j-1}+e_{ij}T_{i,j}-f_{ij}=0 $$

Esto, mediante la iteración de los terminos variables $T$ junto a sus vecinos espaciales. Es decir, se cambian los valores $T_{i,j}$ junto a los $T_{i+1,j}$, $T_{i-1,j}$ y demás, hasta que el sistema deja de presentar variaciones de estos términos (cuando alcanza una variación menor a un valor especificado). La iteración se hace prediciendo los valores siguientes tal que $T^{k+1} = T^k + \omega \cdot (T^{\text{GS}} - T^k)$.

Este punto se alcanza más rápido o lento dependiendo de la combinación de tolerancia a variaciones y al parámetro de sobre relajación &\omega&. Una tolerancia menor significa que el sistema alcanzará la convergencia más rápido pero con menos iteraciones, por lo que parte del sistema puede quedar sin resolverse correctamente. El término $\omega$ por otro lado controla que tan alejada del valor actual se encontrará la predicción, acelerando la convergencia a costa de menos estabilidad del modelo. 

La estabilidad se refiere a cómo se comporta la iteración respecto a los errores de cálculo que se pueden producir. Un modelo estable hace que el error al calcular valores disminuyan conforme se itera, permitiendo llegar a una convergencia. Sin embargo, una iteración inestable implica que este error crece y en cada iteración se aleja más de una solución, por lo que no se alcanza la convergencia.

### 2.2 Definición de parámetros
Algunos de los parámetros más relevantes entonces para la convergencia del modelo planteado son la tolerancia (expresada como la variable tol en el código del repositorio), el parámetro de sobre relajación (llamado w en el código) y el número máximo de iteraciones. Este último es implementado por motivos prácticos para parar la iteración en caso no se alcance la convergencia antes de cierto umbral. 

En este caso, se determinó que los valores $tol = 1e-4$, $w = 1.6$, y max_iter = 25000 lograron convergencia en buen tiempo y calidad de resultados.

### 2.3 Descripción de la iteración
La iteración en sí se divide en bloques dentro de un ciclo for que itera sobre todos los valores i,j. Estos bloques atacan con distintos métodos la predicción dependiendo de la ubicación dentro del tubo; el primer bloque trata los nodos centrales del tubo (i=0, j>0); el segundo la columna siguiente a la entrada del tubo (i>0, j=1); el tercero implementa la condición de borde de flux a la entrada; el cuarto itera sobre todos los nodos centrales del tubo; y el quinto actualiza los valores de temperatura de salida, es decir, en el último nodo.

El código además actualiza en directo el error (residuo) entre los valores obtenidos en la iteración pasada y la actual, manteniendo la iteración hasta alcanzar el máximo definido o que este residuo sea menor a la tolerancia especificada. 

## 3. Resultados
### 3.1 Ejecución del código
El código adjunto al repositorio se debe correr de forma lineal, procurando que los coeficientes y parámetros sean definidos antes de iniciar la iteración. Es posible cambiar los valores de tol, w y max_iter aún para experimentar con ellos.

### 3.2 Datos y gráficos obtenidos 
A continuación, se muestra el perfil de temperaturas obtenido.
<img width="846" height="547" alt="image" src="https://github.com/user-attachments/assets/6afd2d40-cece-4d63-a21e-21a0a8443cdd" />
