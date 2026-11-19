# Grupo 1 - puntos 1, 2 y 3
## Punto 1: diseno y verificacion
|Tipo|Familia|Orden|Frecuencia de diseno (Hz)|A(fp) (dB)|A(fs) (dB)|
|---|---|---:|---:|---:|---:|
|Pasabajos|Butterworth|3|1252.5764|1.0000|30.2594|
|Pasabajos|Chebyshev I|3|1000.0000|1.0000|41.8798|
|Pasaaltos|Butterworth|3|3193.4180|1.0000|30.2594|
|Pasaaltos|Chebyshev I|3|4000.0000|1.0000|41.8798|

## Punto 2: comparacion
Se comparan magnitud, fase continua, retardo de grupo y respuesta al escalon. Los graficos estan en figuras/. El retardo se calcula como -d fase/d omega, no como fase/frecuencia. La variacion de retardo se mide entre DC y 1 kHz para LP y de 4 a 40 kHz para HP. El establecimiento usa una banda de +/-0.02 respecto de la salida final para un escalon unitario.

## Punto 3: filtro elegido
Se selecciona el pasabajos Butterworth de orden 3. Cumple sin rizado y permite una cascada RC + buffer + Sallen-Key de ganancia unitaria. El factor de calidad de la seccion de segundo orden es Q=1, no 0.707. Este ultimo corresponde a un Butterworth completo de segundo orden.

La frecuencia caracteristica Butterworth es su corte a -3.0103 dB; en Chebyshev I la frecuencia de diseno es el borde de paso de rizado, no necesariamente el corte a -3 dB. 
