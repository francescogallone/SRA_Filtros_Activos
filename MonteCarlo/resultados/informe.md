# Nivel II: analisis estadistico en Python (modelo ideal)

**Estado:** estudio reproducible calculado; no se ejecutaron corridas LTspice. No se ha demostrado una funcion no disponible en esa instalacion. La aceptacion de esta sustitucion debe confirmarse con la catedra.

Se suponen R0=12 kohm, R1=R2=9,1 kohm; C0=10 nF; C realimentacion=22 nF || 5,6 nF; C a masa=6,8 nF. Confirmar estos nominales con el montaje antes de usar resultados.

Se usa H(s)=1/[(1+s R0 C0)(1+s Cg(R1+R2)+s^2 R1 R2(C22+C56)Cg)]. Los seguidores son ideales: no incluye GBW, slew rate, offset, ruido, carga de salida ni dispersion de los TL082. No explica la ganancia fisica mayor que uno a 100 Hz.

Se generan 5000 vectores distintos con PCG64, semilla 20261118. Cada coordenada es uniforme independiente en [-1,1]. p=p_nom(1+t*e). Las mismas muestras normalizadas se reutilizan para comparar las tolerancias. Los CSV documentan todos los valores.

Cumplimiento: minimo de ganancia entre DC y 1 kHz >= -1 dB, y ganancia a 4 kHz <= -24 dB. El minimo en paso se calcula incluyendo extremos interiores del polinomio de |H|^-2. Corte: -3,0103 dB respecto de H(0)=1. No se usa la ganancia fisica medida como nueva referencia.

Sensibilidad lineal: Delta G_dB = suma a_i Delta p_i/p_i; a_i=dG_dB/dln(p_i). Derivadas centrales verificadas con dos pasos. sigma^2=suma a_i^2 t_i^2/3. Se comparan tambien predicciones muestra a muestra. La aproximacion lineal puede perder precision con tolerancias amplias.

| Escenario | Paso (%) | Rechazo (%) | Ambos (%) | Corte medio (Hz) | sigma G1k (dB) | sigma lineal G1k (dB) |

|---|---:|---:|---:|---:|---:|---:|

|A_R5_C10|75.80|100.00|75.80|1299.96|0.3450|0.3447|

|B_R1_C2|100.00|100.00|100.00|1298.74|0.0690|0.0689|

Los porcentajes describen el modelo y las distribuciones elegidas, no garantizan todas las unidades ni el peor caso. El JSON incluye intervalos de confianza de Wilson del 95%.

Para repetir: python -m pip install -r requirements.txt ; python montecarlo_filtro.py --n 5000 --seed 20261118 --out resultados. Para otro conjunto cambiar --seed; para R0 distinto, usar --r0 en ohmios. Mantener identica la configuracion y guardar el resumen.json para trazabilidad.

Antes de presentar: confirmar nominales y topologia; validar la respuesta nominal frente a LTspice y justificar la alternativa segun la consigna. Si LTspice admite las funciones y no se autoriza la sustitucion, estas muestras pueden emplearse para automatizar las corridas exigidas en LTspice.