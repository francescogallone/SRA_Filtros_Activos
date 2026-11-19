# Grupo 1 - Nivel I: puntos 1, 2 y 3

- Punto 1: orden minimo, frecuencia de diseno, polos, ceros, ganancia, TF y secciones SOS de LP/HP Butterworth y Chebyshev I.
- Punto 2: comparaciones de magnitud, fase continua, retardo de grupo y escalon. Tablas de establecimiento, sobreimpulso y variacion de retardo.
- Punto 3: seleccion del LP Butterworth de orden 3 y propuesta RC + buffer + Sallen-Key (Q=1).

## Ejecutar

```bash
python3 -m pip install -r requirements.txt
python3 puntos_1_2_3.py
```

Se genera resultados_puntos_1_2_3/ con PNG, PDF, CSV, JSON y tabla_temporal.tex. Para elegir otra carpeta: `python3 puntos_1_2_3.py --out resultados`.


Especificaciones: Ap=1 dB, As=24 dB; LP fp=1 kHz, fs=4 kHz; HP fp=4 kHz, fs=1 kHz. El script usa Hz como entrada y convierte a rad/s para SciPy. Se verifican estabilidad, cumplimiento y equivalencia de ZPK/TF/SOS.
