import pandas as pd
import sympy as sp

# 1. Definición simbólica de los componentes (independientes, sin igualdades)
R1, R2, C1, C2 = sp.symbols("R1 R2 C1 C2", positive=True, real=True)

w0 = 1 / sp.sqrt(R1 * R2 * C1 * C2)
Q = sp.sqrt(R1 * R2 * C1 * C2) / (C2 * (R1 + R2))


def sensibilidad(parametro, var):
    return sp.simplify(sp.diff(parametro, var) * var / parametro)


componentes = [R1, R2, C1, C2]
S_w0 = {c: sensibilidad(w0, c) for c in componentes}
S_Q = {c: sensibilidad(Q, c) for c in componentes}

print("--- 1. SENSIBILIDADES SIMBÓLICAS (ANTES DE IMPONER IGUALDADES) ---")
for c in componentes:
    print(f"{c}:  S_w0 = {S_w0[c]}   |   S_Q = {S_Q[c]}")

# Sensibilidades de la etapa RC de primer orden (sin Q definido)
R0, C0 = sp.symbols("R0 C0", positive=True, real=True)
w_rc = 1 / (R0 * C0)
print("\n--- ETAPA DE PRIMER ORDEN (RC) ---")
print(f"R0: S_wRC = {sensibilidad(w_rc, R0)}")
print(f"C0: S_wRC = {sensibilidad(w_rc, C0)}")

# 2. Valores comerciales E12 elegidos por el Grupo 1
nominales = {R1: 9.1e3, R2: 9.1e3, C1: 27e-9, C2: 6.8e-9}
w0_nom = float(w0.subs(nominales))
Q_nom = float(Q.subs(nominales))

# 3. Validación con perturbaciones finitas del ±1%
filas = []
for c in componentes:
    sw0_num = float(S_w0[c].subs(nominales))
    sQ_num = float(S_Q[c].subs(nominales))

    for signo, etiqueta in [(1, "+1%"), (-1, "-1%")]:
        pert = nominales.copy()
        pert[c] = nominales[c] * (1 + signo * 0.01)

        w0_nuevo = float(w0.subs(pert))
        Q_nuevo = float(Q.subs(pert))

        # Variación relativa exacta (%)
        var_w0_exacta = ((w0_nuevo - w0_nom) / w0_nom) * 100
        var_Q_exacta = ((Q_nuevo - Q_nom) / Q_nom) * 100

        # Variación estimada mediante sensibilidad lineal (%)
        var_w0_lin = sw0_num * (signo * 1.0)
        var_Q_lin = sQ_num * (signo * 1.0)

        filas.append({
            "Componente": str(c),
            "Perturbacion": etiqueta,
            "S_w0": sw0_num,
            "Delta_w0_lin_%": round(var_w0_lin, 4),
            "Delta_w0_exacta_%": round(var_w0_exacta, 4),
            "S_Q": sQ_num,
            "Delta_Q_lin_%": round(var_Q_lin, 4),
            "Delta_Q_exacta_%": round(var_Q_exacta, 4),
        })

df = pd.DataFrame(filas)
print("\n--- 2. COMPARACIÓN LINEAL VS PERTURBACIÓN FINITA (±1%) ---")
print(df.to_string(index=False))

# Exportación de resultados a CSV para el repositorio Git
df.to_csv("sensibilidades_punto6.csv", index=False)
print("\nArchivo 'sensibilidades_punto6.csv' generado correctamente.")