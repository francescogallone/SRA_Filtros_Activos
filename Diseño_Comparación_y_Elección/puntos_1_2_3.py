"""Grupo 1 - Nivel I, puntos 1, 2 y 3.
1. Diseno de LP/HP Butterworth y Chebyshev I.
2. Magnitud, fase, retardo de grupo y respuesta al escalon.
3. Comparacion y seleccion del LP Butterworth de orden 3.

Ejecutar: python3 puntos_1_2_3.py
Dependencias: numpy, scipy y matplotlib.
Modelos analogicos ideales: no simula componentes ni operacionales reales.
"""
import json
import platform
import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

import numpy as np
import scipy
from scipy import signal

AP, AS = 1.0, 24.0
ESPECIFICACIONES = {
    'Pasabajos': (1000.0, 4000.0, 'lowpass'),
    'Pasaaltos': (4000.0, 1000.0, 'highpass'),
}

def disenar(tipo, familia):
    fp, fs, btype = ESPECIFICACIONES[tipo]
    wp, ws = 2*np.pi*fp, 2*np.pi*fs
    if familia == 'Butterworth':
        n, wn = signal.buttord(wp, ws, AP, AS, analog=True)
        z, p, k = signal.butter(n, wn, btype=btype, analog=True, output='zpk')
    else:
        n, wn = signal.cheb1ord(wp, ws, AP, AS, analog=True)
        z, p, k = signal.cheby1(n, AP, wn, btype=btype, analog=True, output='zpk')
    b, a = signal.zpk2tf(z, p, k)
    sos = signal.zpk2sos(z, p, k, analog=True)
    _, h = signal.freqs_zpk(z, p, k, worN=2*np.pi*np.array([fp, fs]))
    atenuacion = -20*np.log10(np.abs(h))
    # Verificación numérica adicional dentro de las bandas, sin usar
    # polinomios expandidos para evaluar la respuesta.
    if btype == 'lowpass':
        paso = np.r_[0.0, np.geomspace(fp*1e-8, fp, 20001)]
        rechazo = np.geomspace(fs, fs*1e6, 20001)
    else:
        paso = np.geomspace(fp, fp*1e8, 20001)
        rechazo = np.geomspace(fs*1e-6, fs, 20001)
    def atenuar(f):
        _, hf = signal.freqs_zpk(z, p, k, worN=2*np.pi*f)
        return -20*np.log10(np.maximum(np.abs(hf), np.finfo(float).tiny))
    max_paso = float(np.max(atenuar(paso)))
    min_rechazo = float(np.min(atenuar(rechazo)))
    estable = bool(np.all(np.real(p) < 0))
    assert estable and max_paso <= AP+1e-8 and min_rechazo >= AS-1e-8
    return dict(tipo=tipo, familia=familia, orden=int(n),
                wn_rad_s=float(wn), fn_Hz=float(wn/(2*np.pi)),
                ceros=[[float(v.real), float(v.imag)] for v in z],
                polos=[[float(v.real), float(v.imag)] for v in p],
                k=float(k), numerador=b.tolist(), denominador=a.tolist(),
                sos=sos.tolist(), estable=estable,
                atenuacion_fp_dB=float(atenuacion[0]),
                atenuacion_fs_dB=float(atenuacion[1]),
                max_atenuacion_paso_malla_dB=max_paso,
                min_atenuacion_rechazo_malla_dB=min_rechazo)


def recuperar_zpk(r):
    z=np.array([complex(*v) for v in r['ceros']],dtype=complex)
    p=np.array([complex(*v) for v in r['polos']],dtype=complex)
    return z,p,r['k']


def respuesta(r,f):
    z,p,k=recuperar_zpk(r)
    return signal.freqs_zpk(z,p,k,worN=2*np.pi*np.asarray(f))[1]


def retardo_grupo(r,f):
    """-d arg(H)/dw, en segundos. Derivada analitica por polos y ceros.
    Los ceros en el origen del HP no aportan pendiente para f>0.
    """
    z,p,_=recuperar_zpk(r)
    jw=2j*np.pi*np.asarray(f)
    tau=np.zeros(jw.shape,dtype=float)
    for pole in p:
        tau+=np.real(1/(jw-pole))
    for zero in z:
        if zero==0:
            continue
        tau-=np.real(1/(jw-zero))
    return tau


def guardar_tabla(path,header,rows):
    with path.open('w',newline='',encoding='utf-8') as stream:
        wr=csv.writer(stream);wr.writerow(header);wr.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',default='resultados_puntos_1_2_3',help='Carpeta de resultados')
    args=parser.parse_args()
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    figuras=out/'figuras';figuras.mkdir(exist_ok=True)
    resultados=[disenar(t,f) for t in ESPECIFICACIONES for f in ('Butterworth','Chebyshev I')]
    versiones=dict(Python=platform.python_version(),NumPy=np.__version__,SciPy=scipy.__version__,Matplotlib=matplotlib.__version__)
    f=np.geomspace(10,1e6,6000)

    # Punto 1: estabilidad y equivalencia numerica de ZPK, TF y cascada SOS.
    for r in resultados:
        assert r['orden']==3
        href=respuesta(r,f)
        htf=signal.freqs(r['numerador'],r['denominador'],worN=2*np.pi*f)[1]
        hsos=np.ones(f.shape,dtype=complex)
        for row in r['sos']:
            hsos*=np.polyval(row[:3],2j*np.pi*f)/np.polyval(row[3:],2j*np.pi*f)
        assert np.allclose(href,htf,rtol=1e-8,atol=1e-12)
        assert np.allclose(href,hsos,rtol=1e-8,atol=1e-12)
        dphase=-np.gradient(np.unwrap(np.angle(href)),2*np.pi*f)
        assert np.max(np.abs(dphase[5:-5]-retardo_grupo(r,f[5:-5])))<1e-7
        poles=recuperar_zpk(r)[1]
        r['Q_par_complejo']=float(abs(poles[np.argmax(poles.imag)])/(-2*poles[np.argmax(poles.imag)].real))
        name=('LP' if r['tipo']=='Pasabajos' else 'HP')+('_Butterworth' if r['familia']=='Butterworth' else '_Chebyshev_I')
        guardar_tabla(out/(name+'_frecuencia.csv'),['f_Hz','ganancia_dB','fase_deg','retardo_grupo_s'],
                      zip(f,20*np.log10(abs(href)),np.rad2deg(np.unwrap(np.angle(href))),retardo_grupo(r,f)))

    guardar_tabla(out/'verificacion.csv',
        ['tipo','familia','orden','fn_Hz','atenuacion_fp_dB','atenuacion_fs_dB','estable','Q_par_complejo'],
        ([r[k] for k in ['tipo','familia','orden','fn_Hz','atenuacion_fp_dB','atenuacion_fs_dB','estable','Q_par_complejo']] for r in resultados))

    fig,axs=plt.subplots(2,2,figsize=(9,7))
    for ax,r in zip(axs.flat,resultados):
        z,p,_=recuperar_zpk(r)
        ax.scatter(p.real/(2*np.pi),p.imag/(2*np.pi),marker='x',s=70,label='Polos')
        if len(z):ax.scatter(z.real,z.imag,facecolors='none',edgecolors='tab:orange',s=75,label='3 ceros en origen')
        ax.axhline(0,c='gray',lw=.6);ax.axvline(0,c='gray',lw=.6);ax.grid(alpha=.3)
        ax.set(title=r['tipo']+' '+r['familia'],xlabel='Re(s)/(2 pi) [Hz]',ylabel='Im(s)/(2 pi) [Hz]');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(figuras/'polos_ceros.png',dpi=170);plt.close(fig)

    # Punto 2: comparaciones en frecuencia. Hz en el eje; rad/s en el calculo.
    for tipo,suffix in [('Pasabajos','LP'),('Pasaaltos','HP')]:
        fig,axs=plt.subplots(3,1,figsize=(8,9),sharex=True)
        for r in [r for r in resultados if r['tipo']==tipo]:
            hf=respuesta(r,f)
            axs[0].semilogx(f,20*np.log10(abs(hf)),label=r['familia'])
            axs[1].semilogx(f,np.rad2deg(np.unwrap(np.angle(hf))),label=r['familia'])
            axs[2].semilogx(f,retardo_grupo(r,f)*1000,label=r['familia'])
        for ax in axs:
            for borde in [1000,4000]:ax.axvline(borde,color='gray',ls=':',lw=.8)
            ax.grid(True,which='both',alpha=.25);ax.legend()
        axs[0].axhline(-AP,color='red',ls='--',lw=.8);axs[0].axhline(-AS,color='red',ls='--',lw=.8)
        axs[0].set(ylabel='Ganancia (dB)',ylim=(-100,3),title=tipo+': Butterworth y Chebyshev I')
        axs[1].set(ylabel='Fase (grados)');axs[2].set(ylabel='Retardo de grupo (ms)',xlabel='Frecuencia (Hz)')
        fig.tight_layout()
        for extension in ['png','pdf']:fig.savefig(figuras/('comparacion_'+suffix+'.'+extension),dpi=170)
        plt.close(fig)

    # Punto 2: escalon unitario, con paso temporal de 1 us durante 12 ms.
    t=np.arange(12001)*1e-6
    fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True)
    temporales=[]
    for r in resultados:
        sys=signal.TransferFunction(r['numerador'],r['denominador'])
        tt,y=signal.step(sys,T=t)
        final=1.0 if r['tipo']=='Pasabajos' else 0.0
        # Banda absoluta +/-0.02 respecto del valor final; escalon de amplitud uno.
        fuera=np.flatnonzero(abs(y-final)>.02)
        settling=float(tt[fuera[-1]+1]) if len(fuera) and fuera[-1]+1<len(tt) else (0. if not len(fuera) else None)
        banda=np.linspace(0,1000,10001) if r['tipo']=='Pasabajos' else np.linspace(4000,40000,10001)
        fila=dict(tipo=r['tipo'],familia=r['familia'],ymax=float(y.max()),ymin=float(y.min()),
                  t2pct_s=settling,variacion_tau_s=float(np.ptp(retardo_grupo(r,banda))),
                  sobreimpulso_pct=float(max(0,(y.max()-1)*100)) if final==1 else None)
        temporales.append(fila)
        ax=axs[0 if r['tipo']=='Pasabajos' else 1];ax.plot(tt*1000,y,label=r['familia'])
        key=('LP' if final else 'HP')+('_Butterworth' if r['familia']=='Butterworth' else '_Chebyshev_I')
        guardar_tabla(out/(key+'_escalon.csv'),['t_s','salida'],zip(tt,y))
    for ax,name in zip(axs,['LP','HP']):
        ax.set(title='Escalon unitario: '+name,ylabel='Salida / amplitud del escalon',xlim=(0,5));ax.legend();ax.grid(alpha=.3)
    axs[-1].set_xlabel('Tiempo (ms)');fig.tight_layout()
    for extension in ['png','pdf']:fig.savefig(figuras/('escalon.'+extension),dpi=170)
    plt.close(fig)
    guardar_tabla(out/'tabla_temporal.csv',list(temporales[0]),([row[k] for k in temporales[0]] for row in temporales))
    tabla=[r'\begin{tabular}{llrrrr}\toprule',
           r'Tipo & Familia & $y_{\max}$ & $y_{\min}$ & $t_{2\%}$ (ms) & $\Delta\tau_g$ (ms)\\\midrule']
    for v in temporales:
        ts='--' if v['t2pct_s'] is None else f"{v['t2pct_s']*1e3:.3f}"
        tabla.append(f"{'LP' if v['tipo']=='Pasabajos' else 'HP'} & {v['familia']} & {v['ymax']:.3f} & {v['ymin']:.3f} & {ts} & {v['variacion_tau_s']*1e3:.3f} "+r'\\')
    tabla+=[r'\bottomrule\end{tabular}'];(out/'tabla_temporal.tex').write_text('\n'.join(tabla))

    # Punto 3: seleccion fundada en las respuestas obtenidas, no optimizacion de componentes.
    elegido=next(r for r in resultados if r['tipo']=='Pasabajos' and r['familia']=='Butterworth')
    seleccion={
        'tipo':'Pasabajos','familia':'Butterworth','orden':3,
        'motivos':['Cumple Ap=1 dB y As=24 dB.','No presenta rizado en banda de paso.',
                   'Se implementa como RC + buffer + Sallen-Key de ganancia unitaria.',
                   'El par complejo tiene Q=1; menor que el del Chebyshev estudiado.'],
        'omega_c_rad_s':elegido['wn_rad_s'],'Q_segundo_orden':elegido['Q_par_complejo'],
        'etapa_1':'wc/(s+wc)','etapa_2':'wc^2/(s^2+wc*s+wc^2)',
        'alcance':'Seleccion y propuesta de topologia. No dimensiona ni simula el circuito real.'}
    contenido=dict(versiones=versiones,Ap_dB=AP,As_dB=AS,filtros=resultados,temporales=temporales,seleccion=seleccion)
    (out/'resultados.json').write_text(json.dumps(contenido,indent=2,ensure_ascii=False),encoding='utf-8')
    md=['# Grupo 1 - puntos 1, 2 y 3','## Punto 1: diseno y verificacion',
        '|Tipo|Familia|Orden|Frecuencia de diseno (Hz)|A(fp) (dB)|A(fs) (dB)|',
        '|---|---|---:|---:|---:|---:|']
    for r in resultados:md.append(f"|{r['tipo']}|{r['familia']}|{r['orden']}|{r['fn_Hz']:.4f}|{r['atenuacion_fp_dB']:.4f}|{r['atenuacion_fs_dB']:.4f}|")
    md+=['','## Punto 2: comparacion','Se comparan magnitud, fase continua, retardo de grupo y respuesta al escalon. Los graficos estan en figuras/. El retardo se calcula como -d fase/d omega, no como fase/frecuencia. La variacion de retardo se mide entre DC y 1 kHz para LP y de 4 a 40 kHz para HP. El establecimiento usa una banda de +/-0.02 respecto de la salida final para un escalon unitario.',
        '','## Punto 3: filtro elegido','Se selecciona el pasabajos Butterworth de orden 3. Cumple sin rizado y permite una cascada RC + buffer + Sallen-Key de ganancia unitaria. El factor de calidad de la seccion de segundo orden es Q=1, no 0.707. Este ultimo corresponde a un Butterworth completo de segundo orden.',
        '','La frecuencia caracteristica Butterworth es su corte a -3.0103 dB; en Chebyshev I la frecuencia de diseno es el borde de paso de rizado, no necesariamente el corte a -3 dB. Todos los resultados corresponden a filtros analogicos ideales.']
    (out/'resumen_puntos_1_2_3.md').write_text('\n'.join(md),encoding='utf-8')
    print('Puntos 1, 2 y 3 calculados. Resultados:',out)
    for r in resultados:print(f"{r['tipo']} {r['familia']}: orden {r['orden']}, fn={r['fn_Hz']:.4f} Hz, A(fp)={r['atenuacion_fp_dB']:.4f} dB, A(fs)={r['atenuacion_fs_dB']:.4f} dB")
    print('Seleccion: LP Butterworth, orden 3; RC + buffer + Sallen-Key, Q=1.')


if __name__=='__main__':
    main()

