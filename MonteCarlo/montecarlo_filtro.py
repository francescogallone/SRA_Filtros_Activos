"""Monte Carlo del RC + Sallen-Key con seguidores ideales. No ejecuta LTspice."""
from pathlib import Path
import argparse, csv, json, hashlib, platform
import numpy as np
import scipy
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES=['R0','R1','R2','C0','C22','C56','Cg']
NOM=np.array([12e3,9.1e3,9.1e3,10e-9,22e-9,5.6e-9,6.8e-9])

def transfer(p,f):
    r0,r1,r2,c0,ca,cb,cg=p
    s=2j*np.pi*np.asarray(f)
    return 1/((1+s*r0*c0)*(1+s*cg*(r1+r2)+s*s*r1*r2*(ca+cb)*cg))

def gain(p,f):
    return 20*np.log10(np.abs(transfer(p,f)))

def metrics(p):
    # Extremos exactos en paso: |H|^-2 es un polinomio cubico en z=(f/1000)^2.
    r0,r1,r2,c0,ca,cb,cg=p
    w=2*np.pi*1000
    a=w*w*r1*r2*(ca+cb)*cg
    b=w*cg*(r1+r2)
    t=(w*r0*c0)**2
    u=b*b-2*a
    v=a*a
    roots=np.roots([3*t*v,2*(v+t*u),t+u])
    freqs=[0.,1000.]
    freqs += [1000*np.sqrt(z.real) for z in roots if abs(z.imag)<1e-10 and 0<z.real<1]
    gs=gain(p,np.array(freqs))
    fc=brentq(lambda f: float(gain(p,f))+10*np.log10(2),1.,1e6)
    return np.array([gain(p,1000),gain(p,4000),fc,min(gs),max(gs)])

MET=['G_1k_dB','G_4k_dB','fc_3dB_Hz','G_paso_min_dB','G_paso_max_dB']

def writecsv(path,header,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        wr=csv.writer(f);wr.writerow(header);wr.writerows(rows)

def wilson(k,n):
    z=1.95996398454;p=k/n;den=1+z*z/n
    c=(p+z*z/(2*n))/den
    h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [100*(c-h),100*(c+h)]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n',type=int,default=5000)
    ap.add_argument('--seed',type=int,default=20261118)
    ap.add_argument('--out',default='resultados')
    ap.add_argument('--r0',type=float,default=12000,help='Ohmios; confirmar contra el construido')
    args=ap.parse_args()
    if args.n<100:ap.error('Se requieren al menos 100 muestras')
    NOM[0]=args.r0
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    rng=np.random.Generator(np.random.PCG64(args.seed))
    eps=rng.uniform(-1,1,(args.n,7))
    assert len(np.unique(eps,axis=0))==args.n
    writecsv(out/'muestras_uniformes.csv',['corrida']+NAMES,([i+1,*e] for i,e in enumerate(eps)))
    nominal=metrics(NOM)
    # Derivadas de G(dB) respecto de ln(p), calculadas por diferencia central.
    def sensitivity(h):
        arr=[]
        for i in range(7):
            pp=NOM.copy();pm=NOM.copy();pp[i]*=np.exp(h);pm[i]*=np.exp(-h)
            arr.append((gain(pp,np.array([1000,4000]))-gain(pm,np.array([1000,4000])))/(2*h))
        return np.array(arr)
    sens=sensitivity(1e-5)
    assert np.allclose(sens,sensitivity(5e-6),rtol=1e-6,atol=1e-7)
    assert np.isclose(transfer(NOM,0),1)
    assert np.allclose(nominal[:2],[-.7285561491578568,-29.457275440502148],atol=1e-8) if args.r0==12000 else True
    writecsv(out/'sensibilidades.csv',['componente','dG1k_dlnp','dG4k_dlnp'],([name,*row] for name,row in zip(NAMES,sens)))
    summary={};fig,axs=plt.subplots(2,3,figsize=(13,7))
    f=np.geomspace(10,1e5,600);fb,ab=plt.subplots(figsize=(8,5))
    for index,(label,tr,tc) in enumerate([('A_R5_C10',.05,.10),('B_R1_C2',.01,.02)]):
        tol=np.array([tr]*3+[tc]*4)
        samples=NOM*(1+eps*tol)
        assert np.all(abs(samples/NOM-1)<=tol+1e-12)
        vals=np.array([metrics(p) for p in samples])
        goodp=vals[:,3]>=-1
        goods=vals[:,1]<=-24
        good=goodp&goods
        rows=([i+1,*p,*m,int(gp),int(gr),int(gt)] for i,(p,m,gp,gr,gt) in enumerate(zip(samples,vals,goodp,goods,good)))
        writecsv(out/(label+'.csv'),['corrida']+NAMES+MET+['cumple_paso','cumple_rechazo','cumple_total'],rows)
        lin=nominal[:2]+(eps*tol)@sens
        stats={name:{'media':float(v.mean()),'sigma_muestral':float(v.std(ddof=1)),
               'min':float(v.min()),'max':float(v.max()),'p5':float(np.percentile(v,5)),
               'p95':float(np.percentile(v,95))} for name,v in zip(MET,vals.T)}
        stats['cumplimiento']={'paso_pct':100*float(goodp.mean()),'rechazo_pct':100*float(goods.mean()),
            'total_pct':100*float(good.mean()),'total_aprobados':int(good.sum()),'IC95_Wilson_pct':wilson(good.sum(),args.n)}
        stats['lineal']={'sigma_teorica_dB':np.sqrt(np.sum((sens*tol[:,None])**2/3,axis=0)).tolist(),
            'sigma_muestreada_dB':lin.std(axis=0,ddof=1).tolist(),
            'RMSE_no_lineal_vs_lineal_dB':np.sqrt(np.mean((vals[:,:2]-lin)**2,axis=0)).tolist()}
        summary[label]=stats
        for j,(metric,title) in enumerate([(0,'Ganancia a 1 kHz (dB)'),(1,'Ganancia a 4 kHz (dB)'),(2,'Corte a -3,0103 dB (Hz)')]):
            ax=axs[index,j];ax.hist(vals[:,metric],bins=40,alpha=.7,label='Transferencia completa')
            if j<2:
                ax.hist(lin[:,j],bins=40,histtype='step',color='black',label='Aproximacion lineal')
                ax.axvline([-1,-24][j],color='red',ls='--',label='Limite')
            ax.set(xlabel=title,ylabel='Cantidad',title=label);ax.legend(fontsize=7)
        for p in samples[:100]:ab.semilogx(f,gain(p,f),alpha=.04,color=['tab:orange','tab:blue'][index])
    fig.suptitle('Monte Carlo: 5000 muestras por escenario; seguidores ideales' if args.n==5000 else f'Monte Carlo: {args.n} muestras por escenario; seguidores ideales')
    fig.tight_layout();fig.savefig(out/'histogramas.png',dpi=160);plt.close(fig)
    ab.semilogx(f,gain(NOM,f),'k',lw=2,label='Nominal ideal')
    ab.set(xlabel='Frecuencia (Hz)',ylabel='Ganancia (dB)',ylim=(-65,3),title='100 respuestas de cada escenario: naranja A, azul B')
    ab.grid(True,which='both',alpha=.3);ab.legend();fb.tight_layout();fb.savefig(out/'dispersion_bode.png',dpi=160);plt.close(fb)
    payload={'modelo':'RC + Sallen-Key, seguidores ideales; no es una ejecucion LTspice',
        'n':args.n,'seed':args.seed,'generador':'PCG64','nominales_SI':dict(zip(NAMES,NOM.tolist())),
        'nominal':dict(zip(MET,nominal.tolist())), 'escenarios':summary,
        'versiones':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'sha256_muestras':hashlib.sha256((out/'muestras_uniformes.csv').read_bytes()).hexdigest()}
    (out/'resumen.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8')
    paragraphs=['# Nivel II: analisis estadistico en Python (modelo ideal)',
        '**Estado:** estudio reproducible calculado; no se ejecutaron corridas LTspice. No se ha demostrado una funcion no disponible en esa instalacion. La aceptacion de esta sustitucion debe confirmarse con la catedra.',
        'Se suponen R0=12 kohm, R1=R2=9,1 kohm; C0=10 nF; C realimentacion=22 nF || 5,6 nF; C a masa=6,8 nF. Confirmar estos nominales con el montaje antes de usar resultados.',
        'Se usa H(s)=1/[(1+s R0 C0)(1+s Cg(R1+R2)+s^2 R1 R2(C22+C56)Cg)]. Los seguidores son ideales: no incluye GBW, slew rate, offset, ruido, carga de salida ni dispersion de los TL082. No explica la ganancia fisica mayor que uno a 100 Hz.',
        f'Se generan {args.n} vectores distintos con PCG64, semilla {args.seed}. Cada coordenada es uniforme independiente en [-1,1]. p=p_nom(1+t*e). Las mismas muestras normalizadas se reutilizan para comparar las tolerancias. Los CSV documentan todos los valores.',
        'Cumplimiento: minimo de ganancia entre DC y 1 kHz >= -1 dB, y ganancia a 4 kHz <= -24 dB. El minimo en paso se calcula incluyendo extremos interiores del polinomio de |H|^-2. Corte: -3,0103 dB respecto de H(0)=1. No se usa la ganancia fisica medida como nueva referencia.',
        'Sensibilidad lineal: Delta G_dB = suma a_i Delta p_i/p_i; a_i=dG_dB/dln(p_i). Derivadas centrales verificadas con dos pasos. sigma^2=suma a_i^2 t_i^2/3. Se comparan tambien predicciones muestra a muestra. La aproximacion lineal puede perder precision con tolerancias amplias.',
        '| Escenario | Paso (%) | Rechazo (%) | Ambos (%) | Corte medio (Hz) | sigma G1k (dB) | sigma lineal G1k (dB) |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for label,st in summary.items():
        c=st['cumplimiento'];paragraphs.append(f"|{label}|{c['paso_pct']:.2f}|{c['rechazo_pct']:.2f}|{c['total_pct']:.2f}|{st['fc_3dB_Hz']['media']:.2f}|{st['G_1k_dB']['sigma_muestral']:.4f}|{st['lineal']['sigma_teorica_dB'][0]:.4f}|")
    paragraphs+=['Los porcentajes describen el modelo y las distribuciones elegidas, no garantizan todas las unidades ni el peor caso. El JSON incluye intervalos de confianza de Wilson del 95%.',
        'Para repetir: python -m pip install -r requirements.txt ; python montecarlo_filtro.py --n 5000 --seed 20261118 --out resultados. Para otro conjunto cambiar --seed; para R0 distinto, usar --r0 en ohmios. Mantener identica la configuracion y guardar el resumen.json para trazabilidad.',
        'Antes de presentar: confirmar nominales y topologia; validar la respuesta nominal frente a LTspice y justificar la alternativa segun la consigna. Si LTspice admite las funciones y no se autoriza la sustitucion, estas muestras pueden emplearse para automatizar las corridas exigidas en LTspice.']
    (out/'informe.md').write_text('\n\n'.join(paragraphs),encoding='utf-8')
    print(json.dumps({'nominal':payload['nominal'],'resultados':{k:v['cumplimiento'] for k,v in summary.items()}},indent=2))

if __name__=='__main__':main()
