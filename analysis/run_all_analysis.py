"""
=============================================================================
Binary Bio-adsorbent System for Azo Dye Removal: RH + STL
ICESSM-2026 — Full Reproducible Analysis Script
=============================================================================
Runs all modelling steps:
  1. Kinetic analysis  (PFO, PSO, Weber-Morris)
  2. Isotherm analysis (Langmuir, Freundlich, Temkin, Dubinin-Radushkevich)
  3. Binary synergy    (Synergy Factor)
  4. ANN optimisation  (dose minimisation for RE >= 95%)

Requirements: numpy, scipy, scikit-learn, matplotlib, pandas
Install     : pip install -r requirements.txt
Usage       : python run_all_analysis.py
=============================================================================
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, brentq
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import r2_score
import warnings, os
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, '..', 'data')
OUT  = os.path.join(BASE, '..', 'results')
os.makedirs(OUT, exist_ok=True)

rh_ct  = pd.read_csv(f'{DATA}/01_RH_contact_time.csv')
stl_ct = pd.read_csv(f'{DATA}/02_STL_contact_time.csv')
rh_d   = pd.read_csv(f'{DATA}/03_RH_dose_variation.csv')
stl_d  = pd.read_csv(f'{DATA}/04_STL_dose_variation.csv')
bin_rh = pd.read_csv(f'{DATA}/05_binary_RH_constant.csv')
bin_stl= pd.read_csv(f'{DATA}/06_binary_STL_constant.csv')

def pfo_model(t, qe, k1):  return qe * (1 - np.exp(-k1 * t))
def pso_model(t, qe, k2):  return k2 * qe**2 * t / (1 + k2 * qe * t)

results = {}

for label, df in [('RH', rh_ct), ('STL', stl_ct)]:
    t  = df['t_min'].values[1:].astype(float)
    qt = df['qt_mg_g'].values[1:].astype(float)
    qe_exp = qt[-1]
    p_pfo, _ = curve_fit(pfo_model, t, qt, p0=[12, 0.05], maxfev=10000)
    p_pso, _ = curve_fit(pso_model, t, qt, p0=[12, 0.008], maxfev=10000)
    r2_pfo = r2_score(qt, pfo_model(t, *p_pfo))
    r2_pso = r2_score(qt, pso_model(t, *p_pso))
    tqt = t / qt
    m, b = np.polyfit(t, tqt, 1)
    qe_lin = 1 / m
    k2_lin = m**2 / b
    r2_lin = r2_score(tqt, np.polyval([m, b], t))
    sq_t = np.sqrt(t)
    m1, c1 = np.polyfit(sq_t[t <= 60], qt[t <= 60], 1)
    m2, c2 = np.polyfit(sq_t[t >= 60], qt[t >= 60], 1)
    results[f'{label}_kinetics'] = {
        'qe_exp': qe_exp,
        'PFO_qe': p_pfo[0], 'PFO_k1': p_pfo[1], 'PFO_R2': r2_pfo,
        'PSO_qe': p_pso[0], 'PSO_k2': p_pso[1], 'PSO_R2': r2_pso,
        'PSO_lin_qe': qe_lin, 'PSO_lin_k2': k2_lin, 'PSO_lin_R2': r2_lin,
        'WM_kid1': m1, 'WM_C1': c1, 'WM_kid2': m2, 'WM_C2': c2,
    }
    print(f"\n{label} Kinetics: PFO R2={r2_pfo:.4f}, PSO R2={r2_pso:.4f}")

def langmuir(Ce, qm, KL): return qm * KL * Ce / (1 + KL * Ce)
def freundlich(Ce, KF, n): return KF * Ce**(1/n)
def temkin(Ce, AT, bT):
    B = 8.314 * 302 / bT
    return B * np.log(AT * Ce)
def dr_epsilon(Ce_mg_L, MW=327.33):
    Ce_mol = Ce_mg_L / (MW * 1000)
    return 8.314 * 302 * np.log(1 + 1 / Ce_mol)

for label, df in [('RH', rh_d), ('STL', stl_d)]:
    Ce = df['Ce_mg_L'].values.astype(float)
    qe = df['qe_mg_g'].values.astype(float)
    idx = np.argsort(Ce); Ce, qe = Ce[idx], qe[idx]
    pL, _ = curve_fit(langmuir,  Ce, qe, p0=[100, 0.01], maxfev=20000, bounds=([0,0],[100000,100]))
    pF, _ = curve_fit(freundlich, Ce, qe, p0=[3, 1.5],   maxfev=20000, bounds=([0,0.01],[500,20]))
    pT, _ = curve_fit(temkin,     Ce, qe, p0=[1, 200],   maxfev=20000, bounds=([1e-6,1],[1000,100000]))
    eps = dr_epsilon(Ce)
    ln_qe = np.log(qe / (327.33 * 1000))
    b_dr, a_dr = np.polyfit(eps**2, ln_qe, 1)
    E_kJ = 1 / np.sqrt(-2 * b_dr) / 1000
    RL = 1 / (1 + pL[1] * 50)
    nature = ('Physisorption' if E_kJ < 8 else 'Ion exchange' if E_kJ < 16 else 'Chemisorption')
    results[f'{label}_isotherms'] = {
        'Langmuir_qm': pL[0], 'Langmuir_KL': pL[1], 'RL': RL, 'Langmuir_R2': r2_score(qe, langmuir(Ce, *pL)),
        'Freundlich_KF': pF[0], 'Freundlich_n': pF[1], 'Freundlich_R2': r2_score(qe, freundlich(Ce, *pF)),
        'Temkin_AT': pT[0], 'Temkin_bT': pT[1], 'Temkin_R2': r2_score(qe, temkin(Ce, *pT)),
        'DR_E_kJ_mol': E_kJ, 'DR_nature': nature,
    }

RE_RH  = rh_d[rh_d['dose_g'] == 1.5]['RE_pct'].values[0]
RE_STL = stl_d[stl_d['dose_g'] == 1.5]['RE_pct'].values[0]
RE_BIN = 95.1
Sf = RE_BIN / (0.5 * RE_RH + 0.5 * RE_STL)
results['synergy'] = {'RE_RH_1p5g': RE_RH, 'RE_STL_1p5g': RE_STL, 'RE_binary_equal': RE_BIN, 'Sf': Sf}
print(f"\nSynergy Factor: Sf = {Sf:.4f}")

X_all = np.vstack([
    np.c_[rh_d['dose_g'].values, np.zeros(16)],
    np.c_[np.zeros(16), stl_d['dose_g'].values],
    np.c_[np.full(15, 1.5), bin_rh['STL_g'].values],
    np.c_[bin_stl['RH_g'].values, np.full(15, 1.5)],
])
y_all = np.concatenate([rh_d['RE_pct'].values, stl_d['RE_pct'].values,
                         bin_rh['RE_pct'].values, bin_stl['RE_pct'].values])
sX, sy = MinMaxScaler(), MinMaxScaler()
Xsc = sX.fit_transform(X_all)
ysc = sy.fit_transform(y_all.reshape(-1,1)).ravel()
ann = MLPRegressor(hidden_layer_sizes=(16,16), activation='relu', solver='lbfgs',
                   max_iter=5000, random_state=42, alpha=0.01)
ann.fit(Xsc, ysc)
y_pred = sy.inverse_transform(ann.predict(Xsc).reshape(-1,1)).ravel()
r2_ann = r2_score(y_all, y_pred)
rmse_ann = np.sqrt(np.mean((y_all - y_pred)**2))
rg = np.linspace(0, 1.8, 100)
sg = np.linspace(0, 1.8, 100)
RG, SG = np.meshgrid(rg, sg)
RE_G = np.clip(sy.inverse_transform(ann.predict(sX.transform(np.c_[RG.ravel(), SG.ravel()])).reshape(-1,1)).ravel(), 0, 100).reshape(RG.shape)
opt_pts = [(r, s, RE_G[j,i]) for i,r in enumerate(rg) for j,s in enumerate(sg) if RE_G[j,i] >= 95.0]
opt_pts.sort(key=lambda x: x[0]+x[1])
opt_rh, opt_stl, opt_re = opt_pts[0]
results['ANN'] = {'R2': r2_ann, 'RMSE': rmse_ann, 'opt_RH_g': opt_rh, 'opt_STL_g': opt_stl, 'opt_RE_pct': opt_re}
print(f"ANN: R2={r2_ann:.4f}, Optimum RH={opt_rh:.2f}g + STL={opt_stl:.2f}g -> RE={opt_re:.1f}%")

rows = [{'Section': s, 'Parameter': k, 'Value': v} for s, vals in results.items() for k, v in vals.items()]
pd.DataFrame(rows).to_csv(f'{OUT}/parameters_summary.csv', index=False)
print(f"\nResults saved to {OUT}/parameters_summary.csv")
print("All analysis complete.")
