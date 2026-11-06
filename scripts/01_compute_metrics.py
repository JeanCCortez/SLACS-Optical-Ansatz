import pandas as pd
import numpy as np
from scipy.optimize import root_scalar

beta = 0.028
c = 3e5  # km/s
arcsec_to_rad = np.pi / (180 * 3600)  # 1 arcsec = π/648000 rad

# Leitura correta do arquivo
with open('SLACS_table.cat', 'r') as f:
    lines = f.readlines()
header = lines[0].strip().lstrip('#').split()
data_lines = [line.strip().split() for line in lines[1:] if line.strip()]
df = pd.DataFrame(data_lines, columns=header)
for col in header:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Renomeia
df.rename(columns={
    'zd': 'z_L',
    'zs': 'z_S_obs',
    'theta_Ein': 'theta_E_arcsec',
    'veldisp(km/s)': 'sigma_v'
}, inplace=True)

# Converte theta_E para radianos
df['theta_E_rad'] = df['theta_E_arcsec'] * arcsec_to_rad

# Remove NaN
df_clean = df.dropna(subset=['z_L', 'z_S_obs', 'theta_E_rad', 'sigma_v'])
print(f"Sistemas com dados completos: {len(df_clean)}")
print(f"Exemplo: z_L={df_clean.iloc[0]['z_L']}, sigma_v={df_clean.iloc[0]['sigma_v']}, theta_E_rad={df_clean.iloc[0]['theta_E_rad']:.6f}")

# Função residual correta
def residual(zS, zL, theta_E_rad, sigma_v):
    eta_C = 1 + beta * np.log(c / sigma_v)
    eta_thick = eta_C * (1 + beta * (zS - zL))
    RHS = theta_E_rad / (4 * np.pi * (sigma_v / c)**2 * eta_thick)
    LHS = (zS - zL) / zS
    return LHS - RHS

# Testa o primeiro sistema
row = df_clean.iloc[0]
zL, zS_obs, theta_E_rad, sigma_v = row['z_L'], row['z_S_obs'], row['theta_E_rad'], row['sigma_v']
print(f"\nTeste com zL={zL}, sigma_v={sigma_v}, theta_E_rad={theta_E_rad:.6f}")
for zS_teste in np.linspace(zL+0.01, 5.0, 10):
    res = residual(zS_teste, zL, theta_E_rad, sigma_v)
    print(f"zS={zS_teste:.2f}, residual={res:.6e}")

# Resolve para todos
predicoes = []
for idx, row in df_clean.iterrows():
    zL = row['z_L']
    theta_E_rad = row['theta_E_rad']
    sigma_v = row['sigma_v']
    try:
        sol = root_scalar(residual, args=(zL, theta_E_rad, sigma_v),
                          bracket=[zL+0.01, 10.0], method='bisect', xtol=1e-6)
        predicoes.append(sol.root)
    except:
        predicoes.append(np.nan)

df_clean['z_S_pred'] = predicoes
df_valid = df_clean.dropna(subset=['z_S_pred'])
print(f"\nSistemas com convergência: {len(df_valid)}")

if len(df_valid) > 0:
    mae_modelo = np.mean(np.abs(df_valid['z_S_pred'] - df_valid['z_S_obs']))
    mae_nulo = np.mean(np.abs(2*df_valid['z_L'] - df_valid['z_S_obs']))
    corr = df_valid['z_L'].corr(df_valid['z_S_pred'] - df_valid['z_S_obs'])
    print(f"MAE modelo nulo: {mae_nulo:.4f}")
    print(f"MAE modelo proposto: {mae_modelo:.4f}")
    print(f"Correlação resíduos vs z_L: {corr:.4f}")
    df_valid.to_csv('resultados_finais.csv', index=False)
else:
    print("Ainda sem convergência. Verifique se os resíduos mudam de sinal.")