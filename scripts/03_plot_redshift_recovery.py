import pandas as pd
import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# Constantes do modelo
# ==============================================================================
beta = 0.028          # impedância de fase do vácuo (fixa, sem ajuste)
c = 3e5               # velocidade da luz em km/s
arcsec_to_rad = np.pi / (180 * 3600)   # 1 arcsec = π/648000 rad

# ==============================================================================
# 1. Leitura do catálogo SLACS_table.cat (formato oficial)
# ==============================================================================
with open('SLACS_table.cat', 'r') as f:
    lines = f.readlines()

# Primeira linha contém o cabeçalho com '#'
header = lines[0].strip().lstrip('#').split()
data_lines = [line.strip().split() for line in lines[1:] if line.strip()]

df = pd.DataFrame(data_lines, columns=header)
for col in header:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Renomeia colunas para facilitar
df.rename(columns={
    'zd': 'z_L',
    'zs': 'z_S_obs',
    'theta_Ein': 'theta_E_arcsec',
    'veldisp(km/s)': 'sigma_v'
}, inplace=True)

# Converte theta_E de arcsec para radianos
df['theta_E_rad'] = df['theta_E_arcsec'] * arcsec_to_rad

# Remove linhas com dados faltantes nas colunas essenciais
df_clean = df.dropna(subset=['z_L', 'z_S_obs', 'theta_E_rad', 'sigma_v'])
print(f"Sistemas com dados completos: {len(df_clean)}")

# ==============================================================================
# 2. Função residual e solução da equação (4)
# ==============================================================================
def residual(zS, zL, theta_E_rad, sigma_v):
    """Equação (4) reescrita como LHS - RHS = 0"""
    eta_C = 1 + beta * np.log(c / sigma_v)
    eta_thick = eta_C * (1 + beta * (zS - zL))
    RHS = theta_E_rad / (4 * np.pi * (sigma_v / c)**2 * eta_thick)
    LHS = (zS - zL) / zS
    return LHS - RHS

def predict_zS(zL, theta_E_rad, sigma_v):
    """Encontra a raiz da equação no intervalo (zL+0.01, 10.0]"""
    try:
        sol = root_scalar(residual, args=(zL, theta_E_rad, sigma_v),
                          bracket=[zL + 0.01, 10.0], method='bisect', xtol=1e-6)
        return sol.root
    except:
        return np.nan

# ==============================================================================
# 3. Calcular predições para todos os sistemas
# ==============================================================================
predicoes = []
for _, row in df_clean.iterrows():
    z_pred = predict_zS(row['z_L'], row['theta_E_rad'], row['sigma_v'])
    predicoes.append(z_pred)

df_clean['z_S_pred'] = predicoes
df_valid = df_clean.dropna(subset=['z_S_pred'])
print(f"Sistemas com solução convergente: {len(df_valid)}")

# ==============================================================================
# 4. Métricas estatísticas
# ==============================================================================
residuos = df_valid['z_S_pred'] - df_valid['z_S_obs']
mae = np.mean(np.abs(residuos))
rmse = np.sqrt(np.mean(residuos**2))
corr = np.corrcoef(df_valid['z_L'], residuos)[0, 1]

print(f"\n--- Resultados finais (β = {beta}) ---")
print(f"MAE  = {mae:.4f}")
print(f"RMSE = {rmse:.4f}")
print(f"Correlação resíduos vs z_L = {corr:.4f}")

# ==============================================================================
# 5. Geração da figura (slacs_rigorous_physics_audit.png)
# ==============================================================================
plt.style.use('default')   # estilo limpo
fig, ax = plt.subplots(figsize=(8, 7))

# Pontos: predito vs observado
ax.scatter(df_valid['z_S_obs'], df_valid['z_S_pred'], 
           c='darkblue', s=40, alpha=0.7, edgecolors='k', linewidth=0.5, label='SLACS systems (N=59)')

# Linha de identidade (perfeita predição)
lim_max = max(df_valid['z_S_obs'].max(), df_valid['z_S_pred'].max()) + 0.1
ax.plot([0, lim_max], [0, lim_max], 'k--', lw=1.5, label='Perfect recovery (y=x)')

# Ajustes dos eixos
ax.set_xlim(0, lim_max)
ax.set_ylim(0, lim_max)
ax.set_xlabel(r'Observed Source Redshift $z_{S}^{\rm obs}$', fontsize=12)
ax.set_ylabel(r'Predicted Source Redshift $z_{S}^{\rm pred}$', fontsize=12)
ax.set_title(r'Blind Source Redshift Recovery (Parameter-Free, $\beta=0.028$, $N=0$)', fontsize=13)

# Caixa de texto com as estatísticas
stats_text = (f'N = {len(df_valid)}\n'
              f'MAE = {mae:.4f}\n'
              f'RMSE = {rmse:.4f}\n'
              f'Corr(res, z_L) = {corr:.4f}')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.70, 0.05, stats_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='bottom', bbox=props)

# Legenda
ax.legend(loc='upper left', fontsize=10)

# Grade leve
ax.grid(True, linestyle=':', alpha=0.5)

plt.tight_layout()
plt.savefig('slacs_rigorous_physics_audit.png', dpi=300, bbox_inches='tight')
print("\nFigura salva como 'slacs_rigorous_physics_audit.png'")
plt.show()