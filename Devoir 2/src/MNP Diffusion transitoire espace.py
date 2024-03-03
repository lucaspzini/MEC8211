# MEC8211 - Devoir 2
# Lucas P. Zini

# Ce code vise à calculer l'équation de diffusion en régime transitoire

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.interpolate import interp1d

# Définition des variables
R = 0.5         # Rayon du cylindre
Deff = 1e-10    # Coefficient de diffusion
k = 4e-9        # Constante de vitesse
Ce = 12         # Concentration à r = R
dt = 1e6        # Pas de temps
print("dt =", dt)
print()

# Fonction pour résoudre le système linéaire : approximation centrée de la dérivée
def solve_diffusion_cylindrical_central(r, h, N, Deff, k, Ce, dt):
    C = np.zeros(N)   # Concentration à t = 0
    C_aux = C
    
    residuMax = 1e-10
    residu = 10.*residuMax

    while residu > residuMax:
        A = np.zeros((N, N))
        b = np.zeros(N)

        for i in range(1, N-1):
            A[i, i-1] = -Deff*dt*(1 - h/(2*r[i]))
            A[i, i] = 2*Deff*dt + h**2 + h**2*dt*k
            A[i, i+1] = -Deff*dt*(1 + h/(2*r[i]))
            b[i] = C[i]*h**2
        A[0, 0] = -3
        A[0, 1] = 4
        A[0, 2] = -1
        b[0] = 0
        A[N-1, N-1] = 1
        b[N-1] = Ce
        C = np.linalg.solve(A, b)

        residu = abs(sum(C - C_aux))
        C_aux = C
    return C

# Obtenir une spline d'interpolation
def interpolation(r, C, kind='cubic'):
    return interp1d(r, C, kind)

# Résoudre l'équation de diffusion

# Initialiser les variables de stockage pour le traçage
all_r = []
all_C = []
all_h = []
all_l1 = []
all_l2 = []
all_l_inf = []

# Définition du maillage pour la solution analytique
N_anly = 100
h = R/(N_anly-1)     # Pas du schéma
r_anly = np.linspace(0, R, N_anly)    # Discrétisation du domaine

# Calculer la concentration
C_anly = solve_diffusion_cylindrical_central(r_anly, h, N_anly, Deff, k, Ce, dt)
print('C_analytique =', C_anly)

# Fonction d'interpolation analytique - spline
spline = interpolation(r_anly, C_anly)

# Définition du maillage
noeuds = [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50]

for N in noeuds:
    h = R/(N-1)     # Pas du schéma
    r = np.linspace(0, R, N)    # Discrétisation du domaine
    print("N =", N)
    print()
    print("h =", h)
    
    # Solution analytique
    sol_anly = []
    for i in range(N):
        sol_anly.append(spline(r[i]))

    # Calculer la concentration avec la MNP
    C = solve_diffusion_cylindrical_central(r, h, N, Deff, k, Ce, dt)
    
    # Stocker les données pour le traçage
    all_r.append(r)
    all_C.append(C)
    all_h.append(h)

    # Calculer les erreurs (normes standardisées)
    erreur = abs(C - sol_anly)
    
    l1 = 1/N * sum(erreur)
    print()
    print("L1 =", l1)
    
    l2 = np.sqrt(1/N * sum(erreur**2))
    print()
    print("L2 =", l2)
    
    l_inf = max(erreur)
    print()
    print("L∞ =", l_inf)

    # Stocker les données pour le traçage
    all_l1.append(l1)
    all_l2.append(l2)
    all_l_inf.append(l_inf)

    print("\n#######\n")

# Tracer le résultat
plt.figure()
plt.plot(r_anly, spline(r_anly),label='Solution proche')
for i in range(len(all_r)):
    plt.plot(all_r[i], all_C[i], label=f'Solution avec {noeuds[i]} nœuds')

# Rendre les axes plus gras
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['left'].set_linewidth(2)
plt.gca().spines['right'].set_linewidth(2)
plt.gca().spines['top'].set_linewidth(2)

# Marques de coche plus longues
plt.tick_params(width=2, length=6)

plt.xlabel('r (m)', fontsize=12, fontweight='bold')
plt.ylabel('Concentration (mol/m³)', fontsize=12, fontweight='bold')
plt.xlim(0, 0.5)
plt.legend()
plt.show()


# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = []
coefficients = np.polyfit(np.log(all_h[6:11]), np.log(all_l2[6:11]), 1)
exponent = coefficients[0]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent * x + coefficients[1]

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Tracer L1, L2 et L∞ en fonction de h
plt.figure(figsize=(8, 6))

plt.loglog(all_h, all_l1, marker='o', linestyle='', color="b", label='L1')
plt.loglog(all_h, all_l2, marker='s', linestyle='', color="k", label='L2')
plt.loglog(all_h, all_l_inf, marker='^', linestyle='', color="g", label='L∞')
plt.loglog(all_h, fit_function(all_h), linestyle='--', color='r', label='Régression en loi de puissance')

# Rendre les axes plus gras
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['left'].set_linewidth(2)
plt.gca().spines['right'].set_linewidth(2)
plt.gca().spines['top'].set_linewidth(2)

# Placer les marques de coche à l'intérieur et les rendre un peu plus longues
plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)

# Afficher l'équation de la régression en loi de puissance
equation_text = f'$L_2 = {np.exp(coefficients[1]):.4f} \\times Δr^{{{exponent:.4f}}}$'
equation_text_obj = plt.text(0.2, 0.2, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel('h ou Δr (m)', fontsize=12, fontweight='bold')
plt.ylabel('Erreur (mol/m³)', fontsize=12, fontweight='bold')
plt.legend()
plt.grid(True)
plt.show()