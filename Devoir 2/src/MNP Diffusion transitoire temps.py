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
T = 1e6         # Temps stationnaire

# Fonction pour résoudre le système linéaire : approximation centrée de la dérivée
def solve_diffusion_cylindrical_central(r, h, N, Deff, k, Ce, dt):
    C = np.zeros(N)   # Concentration à t = 0
    C_aux = C
    
    residuMax = 1e-8
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
all_dt = []
all_l1 = []
all_l2 = []
all_l_inf = []

# Définition du maillage pour la solution analytique
N = 50          # Nombre de nœuds en espace
h = R/(N-1)     # Pas du schéma en espace
r = np.linspace(0, R, N)    # Discrétisation du domaine spatial
Nt = 1001       # Nombre de nœuds en temps
dt = T/(Nt-1)   # Pas du schéma en temps
t = np.linspace(0, T, Nt)   # Discrétisation du domaine temporel

# Calculer la concentration sur le domaine 2D
C_anly = np.zeros((Nt, N))
for i in range(Nt):
    C_anly[i] = solve_diffusion_cylindrical_central(r, h, N, Deff, k, Ce, (i*dt))
    # print(i)

# Fonction d'interpolation analytique - une spline pour chaque collone dans l'espace
splines = []
for i in range(N):
    spline = interpolation(t, C_anly[:,i])
    splines.append(spline)

# Définition du maillage
noeuds_t = [11, 51, 101, 501]

for Nt in noeuds_t:
    dt = T/(Nt-1)     # Pas du schéma en temps
    t = np.linspace(0, T, Nt)
    print("Nœuds en temps =", Nt)
    print()
    print("dt =", dt)
    
    # Stocker les données pour le traçage
    all_dt.append(dt)

    # Calculer la concentration sur le domaine 2D
    C = np.zeros((Nt, N))
    for i in range(Nt):
        C[i] = solve_diffusion_cylindrical_central(r, h, N, Deff, k, Ce, (i*dt))
        # print(i)

    # Calculer l'erreur pour chaque collone de la matrice C
    erreurs = []
    for i in range(N):
        C_spline = splines[i](t)
        erreur = abs(C[:,i] - C_spline)
        erreurs.append(erreur)
    erreurs_array = np.array(erreurs) # Convertir la variable erreurs en NumPy array
    
    # Calculer les normes L1, L2 et L∞ pour chaque pas de temps
    l1 = 1/Nt * 1/N * np.sum(erreurs_array)
    print()
    print("L1 =", l1)

    l2 = np.sqrt(1/Nt * 1/N * np.sum(erreurs_array**2))
    print()
    print("L2 =", l2)

    l_inf = np.max(erreurs_array)
    print()
    print("L∞ =", l_inf)

    # Stocker les données pour le traçage
    all_l1.append(l1)
    all_l2.append(l2)
    all_l_inf.append(l_inf)

    print("\n#######\n")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(all_dt), np.log(all_l2), 1)
exponent = coefficients[0]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent * x + coefficients[1]

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Tracer L1, L2 et L∞ en fonction de h
plt.figure(figsize=(8, 6))

plt.loglog(all_dt, all_l1, marker='o', linestyle='', color="b", label='L1')
plt.loglog(all_dt, all_l2, marker='s', linestyle='', color="k", label='L2')
plt.loglog(all_dt, all_l_inf, marker='^', linestyle='', color="g", label='L∞')
# plt.loglog(all_dt, fit_function(all_dt), linestyle='--', color='r', label='Régression en loi de puissance')

# Rendre les axes plus gras
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['left'].set_linewidth(2)
plt.gca().spines['right'].set_linewidth(2)
plt.gca().spines['top'].set_linewidth(2)

# Placer les marques de coche à l'intérieur et les rendre un peu plus longues
plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)

# Afficher l'équation de la régression en loi de puissance
# equation_text = f'$L_2 = {np.exp(coefficients[1]):.4f} \\times Δt^{{{exponent:.4f}}}$'
# equation_text_obj = plt.text(0.2, 0.2, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel('dt (s)', fontsize=12, fontweight='bold')
plt.ylabel('Erreur (mol/m³)', fontsize=12, fontweight='bold')
plt.legend()
plt.grid(True)
plt.show()