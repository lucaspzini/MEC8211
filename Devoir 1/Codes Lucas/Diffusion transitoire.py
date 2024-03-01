# MEC8211 - Devoir 1
# Bradley, James et Lucas

# Ce code vise à calculer l'équation de diffusion en régime transitoire

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm

# Définition des variables
R = 0.5         # Rayon du cylindre
Deff = 1e-10    # Coefficient de diffusion
S = 8e-9        # Consommation
Ce = 12         # Concentration à r = R
steps = 10000   # Nombre de pas de temps
dt = 1e6        # Pas de temps
print("dt =", dt)
print()

# Fonction pour résoudre le système linéaire : approximation avant de la dérivée
def solve_diffusion_cylindrical_forward(r, h, N, Deff, S, Ce, step, dt):
    # Construction des matrices (A*c = b)
    A = []
    b = []
    C = [0 for i in range(N)]   # Concentration à t = 0
    
    for _ in range(step):
        for i in range(N):
            vecteur = []
            if i == 0:
                for j in range(N):
                    if j == 0:
                        vecteur.append(-3)
                    elif j == 1:
                        vecteur.append(4)
                    elif j == 2:
                        vecteur.append(-1)
                    else:
                        vecteur.append(0)
                A.append(vecteur)
                b.append(0)
            elif i == N-1:
                for j in range(N):
                    if j == N-1:
                        vecteur.append(h**2 + Deff*dt*(2 + h/r[i]))
                    elif j == N-2:
                        vecteur.append(-Deff*dt)
                    else:
                        vecteur.append(0)
                A.append(vecteur)
                b.append(Ce*h**2 - S*h**2*dt - (-Deff*dt*(1 + h/r[i]))*Ce)
            else:
                for j in range(N):
                    if j == i-1:
                        vecteur.append(-Deff*dt)
                    elif j == i:
                        vecteur.append(h**2 + Deff*dt*(2 + h/r[i]))
                    elif j == i+1:
                        vecteur.append(-Deff*dt*(1 + h/r[i]))
                    else:
                        vecteur.append(0)
                A.append(vecteur)
                b.append(C[i]*h**2 - S*dt*h**2)
    
        C = np.linalg.solve(A,b)
        A = []
        b = []
    
    return C

# Fonction pour résoudre le système linéaire : approximation centrée de la dérivée
def solve_diffusion_cylindrical_central(r, h, N, Deff, S, Ce, step, dt):
    # Construction des matrices (A*c = b)
    A = []
    b = []
    C = [0 for i in range(N)]   # Concentration à t = 0
    
    for _ in range(step):
        for i in range(N):
            vecteur = []
            if i == 0:
                for j in range(N):
                    if j == 0:
                        vecteur.append(-3)
                    elif j == 1:
                        vecteur.append(4)
                    elif j == 2:
                        vecteur.append(-1)
                    else:
                        vecteur.append(0)
                A.append(vecteur)
                b.append(0)
            elif i == N-1:
                for j in range(N):
                    if j == N-1:
                        vecteur.append(2*Deff*dt + h**2)
                    elif j == N-2:
                        vecteur.append(-Deff*dt*(1 - h/(2*r[i])))
                    else:
                        vecteur.append(0)
                A.append(vecteur)
                b.append(Ce*h**2 - S*h**2*dt - (-Deff*dt*(1 + h/(2*r[i]))))
            else:
                for j in range(N):
                    if j == i-1:
                        vecteur.append(-Deff*dt*(1 - h/(2*r[i])))
                    elif j == i:
                        vecteur.append(2*Deff*dt + h**2)
                    elif j == i+1:
                        vecteur.append(-Deff*dt*(1 + h/(2*r[i])))
                    else:
                        vecteur.append(0)
                A.append(vecteur)
                b.append(C[i]*h**2 - S*dt*h**2)
    
        C = np.linalg.solve(A,b)
        A = []
        b = []
    
    return C

# Résoudre l'équation de diffusion

# Initialiser les variables de stockage pour le traçage
all_r = []
all_C = []
all_h = []
all_l1 = []
all_l2 = []
all_l_inf = []

# Définition du maillage
noeuds = [5, 10, 20, 50, 100]

for N in noeuds:
    h = R/(N-1)     # Pas du schéma
    r = np.linspace(0, R, N)    # Discrétisation du domaine
    print("N =", N)
    print()
    print("h =", h)

    # Solution analytique
    sol_anly = []

    for i in range(N):
        sol_anly.append(0.25*(S/Deff)*(R**2)*((r[i]/R)**2 - 1) + Ce)
    
    # Choisissez l'une des méthodes pour calculer la concentration
    # C = solve_diffusion_cylindrical_forward(r, h, N, Deff, S, Ce, steps, dt)
    C = solve_diffusion_cylindrical_central(r, h, N, Deff, S, Ce, steps, dt)

    print()
    print("C =", C)
    
    # Calculer les erreurs (normes standardisées)
    # Il y a deux façons de calculer les normes
    erreur = abs(C - sol_anly)
    
    # l1 = 1/N * norm(erreur, 1)
    l1 = 1/N * sum(erreur)
    print()
    print("L1 =", l1)
    
    # l2 = 1/np.sqrt(N) * norm(erreur)
    l2 = np.sqrt(1/N * sum(erreur**2))
    print()
    print("L2 =", l2)
    
    # l_inf = norm(erreur, np.inf)
    l_inf = max(erreur)
    print()
    print("L∞ =", l_inf)
    
    # Stocker les données pour le traçage
    all_r.append(r)
    all_C.append(C)
    all_h.append(h)
    all_l1.append(l1)
    all_l2.append(l2)
    all_l_inf.append(l_inf)
    
    print("\n#######\n")

# Tracer le résultat
plt.figure()
for i in range(len(all_r)):
    plt.plot(all_r[i], all_C[i], label=f'N = {noeuds[i]}')

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
coefficients = np.polyfit(np.log(all_h), np.log(all_l2), 1)
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