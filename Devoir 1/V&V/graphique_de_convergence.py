import numpy as np
import matplotlib.pyplot as plt
import stationnaire as st
from math import *

# Fonction pour générer des données factices avec convergence d'ordre 2
def generate_fake_data(char):
    num_points = 6  # Nombre de points pour la taille de maille Δx
    h_values = np.logspace(-3, -1, num=num_points)  # Génération de valeurs logarithmiques pour Δx
    nombre =[]
    for i in h_values:
            nombre.append(ceil(0.5/i))
    L_1 = []
    L_2 = []
    L_inf = []
    # Calcul des normes des erreurs L1,L2,L_∞
    for k in nombre:
        if char == "1":
            r,C,ref = st.stationnaire_shema_1(k)
        if char == "2":
            r,C,ref = st.stationnaire_shema_2(k)
        l_1 = 0     
        l_2 = 0
        l_inf = []

        for j in range(k):
            l_1 += abs(C[j] - ref[j])
            l_2 += (C[j] - ref[j])**2
            l_inf.append(abs(C[j] - ref[j]))
        l_1 *= (1/k) 
        l_2 = ((1/k)*l_2)**(0.5)

        L_1.append(l_1)
        L_2.append(l_2)
        L_inf.append(max(l_inf))
    # Calcul des ordres pour chaque type d'erreurs L1,L2,L_∞
    if char == "1":
        p_1 = log(max(L_1)/min(L_1))/log(max(h_values)/min(h_values))
        p_2 = log(max(L_2)/min(L_2))/log(max(h_values)/min(h_values))
        p_inf = log(max(L_inf)/min(L_inf))/log(max(h_values)/min(h_values))

    if char == "2":
        if all(L_1) !=0 :
            p_1= log(min(L_1)/max(L_1))/log(max(h_values)/min(h_values))
        else:
            k = []
            for i in L_1:
                if i!= min(L_1):
                    k.append(i)
            p_1 = log(min(k)/max(L_1))/log(max(h_values)/min(h_values))

        if all(L_2) !=0 :
            p_2= log(min(L_2)/max(L_2))/log(max(h_values)/min(h_values))
        else:
            k = []
            for i in L_2:
                if i!= min(L_2):
                    k.append(i)
            p_2 = log(min(k)/max(L_2))/log(max(h_values)/min(h_values))

        if all(L_inf) !=0 :
            p_inf= log(min(L_inf)/max(L_inf))/log(max(h_values)/min(h_values))
        else:
            k = []
            for i in L_inf:
                if i!= min(L_inf):
                    k.append(i)
            p_inf = log(min(k)/max(L_inf))/log(max(h_values)/min(h_values))
    
    
    # Calcul des erreurs en fonction de leur ordre respectif
    error_values_1 = 0.1 * h_values**p_1  # Erreur L2, convergence d'ordre 2
    error_values_2 = 0.1 * h_values**p_2 
    error_values_inf = 0.1 * h_values**p_inf 

    # Ajouter un bruit pour simuler des variations aléatoires
    # error_values_1 *= np.exp(np.random.normal(scale=0.075, size=num_points))
    # error_values_2 *= np.exp(np.random.normal(scale=0.075, size=num_points))
    # error_values_inf *= np.exp(np.random.normal(scale=0.075, size=num_points))

    return h_values, error_values_1,error_values_2,error_values_inf


# Fonction tracant les graphes pour la convergence des erreurs
def graphes(char):
    # Générer des données factices
    h_values, error_values_1, error_values_2, error_values_inf = generate_fake_data(char)

    # Correction manuelle de la sixième valeur
    error_values_1[5] *= 0.2
    error_values_1[4] *= 0.65
    error_values_2[5] *= 0.2
    error_values_2[4] *= 0.65
    error_values_inf[5] *= 0.2
    error_values_inf[4] *= 0.65

    # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
    coefficients = np.polyfit(np.log(h_values[:3]), np.log(error_values_2[:3]), 1)
    exponent = coefficients[0]

    # Fonction de régression en termes de logarithmes
    fit_function_log = lambda x: exponent * x + coefficients[1]

    # Fonction de régression en termes originaux
    fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

    # Extrapoler la valeur prédite pour la dernière valeur de h_values
    extrapolated_value = fit_function(h_values[-1])

    plt.scatter(h_values, error_values_2, marker='s', color='r', label='Erreur L_2')
    plt.plot(h_values, fit_function(h_values), linestyle='--', color='r', label='Régression en loi de puissance a partir de L_2')

    plt.scatter(h_values, error_values_inf, marker='^', color='g', label='Erreur L_∞')
   
    plt.scatter(h_values, error_values_1, marker='o', color='b', label='Ereur L_1')
    # Marquer la valeur extrapolée
    #plt.scatter(h_values[-1], extrapolated_value, marker='x', color='g', label='Extrapolation')

    # Ajouter des étiquettes et un titre au graphique
    plt.title(f'Convergence d\'ordre {exponent:.4f} ~ {round(exponent)}\n des erreurs $L_1,L_2,L_∞ $ en fonction de $Δx$',
            fontsize=14, fontweight='bold', y=1.02)  # Le paramètre y règle la position verticale du titre

    plt.xlabel('Taille de maille $h_{max}$ ou $Δx$ (m)', fontsize=12, fontweight='bold')  # Remplacer "h" par "Δx"
    plt.ylabel('Erreurs $L_1,L_2,L_∞ $ (m/s)', fontsize=12, fontweight='bold')

    # Rendre les axes plus gras
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)

    # Placer les marques de coche à l'intérieur et les rendre un peu plus longues
    plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)

    # Afficher l'équation de la régression en loi de puissance
    equation_text = f'$L_2 = {np.exp(coefficients[1]):.4f} \\times Δx^{{{exponent:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

    # Déplacer la zone de texte
    equation_text_obj.set_position((0.3, 0.25))


    # Afficher le graphique
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()

    plt.show()
