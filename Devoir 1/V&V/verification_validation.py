# MEC8211 - Devoir 1
# Bradley, James et Lucas

# Ce code vise à determiner la diffusion du sel a l'interieur d'un pilier d'un pont 
# plonge dans un etang saumatre
# La plupart des variables sont definies dans des fonctions appelees par ce module principal

import numpy as np
import matplotlib.pyplot as plt
import stationnaire as st
from math import *
import graphique_de_convergence as g

# Demande a l'utilisateur comment il veut discretiser son systeme en lui demandant de fournir un nombre de noeuds,
# et quel schema il souhaite utiliser pour approximer son systeme
char = input("voulez-vous avoir la reponse directement en regime stationnaire? Y/N: ")   
if char == "Y":
    rep = input("schema_1 (1) ou schema_2 (2) : ")
    if rep == "1":
        nombre = input("Choisissez le nombre de noeuds: ")
        nombre_de_noeuds = int(nombre)
        grille_de_sol = []
        inter = []
        r,C,ref = st.stationnaire_shema_1(nombre_de_noeuds) # reception des solutions numerique et analytique
        # Graphes des solutions
        fig,axe = plt.subplots()
        plt.plot(r, ref, linestyle='--', color='r', label='Solution analytique')
        plt.plot(r, C, color='b', label='Solution numerique')
        plt.title(f'Solutions en fonction de la taille du maillage $Δx$',
            fontsize=14, fontweight='bold', y=1.02)  # Le paramètre y règle la position verticale du titre
        plt.xlabel('Taille de maille $Δx$ (m)', fontsize=12, fontweight='bold')  # Remplacer "h" par "Δx"
        plt.ylabel('concentration de sel (mol/m3) ', fontsize=12, fontweight='bold')
        plt.legend()
        plt.show()

    #-------------------------------------------Erreur----------------------------------------------------
        # Graphes des erreurs L1,L2,L_∞
        g.graphes(rep)

    #------------------------------------------------------------------------------------------

    else:
        nombre = input("Choisissez le nombre de noeuds: ")
        nombre_de_noeuds = int(nombre)
        grille_de_sol = []
        inter = []
        r,C,ref = st.stationnaire_shema_2(nombre_de_noeuds)   # reception des solutions numerique et analytique
        # Graphes des solutions
        fig,axe = plt.subplots()
        plt.plot(r, ref, linestyle='--', marker='o', color='r', label='Solution analytique')
        plt.plot(r, C, color='b', label='Solution numerique')
        plt.title(f'Solutions en fonction de la taille du maillage $Δx$',
            fontsize=14, fontweight='bold', y=1.02)  # Le paramètre y règle la position verticale du titre
        plt.xlabel('Taille de maille $Δx$ (m)', fontsize=12, fontweight='bold')  # Remplacer "h" par "Δx"
        plt.ylabel('concentration de sel (mol/m3) ', fontsize=12, fontweight='bold')
        plt.legend()
        plt.show()

    #erreurs--------------------------------------------------------------------------------------------------
        # Graphes des erreurs L1,L2,L_∞
        g.graphes(rep)
        
    #---------------------------------------codes pour l'etat transitoire--------------------------------------------------------------------------- 
else:   
    rep = input("schema_1 (1) ou schema_2 (2) : ")
    if rep == "1":
        nombre = input("donnez le nombre de noeuds: ")
        nombre_de_noeuds = int(nombre)
        r,C,ref = st.transit_shema_1(nombre_de_noeuds)   # reception des solutions numerique et analytique
        # Graphes des solutions
        fig,axe = plt.subplots()
        plt.plot(r, ref, linestyle='--', color='r', label='Solution analytique')
        plt.plot(r, C, color='b', label='Solution numerique')
        plt.title(f'Solutions en fonction de la taille du maillage $Δx$',
            fontsize=14, fontweight='bold', y=1.02)  # Le paramètre y règle la position verticale du titre

        plt.xlabel('Taille de maille $Δx$ (m)', fontsize=12, fontweight='bold')  # Remplacer "h" par "Δx"
        plt.ylabel('concentration de sel (mol/m3) ', fontsize=12, fontweight='bold')
        plt.legend()
        plt.show()
    else:
        nombre = input("donnez le nombre de noeuds: ")
        nombre_de_noeuds = int(nombre)
        r,C,ref = st.transit_shema_2(nombre_de_noeuds)  # reception des solutions numerique et analytique
        fig,axe = plt.subplots()
        plt.plot(r, ref, linestyle='--',marker='o', color='r', label='Solution analytique')
        plt.plot(r, C, color='b', label='Solution numerique')
        plt.title(f'Solutions en fonction de la taille du maillage $Δx$',
            fontsize=14, fontweight='bold', y=1.02)  # Le paramètre y règle la position verticale du titre
        plt.xlabel('Taille de maille $Δx$ (m)', fontsize=12, fontweight='bold')  # Remplacer "h" par "Δx"
        plt.ylabel('concentration de sel (mol/m3) ', fontsize=12, fontweight='bold')
        plt.legend()
        plt.show()

    
         



