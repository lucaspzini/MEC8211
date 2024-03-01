# Module contenant les fonctions qui retournent les solutions analytique et
# numerique soit directement en stationnaire ou en transitoire

import numpy as np
import matplotlib.pyplot as plt


R = 1/2          # Rayon du pilier
Deff = 10**(-10)  # coefficient de diffusion effectif du sel dans le béton poreux
S = 8*10**(-9)    # terme source
Ce = 12           # concentration en sel constante à la surface du pilier

#----------------------------------schema_1-----------------------------------------------

# Fonction qui retourne les solutions analytique et
# numerique directement a partir d'un solveur stationnaire en schema 1

def stationnaire_shema_1(nombre_de_noeuds):

    delta_r = R/(nombre_de_noeuds-1)
    r = np.linspace(0,0.5,nombre_de_noeuds)
    matrix = []
    b = []
    C=[]
    ref = []
    # solution analytique
    for i in range(nombre_de_noeuds):
        ref.append(0.25*(S/Deff)*(R**2)*(((r[i]/R)**2)-1) + 12 )
    # solution numerique
    for i in range(nombre_de_noeuds):
        sous_liste = []
        if i == 0:
            for j in range(nombre_de_noeuds):
                
                if j == 0:
                    sous_liste.append(-3)
                elif j == 1:
                    sous_liste.append(4)
                elif j == 2:
                    sous_liste.append(-1)
                else:
                    sous_liste.append(0)
            matrix.append(sous_liste)
            b.append(0)
        elif i == nombre_de_noeuds-1:
            sous_liste = [1 if j == nombre_de_noeuds-1 else 0 for j in range(nombre_de_noeuds)]
            matrix.append(sous_liste)
            b.append(12)
        else:
            for j in range(nombre_de_noeuds):
                if j == i-1:
                    sous_liste.append(-Deff)
                elif j == i:
                    sous_liste.append(2*Deff + delta_r*Deff/(r[i]))
                elif j == i+1:
                    sous_liste.append(-(Deff + delta_r*Deff/(r[i])))
                else:
                    sous_liste.append(0)
            matrix.append(sous_liste)
            b.append(-(delta_r**2)*S)
    C = np.linalg.solve(matrix,b)
    return r,C,ref

#----------------------------------schema_2-----------------------------------------------

# Fonction qui retourne les solutions analytique et
# numerique a partir d'un solveur stationnaire en schema 2

def stationnaire_shema_2(nombre_de_noeuds):
    delta_r = R/(nombre_de_noeuds-1)
    r = np.linspace(0,0.5,nombre_de_noeuds)
    matrix = []
    b = []
    C=[]
    ref = []
    # Solution analytique
    for i in range(nombre_de_noeuds):
        ref.append(0.25*(S/Deff)*(R**2)*(((r[i]/R)**2)-1) + 12 )
    # Solution numerique
    for i in range(nombre_de_noeuds):
        sous_liste = []
        if i == 0:
            for j in range(nombre_de_noeuds):
                
                if j == 0:
                    sous_liste.append(-3)
                elif j == 1:
                    sous_liste.append(4)
                elif j == 2:
                    sous_liste.append(-1)
                else:
                    sous_liste.append(0)
            matrix.append(sous_liste)
            b.append(0)
        elif i == nombre_de_noeuds-1:
            sous_liste = [1 if j == nombre_de_noeuds-1 else 0 for j in range(nombre_de_noeuds)]
            matrix.append(sous_liste)
            b.append(12)
        else:
            for j in range(nombre_de_noeuds):
                if j == i-1:
                    sous_liste.append((delta_r*Deff)/(2*r[i]) -Deff)
                elif j == i:
                    sous_liste.append(2*Deff)
                elif j == i+1:
                    sous_liste.append(-(Deff + (delta_r*Deff)/(2*r[i])))
                else:
                    sous_liste.append(0)
            matrix.append(sous_liste)
            b.append(-(delta_r**2)*S)
    C = np.linalg.solve(matrix,b)
    return r,C,ref

#------------------------------Code transitoire--------------------------------------------------------------
#------------------------------Schema 1--------------------------------------------------------------

# Fonction qui retourne les solutions analytique et
# numerique a partir d'un solveur transitoire en schema 1

def transit_shema_1(nombre_de_noeuds):

    delta_r = R/(nombre_de_noeuds-1)
    r = np.linspace(0,0.5,nombre_de_noeuds)
    matrix = []
    b = []
    C=[0 for i in range(nombre_de_noeuds)]
    ref = []
    # Solution analytique
    for i in range(nombre_de_noeuds):
        ref.append(0.25*(S/Deff)*(R**2)*(((r[i]/R)**2)-1) + 12 )
    # Solution numerique
    t= 0
    t_station = 20000000000
    delta_t = 2000
    for i in range(nombre_de_noeuds):
            sous_liste = []
            if i == 0:
                for j in range(nombre_de_noeuds):
                    
                    if j == 0:
                        sous_liste.append(-3)
                    elif j == 1:
                        sous_liste.append(4)
                    elif j == 2:
                        sous_liste.append(-1)
                    else:
                        sous_liste.append(0)
                matrix.append(sous_liste)
            elif i == nombre_de_noeuds-1:
                sous_liste = [1 if j == nombre_de_noeuds-1 else 0 for j in range(nombre_de_noeuds)]
                matrix.append(sous_liste)      
            else:
                for j in range(nombre_de_noeuds):
                    if j == i-1:
                        sous_liste.append(-delta_t*Deff)
                    elif j == i:
                        sous_liste.append((delta_r)**2 + 2*delta_t*Deff + delta_t*delta_r*Deff/(r[i]))
                    elif j == i+1:
                        sous_liste.append(-(delta_t*Deff + delta_t*delta_r*Deff/(r[i])))
                    else:
                        sous_liste.append(0)
                matrix.append(sous_liste)
                
    while(t < t_station):    
        if t==0:
            t+=delta_t
            continue
        b.append(0)
        for i in range(nombre_de_noeuds-1):
            if i ==0:
                continue
            b.append(C[i]*(delta_r)**2 - delta_t*(delta_r**2)*S)
        b.append(12)
        C = np.linalg.solve(matrix,b) 
        t= t + delta_t
        b = []

    return r,C,ref

#------------------------------Schema 2--------------------------------------------------------------

# Fonction qui retourne les solutions analytique et
# numerique a partir d'un solveur transitoire en schema 2

def transit_shema_2(nombre_de_noeuds):

    delta_r = R/(nombre_de_noeuds-1)
    r = np.linspace(0,0.5,nombre_de_noeuds)
    matrix = []
    b = []
    C=[0 for i in range(nombre_de_noeuds)]
    ref = []
    # Solution analytique
    for i in range(nombre_de_noeuds):
        ref.append(0.25*(S/Deff)*(R**2)*(((r[i]/R)**2)-1) + 12 )
    # Solution numerique
    t= 0
    t_station = 20000000000
    delta_t = 2000
    for i in range(nombre_de_noeuds):
            sous_liste = []
            if i == 0:
                for j in range(nombre_de_noeuds):
                    
                    if j == 0:
                        sous_liste.append(-3)
                    elif j == 1:
                        sous_liste.append(4)
                    elif j == 2:
                        sous_liste.append(-1)
                    else:
                        sous_liste.append(0)
                matrix.append(sous_liste)
            elif i == nombre_de_noeuds-1:
                sous_liste = [1 if j == nombre_de_noeuds-1 else 0 for j in range(nombre_de_noeuds)]
                matrix.append(sous_liste)      
            else:
                for j in range(nombre_de_noeuds):
                    if j == i-1:
                        sous_liste.append(delta_t*delta_r*Deff/(2*r[i]) - delta_t*Deff)
                    elif j == i:
                        sous_liste.append((delta_r)**2 + 2*delta_t*Deff)
                    elif j == i+1:
                        sous_liste.append(-(delta_t*Deff + delta_t*delta_r*Deff/(2*r[i])))
                    else:
                        sous_liste.append(0)
                matrix.append(sous_liste)
                
    while(t < t_station):
        
        if t==0:
            t+=delta_t
            continue
        b.append(0)
        for i in range(nombre_de_noeuds-1):
            if i ==0:
                continue
            b.append(C[i]*(delta_r)**2 - delta_t*(delta_r**2)*S)
        b.append(12)
     
        C = np.linalg.solve(matrix,b) 
        t= t + delta_t
        b = []

    return r,C,ref