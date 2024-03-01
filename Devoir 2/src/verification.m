% Demande à l'utilisateur le choix du schéma de simulation
schema = input('Choisissez le schema de simulation: 1 pour schema_1, 2 pour schema_2: ');

% Paramètres de simulation
nbr = input('Entrez le nombre de noeuds: ');
R = 0.5; % Rayon du pilier
Deff = 1e-10; % Coefficient de diffusion effectif
S = 8e-9; % Terme source
Ce = 12; % Concentration à la surface

% Vérification du choix de l'utilisateur et appel de la fonction correspondante
if schema == 1
    [r, C, ref] = transit_shema_1(nbr, R, Deff, S, Ce);
elseif schema == 2
    [r, C, ref] = transit_shema_2(nbr, R, Deff, S, Ce);
else
    disp('Choix invalide. Veuillez choisir 1 ou 2.');
end

% Affichage des résultats
figure; % Crée une nouvelle figure
plot(r, C, '-b', 'LineWidth', 2, 'DisplayName', 'Solution Numérique'); % Trace la solution numérique
hold on; % Garde la figure active pour le prochain tracé
plot(r, ref, '--r', 'LineWidth', 2, 'DisplayName', 'Solution Analytique'); % Trace la solution analytique de référence
xlabel('Distance (m)'); % Étiquette de l'axe des x
ylabel('Concentration de sel (mol/m^3)'); % Étiquette de l'axe des y
title('Comparaison entre la solution numérique et la solution analytique'); % Titre du graphique
legend('show'); % Affiche la légende
grid on; % Active la grille pour une meilleure lisibilité
