function [r, C, ref] = transit_shema_2(nbr, R, Deff, S, Ce)

    % Initialisation
    delta_r = R / (nbr - 1);
    r = linspace(0, R, nbr); % points du domaine
    matrix = zeros(nbr);
    b = zeros(1, nbr);
    C = zeros(1, nbr);
    
    % Solution analytique
    ref = 0.25 * (S / Deff) * (R^2) * ((r / R).^2 - 1) + Ce;

    % paramètre de temps
    t = 0;
    t_station = 2e10;
    delta_t = 2000;
    
    % Construction de la matrice pour le schéma transitoire 2
    matrix(1,1) = -3;
    matrix(1,2) = 4;
    matrix(1,3) = -1;
    for i = 2:nbr-1
        matrix(i, i-1) = delta_t * delta_r * Deff / (2 * r(i)) - delta_t * Deff;
        matrix(i, i) = delta_r^2 + 2 * delta_t * Deff;
        matrix(i, i+1) = -(delta_t * Deff + delta_t * delta_r * Deff / (2 * r(i)));
    end
    matrix(nbr, nbr) = 1;

    % Boucle temporelle
    while t < t_station
        b(1) = 0;
        for i = 2:nbr-1
            b(i) = C(i) * delta_r^2 - delta_t * S * delta_r^2;
        end
        b(nbr) = Ce;

        % Résolution du système linéaire
        C = matrix \ b';

        t = t + delta_t;
        b = zeros(1, nbr); % Réinitialise b pour le prochain pas de temps
    end
end
