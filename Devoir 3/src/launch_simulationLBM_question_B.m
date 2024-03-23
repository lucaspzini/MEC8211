
% MATLAB script to launch a fiber structure generation and the corresponding LBM simulation
%
%INPUT VARIABLES:
%
% SEED: integer representing the seed for initializing the random
% generator. If seed=0, automatic seed generation. If you want to reproduce
% the same fiber structure, use the same seed (fibers will be located at the same place). 
%
% MEAN_D: contains the mean fiber to be used
%
% STD_D: contains the standard deviation of the fiber diameters
%
% PORO: estimated porosity of the fiber structure to be generated
% 
% NX: domain lateral size in grid cell

seed=0;
rng('shuffle')
deltaP= 0.1 ; % pressure drop in Pa
NX= 200 ;
% poro= 0.9 ;
mean_fiber_d= 12.5 ; % in microns
std_d= 2.85 ; % in microns
dx= 1e-6 ; % grid size in m
filename= 'fiber_mat.tiff' ;

% Normal random distribution for porosity
mean_poro = 0.9 ;
std_poro = 0.0075 ;
poro = normrnd(mean_poro,std_poro,[1,1000]);

% Vector to store permeability
k = zeros(1,length(poro)) ;

for i = 1:length(poro)
    % generation of the fiber structure
    [d_equivalent]=Generate_sample(seed,filename,mean_fiber_d,std_d,poro(i),NX,dx);
    
    % calculation of the flow field and the permeability from Darcy Law
    k(i) = LBM_B(filename,NX,deltaP,dx,d_equivalent);
    i
end

%% PDF of permeability
%% Normal distribution
figure()
histfit(k)
pd = fitdist(k','Normal');
xlabel('k (\mum^2)');
ylabel('Fréquence');
title('Histogramme et PDF des perméabilités (k)', 'FontWeight', 'bold');
str = {['Moyenne = ', num2str(pd.mu)], ['Écart type = ', num2str(pd.sigma)]};
dim = [0.55, 0.65, 0.1, 0.1]; % Position of annotation box [x, y, width, height]
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

%% Lognormal distribution
figure()
histfit(k,40,'Lognormal')
pd = fitdist(k','Lognormal');
xlabel('k (\mum^2)');
ylabel('Fréquence');
title('Histogramme et PDF des perméabilités (k)', 'FontWeight', 'bold');
str = {['Moyenne = ', num2str(pd.mu)], ['Écart type = ', num2str(pd.sigma)]};
dim = [0.55, 0.65, 0.1, 0.1]; % Position of annotation box [x, y, width, height]
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

%% CDF of permeability
figure()
cdfplot(k)
xlabel("k [\mum^2]")
ylabel("Probabilité cumulative")
title('CDF des perméabilités (k)', 'FontWeight', 'bold');

%% Porosité
figure()
histfit(poro)
pd = fitdist(poro','Normal');
xlabel('\epsilon');
ylabel('Fréquence');
title('Histogramme des porosités (\epsilon)', 'FontWeight', 'bold');
str = {['Moyenne = ', num2str(pd.mu)], ['Écart type = ', num2str(pd.sigma)]};
dim = [0.15, 0.75, 0.1, 0.1]; % Position of annotation box [x, y, width, height]
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

%% Médiane
median = median(k);
fprintf('Médiane = S = %f\n', median);