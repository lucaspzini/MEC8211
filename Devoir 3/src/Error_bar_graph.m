%% Error bar graph
E = 24.7 - 80.6;
error_rel = 100 * E/80.6; % percentage error
u_val = 2*17.8;
u_val_rel = 100 * u_val/80.6; % percentage uncertainty

figure(1)
fig = gcf;
fig.Position(3:4) = [300 400]; % Set the width and height of the figure
hold on
grid on
plot([0 2],[0 0],'--','Color','r','LineWidth',2)
errorbar(1,error_rel,u_val_rel,'k','LineWidth',2);
plot(1,error_rel,'o','MarkerFaceColor','black', "MarkerEdgeColor","black")
set(gca,'xtick',[],'YMinorTick','on','TickLength', [0.02 0.02]);
axis([0.9 1.1 -120 10])
ylabel("\delta_{model} = E/k_{exp} = (S-D)/D (%)")
title('Erreur relative du modèle', 'FontWeight', 'bold');

% Draw a rectangle around the plot
x = 0.9; % x-coordinate of the lower-left corner of the rectangle
y = -120; % y-coordinate of the lower-left corner of the rectangle
width = 0.2; % Width of the rectangle
height = 130; % Height of the rectangle
rectangle('Position',[x y width height],'EdgeColor','black');

hold off