% Fitting of the uniaxial data using Diani-Rey model
% Authors: Serkan Can, Onur Ata, Atakan Aygun



% Reading Data from utdata.csv file
lambda = readtable("utdata.csv",Range="A1:A101",ReadVariableNames=false);
P11 = readtable("utdata.csv",Range="B1:B101",ReadVariableNames=false);
lambda = table2array(lambda);
P11 = table2array(P11);

% Derivative functions
c_1 = @(x)exp(x(1)+x(2)*(lambda.^2+(2./lambda)-3)+x(3)*(lambda.^2+(2./lambda)-3).^2);
c_2 = @(x)exp(x(4)+x(5)*log(2.*lambda+(1./lambda.^2)));

% Objective function
DR_biaxial = @(x)2*(c_1(x)+c_2(x)./lambda).*(lambda-1./lambda.^2);
ErrorFunc =@(x)sum((DR_biaxial(x)-P11).^2);

% Initialize and start optimization
x0 = 0.0000000000001*ones(5,1);
options = optimoptions('fmincon','Display','off');

problem.options = options;
problem.solver = 'fmincon';
problem.objective = ErrorFunc;
problem.x0 = x0;

x = fmincon(problem);
fprintf('Error: %.4g\n', ErrorFunc(x))

% x = [log(0.0338) 0.0221 0.000315 log(0.0552) 0.9985 0 1.0];

% Quality of fit
quality = sum((DR_biaxial(x) - P11).^2 ./ P11);
fprintf('Quality of the fit: %f\n', quality)

fit = DR_biaxial(x);

% Plotting
fs = 10;

fgh =figure(1); h1 = axes;
set(fgh,'Units','centimeters ') 
set(fgh, 'Position', [4.0,4.0,32.0,16.0]); 

plot(lambda, fit, 'k', LineWidth=2)
hold on
scatter(lambda, P11,'filled')

set(gca,'FontSize',fs+2)

leg1 = legend('Fit', 'Data');

set(leg1,'interpreter','latex','fontsize',fs+2,'Location','NorthWest','box','off');
xlabel('$$\lambda$$','interpreter','latex','fontsize',fs+4);
ylabel('$$P_{11}$$','interpreter','latex','fontsize',fs+4);
title('Uniaxial Test', 'interpreter','latex','fontsize',fs+4)


% axis square ;
grid on; grid minor; 
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
set(gca,'FontName','Times New Roman');


print -dpdf -painters ut_fit
hold off;
 
