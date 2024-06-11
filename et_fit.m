% Fitting of the equiaxial and uniaxial data using Diani-Rey model
% Authors: Serkan Can, Onur Ata, Atakan Aygun

% Reading Uniaxial Tension (UT) Data from utdata.csv file
lambda_ut = readtable("utdata.csv",Range="A1:A101",ReadVariableNames=false);
P11_ut = readtable("utdata.csv",Range="B1:B101",ReadVariableNames=false);
lambda_ut = table2array(lambda_ut);
P11_ut = table2array(P11_ut);

% Reading Equiaxial Tension (ET) Data from d10.csv file
lambda_et = readtable("d10.csv",Range="A1:A101",ReadVariableNames=false);
P11_et = readtable("d10.csv",Range="B1:B101",ReadVariableNames=false);
lambda_et = table2array(lambda_et);
P11_et = table2array(P11_et);

% The first derivatives of strain energy function with respect to I1 and I2 
% For UT
c_1 = @(x)exp(x(1)+x(2)*(lambda_ut.^2+(2./lambda_ut)-3)+x(3)*(lambda_ut.^2+(2./lambda_ut)-3).^2);
c_2 = @(x)exp(x(4)+x(5)*log(2.*lambda_ut+(1./lambda_ut.^2)));

% The first derivatives of strain energy function with respect to I1 and I2 
% For ET
d_1 = @(x)exp(x(1)+x(2)*(2.*lambda_et.^2+(1./lambda_et.^4)-3)+x(3)*(2.*lambda_et.^2+(1./lambda_et.^4)-3).^2);
d_2 = @(x)exp(x(4)+x(5)*log(lambda_et.^4+(2./lambda_et.^2)));

% The first Piola–Kirchhoff stress tensor and Objective function for UT
DR_ut = @(x)2*(c_1(x)+c_2(x)./lambda_ut).*(lambda_ut-1./lambda_ut.^2);
ErrorFunc_ut =@(x)sum((DR_ut(x)-P11_ut).^2);

% The first Piola–Kirchhoff stress tensor and Objective function for ET
DR_et = @(x)2*(d_1(x)+d_2(x).*lambda_et.^2).*(lambda_et-(1./lambda_et.^5));
ErrorFunc_et =@(x)sum((DR_et(x)-P11_et).^2);

% Total Error with weights x(6), x(7)
TotalError =@(x)(x(6)*ErrorFunc_et(x)+x(7)*ErrorFunc_ut(x));

% Optimize Equiaxial Only
x0_0 = 0.0000000000001*ones(5,1);
x1 = fmincon(ErrorFunc_et, x0_0);
ErrorFunc_et(x1);
% Initialize and start optimization

% Initial Conditions
% x0= [log(0.135) 0.0018 0.000326 log(0.102) -0.625 0.8 0.2]; % From Reference Paper
x0= [-2.12958799550019 0.00911420261928685 0.000262871856779523 x1(4) x1(5) 0.6 0.4];
% x0 = [0 0 0 0 0 0.2 0.8]; % Zero initial guess

options = optimoptions('fmincon','Display','off');
problem.options = options;
problem.solver = 'fmincon';
problem.objective = TotalError;
problem.x0 = x0;   %IC

% Inequality Conditions (did not used)
problem.A = [];          
problem.b = [];

% Equality Condition for weights (W1+W2=1)
problem.Aeq=[0 0 0 0 0 1 1];
problem.beq=1.0;

% Lower and Upper bounds for parameters and weights
problem.lb = [-Inf,-Inf,-Inf,-Inf,-Inf,0.0,0.0];
problem.ub = [Inf,Inf,Inf,Inf,Inf,1,1];

% Optimization with FMINCON function
x = fmincon(problem);

fprintf('Error: %.4g\n', TotalError(x)) % Total Error

% Quality of fit
quality_ut = sum((DR_ut(x) - P11_ut).^2 ./ P11_ut, "omitnan");
quality_et = sum((DR_et(x) - P11_et).^2 ./ P11_et, "omitnan");

quality = quality_ut + quality_et;

fprintf('Quality of the UT fit: %f\n', quality_ut)
fprintf('Quality of the ET fit: %f\n', quality_et)
fprintf('Quality of the overall fit: %f\n', quality)


% Plotting

fit1 = DR_et(x);
fit2 = DR_ut(x);

fs = 10;

fgh =figure(1); h1 = axes;
set(fgh,'Units','centimeters ') 
set(fgh, 'Position', [4.0,4.0,32.0,16.0]); 

hold on

plot(lambda_et, fit1, 'k', LineWidth=2)
plot(lambda_ut, fit2, 'r', LineWidth=2)
scatter(lambda_et, P11_et,'filled')
scatter(lambda_ut, P11_ut,'filled')
% scatter(lambda_et, P11_ut,'filled')

set(gca,'FontSize',fs+2)

leg1 = legend('ET Model Fit', 'UT Model Fit', 'ET Data', 'UT Data');

set(leg1,'interpreter','latex','fontsize',fs+2,'Location','NorthWest','box','off');
xlabel('$$\lambda$$','interpreter','latex','fontsize',fs+4);
ylabel('$$P_{11}$$','interpreter','latex','fontsize',fs+4);
title('Simultaneous test (UT+ET)', 'interpreter','latex','fontsize',fs+4)


% axis square ;
grid on; grid minor; 
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
set(gca,'FontName','Times New Roman');


% print -dpdf -painters et_fit_failed
hold off;
 
