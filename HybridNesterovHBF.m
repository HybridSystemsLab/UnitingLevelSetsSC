%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HybridNesterovHBF.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Strongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2015 1:34:00

clear all

% global variables
global delta d beta M gamma lambda a c_0 c_10 tauMin tauMax c

%%%%%%%%%%%%%%%%%%%%%
% setting the globals
%%%%%%%%%%%%%%%%%%%%%
setMinima();

% Nesterov constants
kappa = 1; % kappa >= 1, kappa = M/mu = 2/2
M = 2;
d = 1/(sqrt(kappa) + 1);
beta = (sqrt(kappa) - 1)/(sqrt(kappa) + 1);
a = d + beta/(2*kappa);

% Heavy Ball constants
gamma = 2/3; 
lambda = 40; 

% Uniting level sets
c_0 = 400; % \mathcal{U}_0 
c_10 = 111; % \mathcal{T}_{1,0} 

% Poveda-Na constants
tauMin = 3; 
tauMax = 5.63; 
c = 1/4; 

% Convergence within one percent from min
delta = 0.5;

%%%%%%%%%%%%%%%%%%%%%
% setting the locals
%%%%%%%%%%%%%%%%%%%%%

deltaVecNest = [0,0,0];
deltaVecHBF = [0,0,0];
deltaVecUniting = [0,0,0];
deltaVecPN = [0,0,0];

% initial conditions
z1_0 = 50; 
z2_0 = 0;
z2_00 = 50;
q_0 = 0;
tau_0 = 3; 

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];
x00 = [z1_0;z2_00;tau_0];

% simulation horizon
TSPAN=[0 400];
JSPAN = [0 5000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

[tNest,jNest,xNest] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options);

deltaVecNest = timeToConv(xNest,tNest);

[tHBF,jHBF,xHBF] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x0,TSPAN,JSPAN,rule,options);

deltaVecHBF = timeToConv(xHBF,tHBF);

[tPN,jPN,xPN] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options);

deltaVecPN = timeToConv(xPN,tPN);

% Finally, simulate the hybrid closed-loop heavy ball system
[tUniting,jUniting,xUniting] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options);

deltaVecUniting = timeToConv(xUniting,tUniting);

minarc = min([length(xUniting),length(xPN),length(xHBF),length(xNest)]);
ta = [tUniting(1:minarc),tPN(1:minarc),tHBF(1:minarc),tNest(1:minarc)];
ja = [jUniting(1:minarc),jPN(1:minarc),jHBF(1:minarc),jNest(1:minarc)];
xa = [xUniting(1:minarc,1),xPN(1:minarc,1),xHBF(1:minarc,1),xNest(1:minarc,1)];
xb = [xUniting(1:minarc,2),xPN(1:minarc,2),xHBF(1:minarc,2),xNest(1:minarc,2)];

figure(1)
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 3;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 3;
subplot(2,1,1), plotHarc(ta,ja,xa,[],modificatorF,modificatorJ);
hold on
plot(deltaVecPN(3),deltaVecPN(1),'k.','MarkerSize', 20)
strDeltaPN = [num2str(deltaVecPN(3)), 's'];
text(deltaVecPN(3),deltaVecPN(1),strDeltaPN,'HorizontalAlignment','left','VerticalAlignment','bottom');
plot(deltaVecUniting(3),deltaVecUniting(1),'k.','MarkerSize', 20)
strDeltaUniting = [num2str(deltaVecUniting(3)), 's'];
text(deltaVecUniting(3),deltaVecUniting(1),strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom');
plot(deltaVecNest(3),deltaVecNest(1),'k.','MarkerSize', 20)
strDeltaNest = [num2str(deltaVecNest(3)), 's'];
text(deltaVecNest(3),deltaVecNest(1),strDeltaNest,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([0 15 -20 70])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
axes('Position',[0.7 0.78 0.15 0.08])
box on
hold on
plot(tHBF,xHBF(:,1),'LineWidth',3)
plot(deltaVecHBF(3),deltaVecHBF(1),'k.','MarkerSize', 20)
strDeltaHBF = [num2str(deltaVecHBF(3)), 's'];
text(deltaVecHBF(3),deltaVecHBF(1),strDeltaHBF,'HorizontalAlignment','left','VerticalAlignment','bottom');
hold off
set(gca,'xtick',[])
set(gca,'ytick',[])
axis([0 200 -20 70])
hold off
subplot(2,1,2), plotHarc(ta,ja,xb,[],modificatorF,modificatorJ);
hold on
plot(deltaVecPN(3),deltaVecPN(2),'k.','MarkerSize', 20)
plot(deltaVecHBF(3),deltaVecHBF(2),'k.','MarkerSize', 20)
plot(deltaVecNest(3),deltaVecNest(2),'k.','MarkerSize', 20)
plot(deltaVecUniting(3),deltaVecUniting(2),'k.','MarkerSize', 20)
axis([0 15 -50 70])
grid on
ylabel('z_2','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
saveas(gcf,'Plots\ComparisonPlots2','png')
saveas(gcf,'Plots\ComparisonPlots2','epsc')

minarc = min([length(xUniting),length(xPN)]);
te = [tUniting(1:minarc),tPN(1:minarc)];
je = [jUniting(1:minarc),jPN(1:minarc)];
xe = [xUniting(1:minarc,1),xPN(1:minarc,1)];
xf = [xUniting(1:minarc,2),xPN(1:minarc,2)];

figure(2) 
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 3;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 3;
subplot(2,1,1), [x1,t1] = plotHarc(te,je,xe,[],modificatorF,modificatorJ);
hold on
plot(deltaVecPN(3),deltaVecPN(1),'k.','MarkerSize', 20)
strDeltaPN = [num2str(deltaVecPN(3)), 's'];
text(deltaVecPN(3),deltaVecPN(1),strDeltaPN,'HorizontalAlignment','left','VerticalAlignment','top');
plot(deltaVecUniting(3),deltaVecUniting(1),'k.','MarkerSize', 20)
strDeltaUniting = [num2str(deltaVecUniting(3)), 's'];
text(deltaVecUniting(3),deltaVecUniting(1),strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([0 5 -20 70])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
subplot(2,1,2), plotHarc(te,je,xf,[],modificatorF,modificatorJ);
hold on
plot(deltaVecPN(3),deltaVecPN(2),'k.','MarkerSize', 20)
plot(deltaVecUniting(3),deltaVecUniting(2),'k.','MarkerSize', 20)
axis([0 5 -50 70])
grid on
ylabel('z_2','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
saveas(gcf,'Plots\ComparisonPlots4','png')
saveas(gcf,'Plots\ComparisonPlots4','epsc')