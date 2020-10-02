function xdot = fU(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: fU.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Strongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2015 1:34:00
   
% The global variables
global lambda gamma d beta M

% state
z1 = x(1);
z2 = x(2);
q = x(3);

% Black box: L could be anything. 
if (q == 0) 
    u = - lambda*z2 - gamma*GradientL(z1); 
elseif (q == 1)
    u = - 2*d*z2 - (1/M)*GradientL(z1 + beta*z2); 
end

xdot = [z2;u;0]; 
end