function [value] = CU(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: CU.m
%--------------------------------------------------------------------------
% Description: Flow set
% Return 0 if outside of C, and 1 if inside C
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2015 1:34:00

global gamma c_0 c_10 M a

% state
z1 = x(1);
z2 = x(2);
q = x(3);

CalcL = CalculateL(z1);

if(q == 0)
    V0 = gamma*(CalcL) + (1/2)*z2^2;
elseif(q == 1)
    V1 = (1/2)*norm(a*z1 + z2)^2 + (1/M)*(CalcL);
end

if (q == 0 && V0 <= c_0)||(q == 1 && V1 >= c_10) 
    value = 1;
else 
    value = 0;
end

end