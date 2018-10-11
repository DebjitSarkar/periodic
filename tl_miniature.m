% Debjit Sarkar
% October 2018

clear; clc;
%% Structure
%
% ______Bd/2_______Bd/2______
%   |  ------  |  ------  |
%   Y1         Y2         Y1
%   |          |          |
%  GND        GND        GND
%

%% Solution
fun = @optimABCD;
y = fsolve(fun, [1j,1j]);
c = y / (1j * 2 * pi * 10e9);
Y1 = y(1);
Y2 = y(2);

%% Verification
Z0 = 50; Y0 = 1/Z0;
betad = pi/10;
tl = [cos(betad/2), 1j*Z0*sin(betad/2);
    1j*Y0*sin(betad/2), cos(betad/2)];
ABCD_tot = [1,0;Y1,1]*tl*[1,0;Y2,1]*tl*[1,0;Y1,1];

A = ABCD_tot(1,1);
B = ABCD_tot(1,2);
C = ABCD_tot(2,1);
D = A;

S21 = 2 / (A + B / Z0 + C * Z0 + D);
assert(abs(abs(rad2deg(angle(S21))) - 90) < 1);
assert(abs(S21) - 1 < 0.01);

%% Printing
fprintf('Solution:\n');
fprintf('Y1 = %fj\n',imag(Y1));
fprintf('Y2 = %fj\n',imag(Y2));

%% Function for fsolve
function F = optimABCD(Y)
Z0 = 50; Y0 = 1/Z0;
betad = pi/10;
tl = [cos(betad/2), 1j*Z0*sin(betad/2);
    1j*Y0*sin(betad/2), cos(betad/2)];
ABCD_tot = [1,0;Y(1),1]*tl*[1,0;Y(2),1]*tl*[1,0;Y(1),1];

F(1) = ABCD_tot(1,1);
F(2) = ABCD_tot(1,2) - 1j*Z0;
F(3) = ABCD_tot(2,1) - 1j*Y0;
end


