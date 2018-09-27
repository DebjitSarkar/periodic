% Debjit Sarkar
% September 2018
% stages2.m - Optimizing the number of stages necessary for
% phase shifting from [0,360] above a power threshold

%% Structure
%
%  ____45deg_________________45deg____
%  ------------     |     ------------
%           |     -----     |
%          GND    |   |    GND
%                 C   L
%                 |   |
%                 -----
%                   |
%                  GND
%

%% Inputs

Ptol = 0.5; % Power tolerance
Dm = 0.2; % Modulation depth as %-tage of C
Nmax = 30; % Maximum number of stages considered
Clen = 2000; % Number of values used for C

%% Setup
BL = pi / 4;
Z0 = 50;
Y0 = 1 / Z0;
f = 10e9;
w = 2 * pi * f;

syms L C;

%% One segment ABCD and S-matrix
ABCD = [cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L),...       %A
    1j*Z0*sin(2*BL)+1j*Z0^2*sin(BL)^2*(1-w^2*L*C)/(w*L);...  %B
    1j*Y0*sin(2*BL)-1j*cos(BL)^2*(1-w^2*L*C)/(w*L),...       %C
    cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L)];             %D

S21 = 2/(ABCD(1,1)+ABCD(1,2)/Z0+ABCD(2,1)*Z0+ABCD(2,2));

%% Check det(ABCD) = 1
q = subs(det(ABCD), 'C', 1e-15);
qq = subs(q, 'L', 1e-7);
assert(double(vpa(abs(qq - 1))) < 1e-6, 'Determinant =/= 1');

%%
Zrange = linspace(0, 1/10, Clen); % 2000 pts for Z_C = [0,10]
Crange = Zrange / w;
Cinit = 1e-12; % Start
Cincr = 1e-12; % Increment by
Ccurr = 1e-12; % Current value
L0 = 1 / (w^2 * Ccurr);

S21new = subs(S21, 'L', L0);

S21a = subs(S21new, 'C', Ccurr * (1 + Dm));
S21b = subs(S21new, 'C', Ccurr * (1 - Dm));

phia = angle(S21a);
phib = angle(S21b);

diff = double(abs(rad2deg(vpa(phib - phia))));
diff

% for each n, find optimal c for 360 degrees
% increase n until you reach the power tolerance

%% Main loop

for N = 1:Nmax
    for i = 1:Clen
        % Find vector of diffs
    end
    [~, idx] = min(abs(diffs - 2 * pi));
    
end

