% Debjit Sarkar
% stages2.m

Ptol = 0.5; % Power tolerance
Nmax = 30; % Maximum number of stages considered
Dm = 0.2; % Modulation depth as %-tage of C

BL = pi / 4;
Z0 = 50;
Y0 = 1 / Z0;
f0 = 10e9;
w0 = 2 * pi * f0;

syms L C %w
w = w0;

% Check either B or C, or both A and D are wrong (probably C)
% Check det(ABCD) = 1
ABCD = [cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L),... %A
    1j*Z0*sin(2*BL)+1j*Z0^2*sin(BL)^2*(1-w^2*L*C)/(w*L);...  %B
    1j*Y0*sin(2*BL)-1j*cos(BL)^2*(1-w^2*L*C)/(w*L),... %C
    cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L)];       %D

S21 = 2/(ABCD(1,1)+ABCD(1,2)/Z0+ABCD(2,1)*Z0+ABCD(2,2));

%% det(ABCD) = 1 check
q = subs(det(ABCD), 'C', 1e-15);
qq = subs(q, 'L', 1e-7);
vpa(qq)

%%
%TL_ABCD = [cos(BL), 1j*Z0*sin(BL);
 %          1j*sin(BL)/Z0, cos(BL)];
syms a1 b1 c1 d1 a2 b2 c2 d2;
syms BL Y;
q1 = [a1 b1; c1 d1];
q2 = [a2 b2; c2 d2];
q3 = q1 * q2 * q1;
q4 = subs(q3, 'a1', cos(BL));
q4 = subs(q4, 'b1', 1j*Z0*sin(BL));
q4 = subs(q4, 'c1', 1j*sin(BL)/Z0);
q4 = subs(q4, 'd1', cos(BL));
q4 = subs(q4, 'a2', 1);
q4 = subs(q4, 'b2', 0);
q4 = subs(q4, 'c2', Y); % c2 = Y
q4 = subs(q4, 'd2', 1);

q5 = subs(det(q4), 'BL', pi/4); % This works!

%%
Cinit = 1e-12; % Start
Cincr = 1e-12; % Increment by
Ccurr = 1e-12; % Current value
L0 = 1 / (w0^2 * Ccurr);

S21new = subs(S21, 'L', L0);

Ca = Ccurr * (1 + Dm);
Cb = Ccurr * (1 - Dm);

phia = angle(subs(S21new, 'C', Ca));
phib = angle(subs(S21new, 'C', Cb));

diff = abs(rad2deg(vpa(phib - phia)));
%diff % Looks wrong

% for each n, find optimal c for 360 degrees
% increase n until you reach the power tolerance
