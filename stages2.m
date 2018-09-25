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

ABCD = [cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L),... %A
    1j*Z0*sin(2*BL)+1j*Z0^2*sin(BL)^2*(1-w^2*L*C);... %B
    1j*Y0*sin(2*BL)-1j*cos(BL)^2*(1-w^2*L*C)/(w*L),... %C
    cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L)]; %D

S21 = 2/(ABCD(1,1)+ABCD(1,2)/Z0+ABCD(2,1)*Z0+ABCD(2,2));

Cinit = 1e-12; % Start
Cinc = 1e-12; % Increment by
Ccurr = 1e-12; % Current value
L0 = 1 / (w0^2 * Ccurr);

S21new = subs(S21, 'L', L0);

Ca = Ccurr * (1 + Dm);
Cb = Ccurr * (1 - Dm);

phia = angle(subs(S21new, 'C', Ca));
phib = angle(subs(S21new, 'C', Cb));
diff = abs(rad2deg(vpa(phib - phia)));
diff % Looks wrong

% for each n, find optimal c for 360 degrees
% increase n until you reach the power tolerance
