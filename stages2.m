% Debjit Sarkar
% stages_v2.m

Ptol = 0.5; % Power tolerance
Nmax = 30; % Maximum number of stages considered
Dm = 0.2; % Modulation depth as %-tage of C
BL = pi / 4;
Z0 = 50;
Y0 = 1 / Z0;
syms w L C
ABCD = [cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L),...
    1j*Z0*sin(2*BL)+1j*Z0^2*sin(BL)^2*(1-w^2*L*C);...
    1j*Y0*sin(2*BL)-1j*cos(BL)^2*(1-w^2*L*C)/(w*L),...
    cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L)];
S21 = 2/(ABCD(1,1)+ABCD(1,2)/Z0+ABCD(2,1)*Z0+ABCD(2,2));


for n = 1:Nmax
    
    
    
end
