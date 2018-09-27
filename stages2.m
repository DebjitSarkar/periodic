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
Nmax = 5; %30 % Maximum number of stages considered
Clen = 250; %2000 % Number of values used for C

%% Setup
digits(15); % Used to speed up symbolic calculations, default = 32

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
% Ccurr = 1e-12; % Current value
%
% S21new = subs(S21, 'L', 1 / (w^2 * Ccurr));
%
% S21a = subs(S21new, 'C', Ccurr * (1 + Dm));
% S21b = subs(S21new, 'C', Ccurr * (1 - Dm));
%
% phia = angle(S21a);
% phib = angle(S21b);
%
% diff = double(abs(rad2deg(vpa(phib - phia))));
% diff

% for each n, find optimal c for 360 degrees
% increase n until you reach the power tolerance

%% Main loop
diffs = zeros(1,Clen);
phias = diffs;
phibs = phias;
for N = 1:Nmax
    ABCD_new = ABCD ^ N;
    S21 = 2/(ABCD_new(1,1)+ABCD_new(1,2)/Z0+ABCD_new(2,1)*Z0+ABCD_new(2,2));
    for i = 1:Clen
        S21new = subs(S21, 'L', 1 / (w^2 * Crange(i)));
        
        S21a = subs(S21new, 'C', Crange(i) * (1 + Dm));
        S21b = subs(S21new, 'C', Crange(i) * (1 - Dm));
        
        phias(i) = rad2deg(vpa(angle(S21a)));
        phibs(i) = rad2deg(vpa(angle(S21b)));
        
        diff = double(abs(rad2deg(vpa(angle(S21b) - angle(S21a)))));
        % Find vector of diffs
        diffs(i) = diff;
        if(mod(i, 75) == 0)
            fprintf('On i = %i\n', i);
        end
    end
    [~, idx] = min(abs(diffs - 360));
    
    % Calculate S21^2 for Crange(idx)
    S21new = subs(S21, 'L', 1 / (w^2 * Crange(idx)));
    S21a = subs(S21new, 'C', Crange(idx) * (1 + Dm));
    S21b = subs(S21new, 'C', Crange(idx) * (1 - Dm));
    fprintf('N = %i Capacitance: C = %e F\n', N, Crange(idx));
    fprintf('Phase difference (deg): %i\n', diffs(idx));
    fprintf('Power:\n\tP_A (+Dm) = %4.2f\n\tP_B (+Dm) = %4.2f\n', abs(S21a)^2, abs(S21b)^2);
    
    %if(min(abs(S21a)^2, abs(S21b)^2) > Ptol)
    %    fprintf('Done\n');
    %    break;
    %end
    figure;
    hold on;
    plot(1:Clen, phias);
    plot(1:Clen, phibs);
    plot(1:Clen, diffs);
    legend('PHI_A', 'PHI_B', 'DIFF');
end



