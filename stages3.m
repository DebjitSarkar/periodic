% Debjit Sarkar
% October 2018
% stages3.m - Optimizing the number of stages necessary for
% phase shifting from [0,360] above a power threshold
% Difference from stages2.m: Different structure

%% Structure
%
%  ____________45deg____________
%    |     -------------     |
%  -----     |       |     -----
%  |   |    GND     GND    |   |
%  C   L                   C   L
%  |   |                   |   |
%  -----                   -----
%    |                       |
%   GND                     GND
%

%% Inputs
Ptol = 0.5; % Power tolerance
Dm = 0.2; % Modulation depth as %-tage of C

Nmax = 8; %30 % Maximum number of stages considered
Clen = 250; %2000 % Number of values used for C
phase_err = 0.5; % Degrees away from 360 for total phase shift
verbose = true; % Printing of statements during code execution
verbose_check = 75; % How often to print statements

%% Setup
digits(15); % Used to speed up symbolic calculations, default = 32

BL = pi / 4;
Z0 = 50;
Y0 = 1 / Z0;
f = 10e9;
w = 2 * pi * f;

syms L C;

%% Storage
Zrange = linspace(0, 1/10, Clen); % 2000 pts for Z_C = [0,10]
Crange = Zrange / w; % Capacitance values

phias = zeros(1,Clen); % PHI_A
phibs = zeros(1,Clen); % PHI_B
diffs = zeros(1,Clen); % |PHI_B - PHI_A|

%% ABCD Calculation
syms Y BL Z0;
TL_ABCD = [cos(BL), 1j*Z0*sin(BL);
    1j*sin(BL)/Z0, cos(BL)];
shunt_ABCD = [1 0; Y 1];
ABCD = shunt_ABCD * TL_ABCD * shunt_ABCD;

syms b;
% effective ABCD
ABCDef = [sqrt(2)*(0.5-25*b), sqrt(2)*25j;
    sqrt(2)*1j*(b/2+1/100+b.*(0.5-25*b)),sqrt(2)*(0.5-25*b)];

%% Unit cell ABCD and S-matrix
% ABCD = [cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L),...       %A
%     1j*Z0*sin(2*BL)+1j*Z0^2*sin(BL)^2*(1-w^2*L*C)/(w*L);...  %B
%     1j*Y0*sin(2*BL)-1j*cos(BL)^2*(1-w^2*L*C)/(w*L),...       %C
%     cos(2*BL)+Z0*sin(2*BL)*(1-w^2*L*C)/(2*w*L)];             %D

S21 = 2/(ABCD(1,1)+ABCD(1,2)/Z0+ABCD(2,1)*Z0+ABCD(2,2));

%% Check det(ABCD) = 1
q = subs(det(ABCD), 'C', 1e-15);
qq = subs(q, 'L', 1e-7);
assert(double(vpa(abs(qq - 1))) < 1e-6, 'Determinant =/= 1');

%% N = 1 Case
% b = w^2LC-1 / wL
L = 1 ./ (w^2 * Crange);
b0 = (w^2 .* L .* Crange * (1 - Dm) - 1) ./ (w * L);
b1 = (w^2 .* L .* Crange - 1) ./ (w * L);
b2 = (w^2 .* L .* Crange * (1 + Dm) - 1) ./ (w * L);

ABCD0 = zeros(2,2,Clen);
ABCD1 = ABCD0;
ABCD2 = ABCD0;

for i = 1:Clen
    ABCD0(:,:,i) = [cos(BL)-sin(BL)*b0(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b0(i)*Z0*2*cos(BL)-1j*(Z0*b0(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b0(i)*Z0];
    ABCD1(:,:,i) = [cos(BL)-sin(BL)*b1(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b1(i)*Z0*2*cos(BL)-1j*(Z0*b1(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b1(i)*Z0];
    ABCD2(:,:,i) = [cos(BL)-sin(BL)*b2(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b2(i)*Z0*2*cos(BL)-1j*(Z0*b2(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b2(i)*Z0];
end

%% Main loop
for N = 1:Nmax
    fprintf('\n========== N=%i ==========\n', N);
    ABCD_new = ABCD ^ N;
    S21 = 2/(ABCD_new(1,1)+ABCD_new(1,2)/Z0+ABCD_new(2,1)*Z0+ABCD_new(2,2));
    for i = 1:Clen
        S21new = subs(S21, 'L', 1 / (w^2 * Crange(i)));
        
        S21a = subs(S21new, 'C', Crange(i) * (1 + Dm));
        S21b = subs(S21new, 'C', Crange(i) * (1 - Dm));
        
        phias(i) = vpa(angle(S21a)); % radians
        phibs(i) = vpa(angle(S21b));
        
        if(verbose && mod(i, verbose_check) == 0) % Check that code is running
            fprintf('i = %i\n', i);
        end
    end
    
    phias = rad2deg(unwrap(phias)); % degrees
    phibs = rad2deg(unwrap(phibs));
    
    diffs = abs(phibs - phias);
    
    [~, idx] = min(abs(diffs - 360));
    Copt = Crange(idx);  % Optimal value for C
    
    % S21 calculations for Copt
    S21new = subs(S21, 'L', 1 / (w^2 * Copt));
    S21a = subs(S21new, 'C', Copt * (1 + Dm));
    S21b = subs(S21new, 'C', Copt * (1 - Dm));
    
    fprintf('Capacitance: C = %e F\n', Copt);
    fprintf('Phase difference: %f deg\n', diffs(idx));
    fprintf('Power:\n\tP_A (+Dm) = %4.2f dB\n\tP_B (-Dm) = %4.2f dB\n', db(abs(S21a)^2,'power'), db(abs(S21b)^2,'power'));
    
    figure;
    hold on;
    plot(1:Clen, phias);
    plot(1:Clen, phibs);
    plot(1:Clen, diffs);
    legend('PHI_A', 'PHI_B', 'DIFF');
    
    % Fix the abs(diffs(idx)... statement
    if(abs(diffs(idx) - 360) < phase_err && min(abs(S21a)^2, abs(S21b)^2) > Ptol)
        fprintf('Done\n');
        break;
    end
end


%% TODO
% Change the unit cell to TANK + lambda/4 + TANK

