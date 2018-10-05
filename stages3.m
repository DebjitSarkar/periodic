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

%% System Parameters
BL = pi / 4;
Z0 = 50;
Y0 = 1 / Z0;
f = 10e9;
w = 2 * pi * f;

%% Storage
Zrange = linspace(0, 1/10, Clen); % 2000 pts for Z_C = [0,10]
Crange = Zrange / w; % Capacitance values

ABCD_a = zeros(2,2,Clen);
ABCD0 = ABCD_a;
ABCD_b = ABCD_a;
S21_a = zeros(1,Clen);
S210 = S21_a;
S21_b = S21_a;

phias = zeros(1,Clen); % PHI_A
phibs = zeros(1,Clen); % PHI_B
diffs = zeros(1,Clen); % |PHI_B - PHI_A|

%S21 = 2/(ABCD(1,1)+ABCD(1,2)/Z0+ABCD(2,1)*Z0+ABCD(2,2));

%% N = 1 Case
L = 1 ./ (w^2 * Crange);
b_a = (w^2 .* L .* Crange * (1 - Dm) - 1) ./ (w * L);
b0 = (w^2 .* L .* Crange - 1) ./ (w * L);
b_b = (w^2 .* L .* Crange * (1 + Dm) - 1) ./ (w * L);

for i = 1:Clen
    ABCD_a(:,:,i) = [cos(BL)-sin(BL)*b_a(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b_a(i)*Z0*2*cos(BL)-1j*(Z0*b_a(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b_a(i)*Z0];
    ABCD0(:,:,i) = [cos(BL)-sin(BL)*b0(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b0(i)*Z0*2*cos(BL)-1j*(Z0*b0(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b0(i)*Z0];
    ABCD_b(:,:,i) = [cos(BL)-sin(BL)*b_b(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b_b(i)*Z0*2*cos(BL)-1j*(Z0*b_b(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b_b(i)*Z0];
    
    S21_a(i) = 2/(ABCD_a(1,1,i)+ABCD_a(1,2,i)/Z0+ABCD_a(2,1,i)*Z0+ABCD_a(2,2,i));
    S210(i) = 2/(ABCD0(1,1,i)+ABCD0(1,2,i)/Z0+ABCD0(2,1,i)*Z0+ABCD0(2,2,i));
    S21_b(i) = 2/(ABCD_b(1,1,i)+ABCD_b(1,2,i)/Z0+ABCD_b(2,1,i)*Z0+ABCD_b(2,2,i));
    
    phias(i) = angle(S21_a(i));
    phibs(i) = angle(S21_b(i));
end

phias = rad2deg(unwrap(phias)); % degrees
phibs = rad2deg(unwrap(phibs));

diffs = abs(phibs - phias);

[~, idx] = min(abs(diffs - 360));
Copt = Crange(idx);  % Optimal value for C

%% Printing
fprintf('Capacitance: C = %e F\n', Copt);
fprintf('Phase difference: %f deg\n', diffs(idx));
fprintf('Power:\n\tP_A (+Dm) = %4.2f dB\n\tP_B (-Dm) = %4.2f dB\n', db(abs(S21_a(idx))^2,'power'), db(abs(S21_b(idx))^2,'power'));

figure;
hold on;
plot(1:Clen, phias);
plot(1:Clen, phibs);
plot(1:Clen, diffs);
legend('PHI_A', 'PHI_B', 'DIFF');


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

