% Debjit Sarkar
% October 2018
% stages3.m - Optimizing the number of stages necessary for
% phase shifting from [0,360] above a power threshold
% Difference from stages2.m: Different structure

%% Structure
%
%  ____________90deg____________
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
Ptol = 0.8; % Power tolerance
Dm = 0.2; % Modulation depth as %-tage of C

Nmax = 5; %30 % Maximum number of stages considered
Clen = 500; %2000 % Number of values used for C
phase_err = 2; % Degrees away from 360 for total phase shift

%% System Parameters
BL = pi / 2;
Z0 = 50;
Y0 = 1 / Z0;
f = 10e9;
w = 2 * pi * f;

%% Storage
Zrange = linspace(0, 1/10, Clen+1); % 2000 pts for Z_C = [0,10]
Zrange = Zrange(2:end);
Crange = Zrange / w; % Capacitance values

ABCD_a = zeros(2,2,Clen);
ABCD_0 = ABCD_a;
ABCD_b = ABCD_a;
S21_a = zeros(1,Clen);
S21_0 = S21_a;
S21_b = S21_a;

ABCD_aN = zeros(2,2,Clen);
ABCD_0N = ABCD_aN;
ABCD_bN = ABCD_aN;
S21_aN = zeros(1,Clen);
S21_0N = S21_aN;
S21_bN = S21_aN;

phias = zeros(1,Clen); % PHI_A
phibs = zeros(1,Clen); % PHI_B
diffs = zeros(1,Clen); % |PHI_B - PHI_A|

%S21 = 2/(ABCD(1,1)+ABCD(1,2)/Z0+ABCD(2,1)*Z0+ABCD(2,2));

%% N = 1
L = 1 ./ (w^2 * Crange);
b_a = (w^2 .* L .* Crange * (1 - Dm) - 1) ./ (w * L);
b0 = (w^2 .* L .* Crange - 1) ./ (w * L);
b_b = (w^2 .* L .* Crange * (1 + Dm) - 1) ./ (w * L);

for i = 1:Clen
    ABCD_a(:,:,i) = [cos(BL)-sin(BL)*b_a(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b_a(i)*Z0*2*cos(BL)-1j*(Z0*b_a(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b_a(i)*Z0];
    ABCD_0(:,:,i) = [cos(BL)-sin(BL)*b0(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b0(i)*Z0*2*cos(BL)-1j*(Z0*b0(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b0(i)*Z0];
    ABCD_b(:,:,i) = [cos(BL)-sin(BL)*b_b(i)*Z0, 1j*sin(BL)*Z0;
        (1j.*b_b(i)*Z0*2*cos(BL)-1j*(Z0*b_b(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b_b(i)*Z0];
    
    S21_a(i) = 2/(ABCD_a(1,1,i)+ABCD_a(1,2,i)/Z0+ABCD_a(2,1,i)*Z0+ABCD_a(2,2,i));
    S21_0(i) = 2/(ABCD_0(1,1,i)+ABCD_0(1,2,i)/Z0+ABCD_0(2,1,i)*Z0+ABCD_0(2,2,i));
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
fprintf('\n========== N=1 ==========\n');
fprintf('Capacitance: C = %e F\n', Copt);
fprintf('Phase difference: %f deg\n', diffs(idx));
fprintf('Power:\n\tP_A (-Dm) = %4.2f dB\n\tP_B (+Dm) = %4.2f dB\n', db(abs(S21_a(idx))^2,'power'), db(abs(S21_b(idx))^2,'power'));

figure;
hold on;
plot(Crange, phias);
plot(Crange, phibs);
plot(Crange, diffs);
xlabel('Capacitance [F]');
legend('PHI_A', 'PHI_B', 'DIFF');

%% N > 1
for n = 1:Nmax
%     if(abs(diffs(idx) - 360) < phase_err && min(abs(S21_a(idx))^2, abs(S21_b(idx))^2) > Ptol)
%         fprintf('Done\n');
%         break;
%     end
    for i = 1:Clen
        ABCD_aN(:,:,i) = ABCD_a(:,:,i)^n;
        ABCD_0N(:,:,i) = ABCD_0(:,:,i)^n;
        ABCD_bN(:,:,i) = ABCD_b(:,:,i)^n;
        
        S21_aN(i) = 2/(ABCD_aN(1,1,i)+ABCD_aN(1,2,i)/Z0+ABCD_aN(2,1,i)*Z0+ABCD_aN(2,2,i));
        S21_0N(i) = 2/(ABCD_0N(1,1,i)+ABCD_0N(1,2,i)/Z0+ABCD_0N(2,1,i)*Z0+ABCD_0N(2,2,i));
        S21_bN(i) = 2/(ABCD_bN(1,1,i)+ABCD_bN(1,2,i)/Z0+ABCD_bN(2,1,i)*Z0+ABCD_bN(2,2,i));
        
        phias(i) = angle(S21_aN(i));
        phibs(i) = angle(S21_bN(i));
    end
    
    phias = rad2deg(unwrap(phias)); % degrees
    phibs = rad2deg(unwrap(phibs));
    
    if(n == 4)
        phias = phias + 360; % FIX THIS CASE
    end
    
    diffs = abs(phibs - phias);
    
    [~, idx] = min(abs(diffs - 360));
    Copt = Crange(idx);  % Optimal value for C
    
    fprintf('\n========== N=%i ==========\n', n);
    fprintf('Capacitance: C = %e F\n', Copt);
    fprintf('Phase difference: %f deg\n', diffs(idx));
    fprintf('Power:\n\tP_A (-Dm) = %4.2f dB\n\tP_B (+Dm) = %4.2f dB\n', db(abs(S21_a(idx))^2,'power'), db(abs(S21_b(idx))^2,'power'));
    
    figure;
    hold on;
    plot(Crange, phias);
    plot(Crange, phibs);
    plot(Crange, diffs);
    xlabel('Capacitance [F]');
    legend('PHI_A', 'PHI_B', 'DIFF');
    
    %if(abs(diffs(idx) - 360) < phase_err && min(abs(S21_aN(idx))^2, abs(S21_bN(idx))^2) > Ptol)
    %    fprintf('Done\n');
    %    break;
    %end
end

