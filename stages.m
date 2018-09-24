% Debjit Sarkar
% September 2018
% Optimizing the number of stages need for 2PI phase shifting

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
Ptol = 1; % Input power tolerance
Dm = 100e-12; % Depth of modulation for saw, +/- to capacitance
C0 = 101e-12; % For example

%Y_C = 1j * 50; % Admittance of tunable capacitor
%Y_L = -1j * 50; % Admittance of fixed inductor

%% Setup
f = 10e9;
w = 2 * pi * f;
L0 = 1 / (w^2 * C0);

T_saw = 50; % Number of values in one saw period

Call = linspace(C0 - Dm, C0 + Dm, T_saw);
Yall = 1./(1/(1j*w*L0)+1j*w*Call);

%Y_Cmin = Y_C - 1j * Dm;
%Y_Cmax = Y_C + 1j * Dm;

%Y_Call = linspace(Y_Cmin, Y_Cmax, T_saw);

toggle = true; % Used to initialize phase_min/max
phase_min = 0; % in degrees
phase_max = 0;
phase_vec = zeros(1,T_saw); % Stores phases for each cap val
S_vec(:,:,1:T_saw) = zeros(2,2,T_saw); % Stores S-params

%% Finding phases for single segment

for foo = 1:T_saw
    %Y = Y_Call(foo) + Y_L;
    Y = Yall(foo);
    
    BL = pi / 4;
    Z0 = 50;
    TL_ABCD = [cos(BL), 1j*Z0*sin(BL);
        1j*sin(BL)/Z0, cos(BL)];
    shunt_ABCD = [1 0; Y 1];
    
    full_ABCD = TL_ABCD * shunt_ABCD * TL_ABCD; % ABCD for 1 segment
    
    %full_S = abcd2s(full_ABCD, Z0); doesn't work for symbolic vars
    A = full_ABCD(1,1);
    B = full_ABCD(1,2);
    C = full_ABCD(2,1);
    D = full_ABCD(2,2);
    
    denom = A + B / Z0 + C * Z0 + D; % ABCD to S
    S_vec(1,1,foo) = (A + B / Z0 - C * Z0 - D) / denom;
    S_vec(1,2,foo) = 2 * (A * D - B * C) / denom;
    S_vec(2,1,foo) = 2 / denom;
    S_vec(2,2,foo) = (-A + B / Z0 - C * Z0 + D) / denom;
    %S11 = (A + B / Z0 - C * Z0 - D) / denom;
    %S12 = 2 * (A * D - B * C) / denom;
    %S21 = 2 / denom;
    %S22 = (-A + B / Z0 - C * Z0 + D) / denom;
    
    %phase = rad2deg(wrapTo2Pi(angle(S21)));
    %phase = rad2deg(angle(S21));
    phase = rad2deg(angle(S_vec(2,1,foo)));
    
    phase_vec(foo) = phase; % Write phase for plotting
    
    % Update phase min/max
    if(toggle) % Initialize
        phase_min = phase;
        phase_max = phase;
        toggle = false;
    elseif(phase < phase_min)
        phase_min = phase;
    elseif(phase > phase_max)
        phase_max = phase;
    end
    
end

%% Display one segment
figure;
hold on;
scatter(Call, phase_vec);
title('Phase range vs capacitance for one segment');
xlabel('C [F]');
ylabel('Phase [\circ]');

%scatter(imag(Y_Call), phase_vec);
%title('Phase range vs Y_C');
%xlabel('Y_C [1/\omega]');

%% Determine minimum number of segments


%% Outputs
% C = 1; % Center capacitance
% N = 1; % # of segments
