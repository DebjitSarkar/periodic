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

%% Setup
f = 10e9;
w = 2 * pi * f;
L0 = 1 / (w^2 * C0);

T_saw = 50; % Number of values in one saw period

Call = linspace(C0 - Dm, C0 + Dm, T_saw);
Yall = 1./(1/(1j*w*L0)+1j*w*Call);

toggle = true; % Used to initialize phase_min/max
phase_min = 0; % in degrees
phase_max = 0;
phase_vec = zeros(1,T_saw); % Stores phases for each cap val
S_vec(:,:,1:T_saw) = zeros(2,2,T_saw); % Stores S-params
ABCD_vec(:,:,1:T_saw) = zeros(2,2,T_saw); % Store ABCD matrices

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
    ABCD_vec(:,:,foo) = full_ABCD;
    
    S_vec(:,:,foo) = abcd2s(full_ABCD, Z0);
    
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

%% Determine minimum number of segments
% Will brute force for several stages first

%% Outputs
% C = 1; % Center capacitance
% N = 1; % # of segments
