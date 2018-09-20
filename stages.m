% Debjit Sarkar
% September 2018
% Finding the minimum number of stages to achieve phase shifting

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

%% Setup
syms Y;
BL = pi / 4;
Z0 = 50;
TL_ABCD = [cos(BL), 1j*Z0*sin(BL);
           1j*sin(BL)/Z0, cos(BL)];
shunt_ABCD = [1 0;
              Y 1];

full_ABCD = TL_ABCD * shunt_ABCD * TL_ABCD;
full_S = abcd2s(full_ABCD, Z0);

%% Inputs
Ptol = 1; % Input power tolerance
Dm = 1; % Depth of modulation

% f = 10e9;
% w = 2 * pi * f;




ZC0 = 1;





%% Outputs
