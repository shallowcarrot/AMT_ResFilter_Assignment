% AP_Filtered_Saw.m
%
% Date of creation: 22nd March 2026
%
% Author: Abhiroop Paul
%
% Preamble:
% This script uses a time-varying version of trapezoidal integration
% technique to create a sweep of a resonant filter, applied on a sawtooth
% wave. The pitch of the tone must be specified in Hz as "f0" in the "INPUT
% PARAMETERS" section at the top of the script. The simulation linearly
% ramps the cutoff frequency from 8,000Hz down to f0. It also ramps the
% resonance (feedback coefficient) from 0.99 down to 0 over the course of
% the simulation. This results in a "ZWEEEP!" sound that is characteristic
% of many synth tones that use subtractive synthesis.
% 
%=========================================================================
clear all; close all; clc;


%-------------------------------------------------------
% INPUT PARAMETERS
%-------------------------------------------------------
Fs  = 48000;                        % Sample Rate [Hz]
Tf  = 0.8;                          % Simulation Duration [s]
f0  = 110;                          % Pitch of Sawtooth Wave [Hz]
gain    = 0.8;                      % Amplitude Scaling factor for Saw Wave (gain <1)               


%-------------------------------------------------------
% DERIVED PARAMETERS
%-------------------------------------------------------
T   = 1/Fs;                                         % Time Step / Sample [s]
Nf  = floor(Fs*Tf);                                 % Number of Samples in Simulation

freq_step   = (8000 - f0)/Nf;                       % Frequency Step / Sample
frequencies = f0 + (flip([0:Nf-1]')*freq_step);     % Frequencies Vector
omega       = 2*pi*frequencies;                     % Angular Frequencies Vector

r_step      = 0.99/Nf;                              % r Step / Sample
r           = flip([0:Nf-1]')*r_step;               % r Vector

PD          = f0/Fs;                                % Phase Delta for Saw Wave


%-------------------------------------------------------
% INITIALISATION
%-------------------------------------------------------
% Creating the required matrices
I   = eye(4);                       % Indentity Matrix of Desired Size
A1  = [-1,0,0,0];
A2  = [1,-1,0,0];
A3  = [0,1,-1,0];
A4  = [0,0,1,-1];
A_base   = [A1;A2;A3;A4];           % Base Matrix to build 'A' upon
b   = [1;0;0;0];
c   = [0,0,0,1];

% Initialising State Vectors
x  = zeros(4,1);                    % x Vector
y  = zeros(Nf,1);                   % Output Vector

u   =   zeros(Nf,1);                % Empty input vector
phase       = 0;                    % Empty Phase Variable for Saw Wave

% Frequency vector for plotting
fvec    = [0:(Nf/2)]'*(Fs/Nf);


%-------------------------------------------------------
% MAIN LOOP
%-------------------------------------------------------
tic;
for n = 2 : Nf
    % Creation of Sawtooth Wave
    phase = phase + PD;
    if phase > 1.0
        phase = phase - 1.0;
    end

    u(n) = gain * ((2*phase) - 1);

    % Updating System Matrices
    A = A_base;
    A(1,4) = -4*r(n-1);
    A = omega(n-1)*A;               % System Matrix 'A'

    Jf = I + (T/2)*A;               % Precomputing I + T/2*A
    Jb = I - (T/2)*A;               % Precomputing I + T/2*A
    K = omega(n-1)*(T/2)*b;

    % Trapezoidal Integration Update Step
    x  =   Jb \ ((Jf * x) + K*(u(n-1) + u(n)));
    y(n) = c*x;

end
simTime = toc


%-------------------------------------------------------
% NORMALISATION & SOUND OUTPUT
%-------------------------------------------------------
% Normalisation
y = y / 0.98*max(abs(y));

% Sound Output
soundsc(y,Fs);

