% AP_ResFilter_Trapezoidal.m
%
% Date of creation: 21st March 2026
%
% Author: Abhiroop Paul
%
% Preamble:   
% This script digitally simulates a resonant filter, such as the ones used
% in the Moog VCF synth containing ladder-circuits. The resonant frequency
% "f0" and the feedback coefficient "r" must be specified in the "INPUT
% PARAMETERS" section of the script. Instead of using forward or backward
% Euler methods, the trapezoidal integration technique is used here, which
% is unconditionally stable. Its output is compared with the theoretical
% transfer function obtained from the exact analytical solution through
% overlapping frequency response plots.
%
%=========================================================================
clear all; close all; clc;


%-------------------------------------------------------
% INPUT PARAMETERS
%-------------------------------------------------------
Fs  = 48000;                        % Sample Rate [Hz]
f0  = 1000;                         % Cutoff / Resonance Frequency [Hz]
r   = 0.71;                         % Feedback Coefficient (0<r<1)

mode = 1;                           % Which approximation will be used for output
                                    % (0 = Trapezoidal; 1 = Exact)

% Audio file input
[x_in] = audioread("TestGuitarPhraseMono.WAV");

% Conversion of audio to mono if stereo audio
if size(x_in,2) == 2
    input = 0.5 * (x_in(:,1) + x_in(:,2));
else
    input = x_in;
end



%-------------------------------------------------------
% DERIVED PARAMETERS
%-------------------------------------------------------
om0 = 2*pi*f0;                      % Resonance Angular Frequency [rad/s]
T   = 1/Fs;                         % Time Step / Sample [s]
Nf  = length(input);                 % Number of Samples in Simulation


%-------------------------------------------------------
% INITIALISATION
%-------------------------------------------------------
% Creating the required matrices
I   = eye(4);                       % Indentity Matrix of Desired Size
A1  = [-1,0,0,(-4*r)];
A2  = [1,-1,0,0];
A3  = [0,1,-1,0];
A4  = [0,0,1,-1];
A   = om0*[A1;A2;A3;A4];            % System Matrix 'A'
b   = om0*[1;0;0;0];
c   = [0,0,0,1];

Jf  = I + ((T/2)*A);                % Precomputing I + T/2*A
Jb  = I - ((T/2)*A);                % Precomputing I - T/2*A

K   = (T/2)*b;                      % Precomputing T*b

% Initialising State Vectors
x  = zeros(4,1);                    % x Vector
y  = zeros(Nf,1);                   % Output Vector

% Input Impulse
u   =   [1;zeros(Nf-1,1)];

% Frequency vector for plotting
fvec    = [0:(Nf/2)]'*(Fs/Nf);

%-------------------------------------------------------
% MAIN LOOP
%-------------------------------------------------------
tic;
for n = 1 : Nf-1
    % Trapezoidal Integration Update Step
    x  =   Jb \ ((Jf * x) + K*(u(n) + u(n+1)));
    y(n) = c*x;

end
simTime = toc


%-------------------------------------------------------
% TRANSFER FUNCTIONS
%-------------------------------------------------------
% Exact Transfer Function
H_exact = zeros(length(fvec),1);
for k = 1:length(fvec)
    w = 2*pi*fvec(k);
    
    H_exact(k) = c * ((1i*w*I - A) \ b);
end

% Numerical Transfer Function of Trapezoidal Simulation
H_Trap = fft(y);
H_Trunc = H_Trap(1:(Nf/2)+1);    % Truncated to Nyquist Frequency


%-------------------------------------------------------
% PLOTTING
%-------------------------------------------------------
figure(1);
semilogx(fvec, 20*log10(abs(H_Trunc)), "r", fvec, 20*log10(abs(H_exact)), "g");
xlim([20,20000]);
title("Comparison of Trapezoidal Scheme, and Exact Transfer Functions");
xlabel("Frequency (Hz)");
ylabel("Amplitude (dB)");
legend("Trapezoidal Solution", "Exact Solution");


%-------------------------------------------------------
% CONVOLUTION & SOUND OUTPUT
%-------------------------------------------------------
if mod(mode,3) == 0
    H = H_Trap;
else
    H = H_exact;
end


% Compute the DFT of the input signal and impulse response
LFFT = length(H);
X = fft(input, LFFT);

% Multiplication (Convolution in frequency domain)
Y = H .* X;

% IFFT
out = ifft(Y, LFFT, "symmetric");

% Normalisation
out = out / max(abs(out));

% Sound Output
soundsc(out,Fs);
