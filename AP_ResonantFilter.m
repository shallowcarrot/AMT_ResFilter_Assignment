% AP_ResonantFilter.m
%
% Date of creation: 17th March 2026
% 
% Author: Abhiroop Paul
%
% Preamble:   
% This script digitally simulates a resonant filter, such as the ones used
% in the Moog VCF synth containing ladder-circuits. The resonant frequency
% "f0" and the feedback coefficient "r" must be specified in the "INPUT
% PARAMETERS" section of the script. The simulation uses the Forward Euler
% (FE)and Backward Euler (BE) numerical methods to create the filter, and compares
% their frequency response with the theoretical transfer function obtained
% from an exact analytical solution. This comparison is visualised by
% plotting the three transfer functions (FE, BE, and Exact) on the same
% graph in different colour traces.
%
%=========================================================================
clear all; close all; clc;


%-------------------------------------------------------
% INPUT PARAMETERS
%-------------------------------------------------------
Fs  = 48000;                        % Sample Rate [Hz]
f0  = 1000;                         % Cutoff / Resonance Frequency [Hz]
r   = 0.71;                         % Feedback Coefficient (0<r<1)

mode = 2;                           % Which approximation will be used for output
                                    % (0 = FE; 1 = BE; 2 = Exact)

% Audio file input
[x_in] = audioread("TestGuitarPhraseMono.WAV");

% Conversion of audio to mono if stereo audio
if size(x_in,2) == 2
    x = 0.5 * (x_in(:,1) + x_in(:,2));
else
    x = x_in;
end


%-------------------------------------------------------
% DERIVED PARAMETERS
%-------------------------------------------------------
om0 = 2*pi*f0;                      % Resonance Angular Frequency [rad/s]
T   = 1/Fs;                         % Time Step / Sample [s]
Nf  = length(x);                    % Number of Samples in Simulation


%-------------------------------------------------------
% INITIALISATION
%-------------------------------------------------------
% Creating the required matrices
I   = eye(4);                       % Identity Matrix of Desired Size
A1  = [-1,0,0,(-4*r)];
A2  = [1,-1,0,0];
A3  = [0,1,-1,0];
A4  = [0,0,1,-1];
A   = om0*[A1;A2;A3;A4];            % System Matrix 'A'
b   = om0*[1;0;0;0];
c   = [0,0,0,1];

Jf  = I + (T*A);                    % Precomputing I + T*A
Jb  = I - (T*A);                    % Precomputing I - T*A

K   = T*b;                          % Precomputing T*b

% Finding Eigenvalues of I + T*A to check stability
z   = eig(Jf);
z_mod = max(abs(z));

% Vectors to hold x values
xf  = zeros(4,1);                   % x vector for FE
xb  = zeros(4,1);                   % x vector for BE

% Output vectors
yf  = zeros(Nf,1);                  % FE output
yb  = zeros(Nf,1);                  % BE output

% Input Impulse
u   =   [1;zeros(Nf-1,1)];

% Frequency vector for plotting
fvec    = [0:(Nf/2)]'*(Fs/Nf);      


%-------------------------------------------------------
% STABILITY CONDITION
%-------------------------------------------------------
assert(z_mod <= 1, "System parameters are unstable. Eigenvalues of I + TA exceed 1");


%-------------------------------------------------------
% MAIN LOOP
%-------------------------------------------------------
tic;
for n = 1 : Nf-1
    % Forward Euler Update
    xf  =   (Jf * xf) + (K*u(n));
    yf(n) = c*xf;

    % Backward Euler Update
    xb  =   Jb \ (xb + K*u(n));
    yb(n) = c*xb;

end
simTime = toc;


%-------------------------------------------------------
% EXACT CONTINUOUS-TIME FREQUENCY RESPONSE
%-------------------------------------------------------
H_exact = zeros(length(fvec),1);

for k = 1:length(fvec)
    w = 2*pi*fvec(k);
    
    H_exact(k) = c * ((1i*w*I - A) \ b);
end


%-------------------------------------------------------
% DISCRETE TRANSFER FUNCTIONS
%-------------------------------------------------------
H_f = fft(yf);                  % FE Transfer Function
H_b = fft(yb);                  % BE Transfer Function

N   = length(H_f);              % Length of FE transfer function

% Transfer Functions Truncated to the Nyquist Frequency
H_f_Nyq = H_f(1:(Nf/2)+1);
H_b_Nyq = H_b(1:(Nf/2)+1);


%-------------------------------------------------------
% PLOTTING
%-------------------------------------------------------
figure(1);
semilogx(fvec, 20*log10(abs(H_f_Nyq)),"r",fvec, 20*log10(abs(H_b_Nyq)),"b",fvec, 20*log10(abs(H_exact)),"g");
xlim([20,20000]);
title("Comparison of Forward Euler, Backward Euler, and Exact Transfer Functions");
xlabel("Frequency (Hz)");
ylabel("Amplitude (dB)");
legend("Forward Euler", "Backward Euler", "Exact Solution");


%-------------------------------------------------------
% CONVOLUTION & SOUND OUTPUT
%-------------------------------------------------------
if mod(mode,3) == 0
    H   =   H_f;
elseif mod(mode,3) == 1
    H   =   H_b;
else
    H = H_exact;
end

LFFT = length(H);

% Compute the DFT of the input signal and impulse response
X = fft(x, LFFT);

% Multiplication (Convolution in frequency domain)
Y = H .* X;

% IFFT
y = ifft(Y, LFFT, "symmetric");

% Sound Output
soundsc(y,Fs);

