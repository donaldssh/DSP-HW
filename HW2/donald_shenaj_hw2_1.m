clearvars
close all
clc

%% Read the audio file
[y, Fp] = audioread('sounds/hw2_vuvuzela.wav');
% so that the samplig period is
T = 1/Fp; % [s]

% show ORIGINAL signal in time and frequency
Ny = length(y); % length
ty = T*(0:Ny-1); % time samples
Y = T*fft(y); % fft
fy = (0:Ny-1)/(T*Ny); % frequency samples
figure
subplot(2,1,1) % show a portion of time
plot(ty,y); grid;  xlabel('time [s]'); 
title('original audio signal in time')
subplot(2,1,2) % show frequency content in dB scale
plot(fy/1e3,20*log10(abs(Y))); grid; 
xlabel('frequency [kHz]'); title('original audio signal in frequency')


%% 1 - Vuvuzela harmonics

f0 = 235; % first harmonic

N = 5; % number of harmonics i'll consider in the first approach
f = f0:f0:f0*5;  % notch frequencies

%% 2 - IIR NOTCH DESIGN M = N = 10

% extract theta for poles and zeros corresponding to the frequency we want to remove
theta = f*2*pi*T; 

% find the poles and zeros by exploiting the theory seen in class
% the values are explained in the report
r = 0.997;   % take a value near one

% first half poles and zeros
p1_5 = r*exp(1i*theta);   
z1_5 = exp(1i*theta);

 % second half poles and zeors
p6_10 = conj(p1_5);  
z6_10 = conj(z1_5);

% build the poles and zeros vector
p = [p1_5 p6_10]';
z = [z1_5 z6_10]';
k = 1;  % constant factor

% show poles and zeros in the zplane
figure
zplane(z,p)

%% 3 - Frequency response

% extract the IIR filter coefficients from poles, zeros and constant
[b,a] = zp2tf(z,p,k);

% show the freq response
fvtool(b,a)

%% 4 - Filter the audio signal

% filter the signal
z = filter(b,a,y);

% write the filtered audio
audiowrite('donald_shenaj_hw2_1_1.wav',z,Fp);

%% Alternative solution - matlab iircomb tool

q = 10;  % example value for Q factor in matlab doc   
harmonics = f0:f0:Fp/2;  % all harmonics of the vuvuzela sound
bw = (f0/(Fp/2))/q;
N = round((Fp)/f0);
[b, a] = iircomb(N,bw,'notch');

% show the freq response
fvtool(b,a);

% filter the audio
z2 = filter(b,a,y);

audiowrite('donald_shenaj_hw2_1_2.wav',z2,Fp);

%% Combined plot of filtered audio

% find some values for the plot
Nz = length(z);
tz = T*(0:Nz-1); % time samples
Z = T*fft(z); % fft
Z2 = T*fft(z2); % fft
fz = (0:Nz-1)/(T*Nz); % frequency samples

figure
subplot(2,1,1)
plot(tz,z); grid; xlabel('time [s]'); 
title('filtered audio signal in time')
subplot(2,1,2)
plot(fz/1e3,20*log10(abs(Z))); grid; 
xlabel('frequency [kHz]'); title('filtered audio signal in frequency')
hold on
subplot(2,1,1)
plot(tz,z2); grid; xlabel('time [s]'); 
title('filtered audio signal in time')
subplot(2,1,2)
plot(fz/1e3,20*log10(abs(Z2))); grid; 
xlabel('frequency [kHz]'); title('filtered audio signal in frequency')
hold off
legend('notch 5 harmonics', 'iircomb: notch 102 harmonics')



