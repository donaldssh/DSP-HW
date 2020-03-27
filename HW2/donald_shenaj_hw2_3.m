clearvars
close all
clc

%% read the audio
% read the corrisponding WAVE file (*.wav)
% x array containing sound samples x(n), n=1:length(x)
% Fp sampling frequency in Hz
[x, F1] = audioread('sounds/hw2_buongiorno.wav');
% so that the samplig period is
T1 = 1/F1; % [s]


% show ORIGINAL signal in time and frequency
Nx = length(x); % length
tx = T1*(0:Nx-1); % time samples
X = T1*fft(x); % fft
fx = (0:Nx-1)/(T1*Nx); % frequency samples
figure
subplot(2,1,1) % show a portion of time
plot(tx,x); grid; xlabel('time [s]'); 
title('original audio signal in time')
subplot(2,1,2) % show frequency content in dB scale
plot(fx/1e3,20*log10(abs(X))); grid;
xlabel('frequency [kHz]'); title('original audio signal in frequency')


%% 1 - Design a rate-conversion algorithm 

% final freq
F2 = 48000;  %[Hz]
T2 = 1/F2;   

L = 2^5 * 5;
M = 7^2 *3;

% intermediate freq
F_p = F1*L;   % 7056000 Hz
T_p = 1/F_p;

%% Interpolation

% number of samples after interpolation
n_samples_interp = length(x)*L - (L-1);

% interpolation, v is the resulting signal
v = zeros(n_samples_interp,1);
j=1;
for i=1:L:n_samples_interp
    v(i) = x(j);
    j = j+1;
end

%% 2 - FIR filter design
% ~ 1 minute time to process

f0 = F_p/(L*2);   % L>M
al = 0.1; % transition bandwidth in percentage
fp = f0*(1-al); % pass band upper limit
fs = f0; % stop band lower limit
err_lim = [1e-3 1e-4];
[N,Fo,Ao,W] = firpmord([fp fs],[1 0],err_lim, F_p);
disp('Filter N = ' +string(N));
h0 = firpm(N,Fo,Ao,W)*F_p;
t = T_p*(-N/2:N/2);

% filter signal
z = T_p*conv(v,h0);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),F_p);
H0 = T_p*H0; % normalization factor

% plots
figure
subplot(2,1,1)
stem(t,h0); grid;
title('time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
yyaxis left;
plot(ff,20*log10(abs(H0))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0)));
ylabel('phase [degrees]')
xlabel('frequency [Hz]')
title('frequency response')

% zoom plot
figure
yyaxis left;
plot(ff,20*log10(abs(H0))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0)));
ylabel('phase [degrees]')
xlabel('frequency [Hz]')
title('frequency response')
xlim([0 5e4]);
hold on
plot([1,1]*fs,ylim,'r--'); plot([1,1]*fp,ylim,'r--'); 
hold off;


%% Decimation

% sample signal
y = z(1:M:end);
y = y/max(y);
Ny = length(y);

% write the rate converted audio
y = y/max(abs(y));
audiowrite('donald_shenaj_hw2_3_1.wav',y,F2);


% show the converted signal in time and frequency
ty = T2*(0:Ny-1); % time samples
Y = T2*fft(y); % fft
fy = (0:Ny-1)/(T2*Ny); % frequency samples
figure
subplot(2,1,1) % show a portion of time
plot(ty,y); grid;  xlabel('time [s]'); 
title('Converted audio signal in time')
subplot(2,1,2) % show frequency content in dB scale
plot(fy/1e3,20*log10(abs(Y))); grid
xlabel('frequency [kHz]'); title('Converted audio signal in frequency')


%% (advanced) cascade of rate-conversion blocks
disp('Multistage solution');

% Rate conversion 1

L = 4;
n_samples_interp = length(x)*L - (L-1);
% intermediate signal
v1 = zeros(n_samples_interp,1);
j=1;
for i=1:L:n_samples_interp
    v1(i) = x(j);
    j = j+1;
end

% filter design
F_p1 = F1*L;
T_p1 = 1/F_p1;
M = 3;

%f0 = F_p1/(L*2);  %L>M  same expression as below
f01 = F1/2;  % because F1 < F_y1 

al = 0.1; % transition bandwidth in percentage
fp1 = f01*(1-al); % pass band upper limit -> for all the filters
fs = f01; % stop band lower limit
err_lim = [1e-3 1e-4];
[N1,Fo,Ao,W] = firpmord([fp1 fs],[1 0],err_lim, F_p1);
h0 = firpm(N1,Fo,Ao,W)*F_p1;
t = T_p1*(-N1/2:N1/2);
disp('Filter 1 N = ' +string(N1));


%filter plot
% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N1+1),F_p1);
H0 = T_p*H0; % normalization factor

% plots
figure
subplot(2,1,1)
stem(t,h0); grid;
title('time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
yyaxis left;
plot(ff,20*log10(abs(H0))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0)));
ylabel('phase [degrees]')
xlabel('frequency [Hz]')
title('frequency response 1')


% filter signal
z1 = T_p*conv(v1,h0);
Nz = length(z1);

% sample signal
y1 = z1(1:M:end);
y1 = y1/max(y1);
Ny1 = length(y1);
F_y1 = F1*L/M;


%% Rate conversion 2

L = 8;
n_samples_interp = length(y1)*L - (L-1);
v2 = zeros(n_samples_interp,1);
j=1;
for i=1:L:n_samples_interp
    v2(i) = y1(j);
    j = j+1;
end

% filter design
F_p2 = F_y1*L;
T_p2 = 1/F_p2;
M = 7;

%f02 = F_p2/(L*2);   % L>M   same expression as below
f02 = F_y1/2; % because F_y1 < F_y2

al = 0.1; % transition bandwidth in percentage
%fp = f0*(1-al); % pass band upper limit
fp = fp1;   % pass band upper limit 
fs = f02; % stop band lower limit
err_lim = [1e-3 1e-4];
[N2,Fo,Ao,W] = firpmord([fp fs],[1 0],err_lim, F_p2);
h0 = firpm(N2,Fo,Ao,W)*F_p2;
t = T_p2*(-N2/2:N2/2);
disp('Filter 2 N = ' +string(N2));

%filter plot
% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N2+1),F_p2);
H0 = T_p*H0; % normalization factor

% plots
figure
subplot(2,1,1)
stem(t,h0); grid;
title('time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
yyaxis left;
plot(ff,20*log10(abs(H0))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0)));
ylabel('phase [degrees]')
xlabel('frequency [Hz]')
title('frequency response 2')

% filter signal
z2 = T_p*conv(v2,h0);
Nz2 = length(z2);

% sample signal
y2 = z2(1:M:end);
y2 = y2/max(y2);
Ny2 = length(y2);
F_y2 = F_y1*L/M;

%% Rate conversion 3

L = 5;
n_samples_interp = length(y2)*L - (L-1);
v3 = zeros(n_samples_interp,1);
j=1;
for i=1:L:n_samples_interp
    v3(i) = y2(j);
    j = j+1;
end

% filter design
F_p3 = F_y2*L;
T_p3 = 1/F_p3;
M = 7;

f03 = F_p3/(M*2);   % F2/2 
al = 0.1; % transition bandwidth in percentage
%fp = f0*(1-al); % pass band upper limit
fp = fp1; 
fs = f03; % stop band lower limit
err_lim = [1e-3 1e-4];
[NN,Fo,Ao,W] = firpmord([fp fs],[1 0],err_lim, F_p3);
h0 = firpm(NN,Fo,Ao,W)*F_p3;
t = T_p3*(-NN/2:NN/2);
N3=NN;
disp('Filter 3 N = ' +string(NN));

%filter plot
% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N3+1),F_p3);
H0 = T_p*H0; % normalization factor

% plots
figure
subplot(2,1,1)
stem(t,h0); grid;
title('time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
yyaxis left;
plot(ff,20*log10(abs(H0))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0)));
ylabel('phase [degrees]')
xlabel('frequency [Hz]')
title('frequency response 3')


% filter signal
z3 = T_p*conv(v3,h0);
Nz2 = length(z3);

% sample signal
y3 = z3(1:M:end);
y3 = y3/max(y3);
Ny3 = length(y3);
F_y3 = F_y2*L/M;  % == F2

% write the audio
y3 = y3/max(abs(y3));
audiowrite('donald_shenaj_hw2_3_2.wav',y3,F2);

% show the converted signal in time and frequency
ty3 = T2*(0:Ny3-1); % time samples
Y3 = T2*fft(y3); % fft
fy3 = (0:Ny3-1)/(T2*Ny3); % frequency samples
figure
subplot(2,1,1) % show a portion of time
plot(ty3,y3); grid;  xlabel('time [s]'); 
title('Multistage converted audio signal  in time')
subplot(2,1,2) % show frequency content in dB scale
plot(fy3/1e3,20*log10(abs(Y3))); grid; 
xlabel('frequency [kHz]'); title('Multistage converted audio signal in frequency')

% disp the sum pf orders of the 3 rate-conv modules
disp('Total N = ' +string(N1+N2+N3));

