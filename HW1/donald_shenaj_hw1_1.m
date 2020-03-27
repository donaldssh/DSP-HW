clearvars
close all
clc

%% 1 - READING THE AUDIO FILE
% read the WAVE file (*.wav)
% y array containing sound samples y(n), n=1:length(y)
% Fp sampling frequency in Hz
[y, Fp] = audioread('sounds/hw1_forest_sound.wav');
% so that the samplig period is
T = 1/Fp; % [s]

% not mandatory: show ORIGINAL signal in time and frequency
Ny = length(y); % length
ty = T*(0:Ny-1); % time samples
Y = T*fft(y); % fft
fy = (0:Ny-1)/(T*Ny); % frequency samples
figure(1)
subplot(2,1,1) % show a portion of time
plot(ty,y); grid; xlim([0.1 0.2]); xlabel('time [s]'); 
title('original audio signal in time')
subplot(2,1,2) % show frequency content in dB scale
plot(fy/1e3,20*log10(abs(Y))); grid; xlim([0 Fp/2e3]); ylim([-150 -20]) 
xlabel('frequency [kHz]'); title('original audio signal in frequency')



%% 2 - DESIGNING OF A SUITABLE TYPE I FILTER

% REMEZ 

% pass band, cuckoo calls freq
fp1 = 640; % [Hz] from
fp2 = 1280; % [Hz] to

% attenuation bands, freq to remove
fs1 = fp1 - 160; % [Hz] from 0 to
fs2 = fp2 + 160; % [Hz] from ... to Fp/2
% required limits, common values used before
err_lim = [1e-3 1e-2 1e-3];

% Filter desing with REMEZ
[N,Fo,Ao,W] = firpmord([fs1 fp1 fp2 fs2],[0 1 0],err_lim,Fp); 
disp(['firpmord suggests a filter of order ' num2str(N)])
N = 100;  %ill use N=100
h0 = firpm(N,Fo,Ao,W)/T; 
t = T*(-N/2:N/2);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor


%% 3 - ILLUSTRATE THE RESULTING FILTER IN TIME AND FREQUENCY

% TIME and FREQUENCY response - REMEZ
figure(2)
subplot(2,1,1)
stem(t,h0); grid; 
title('FIR - band pass - Remez approach - time domain')
xlabel('time [s]')
ylabel('amplitude []')
subplot(2,1,2)
xlabel('frequency [Hz]')
yyaxis left
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]); ylim([-80, 5]); 
title('frequency domain')
hold on; plot([1,1]*fp1,ylim,'r--'); plot([1,1]*fp2,ylim,'r--'); 
plot([1,1]*fs1,ylim,'r--'); plot([1,1]*fs2,ylim,'r--'); hold off;
ylabel('amplitude [dB]')
yyaxis right
plot(ff, unwrap(angle(H0)));   
ylabel('phase [degrees]')


%-----------------------------------------

% LINEAR PROGRAMMING APPROACH (from lab4)

% frequencies samples of interest
F = Fp/(N+1)/32; % min 32 samles per cosine period
f = [0:F:fs1, fp1:fp2, fs2:F:Fp/2].'; % frequency samples, column vector
% build matrices
d = (f>=fp1).*(f<=fp2); % ideal filter shape
w = err_lim(1)*(f<=fs1) + err_lim(2)*(f>=fp1).*(f<=fp2) ...
    + err_lim(3)*(f>=fs2); % weighting function
V = T*ones(size(f)); % cosines matrix
for n = 1:N/2
    V = [V,2*T*cos(2*pi*f*n*T)];
end
% linear programming solution
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-d;d];
x = linprog(g,A,b); % solve the linear program

% define filter
h0 = [x(N/2+1:-1:2);x(1:N/2+1)];
t = T*(-N/2:N/2);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor

% show results
figure(3)
subplot(2,1,1)
stem(t,h0); grid; 
title('FIR - band pass - linear programming approach - time domain')
xlabel('time [s]')
ylabel('amplitude []')
subplot(2,1,2)
xlabel('frequency [Hz]')
yyaxis left
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]); ylim([-80, 5]); 
title('frequency domain')
hold on; plot([1,1]*fp1,ylim,'r--'); plot([1,1]*fp2,ylim,'r--'); 
plot([1,1]*fs1,ylim,'r--'); plot([1,1]*fs2,ylim,'r--'); hold off;
ylabel('amplitude [dB]')
yyaxis right
plot(ff, unwrap(angle(H0)));   
ylabel('phase [degrees]')




%% 4 - FILTER THE AUDIO SIGNAL 

% listen to the original content
%disp('press any key to listen to the ORIGINAL audio signal')
%pause
%sound(y(1:5*Fp), Fp); % play the first 5 seconds 
%pause(5)

% filter the signal through convolution, linear programming coefficients 
z = T*conv(y,h0);

% show the FILTERED signal in time and frequency
Nz = length(z);
tz = T*(0:Nz-1); % time samples
Z = T*fft(z); % fft
fz = (0:Nz-1)/(T*Nz); % frequency samples
figure(4)
subplot(2,1,1)
plot(tz,z); grid; xlim([0.1 0.2]); xlabel('time [s]'); 
title('filtered audio signal in time')
subplot(2,1,2)
plot(fz/1e3,20*log10(abs(Z))); grid; xlim([0 Fp/2e3]); ylim([-150 -20])
xlabel('frequency [kHz]'); title('filtered audio signal in frequency')

% listen to the filtered content
%disp('press any key to listen to the FILTERED audio signal')
%pause
%sound(z(1:5*Fp), Fp); % play the first 5 seconds 

% write the audio filtered 
audiowrite('donald_shenaj_hw1.1.wav',z,Fp);


%% 5 - ADVANCED 1.5 - TYPE II filter

% for a type II filter i've choosen N to be 101 (it must be odd)
N = 101 ;

% frequencies samples of interest
F = Fp/N/2/32; 
f = [0:F:fs1, fp1:fp2, fs2:F:Fp/2].'; % frequency samples, column vector
err_lim = [1e-3 1e-2 1e-3];

% build matrices for the linear programming solution
d = (f>=fp1).*(f<=fp2); % ideal filter shape
w = err_lim(1)*(f<=fs1) + err_lim(2)*(f>=fp1).*(f<=fp2) ...
    + err_lim(3)*(f>=fs2); % weighting function

V = [];
for n=(N-1)/2:-1:0
    V = [V,2*T*cos(2*pi*f*(N/2 -n)*T)];
end

% linear programming solution
g = [zeros((N-1)/2 +1,1);1];
A = [-V, -w; V, -w];
b = [-d;d];
x = linprog(g,A,b); % solve the linear program

% define filter coefficients
h0 = [x((N-1)/2+1:-1:1); x(1:(N-1)/2+1)];
t = T*(-(N-1)/2-0.5:(N-1)/2+0.5);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor

% show results
figure(5)
subplot(2,1,1)
stem(t,h0); grid; 
title('Type II FIR - band pass - time domain')
xlabel('time [s]')
ylabel('amplitude []')
subplot(2,1,2)
title('frequency domain')
xlabel('frequency [Hz]')
yyaxis left
plot(ff,20*log10(abs(H0))); grid; xlim([0 Fp/2]); ylim([-80, 5]); 
ylabel('amplitude [dB]')
hold on; plot([1,1]*fp1,ylim,'r--'); plot([1,1]*fp2,ylim,'r--'); 
plot([1,1]*fs1,ylim,'r--'); plot([1,1]*fs2,ylim,'r--'); hold off;
yyaxis right;
plot(ff, unwrap(angle(H0)));   
ylabel('phase [degrees]')



