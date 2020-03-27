clearvars
close all
clc

%% - 1 - read your audio file

% parameters defined in the assignment
T1 = 75e-6;
T2 = 318e-6;
T3 = 3180e-6;
A = (T2 - T1)/(T3 - T1);

% read the corrisponding WAVE file (*.wav)
% y array containing sound samples y(n), n=1:length(y)
% Fp sampling frequency in Hz
[y, Fp] = audioread('sounds/hw1_peter_and_the_wolf_9.wav');
% so that the samplig period is
T = 1/Fp; % [s]


%% - 2 - filter design by a windowing technique which uses h(t)
N = 100;
h0 = zeros(1,N+1);  %h ha un elemento in piu'

% sampling h(t) with step T, N+1 values
% h0(t) should be 0 for the first half, t<0
t = T*(0:N/2);
h0 = zeros(N+1,1);
h0(N/2+1:N+1) = (A/T1)*exp(-t/T1) + ((1-A)/T3)*exp(-t/T3);
t = T*(-N/2:N/2);  %complete the time grid for the plot


%% - 3 - freq and time response of the filter
% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor
figure
subplot(2,1,1)
stem(t,h0); grid;
title('Windowing filter h(t) - time response')
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

%% - 4 - apply the filter to the audio sample

% filtering the audio sample through the convolution
z = T*conv(y,h0);
audiowrite('donald_shenaj_hw1.3.wav',z,Fp);



%% - 5 - ADVANCED - minimax linprog

%I've determined the num and dem of both Real and Imm components 
% as functions of f, so i can use them easily below (not necessary but
% useful)
num_re = @(f) 1-4*((pi*f).^2)*(T1*T3-T2*(T1+T3));
num_imm = @(f) 2*pi*f.*(T2-T1-T3-4*T1*T2*T3*((pi*f).^2));
den = @(f) 1+16*((pi*f).^4)*((T1*T3)^2)-8*T1*T3*((pi*f).^2)+4*((pi*f*(T1+T3)).^2);
    
% type I filter design ----------------------------------------------------
N = 200;
F = Fp/(N/3)/32; %frequency step
f = [0:F:Fp/2].'; % frequency samples, column vector

dr = num_re(f)./den(f); % ideal filter shape
V = T*ones(size(f)); % cosines matrix
for n = 1:N/2
    V = [V,2*T*cos(2*pi*f*n*T)];
end
w = abs(dr); %weighting factor proportional to the absolute value of the target

%usual steps for linear programming
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-dr;dr];
x = linprog(g,A,b); 
h0r = [x(N/2+1:-1:2);x(1:N/2+1)];
t = T*(-N/2:N/2);
[H0r,ff] = freqz(h0r,1,8*(N+1),Fp);
H0r = T*H0r; % normalization factor

figure
subplot(2,1,1)
stem(t,h0r); grid;
title('Real part of the filter - type I - time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
yyaxis left;
plot(ff,20*log10(abs(H0r))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0r)));
title('frequency response')
ylabel('phase [degrees]')
xlabel('frequency [Hz]')

% type III filter design --------------------------------------------------

di = num_imm(f)./den(f); % ideal filter shape
%built matrix as done before
V = zeros(size(f));
for n = 1:N/2
    V = [V,2*T*sin(2*pi*f*n*T)];
end

w = abs(di); %weighting factor proportional to the absolute value of the target

g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-di;di];
x = linprog(g,A,b); 
h0i = [x(N/2+1:-1:2);-x(1:N/2+1)];
t = T*(-N/2:N/2);
[H0i,ff] = freqz(h0i,1,8*(N+1),Fp);
H0i = T*H0i; % normalization factor


figure
subplot(2,1,1)
stem(t,h0i); grid;
title('Immaginary part of the filter - type III - time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
yyaxis left;
plot(ff,20*log10(abs(H0i))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0i)));  
title('frequency response')
ylabel('phase [degrees]')
xlabel('frequency [Hz]')

h0t = h0r+h0i;
%filtering the audio with the filter h = h0r + h0i
z2 = T*conv(y,h0t);
% WRIING A SECOND AUDIO FOR THIS EXERCISE, USING THE MINIMAX TOTAL FILTER
audiowrite('donald_shenaj_hw1.3_2.wav',z2,Fp);  %minimax filtering

[H0_200,ff] = freqz(h0t,1,8*(N+1),Fp);
H0_200 = T*H0_200; % normalization factor

figure
subplot(2,1,1)
stem(t,h0t); grid;
title('Total filter N=200 - time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
plot(ff,20*log10(abs(H0_200))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff, unwrap(angle(H0_200))); 
title('frequency response')
ylabel('phase [degrees]')
xlabel('frequency [Hz]')


%% ----------- COPY/PAST OF THE CODE ABOVE, USING N=100
N = 100;
V = T*ones(size(f)); % cosines matrix
for n = 1:N/2
    V = [V,2*T*cos(2*pi*f*n*T)];
end
w = abs(dr); %weighting factor proportional to the absolute value of the target
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-dr;dr];
x = linprog(g,A,b); 
h0r = [x(N/2+1:-1:2);x(1:N/2+1)];
t = T*(-N/2:N/2);
[H0r,ff_100] = freqz(h0r,1,8*(N+1),Fp);
H0r = T*H0r; % normalization factor
% type III filter design
V = zeros(size(f));
for n = 1:N/2
    V = [V,2*T*sin(2*pi*f*n*T)];
end
w = abs(di); %weighting factor proportional to the absolute value of the target
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-di;di];
x = linprog(g,A,b); 
h0i = [x(N/2+1:-1:2);-x(1:N/2+1)];
t = T*(-N/2:N/2);
[H0i,ff_100] = freqz(h0i,1,8*(N+1),Fp);
H0i = T*H0i; % normalization factor
[H0_100,ff_100] = freqz(h0t,1,8*(N+1),Fp);
H0_100 = T*H0_100; % normalization factor
h0t = h0r+h0i;
figure
subplot(2,1,1)
stem(t,h0t); grid;
title('Total filter N=100 - time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
plot(ff_100,20*log10(abs(H0_100))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff_100, unwrap(angle(H0_100))); 
title('frequency response')
ylabel('phase [degrees]')
xlabel('frequency [Hz]')

%% ----------- COPY/PAST OF THE CODE ABOVE, USING N=50
N = 50;
V = T*ones(size(f)); % cosines matrix
for n = 1:N/2
    V = [V,2*T*cos(2*pi*f*n*T)];
end
w = abs(dr); %weighting factor proportional to the absolute value of the target
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-dr;dr];
x = linprog(g,A,b); 
h0r = [x(N/2+1:-1:2);x(1:N/2+1)];
t = T*(-N/2:N/2);
[H0r,ff_50] = freqz(h0r,1,8*(N+1),Fp);
H0r = T*H0r; % normalization factor
% type III filter design
V = zeros(size(f));
for n = 1:N/2
    V = [V,2*T*sin(2*pi*f*n*T)];
end
w = abs(di); %weighting factor proportional to the absolute value of the target
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-di;di];
x = linprog(g,A,b); 
h0i = [x(N/2+1:-1:2);-x(1:N/2+1)];
t = T*(-N/2:N/2);
[H0i,ff_50] = freqz(h0i,1,8*(N+1),Fp);
H0i = T*H0i; % normalization factor
[H0_50,ff_50] = freqz(h0t,1,8*(N+1),Fp);
H0_50 = T*H0_50; % normalization factor
h0t = h0r+h0i;
figure
subplot(2,1,1)
stem(t,h0t); grid;
title('Total filter N=50 - time response')
ylabel('amplitude []')
xlabel('time [s]')
subplot(2,1,2)
plot(ff_50,20*log10(abs(H0_50))); grid;
ylabel('amplitude [dB]')
yyaxis right;
plot(ff_50, unwrap(angle(H0_50))); 
title('frequency response')
ylabel('phase [degrees]')
xlabel('frequency [Hz]')



%% COMPARING THE FITTING OF N=50,100,200
figure
plot( f,abs(dr+1i*di), ff_50,abs(H0_50), ff_100,abs(H0_100),ff,abs(H0_200));
title('Approximation of H(f)');
ylabel('amplitude []');   %since is not given the unity of measurement -> []
xlabel('frequency [Hz]');
grid;  
legend('Desired', 'N=50', 'N=100', 'N=200');



