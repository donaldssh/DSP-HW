clearvars
close all
clc

%% 1 - values from HW assignment
Fp = 8e3;  %sampling freq in Hz
B = 3e3;   
alpha = 0.05;
T = 1/Fp;   %sampling period in s

%% 2 - filter design, type III
N = 100;
% ill use type 3 filter, because H has odd simmetry, 
% which is needed for our derivative filter

f1 = (1 - alpha)*B;  %for f<f1 we want a derivative behaviour 
f2 = (1 + alpha)*B;  %for f>f2 we want to have the attenuation band
                     %for f1<f<f2 we have the transition band

% frequencies samples of interest
F = Fp/(N/2)/32; 
f = [0:F:f1, f2:F:Fp/2].'; % frequency samples, column vector
d = zeros(length(f),1);   % desired response 
for i=1:length(f)
    if f(i) <= f1 
        d(i) = 2*pi*f(i);    %pass band response
    else
        d(i) = 0;      %stop band response
    end
end

%building matrices
V = zeros(size(f));  %with this we get a 0 coefficient for the center value in the time response   
for n = 1:N/2
    V = [V,2*T*sin(2*pi*f*n*T)];   
end
    
err_lim = [1e-3 1e-3];  %common values used in labs
w = err_lim(1)*(f<=f1) + err_lim(2)*(f>=f2);  %common weighting function

% linear programming solution
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-d;d];
x = linprog(g,A,b); % solve the linear program

% define filter
h0 = [-x(N/2+1:-1:2);x(1:N/2+1)];
t = T*(-N/2:N/2);

% behaviour in frequency obtained by use of function freqz
[H0,ff] = freqz(h0,1,8*(N+1),Fp);
H0 = T*H0; % normalization factor


%% 3 - plot
% show results
figure(1)
subplot(2,1,1)
stem(t,h0); grid; 
title('FIR - derivative - linear programming approach - time domain');
xlabel('time [s]');
ylabel('amplitude []');
subplot(2,1,2)
xlabel('frequency [Hz]');
title('Frequency domain')
yyaxis left
plot(ff,20*log10(abs(H0)));  
ylabel('amplitude [dB]');
yyaxis right
plot(ff,unwrap(angle(H0)));
ylabel('phase [degrees]');
grid; %xlim([0 Fp/2]); ylim([-80, 5]); 
hold on; plot([1,1]*f1,ylim,'r--'); plot([1,1]*f2,ylim,'r--'); hold off; 

% ALTERNATIVE SOLUTION WITH REMEZ using `derivative` parameter
h01 = firpm(N, [0 f1 f2 Fp/2]/(Fp/2), [0 2*pi*f1 0 0] ,'derivative')/T; 
[H01,ff] = freqz(h01,1,8*(N+1),Fp);
H01 = T*H01; % normalization factor 

figure(2)
subplot(2,1,1)
stem(t,h01); 
grid; title('FIR - derivative - remez approach - time domain')
xlabel('time [s]');
ylabel('amplitude []');
subplot(2,1,2)
title('frequency domain')
yyaxis left
plot(ff,20*log10(abs(H01)))
ylabel('amplitude [dB]');
yyaxis right
plot(ff,unwrap(angle(H01)));
ylabel('phase [degrees]');
grid;  
xlabel('frequency [Hz]');
hold on; plot([1,1]*f1,ylim,'r--'); plot([1,1]*f2,ylim,'r--'); hold off;


%% 4 - test the filter

figure(3)
plot( f,d, ff,abs(H0),ff, abs(H01));
title('Approximation of Href(f)');
ylabel('amplitude []');   %since is not given the unity of measurement -> []
xlabel('frequency [Hz]');
grid;  
hold on; plot([1,1]*f1,ylim,'r--'); plot([1,1]*f2,ylim,'r--'); hold off;
legend('Desired', 'Linear programming', 'Remez');


% creating a sine wave       
t = T*(0:N); % i'll show N samples of the signal in the time grid t
Fc = 300;   % freq of the sine in Hz
x_signal = sin(2*pi*Fc*t);  % input sine 

y_signal = T*conv(x_signal,h0);  %filtering with linear prog filter
y_signal = y_signal(1:length(t));

y1_signal = T*conv(x_signal,h01);   %filtering with remez filter
y1_signal = y1_signal(1:length(t));

% plot of input sine, cosine (which is the derivative) and filtered signal
% with linprog and remez filters (in total 4 signals)
figure(4)
yyaxis left
plot(t,x_signal, t, cos(2*pi*Fc*t));  
xlabel('time [s]');
ylabel('amplitude [.]');
yyaxis right
plot(t,y_signal,t,y1_signal)
ylabel('amplitude [.]');
title('Signals versus Time');
legend('sin(2pi300t)','cos(2pi300t)','filtered (linprog)','filtered (remez)');
zoom xon;



%% 5 - ADVANCED -- smallest order of N 

N = 54;   %smallest order of N for my filter 

% frequencies samples of interest
F = Fp/(N/2)/32; % min 32 samples per cosine period
f = [0:F:f1, f2:F:Fp/2].'; % frequency samples, column vector
d = zeros(length(f),1);
for i=1:length(f)
    if f(i) <= f1 
        d(i) = 2*pi*f(i);   %behaviour in the pass band
    else
        d(i) = 0;     %behaviour in the stop band
    end
end

V = zeros(size(f));
for n = 1:N/2
    V = [V,2*T*sin(2*pi*f*n*T)];
end
err_lim = [1e-3 1e-3];
w = err_lim(1)*(f<=f1) + err_lim(2)*(f>=f2);
% linear programming solution
g = [zeros(N/2+1,1);1];
A = [-V, -w; V, -w];
b = [-d;d];
x = linprog(g,A,b); % solve the linear program
h02 = [-x(N/2+1:-1:2);x(1:N/2+1)];  %impulsive response
t = T*(-N/2:N/2);  %time grid
[H02,ff] = freqz(h02,1,8*(N+1),Fp); %freq response through freqz
H02 = T*H02; % normalization factor

% some calculation for visualizing the stopband attenuation
% not mandatory
i_bp = 0;  % last index in the pass band
i_sb = 0;   % first index in the stop band
flag_bp=0;
flag_sb=0;
% extract useful indexes of frequency
for j=1:length(ff)
    if (ff(j)>= f1 && flag_bp == 0)
        i_bp = j;
        flag_bp = 1;        
    end
    if (ff(j)>= f2 && flag_sb == 0)
        i_sb = j;
        flag_sb = 1;
    end
end

%find the maximum (in dB) in the pass-band and in stop-band
[max_bp, id_bp] = max(20*log10(abs(H02(1:i_bp))));
[max_sb, id_sb] = max(20*log10(abs(H02(i_sb:length(ff)))));
sb_attenuation = max_bp - max_sb;  %find the stop-band attenuation

% show results
figure
title('Smallest order of N')
xlabel('frequency [Hz]');
yyaxis left
plot(ff,20*log10(abs(H02))); 
ylabel('amplitude [dB]');
hold on; plot([1,1]*f1,ylim,'r--'); plot([1,1]*f2,ylim,'r--');hold off 
ylim([0 100])
yline(max_bp);
yline(max_sb);
yyaxis right
plot(ff,unwrap(angle(H02))); 
ylabel('phase [degrees]');
grid; 

message = ['With N=', num2str(N), ' the stop-band attenuation is ', num2str(sb_attenuation), ' dB'];
disp(message);


