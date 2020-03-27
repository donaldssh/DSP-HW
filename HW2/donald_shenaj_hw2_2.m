clearvars
close all
clc

%% 1 - Reading the audio file

% read the corrisponding WAVE file (*.wav)
% y array containing sound samples y(n), n=1:length(y)
% Fp sampling frequency in Hz
[y, Fp] = audioread('sounds/hw1_peter_and_the_wolf_9.wav');
% so that the samplig period is
T = 1/Fp; % [s]

% data
T1 = 75e-6; %[s]
T2 = 318e-6;
T3 = 3180e-6;

%% 2 - IIR transform method filter design

% continous time poles and zero
Z1 = -T/(T2*2);
P1 = -T/(T1*2);
P2 = -T/(T3*2);

% discrete time mapping 
z1 = (1+Z1)/(1-Z1);
z2 = -1;
p1 = (1+P1)/(1-P1);
p2 = (1+P2)/(1-P2);

% constant
const = (T2*T/(2*T1*T3))*(1-Z1)/((1-P1)*(1-P2));

%% 3 - Frequency and time response

% IIR coefficients from poles and zeros
[b,a] = zp2tf([z1; z2], [p2; p1] , const);

%b = [1 -(z1+z2) z1*z2];  %another way to find a and b
%a = [1 -(p1+p2) p1*p2];

% show the freq response
fvtool(b,a)

% show the time response
[h,t] = impz(b,a);   % find the time response
figure
plot(t,abs(h))
grid
title('impulsive response')
title('Transform method filter - time response')
ylabel('amplitude []')
xlabel('time [s]')

%% 4 - Filtering the audio signal

% filter the audio signal
z = filter(b, a, y);

% write the output
z = z/max(abs(z)); % normalize
audiowrite('donald_shenaj_hw2_2_1.wav',z,Fp);

%% 5 - (advanced) DIRECT OPTIMIZATION METHOD
% derived from lab6f_iir_linear_prog.m

% filter order
N = 2; 

% frequencies samples of interest
F = Fp/(N)/64; % min 64 samples per cosine period
f = [0:F:Fp/2].'; % frequency samples, column vector

% build matrices
w = ones(length(f),1); % weight matrix

% D(f) = |H_a(i2pif)|
s = 1i*2*pi*f;
Ha = (1 + s*T2)./((1 + s*T1).*(1 + s*T3));
d = abs(Ha); % ideal filter shape 

% run iterative algorithm
[zk,pk] = zero_pole_identification(d,f,w,N,T);

% extracts poles/zeros of interest (may fail in the general case)
pk = pk(abs(pk)<1);
zk = zk([1:4:end,4:4:end]);
%zk = zk(abs(zk)<1);

% builds and shows the frequency response
ff = 0:Fp/10000:Fp;  % frequency samples
H = ones(size(ff)); % extracts frequency response (normalized at f=0)
for k = 1:N
    num = (1-zk(k)*exp(-2i*pi*ff*T)) / (1-zk(k)); % normalized at f=0
    den = (1-pk(k)*exp(-2i*pi*ff*T)) / (1-pk(k)); % normalized at f=0
    H = H.*num./den;
end

% plots the frequency response
figure 
plot(ff,20*log10(abs(H)))
grid
ylabel('amplitude [dB]')
xlabel('frequency [Hz]')
title('frequency response - amplitude')
hold on
plot(f,20*log10(d))
hold off
xlim([0 Fp/2])
legend('direct optimization filter', 'desired amplitude');


% plots the phase frequency response
figure 
plot(ff,angle(H));
grid
ylabel('phase [degrees]')
xlabel('frequency [Hz]')
xlim([0 Fp/2])
title('frequency response - phase')


[b,a] = zp2tf(zk,pk,1);  % find the iir coefficients

% filter the audio and write the file
z = filter(b,a,y);
z = z/max(z);  % normalize the amplitude
audiowrite('donald_shenaj_hw2_2_2.wav',z,Fp);


% plot the time response
[h,t] = impz(b,a);   % find the time response
figure
plot(t,abs(h))
grid
title('impulsive response')
title('Direct optimization filter - time response')
ylabel('amplitude []')
xlabel('time [s]')

%% find zero and pole - linear programming function  

% from lab6f_iir_linear_prog.m

function [zk,pk] = zero_pole_identification(d,f,w,N,T)
%
% inputs:
% f  : column vector containing the selected frequencies
% d  : column vector containing the target function (squared absolute 
%      value) at frequencies f
% w  : weighting vector (same dimension as f and d)
% N  : ARMA order
% T  : sampling period in the time domain
%
% outputs:
% zk : zeros of the (self-reciprocal) transfer function
% pk : poles of the (self-reciprocal) transfer function
%

V = ones(size(f)); % builds cosines matrix
for n = 1:N
    V = [V,2*cos(2*pi*f*n*T)];
end
wdV = diag(d.*w)*V; % builds product matrices
wV = diag(w)*V;

% starting solution
p_i = poly(0.5*ones(N,1))'; % poles
q_i = [1;zeros(N,1)]; % zeros
Pi = V*p_i; % corresponding numerator
Qi = V*q_i; %           and denominator
delta = max(abs(d-Pi./Qi)); % error

% cycle
for it = 1:100 % set maximum iterations at 100
    
    % vectors for the target function
    g = [zeros(1,2*N+2), 1]; 
    % matrix and vector for inequalities
    A = [ -wV,  wdV-delta*V,  -Qi; ... % absolute bound (upper part)
           wV, -wdV-delta*V,  -Qi; ... % absolute bound (lower part)
           -V,          0*V, 0*Qi; ... % Pi >= 0
          0*V,            V, 0*Qi; ... % Qi <= 1
        ];
    b = [zeros(size(A,1)-size(V,1),1); ones(size(V,1),1)];
    
    % solve the linear program
    options = optimoptions('linprog','Algorithm','interior-point',...
        'Display','off');
    [x,~,exitflag] = linprog(g,A,b,[],[],[],[],[],options); 
    
    % exit criterion
    if (exitflag~=1) % exit if the solver fails
        break
    else
        tmp = max(abs(d-(V*x(1:N+1))./(V*x(N+2:2*N+2)))); % error
        if (tmp/delta>0.999999) % exit also if the improvement is too small
            break
        else % update your solution
            delta = tmp;
            disp(['it#' num2str(it) ... % you may want to comment this !!!
                  ': error = ' num2str(10*log10(delta)) ' dB'])
            p_i = x(1:N+1);
            q_i = x(N+2:2*N+2);    
            Qi = V*q_i;
        end
    end
end

% extract poles and zeros from the final solution
pk = roots([q_i(end:-1:2); q_i]);
zk = roots([p_i(end:-1:2); p_i]);

end
