
function g=BPF(fs,N,f_L,f_H)  
% clc;clear all;close all;
% N=1280;
% fs=32;
% f_L=0.9; 
% f_H=2;

df=fs/N; 
n=-N+1:N-1;
p=ceil(f_L/df);
q=round((f_H-f_L)/df);

H=[zeros(1,p)  ones(1,q)  zeros(1, N-(2*p+2*q-1))  ones(1,q) zeros(1,p-1)]; %频率向量
h1 = ifft(H,N);        % IDFT
h2 = [h1(2:end)  h1]; % 周期延拓
f=hamming(N);
b=boxcar(N)'; 
win=conv(f,b);
C=b*f;
g=h2.*win/C;
dw=2*pi/N;
w=0:dw/20:2*pi;
Hm=g*exp(-1i*n'*w);
subplot(211); stem(n, g, '.')
subplot(212); plot(w/dw,abs(Hm));
grid on;  hold on; plot(0:N-1, H, 'r.')
