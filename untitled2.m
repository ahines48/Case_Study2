Ts = .1; % symbol period (rate 1/Ts)
dt = .01; % sample period
t = -5*Ts:dt:5*Ts; % time vector
x = sinc(t/Ts); % define sinc, note Matlab convention sinc(x) =
sin(pi*x)/(pi*x)
figure
subplot(2,1,1), plot(t,x)
xlabel('time (s)'), ylabel('x(t)'), title('Truncated sinc')
fs = 1/dt; % sample frequency
Nfft = 1024; % length of fft
f = [0:fs/Nfft:fs-fs/Nfft];
subplot(2,1,2), plot(f,abs(fftshift((fft(x,Nfft)))))
xlabel('frequency (Hz)'), ylabel('|X(j\omega)|')

figure, plot(f,angle(fftshift((fft(x,Nfft)))))