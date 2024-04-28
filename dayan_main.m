%% Case study 2
%% Introduction
% * Author:                   Amelia Hines and Dayan Parker
% * Class:                    ESE 351
% * Date:                     Created 4/5/2024, Last Edited 4/5/2024
%% create pulse shapes
Ts=0.1; %define period in miliseconds
bw=5;
dt=Ts/(10*bw); %define sampleing rate
t=-Ts:dt:Ts; %define time vector from 0 to T
fs=1/dt;

triangle=1-abs(t)/Ts;%define the shape of the output vector; %define a triangular pulse shape
triangle_fft = fftshift(fft(triangle)/length(triangle)); %take fft of triangle
sinc_1=sinc(t/Ts); %define truncated sin function
sinc_2=sinc((t-0.1)/Ts); %define shifted truncated sin function
sinc_3=sinc((t+0.1)/Ts); %define shifted sinc

figure 
plot(t,sinc_1), title ("truncated sin shape"),  xlabel('time (seconds)')
hold on
plot(t,sinc_2), plot(t,abs(sinc_1)), 
legend('non-shifted sinc','left shifted sinc', 'right shifted sinc');
hold off
grid on

fftL=1024;
sinc_jw = fftshift(fft(sinc_1,fftL)); %take the fft of p(t)
f = 0:(fs/bw)/fftL:(fs/bw)-(fs/bw)/fftL;
%plot the frequnecy response of p(t)
figure
subplot(3,1,1), plot(f,abs(sinc_jw)), title('magnitude of P(\omega)')
xlabel('frequnecty (rads)')

subplot(3,1,2), plot(f,angle(sinc_jw)), title('phase of P(\omega)')
xlabel('frequnecty (rads)')

subplot(3,1,3), spectrogram(sinc_jw), title('spectrogram of P(\omega)')
%the band width need to for communications with this shape is pi/25 radians

%% create message
p_t=sinc_1;
fb=1/(Ts); %create bit rate vector
upscale=100; %use the length of p(t) to idenitfy the number of zeros 
% needed to get the full convoltuion y(t)
message = 'I love ESE 351! :)';
x_t = str2num(reshape(dec2bin(message)',1,[])')';
x_t(x_t == 0) = -1;

scaled_x = zeros(1, numel(x_t)*upscale);
% Insert bit values into the result vector with upsacle zeros in between each value
scaled_x(1:100:end) = x_t;
y_t=conv(scaled_x,p_t,'same'); %create y[n] vector
noise_time=0:dt:(length(y_t)-1)*dt;%defin a time vector for r(t) and y(t)
figure, plot(noise_time,y_t), title('orignal data')

sigma=0.1; %define noise parameter 
noise=sigma.^2*randn(1,length(y_t)); %create noise vector
r_t=noise+y_t; %create return signals
figure, plot(noise_time,r_t), title('noisy data')

%% decode message
%sign reciver
output=zeros(1,length(x_t)-1); %initailize an output vector
%write a 1 or -1 depending on every 100th entry in r_t
output(r_t(upscale:upscale:end)>0)=1;
output(r_t(upscale:upscale:end)<0)=-1;
message_output=output;
message_output(output==-1)=0;
message_out_sign = char(bin2dec(num2str(reshape(message_output,7,[])')))';
disp(message_out_sign)

%plot the output compared to the input
% Plot x_t
diff_idx=find(output~=x_t); %find indexes were output is different from x_t
same_idx=find(output==x_t); %find indexes were output is the same from x_t
error_rate=length(diff_idx)/length(x_t); %define error rate

figure
subplot(2,1,1)
stem(same_idx,output(same_idx), 'bo'); 
hold on; 
% Plot p_t where p_t is different from x_t
stem(diff_idx, output(diff_idx), 'rx');  % Plot p_t where it is different from x_t with red x marks
title('Output of the sign reciver')
xlabel('bits')
hold off; % Release the hold on the plot
legend('x_t', 'p_t different from x_t');

subplot(2,1,2), stem(x_t), title('orignal signal') %plot original singal
xlabel(sprintf('bits, error rate = %0.3f',error_rate))

%matched filter
p_opposite=p_t(end:-1:1);
z_t = filter(p_opposite,1,r_t);

output_match=zeros(1,length(x_t)-1);
output_match(z_t(upscale:upscale:end)<0)=-1;
output_match(z_t(upscale:upscale:end)>0)=1;
message_output=output_match;
message_output(output_match==-1)=0;
message = char(bin2dec(num2str(reshape(message_output,7,[])')))';
disp(message)

% Plot x_t
diff_idx=find(output_match~=x_t);
same_idx=find(output_match==x_t);
error_rate=length(diff_idx)/length(x_t); %define error rate

figure
subplot(2,1,1)
stem(same_idx,output_match(same_idx), 'bo'); 
hold on; 
% Plot p_t where p_t is different from x_t
stem(diff_idx, output_match(diff_idx), 'rx');  % Plot p_t where it is different from x_t with red x marks
title('Output of the match reciver')
xlabel('bits')
hold off; % Release the hold on the plot
legend('x_t', 'p_t different from x_t');

subplot(2,1,2), stem(x_t), title('orignal signal')%plot orignal signal
xlabel(sprintf('bits, error rate = %0.3f',error_rate))


%calculate the SNR of the system 
signal_power = sum(y_t.^2);
noise_power = sum(noise.^2);
SNR = signal_power / noise_power;

disp(['SNR =',num2str(SNR)])
disp(['bit rate =',num2str(fb)])
disp(['standard deviation of noise =',num2str(sigma)])


