%% Case study 2
%% Introduction
% * Author:                   Amelia Hines, Dayan Parker, Pavly Youssef
% * Class:                    ESE 351
% * Date:                     Created 4/20/2024, Last Edited 4/30/2024
%% create pulse shapes
Ts=0.1; %define period in miliseconds
bw=1;
dt=0.005; %define sampleing rate
t=-Ts*5:dt:Ts*5; %define time vector from 0 to T
fs=1/dt;
N = 1024; %FFT length
beta1 = 1; %Different Beta Values
beta2 = .5;
beta3 = 0;
name1 = 'Truncated Sinc';
name2 = 'Raised Cosine';

raisCos_1B1 = sinc(t/Ts).*(cos(pi*beta1*t/Ts)./(1-(2*beta1*t/Ts).^2)); %define raised cos
raisCos_1B2 = sinc(t/Ts).*(cos(pi*beta2*t/Ts)./(1-(2*beta2*t/Ts).^2)); %define raised cos
raisCos_1B3 = sinc(t/Ts).*(cos(pi*beta3*t/Ts)./(1-(2*beta3*t/Ts).^2)); %define raised cos
raisCos_1B2(t == Ts/(2*beta2) | t == -Ts/(2*beta2)) = pi/4*sinc(1/(2*beta2)); %Deal with 0
raisCos_1B3(t == Ts/(2*beta3) | t == -Ts/(2*beta3)) = pi/4*sinc(1/(2*beta3)); %Deal with 0
raisCos_2 = raisCos_1B1;
raisCos_3 = raisCos_1B1;
raisCos_1B2(t == Ts/(2*beta1) | t == -Ts/(2*beta1)) = pi/4*sinc(1/(2*beta1)); %Deal with 0
raisCos_1 = raisCos_1B1.*cos(2*pi*10*t); %Define shift
raisCos_2(t == Ts/(2*beta1) | t == -Ts/(2*beta1)) = pi/4*sinc(1/(2*beta1));
raisCos_2 = raisCos_2.*cos(2*pi*20*t);
raisCos_3(t == Ts/(2*beta1) | t == -Ts/(2*beta1)) = pi/4*sinc(1/(2*beta1));
raisCos_3 = raisCos_3.*cos(2*pi*30*t);
sinc_1=sinc(t/Ts).*cos(2*pi*10*t); %define truncated sin function
sinc_2=sinc(t/Ts).*cos(2*pi*20*t); %define shifted truncated sin function
sinc_3=sinc(t/Ts).*cos(2*pi*30*t); %define shifted sinc

figure 
plot(t,sinc_1), title ("Truncated Sinc Shape"),  xlabel('time (seconds)')
grid on

figure 
subplot(3,1,1)
plot(t,raisCos_1B1), title ("Raised Cos Shape \beta = 1"),  xlabel('time (seconds)')
grid on

subplot(3,1,2);
plot(t,raisCos_1B2), title ("Raised Cos Shape \beta = .5"),  xlabel('time (seconds)')
grid on

subplot(3,1,3);
plot(t,raisCos_1B3), title ("Raised Cos Shape \beta = 0"),  xlabel('time (seconds)')
grid on

%Unshifted Raised Cos in Freq
rcos1B1_jw = fftshift(fft(raisCos_1B1,N));
rcos1B2_jw = fftshift(fft(raisCos_1B2,N));
rcos1B3_jw = fftshift(fft(raisCos_1B3,N));
f = 0:(fs/bw)/N:(fs/bw)-(fs/bw)/N;

figure
title('P(\omega): Raised Cosine Unshifted \beta = 1')
subplot(3,1,1), plot(f,abs(rcos1B1_jw)), title('Magnitude of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,2), plot(f,angle(rcos1B1_jw)), title('Phase of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,3), spectrogram(rcos1B1_jw), title('Spectrogram of P(\omega)')

figure
title('P(\omega): Raised Cosine Unshifted \beta = .5')
subplot(3,1,1), plot(f,abs(rcos1B2_jw)), title('Magnitude of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,2), plot(f,angle(rcos1B2_jw)), title('Phase of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,3), spectrogram(rcos1B2_jw), title('Spectrogram of P(\omega)')

figure
title('P(\omega): Raised Cosine Unshifted \beta = 0')
subplot(3,1,1), plot(f,abs(rcos1B3_jw)), title('Magnitude of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,2), plot(f,angle(rcos1B3_jw)), title('Phase of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,3), spectrogram(rcos1B3_jw), title('Spectrogram of P(\omega)')

%10Hz
sinc_jw = fftshift(fft(sinc_1,N)); %take the fft of p(t)
rcos_jw1 =  fftshift(fft(raisCos_1,N));
%plot the frequnecy response of p(t)
figure
title('P(\omega): Truncated Sinc (10Hz)')
subplot(3,1,1), plot(f,abs(sinc_jw)), title('Magnitude of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,2), plot(f,angle(sinc_jw)), title('Phase of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,3), spectrogram(sinc_jw), title('Spectrogram of P(\omega)')

figure
title('P(\omega): Raised Cosine (10Hz)')
subplot(3,1,1), plot(f,abs(rcos_jw1)), title('Magnitude of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,2), plot(f,angle(rcos_jw1)), title('Phase of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,3), spectrogram(rcos_jw1), title('Spectrogram of P(\omega)')

%20 Hz
sinc_jw2 = fftshift(fft(sinc_2,N)); %take the fft of p(t)
rcos_jw2 =  fftshift(fft(raisCos_2,N));

figure
title('P(\omega): Truncated Sinc (20Hz)')
subplot(3,1,1), plot(f,abs(sinc_jw2)), title('Magnitude of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,2), plot(f,angle(sinc_jw2)), title('Phase of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,3), spectrogram(sinc_jw2), title('Spectrogram of P(\omega)')

figure
title('P(\omega): Raised Cosine (20Hz)')
subplot(3,1,1), plot(f,abs(rcos_jw2)), title('Magnitude of P(\omega)')
xlabel('frequnecy (Hz)')

subplot(3,1,2), plot(f,angle(rcos_jw2)), title('Phase of P(\omega)')
xlabel('frequnecy (Hz)')

subplot(3,1,3), spectrogram(rcos_jw2), title('Spectrogram of P(\omega)')

%30 Hz
sinc_jw3 = fftshift(fft(sinc_3,N)); %take the fft of p(t)
rcos_jw3 =  fftshift(fft(raisCos_3,N));

figure
subplot(3,1,1), plot(f,abs(sinc_jw3)), title('Magnitude of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,2), plot(f,angle(sinc_jw3)), title('Phase of P(\omega)')
xlabel('frequency (Hz)')

subplot(3,1,3), spectrogram(sinc_jw3), title('Spectrogram of P(\omega)')
%the band width need to for communications with this shape is pi/25 radians
figure
title('P(\omega): Raised Cosine (30Hz)')
subplot(3,1,1), plot(f,abs(rcos_jw3)), title('Magnitude of P(\omega)')
xlabel('frequnecy (Hz)')

subplot(3,1,2), plot(f,angle(rcos_jw3)), title('Phase of P(\omega)')
xlabel('frequnecy (Hz)')

subplot(3,1,3), spectrogram(rcos_jw3), title('Spectrogram of P(\omega)')
%combine plot
figure
plot(f,abs(sinc_jw)), hold on
plot(f,abs(sinc_jw2)),plot(f,abs(sinc_jw3))
hold off
legend('10Hz', '20Hz', '30Hz')
title('magnitude of P(\omega)')
xlabel('frequency (Hz)')

figure
plot(f,abs(rcos_jw1)), hold on
plot(f,abs(rcos_jw2)),plot(f,abs(rcos_jw3))
hold off
legend('10Hz', '20Hz', '30Hz')
title('magnitude of P(\omega)')
xlabel('frequency (Hz)')
%%
message10 = 'I love ESE 351! :)'; 
message20 = 'I love Jason Trobaugh';
message30 = 'I wish I was a LP Filter';
sigmas = [.1 .5 1];
communication_fullSinc = communication(sinc_1,message10,sinc_2,message20,sinc_3,message30,sigmas,name1);
communication_fullRCos = communication(raisCos_1,message10,raisCos_2,message20,raisCos_3,message30,sigmas,name2);
%%
function message_out = communication(pt,message1,pt2,message2,pt3,message3,sigma,name)
    Ts=0.1; %define period in miliseconds
    dt=0.01; %define sampleing rate
    t=-Ts*5:dt:Ts*5; %define time vector from 0 to T
    fs=1/dt;

    fb=1/(Ts); %create bit rate vector
    upscale=100; %use the length of p(t) to idenitfy the number of zeros 
    % needed to get the full convoltuion y(t)
    

    %convert messages to binary
    x_t10 = str2num(reshape(dec2bin(message1)',1,[])')';
    x_t10(x_t10 == 0) = -1; 

    x_t20 = str2num(reshape(dec2bin(message2)',1,[])')';
    x_t20(x_t20 == 0) = -1;

    x_t30 = str2num(reshape(dec2bin(message3)',1,[])')';
    x_t30(x_t30 == 0) = -1;

    scaled_x10 = zeros(1, numel(x_t10)*upscale);
    length_scaled_10=length(scaled_x10);
    scaled_x20 = zeros(1, numel(x_t20)*upscale);
    length_scaled_20=length(scaled_x20);
    scaled_x30 = zeros(1, numel(x_t30)*upscale);
    length_scaled_30=length(scaled_x30);

    % Insert bit values into the result vector with upsacle zeros in between each value
    scaled_x10(1:100:end) = x_t10;
    scaled_x20(1:100:end) = x_t20;
    scaled_x30(1:100:end) = x_t30;

    %convolute eache scaled signal vecotor with the corrisponding pulse
    %shape
    y_t10=conv(scaled_x10,pt,'same');
    y_t20=conv(scaled_x20,pt2,'same');
    y_t30=conv(scaled_x30,pt3,'same');

    %pad the vectors with zeros to be the same legnth
    scaled_vector=[length(y_t10),length(y_t20),length(y_t30)];
    max_scale=max(scaled_vector);

    y_t10=[y_t10,zeros(1,max_scale-length(y_t10))];
    y_t20=[y_t20,zeros(1,max_scale-length(y_t20))];
    y_t30=[y_t30,zeros(1,max_scale-length(y_t30))];

    y_t=y_t10+y_t20+y_t30; %add the seperate signals together

    noise_time=0:dt:(length(y_t)-1)*dt;%defin a time vector for r(t) and y(t)
    figure, plot(noise_time,y_t), title('orignal data')

    for idx = 1:length(sigma)
        noise=sigma(idx).^2*randn(1,length(y_t)); %create noise vector
        r_t=noise+y_t; %create return signals
        figure, plot(noise_time,r_t), title('noisy data at')
    
        time=0:length(r_t)-1;
    
        %matched filter
        %convolute the noisy data with the opposite of the time vector
        p_opposite_10=pt(end:-1:1);
        z_t10=filter(p_opposite_10,1,r_t.*cos(2*pi*time*10));
        z_t10=z_t10(1,1:length_scaled_10);
    
        p_opposite_20=pt2(end:-1:1);
        z_t20=filter(p_opposite_20,1,r_t.*cos(2*pi*20*time));
        z_t20=z_t20(1,1:length_scaled_20);
        p_opposite_30=pt3(end:-1:1);
        z_t30=filter(p_opposite_30,1,r_t.*cos(2*pi*30*time));
        z_t30=z_t30(1,1:length_scaled_30);
    
        %create a binary output
        output_match10=zeros(1,length(x_t10));
        output_match10(z_t10(upscale:upscale:end)<0)=-1;
        output_match10(z_t10(upscale:upscale:end)>0)=1;
        message_output10=output_match10;
        message_output10(output_match10==-1)=0;
        message_out10 = char(bin2dec(num2str(reshape(message_output10,7,[])')))';
    
        output_match20=zeros(1,length(x_t10)-1);
        output_match20(z_t20(upscale:upscale:end)<0)=-1;
        output_match20(z_t20(upscale:upscale:end)>0)=1;
        message_output20=output_match20;
        message_output20(output_match20==-1)=0;
        message_out20 = char(bin2dec(num2str(reshape(message_output20,7,[])')))';
    
        output_match30=zeros(1,length(x_t30)-1);
        output_match30(z_t30(upscale:upscale:end)<0)=-1;
        output_match30(z_t30(upscale:upscale:end)>0)=1;
        message_output30=output_match30;
        message_output30(output_match30==-1)=0;
        message_out30 = char(bin2dec(num2str(reshape(message_output30,7,[])')))';
        
        message_out=cell(3,1);
        message_out{1}=message_out10;
        message_out{2}=message_out20;
        message_out{3}=message_out30;
    
    
        %10 Hz plot
        diff_idx1=find(output_match10~=x_t10);
        same_idx1=find(output_match10==x_t10);
        error_rate1=length(diff_idx1)/length(x_t10); %define error rate
        

        %20 Hz plot
        diff_idx2=find(output_match20~=x_t20);
        same_idx2=find(output_match20==x_t20);
        error_rate2=length(diff_idx2)/length(x_t20); %define error rate
        
        %30 Hz plot
        diff_idx3=find(output_match30~=x_t30);
        same_idx3=find(output_match30==x_t30);
        error_rate3=length(diff_idx3)/length(x_t30); %define error rate

        figure
        title('Sigma = %.3f',sigma(idx))
        subplot(3,2,1)
        stem(same_idx1,output_match10(same_idx1), 'bo'); 
        hold on; 
        % Plot p_t where p_t is different from x_t
        stem(diff_idx1, output_match10(diff_idx1), 'rx');  % Plot p_t where it is different from x_t with red x marks
        title('Output of the match reciver at 10 Hz')
        xlabel('bits')
        hold off; % Release the hold on the plot
        legend('x_t', 'p_t different from x_t');
        
        subplot(3,2,2), stem(x_t10), title('orignal signal')%plot orignal signal
        xlabel(sprintf('bits, error rate = %0.3f, sigma = %0.3f', error_rate1, sigma(idx)))

        subplot(3,2,3)
        stem(same_idx2,output_match20(same_idx2), 'bo'); 
        hold on; 
        % Plot p_t where p_t is different from x_t
        stem(diff_idx2, output_match20(diff_idx2), 'rx');  % Plot p_t where it is different from x_t with red x marks
        title('Output of the match reciver at 20 Hz')
        xlabel('bits')
        hold off; % Release the hold on the plot
        legend('x_t', 'p_t different from x_t');
        
        subplot(3,2,4), stem(x_t20), title('orignal signal')%plot orignal signal
        xlabel(sprintf('bits, error rate = %0.3f, sigma = %0.3f', error_rate2, sigma(idx)))
        
        subplot(3,2,5)
        stem(same_idx3,output_match30(same_idx3), 'bo'); 
        hold on; 
        % Plot p_t where p_t is different from x_t
        stem(diff_idx3, output_match30(diff_idx3), 'rx');  % Plot p_t where it is different from x_t with red x marks
        title('Output of the match reciver at 30 Hz')
        xlabel('bits')
        hold off; % Release the hold on the plot
        legend('x_t', 'p_t different from x_t');
        
        subplot(3,2,6), stem(x_t30), title('orignal signal')%plot orignal signal
        xlabel(sprintf('bits, error rate = %0.3f, sigma = %0.3f', error_rate3, sigma(idx)))
        
        %calculate the SNR of the system 
    
        signal_power = sum(y_t.^2);
        noise_power = sum(noise.^2);
        SNR = signal_power / noise_power;
        SNR_dB = 10 * log10(SNR); % SNR in decibels

        disp(['Carrier Function = ',name])
        disp(['standard deviation of noise = ',num2str(sigma(idx))])
        disp(['10Hz message: ',message_out{1}])
        disp(['20Hz message: ',message_out{2}])
        disp(['30Hz message: ',message_out{3}])
        disp(['Error Rate 10Hz = ',num2str(error_rate1)])
        disp(['Error Rate 20Hz = ',num2str(error_rate2)])
        disp(['Error Rate 30Hz = ',num2str(error_rate3)])
        disp(['SNR 10Hz = ',num2str(SNR_dB)])
        disp(['bit rate 10Hz = ',num2str(fb)])
    end
end
