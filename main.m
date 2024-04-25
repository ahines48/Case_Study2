Tp = .1;
fb = 1/Tp;
dt = Tp/50;
%t = -Tp:dt:Tp;
t = 0:.01:2;
p = 1 - abs(t)/Tp;
P = fftshift(fft(p)/length(p));
x = sinc(fb*pi*(t-1));
x(t == 1) = 1;
X = fftshift(fft(x) / length(x));

disp(['Bit Rate: ', num2str(fb), ' bps']);

figure;
plot(t, p);
title('Signal p(t)');
xlabel('t');
ylabel('x[n]');

figure;
subplot(2,1,1);
plot(t, abs(P));
title('Magnitude of P(jw)');
xlabel('k');
ylabel('P(jw)');

subplot(2,1,2);
plot(t, angle(P));
title('Phase of P(jw)');
xlabel('k');
ylabel('Phase(a_k) [radians]');

figure;
subplot(2,1,1);
plot(t, abs(X));
title('Magnitude of X(jw)');
xlabel('k');
ylabel('P(jw)');

subplot(2,1,2);
plot(t, angle(X));
title('Phase of X(jw)');
xlabel('k');
ylabel('Phase(a_k) [radians]');

%Range of Sigmas analyzed
sigma_values = 0:0.1:1; 
BER_sign_values = zeros(1, length(sigma_values));
BER_matched_values = zeros(1, length(sigma_values));
SNR_values = zeros(1, length(sigma_values));

for i = 1:length(sigma_values)
    sigma = sigma_values(i);

    % Generate binary message
    b = 2*(rand(1, N) > 0.5) - 1; 
    b_up = zeros(1, N*upsample_factor); 
    b_up(1:upsample_factor:end) = b; 
    y1 = conv(b_up, p);
    y2 = conv(b_up, x);

    % Add noise
    noise = sigma * randn(1, length(y));
    r1 = y1 + noise;
    r2 = y2 + noise;

    % Calculate SNR
    Py1 = sum(y1.^2) / length(y1);
    Py2 = sum(y2.^2) / length(y2);

    Pn = sum(noise.^2) / length(noise);
    SNR1 = 10 * log10(Py1 / Pn);
    SNR2 = 10 * log10(Py2 / Pn);
    SNR1_values(i) = SNR1;
    SNR2_values(i) = SNR2;

    %Sign Based Reciever
    midpoints = 1:upsample_factor:length(r)-1;
    downsampled_r = r(midpoints);
    message = double(downsampled_r > 0);
    message(message == 0) = -1;
    message = message(2:N+1);
    original_message = b;

    decoded_message_sign = message;
    errors_sign = sum(decoded_message_sign ~= original_message);
    BER_sign = errors_sign / N;
    
    %Matched Filter
    p_rev = 1 - abs(-t)/Tp;
    z = conv(r, p_rev);
    indices = 1:upsample_factor:length(z)-1;
    downsampled_z = z(indices);
    msg = double(downsampled_z > 0);
    msg(msg == 0) = -1;
    msg = msg(3:N+2);
    decoded_message_matched = msg;
    errors_matched = sum(decoded_message_matched ~= original_message);
    BER_matched = errors_matched / N;

    % Store BER values
    BER_sign_values(i) = BER_sign; 
    BER_matched_values(i) = BER_matched;
end

% Plot BER vs. SNR
figure;
plot(SNR1_values, BER_sign_values, '-o', SNR_values, BER_matched_values, '-x');
legend('Sign-Based Receiver', 'Matched Filter');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR');
grid on;

figure;
plot(SNR2_values, BER_sign_values, '-o', SNR_values, BER_matched_values, '-x');
legend('Sign-Based Receiver', 'Matched Filter');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR');
grid on;
%% Binary Message Bull
b_up = zeros(1, N*upsample_factor); 
b_up(1:upsample_factor:end) = b; 
y1 = conv(b_up, p);
y2 = conv(b_up, x);
t_plot = 0:dt:(length(y)-1)*dt;

figure;
plot(t_plot, y1, 'LineWidth', 2);
title('Y');
xlabel('t');
ylabel('t');

%% Noise Shit
sigma = 1;
noise = sigma * randn(1, length(y1));
r1 = y1 + noise;
r2 = y2 + noise;

figure;
plot(t_plot, r1, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal r(t)');
grid on; 

figure;
plot(t_plot, r2, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal r(t)');
grid on; 

%% Decode that shit
start_index = round(upsample_factor/2);
%generate samples in middle of upsamplped periods

%Sign Based Reciever
midpoints1 = 1:upsample_factor:length(r1)-1;
midpoints2 = 1:upsample_factor:length(r2)-1;
downsampled_r1 = r1(midpoints1);
downsampled_r2 = r2(midpoints2);
message1 = double(downsampled_r1 > 0);
message2 = double(downsampled_r2 > 0);
message1(message1 == 0) = -1;
message1 = message1(2:N+1);
message2(message2 == 0) = -1;
message2 = message2(2:N+1);
disp(message1);
disp(message2);

%Matched Filter
p_rev = 1 - abs(-t)/Tp;
x_rev = sinc(fb*pi*(t-1));
z1 = conv(r, p_rev);
z2 = conv(r, x_rev);
z1_plot = 0:length(z1)-1;
z2_plot = 0:length(z2)-1;

figure;
plot(z1_plot, z1, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Matched Filter Function');
grid on; 

figure;
plot(z2_plot, z2, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Matched Filter Function');
grid on; 

indices1 = 1:upsample_factor:length(z1)-1;
indices2 = 1:upsample_factor:length(z2)-1;

downsampled_z1 = z1(indices1);
downsampled_z2 = z2(indices2);

msg1 = double(downsampled_z1 > 0);
msg2 = double(downsampled_z2 > 0);
msg1(msg1 == 0) = -1;
msg2(msg2 == 0) = -1;
msg1 = msg1(3:N+2);
msg2 = msg2(3:N+2);
disp(msg1);
disp(msg2);
%%Analysis
Py1 = sum(y1.^2) / length(y1);
Py2 = sum(y2.^2) / length(y2);
Pn = sum(noise.^2) / length(noise);
SNR1 = 10 * log10(Py1 / Pn);
SNR2 = 10 * log10(Py2 / Pn);
disp(['SNR (dB) Triangle: ', num2str(SNR1)]);
disp(['SNR (dB) Sinc: ', num2str(SNR2)]);

original_message = b;
disp(['Original Message: ', num2str(original_message)])
decoded_message1_sign = message1;
decoded_message2_sign = message2;
errors_sign1 = sum(decoded_message1_sign ~= original_message);
errors_sign2 = sum(decoded_message2_sign ~= original_message);
BER_sign1 = errors_sign1 / N;
BER_sign2 = errors_sign2 / N;

disp(['BER (Sign-Based Receiver) Triangle: ', num2str(BER_sign1)]);
disp(['BER (Sign-Based Receiver) Sinc: ', num2str(BER_sign2)]);

decoded_message1_matched = msg1;
decoded_message2_matched = msg2;
errors_matched1 = sum(decoded_message1_matched ~= original_message);
errors_matched2 = sum(decoded_message2_matched ~= original_message);
BER_matched1 = errors_matched1 / N;
BER_matched2 = errors_matched2 / N;

disp(['BER (Matched Filter) Triangle: ', num2str(BER_matched1)]);
disp(['BER (Matched Filter) Sinc: ', num2str(BER_matched2)]);

%% NOTE PLEASE READ
%Displayed values are printed along the code

