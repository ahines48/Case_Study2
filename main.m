Tp = .1;
fb = 1/Tp;
dt = Tp/50;
t = -Tp:dt:Tp;
p = 1 - abs(t)/Tp;
P = fftshift(fft(p)/length(p));
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

%Range of Sigmas analyzed
sigma_values = 0:0.1:1; 
BER_sign_values = zeros(1, length(sigma_values));
BER_matched_values = zeros(1, length(sigma_values));
SNR_values = zeros(1, length(sigma_values));

for i = 1:length(sigma_values)
    sigma = sigma_values(i);

    % Generate binary message
    x = 2*(rand(1, N) > 0.5) - 1; % Generates N bits of ±1
    x_up = zeros(1, N*upsample_factor); 
    x_up(1:upsample_factor:end) = x; 
    y = conv(x_up, p);

    % Add noise
    noise = sigma * randn(1, length(y));
    r = y + noise;

    % Calculate SNR
    Py = sum(y.^2) / length(y);
    Pn = sum(noise.^2) / length(noise);
    SNR = 10 * log10(Py / Pn);
    SNR_values(i) = SNR;

    %Sign Based Reciever
    midpoints = 1:upsample_factor:length(r)-1;
    downsampled_r = r(midpoints);
    message = double(downsampled_r > 0);
    message(message == 0) = -1;
    message = message(2:N+1);
    original_message = x;

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
plot(SNR_values, BER_sign_values, '-o', SNR_values, BER_matched_values, '-x');
legend('Sign-Based Receiver', 'Matched Filter');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR');
grid on;

%% Binary Message Bull
N = 5; 
x = 2*(rand(1, N) > 0.5) - 1; % Generates N bits of ±1
disp(x);
upsample_factor = 50; 
x_up = zeros(1, N*upsample_factor); 
x_up(1:upsample_factor:end) = x; 
y = conv(x_up, p);
t_plot = 0:dt:(length(y)-1)*dt;

figure;
plot(t_plot, y, 'LineWidth', 2);
title('Y');
xlabel('t');
ylabel('t');

%% Noise Shit
sigma = 1;
noise = sigma * randn(1, length(y));
r = y + noise;

figure;
plot(t_plot, r, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal r(t)');
grid on; 


%% Decode that shit
start_index = round(upsample_factor/2);
%generate samples in middle of upsamplped periods

%Sign Based Reciever
midpoints = 1:upsample_factor:length(r)-1;
downsampled_r = r(midpoints);
message = double(downsampled_r > 0);
message(message == 0) = -1;
message = message(2:N+1);
disp(message);

%Matched Filter
p_rev = 1 - abs(-t)/Tp;
z = conv(r, p_rev);
z_plot = 0:length(z)-1;
figure;
plot(z_plot, z, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Matched Filter Function');
grid on; 

indices = 1:upsample_factor:length(z)-1;
downsampled_z = z(indices);
msg = double(downsampled_z > 0);
msg(msg == 0) = -1;
msg = msg(3:N+2);
disp(msg);

%%Analysis
Py = sum(y.^2) / length(y);
Pn = sum(noise.^2) / length(noise);
SNR = 10 * log10(Py / Pn);
disp(['SNR (dB): ', num2str(SNR)]);

original_message = x;
disp(['Original Message: ', num2str(original_message)])
decoded_message_sign = message;
errors_sign = sum(decoded_message_sign ~= original_message);
BER_sign = errors_sign / N;
disp(['BER (Sign-Based Receiver): ', num2str(BER_sign)]);

decoded_message_matched = msg;
errors_matched = sum(decoded_message_matched ~= original_message);
BER_matched = errors_matched / N;
disp(['BER (Matched Filter): ', num2str(BER_matched)]);

%% Truncated Sinc Messenger
x = sin(20*pi*(t-1))./(20*pi*(t-1));
x(t == 1) = 1;
X = fft(x) / length(x);

%% NOTE PLEASE READ
%Displayed values are printed along the code

