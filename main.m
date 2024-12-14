close all; clear all; clc;
K_uncertain = ureal('K', 20, 'Percentage', 15);  
T1_uncertain = ureal('T1', 0.2, 'Percentage', 15); 
T2_uncertain = ureal('T2', 0.4, 'Percentage', 15);

s = tf('s');


P0_nominal = 20/((0.2*s+1)*(0.4*s+1))
P0_uncertain = K_uncertain / ((T1_uncertain * s + 1) * (T2_uncertain * s + 1));

S_uncertain = 1 / (1 + P0_uncertain);
figure;
bode(S_uncertain);
grid on;

%% nefunguje
% sirka_pasma = 10;
% amplitude = 6;
% gain = 10^(amplitude/20)/1.5;
% %W1 = (s+0.1)/(s+15) %alright i guess
% W1_inv= (s+2)/(s+sirka_pasma);

%% funguje
% derivacni clanek
K = 0.3;
Td = 0.9;
alpha = 10;
Ti = Td/alpha;

W1_inv = K*(Td*s+1)/(Ti*s+1)

% První graf: Frekvenční charakteristiky W1_inv a S_uncertain
figure;
hold on;
bodemag(W1_inv); % Gain * W1_inv - červená
bodemag(S_uncertain); % S_uncertain - modrá
grid on;
xlabel('Frekvence [rad/s]', 'FontSize', 12);
ylabel('Zisk [dB]', 'FontSize', 12);
title('Frekvenční charakteristika: W1\_inv a S\_uncertain', 'FontSize', 14);
legend({'W1\_inv', 'S\_uncertain'}, 'FontSize', 10, 'Location', 'Best');
set(gca, 'FontSize', 12);
hold off;

% Aktualizace W1_inv a výpočet W1
W1 = inv(W1_inv);

% Výpočet normy H∞
disp('H∞ norma:');
H_inf_norm = norm(W1 * S_uncertain, "inf");
disp(H_inf_norm);


% Čtvrtý graf: Frekvenční charakteristika W1*S_uncertain
figure;
bodemag(W1 * S_uncertain); % W1 * S_uncertain - zelená
grid on;
xlabel('Frekvence [rad/s]', 'FontSize', 12);
ylabel('Zisk [dB]', 'FontSize', 12);
title('Frekvenční charakteristika: W1 * S\_uncertain', 'FontSize', 14);
legend({'W1 * S\_uncertain'}, 'FontSize', 10, 'Location', 'Best');
set(gca, 'FontSize', 12);

%%
% Define uncertain parameters
K_uncertain = ureal('K', 20, 'Percentage', 15);  
T1_uncertain = ureal('T1', 0.2, 'Percentage', 15); 
T2_uncertain = ureal('T2', 0.4, 'Percentage', 15);

% Define Laplace variable
s = tf('s');

% Nominal and uncertain transfer functions
P0_nominal = 20 / ((0.2 * s + 1) * (0.4 * s + 1));
P0_uncertain = K_uncertain / ((T1_uncertain * s + 1) * (T2_uncertain * s + 1));

% Compute the worst-case response using worst-case gain
[~, worstCaseModel] = wcgain(P0_uncertain);
P0_worst = worstCaseModel.K / ((worstCaseModel.T1 * s + 1) * (worstCaseModel.T2* s + 1));

% Plot Bode plots for the nominal and worst-case models
figure;
bode(P0_nominal);
hold on;
bode(P0_worst);
grid on;

% Add legend for clarity
legend('Nominal Model', 'Worst-case Model');
title('Bode Plot: Nominal vs Worst-case Transfer Function');

W2 = P0_worst/P0_nominal -1

Pmin = P0_nominal*(1-W2)
Pmax = P0_nominal*(1+W2)

figure;
bode(P0_uncertain/(1+P0_uncertain))
hold on;
bode(W2*P0_nominal/(1+P0_nominal))
S=1/(1+P0_nominal)
T=P0_nominal/(1+P0_nominal)
figure;
bode(W1*S)
hold on;
bode(W2*T)
bode(W1*S+W2*T)
grid on;
figure;
nyquist(Pmin)
hold on;
nyquist(Pmax)
nyquist(P0_nominal)

figure;
bode(Pmin)
hold on;
bode(Pmax)
bode(P0_nominal)


figure;
bode(W2)
hold on;
bode(W1)

%% Výpočet všeho potřebného pro W1 a W2
f = logspace(-2,4,10000);
FRS = freqresp(S, f); % frekvencni odezva citlivostni funkce
FRW1 = freqresp(W1, f); % frekvencni odezva vahove funkce
FRW1m1 = freqresp(1/W1, f); % frekvencni odezva funkce 1/W1
FRW1S = freqresp(W1*S, f); % frekvencni odezva funkce S*W1
AS = (abs(squeeze(FRS))); % AFCH citlivostni funkce
AW1 = (abs(squeeze(FRW1))); % AFCH vahove funkce
AW1m1 = (abs(squeeze(FRW1m1))); % AFCH funkce 1/W1
AW1S = (abs(squeeze(FRW1S))); % AFCH funkce S*W1

FRT = freqresp(T, f); % frekvencni odezva komplementarni citlivosti
FRW2 = freqresp(W2, f); % frekvencni odezva vahove funkce
FRW2m1 = freqresp(1/W2, f); % frekvencni odezva funkce 1/W2
FRW2T = freqresp(W2*T, f); % frekvencni odezva funkce T*W2
AT = (abs(squeeze(FRT))); % AFCH komplementarni citlivosti
AW2 = (abs(squeeze(FRW2))); % AFCH vahove funkce
AW2m1 = (abs(squeeze(FRW2m1))); % AFCH funkce 1/W2
AW2T = (abs(squeeze(FRW2T))); % AFCH funkce T*W2


%% Testy pro W1 a W2

W1_W2_test(W1, true, W2, true, S, T) % W1, grafy W1, W2, grafy W2, S, T

%%
% Define uncertain parameters
% Define uncertain parameters
% K_uncertain = ureal('K', 20, 'Percentage', 15);  
% T1_uncertain = ureal('T1', 0.2, 'Percentage', 15); 
% T2_uncertain = ureal('T2', 0.4, 'Percentage', 15);
% 
% s = tf('s');
% 
% % Define nominal and uncertain transfer functions
% P0_nominal = 20 / ((0.2 * s + 1) * (0.4 * s + 1));
% P0_uncertain = K_uncertain / ((T1_uncertain * s + 1) * (T2_uncertain * s + 1));
% 
% % Define W2 as given
% W2 = (P0_uncertain - P0_nominal) / P0_nominal;

% Define 5 logarithmically spaced frequencies
frequencies = [0 0.57 1.2 2.5 4 7.58]; % Frequencies: 0.1 to 100 rad/s


% Compute |W2(iω)| for each frequency
W2_radii = zeros(1, length(frequencies)); % Preallocate for radii
nyquist_points = zeros(1, length(frequencies)); % Store nominal points for markers
for f = 1:length(frequencies)
    freq = frequencies(f);
    W2_value = evalfr(W2*P0_nominal, 1j * freq); % Evaluate W2 at s = iω
    W2_radii(f) = abs(W2_value); % Take the absolute value
    
    nominal_response = evalfr(P0_nominal, 1j * freq); % Nominal Nyquist point
    nyquist_points(f) = nominal_response; % Store for plotting markers
end

% Nyquist diagram of the nominal model
figure;
nyquist(P0_nominal, 'b--'); % Nominal model in dashed blue
hold on;

% Adjust the Nyquist plot limits dynamically
% x_min = min(real(nyquist_points)) - max(W2_radii);
% x_max = max(real(nyquist_points)) + max(W2_radii);
% y_min = min(imag(nyquist_points)) - max(W2_radii);
% y_max = max(imag(nyquist_points)) + max(W2_radii);
% xlim([x_min, x_max]);
% ylim([y_min, y_max]);

% Plot circles for uncertainty at each frequency
for f = 1:length(frequencies)
    freq = frequencies(f);
    nominal_response = evalfr(P0_nominal, 1j * freq); % Evaluate P0_nominal at s = iω
    radius = W2_radii(f); % Radius for this frequency

    % Create a circle around the nominal point
    theta = linspace(0, 2 * pi, 100); % Circle perimeter
    x = real(nominal_response) + radius * cos(theta); % X-coordinates
    y = imag(nominal_response) + radius * sin(theta); % Y-coordinates
    plot(x, y, 'r--', 'LineWidth', 1.2); % Plot the circle in red dashed lines

    % Add a marker and frequency annotation
    plot(real(nominal_response), imag(nominal_response), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k'); % Mark nominal point
    text(real(nominal_response) + radius, imag(nominal_response), sprintf('%.2f rad/s', freq), 'FontSize', 8); % Annotate frequency
end

% Finalize the plot
hold on;
nyquist(Pmin)
nyquist(Pmax)
grid on;
legend('P0 (Nominal)', 'Uncertainty Circles');
title('Nyquist Diagram with Uncertainty Circles at Specific Frequencies');
xlabel('Real Axis');
ylabel('Imaginary Axis');
axis equal;


