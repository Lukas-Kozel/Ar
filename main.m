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
K = 0.5;
Td = 0.8;
alpha = 13;
Ti = Td/alpha;

W1_inv = K*(Td*s+1)/(Ti*s+1)
S=1/(1+P0_nominal);
% První graf: Frekvenční charakteristiky W1_inv a S_uncertain
figure;
hold on;
bodemag(W1_inv); % Gain * W1_inv - červená
bodemag(S); % S_uncertain - modrá
grid on;
title('Frekvenční charakteristika: W1\_inv a S', 'FontSize', 14);
legend({'W1\_inv', 'S'}, 'FontSize', 10, 'Location', 'Best');
set(gca, 'FontSize', 12);
hold off;

% Aktualizace W1_inv a výpočet W1
W1 = inv(W1_inv)

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

figure;
bode(P0_uncertain)
grid on;
legend("$\mathcal{P}$",Interpreter="latex")
% Plot Bode plots for the nominal and worst-case models
figure;
bode(P0_nominal);
hold on;
bode(P0_worst);
grid on;

% Add legend for clarity
legend('P0', 'P1');
title('Bodeho diagram');

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


% Nyquist Plot
figure;
nyquist(Pmin, 'b'); % Pmin - modrá
hold on;
nyquist(Pmax, 'g'); % Pmax - zelená
nyquist(P0_nominal, 'r--'); % P0_nominal - červená
grid on;

% Přidání legendy
legend({'P_{min}', 'P_{max}', 'P_{0}'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Best');
title('Nyquistův diagram', 'FontSize', 14);
xlabel('Re', 'FontSize', 12);
ylabel('Im', 'FontSize', 12);
axis equal;
hold off;

% Bode Plot
figure;
bode(Pmin, 'b'); % Pmin - modrá
hold on;
bode(Pmax, 'g'); % Pmax - zelená
bode(P0_nominal, 'r--'); % P0_nominal - červená
grid on;

% Přidání legendy
legend({'P_{min}', 'P_{max}', 'P_{0}'}, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'Best');
title('Bodeho diagram', 'FontSize', 14);
hold off;



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


%random nalezený regulátor, který nesplňuje robustnost ve stabilitě:
%ki=1.5
%kp = 1.311


%random nalezený, který má gm skoro inf a Pm má 102, takže ready
%regulátor:
%ki = 0.001087
%kp = 0.01134
%% PID lab regulátor

% % % % % % ki = 0.001087;
% % % % % % kp = 0.01134;
% % % % % % % ki=1.5;
% % % % % % % kp = 1.311;
% % % % % % 
% % % % % % %funkcni nalezen pidlabem 
% % % % % % ki=0.07513
% % % % % % kp=0.09077

%funkcni - nalezen rucne v hinf regionu
%ki = 0.3
%kp = 0.17

%retardovany, ale potrebuju grafy xd:
ki = 0.6329
kp = 0.439

C = kp + ki/s;
L = P0_nominal * C;
L_uzavr = L/(1 + L);

figure;
step(L_uzavr);
title('Jednotkový skok regulovaného systému');
grid on;

% || W2*T0_reg || < 1
T0_reg = L/(1+L);
W2T0_reg_norm = norm(W2 * T0_reg, inf);
fprintf('|| W2*T0_reg ||_inf = %.4f\n\n', W2T0_reg_norm);


% || W2*S0_reg || < 1
S0_reg = 1/(1+L);
W2S0_reg_norm = norm(W2 * S0_reg, inf);
fprintf('|| W2*S0_reg ||_inf = %.4f\n\n', W2S0_reg_norm);


f = logspace(-2,4,10000);
S0_reg_FRS_reg = freqresp(S0_reg, f); 
AS0_reg = (abs(squeeze(S0_reg_FRS_reg))); 
T0_reg_FRS_reg = freqresp(T0_reg, f); 
AT0_reg = (abs(squeeze(T0_reg_FRS_reg)));


W2T0_reg_FRS_reg = freqresp(W2*T0_reg, f); 
AW2T0_reg = (abs(squeeze(W2T0_reg_FRS_reg)));
W1S0_reg_FRS_reg = freqresp(W1*S0_reg, f); 
AW1S0_reg = (abs(squeeze(W1S0_reg_FRS_reg)));


figure
semilogx(f, squeeze(S0_reg_FRS_reg))
hold on
semilogx(f, AW1m1)
title('Frekvenční charakteristika citlivostní funkce S_0 a W_1^{-1}')
legend('S_0','W_1^{-1}')
grid

% bode diagram S0_reg a 1/W1
figure
hold on
bodemag(S0_reg)
bodemag(1/W1)
legend('S_0','W_1^{-1}')
grid

% bode diagram W1*S0_reg
figure
hold on
bodemag(W1*S0_reg)
legend('S_0 W_1^{-1}')
grid



% bode diagram |W_1S_0| + |W_2T_0|
figure;
semilogx(f, AW1S0_reg + AW2T0_reg, 'b', 'LineWidth', 1.5); % Součet |W_1 S_0| + |W_2 T_0|
hold on;
semilogx(f, AW1S0_reg, 'r--', 'LineWidth', 1.5); % |W_1 S_0|
semilogx(f, AW2T0_reg, 'g-.', 'LineWidth', 1.5); % |W_2 T_0|

% Nastavení grafu
title('Frekvenční charakteristika', 'FontSize', 14);
xlabel('Frekvence [rad/s]', 'FontSize', 12);
ylabel('Hodnota', 'FontSize', 12);
legend({'|W_1 S_0| + |W_2 T_0|', '|W_1 S_0|', '|W_2 T_0|'}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;

normAW1S_AW2_reg = norm( AW1S0_reg + AW2T0_reg, inf);
fprintf('|| |W_1*S| + |W_2*T| ||_inf = %.4f\n\n', normAW1S_AW2_reg);

%%
figure;
nyquist(L)

figure
nyquiststability((L.Numerator{1}),(L.Denominator{1}),0.1)


%% regulátor  matlab
nmeas = 1; % Number of measured outputs (y)
ncont = 1; % Number of control inputs (u)
P = augw(P0_nominal,W1,W2,0);
fprintf('Regulátor matlab: \n\n');

[K_hinf, CL, gamma] = hinfsyn(P, nmeas, ncont);

L_hinf = P0_nominal * K_hinf; % Open-loop transfer function
T_hinf = feedback(L_hinf, 1); % Closed-loop transfer function

% Plot Step Response
figure;
step(T_hinf);
title('Jednotkový skok regulovaného systému');
grid on;
% Nyquist Plot for Stability Check
figure;
nyquist(L_hinf);
title('Nyquist Diagram of Open-Loop Transfer Function');
grid on;


[b,a] = ss2tf(L_hinf.A, L_hinf.B, L_hinf.C, L_hinf.D);
L_Hinf = tf(b, a);


% || W2*T0_reg || < 1
T0_reg_matlab = L_Hinf/(1+L_Hinf);
W2T0_reg_norm_matlab = norm(W2 * T0_reg_matlab, inf);
fprintf('|| W2*T0_reg ||_inf = %.4f\n\n', W2T0_reg_norm_matlab);


% || W2*S0_reg || < 1
S0_reg_matlab = 1/(1+L_Hinf);
W2S0_reg_norm_matlab = norm(W2 * S0_reg_matlab, inf);
fprintf('|| W2*S0_reg ||_inf = %.4f\n\n', W2S0_reg_norm_matlab);


f = logspace(-2,4,10000);
S0_reg_FRS_reg_matlab = freqresp(S0_reg_matlab, f); 
AS0_reg_matlab = (abs(squeeze(S0_reg_FRS_reg_matlab))); 
T0_reg_FRS_reg_matlab = freqresp(T0_reg_matlab, f); 
AT0_reg_matlab = (abs(squeeze(T0_reg_FRS_reg_matlab)));


W2T0_reg_FRS_reg_matlab = freqresp(W2*T0_reg_matlab, f); 
AW2T0_reg_matlab = (abs(squeeze(W2T0_reg_FRS_reg_matlab)));
W1S0_reg_FRS_reg_matlab = freqresp(W1*S0_reg_matlab, f); 
AW1S0_reg_matlab = (abs(squeeze(W1S0_reg_FRS_reg_matlab)));


figure
semilogx(f, squeeze(S0_reg_FRS_reg_matlab))
hold on
semilogx(f, AW1m1)
title('Frekvenční charakteristika citlivostní funkce S_0 a W_1^{-1}')
legend('S_0','W_1^{-1}')
grid

% bode diagram S0_reg a 1/W1
figure
hold on
bodemag(S0_reg_matlab)
bodemag(1/W1)
legend('S_0','W_1^{-1}')
grid

% bode diagram W1*S0_reg
figure
hold on
bodemag(W1*S0_reg_matlab)
legend('S_0 W_1^{-1}')
grid



% bode diagram |W_1S_0| + |W_2T_0|
figure
semilogx(f, AW1S0_reg_matlab + AW2T0_reg_matlab)
title('Frekvenční charakteristika |W_1 S_0| + |W_2 T_0|')
legend('|W_1 S_0| + |W_2 T_0|')
grid

normAW1S_AW2_reg = norm( AW1S0_reg_matlab + AW2T0_reg_matlab, inf);
fprintf('|| |W_1*S| + |W_2*T| ||_inf = %.4f\n\n', normAW1S_AW2_reg);

figure
nyquiststability((L_Hinf.Numerator{1}),(L_Hinf.Denominator{1}),0.1)
