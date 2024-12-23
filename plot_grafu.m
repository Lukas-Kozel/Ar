% Define the transfer function
s = tf('s');
P0_nominal = 20 / ((0.2 * s + 1) * (0.4 * s + 1));

% Generate samples for uncertain system
num_samples = 6; % Number of samples
P0_samples = usample(P0_uncertain, num_samples); % Sampled systems

% Frequency range for the plot
f = logspace(-2, 2, 1000); % Frequency from 10^-2 to 10^2 rad/s

% Plot LAFCH (magnitude only) for P0_nominal and samples
figure;
hold on;

% Plot P0_nominal
bodemag(P0_nominal, f, 'r--'); % Nominal system in blue dashed line

% Plot sampled uncertain systems
bodemag(P0_samples(:,:,1), f, 'b-'); % Sample 1 (Red solid line)
bodemag(P0_samples(:,:,2), f, 'g-'); % Sample 2 (Green solid line)
bodemag(P0_samples(:,:,3), f, 'm-'); % Sample 3 (Magenta solid line)
bodemag(P0_samples(:,:,4), f, 'c-'); % Sample 4 (Cyan solid line)
bodemag(P0_samples(:,:,5), f, 'k-'); % Sample 5 (Black solid line)
bodemag(P0_samples(:,:,6), f, 'y-');
% Add grid and title
grid on;

% Add legend for nominal and sampled systems
legend('P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'FontSize', 10, 'Location', 'Best');

hold off;

% Plot LAFCH (magnitude only) for P0_nominal and samples
figure;
hold on;

% Plot P0_nominal
nyquist(P0_nominal, f, 'r--'); % Nominal system in blue dashed line

% Plot sampled uncertain systems
nyquist(P0_samples(:,:,1), f, 'b-'); % Sample 1 (Red solid line)
nyquist(P0_samples(:,:,2), f, 'g-'); % Sample 2 (Green solid line)
nyquist(P0_samples(:,:,3), f, 'm-'); % Sample 3 (Magenta solid line)
nyquist(P0_samples(:,:,4), f, 'c-'); % Sample 4 (Cyan solid line)
nyquist(P0_samples(:,:,5), f, 'k-'); % Sample 5 (Black solid line)
nyquist(P0_samples(:,:,6), f, 'y-'); % Sample 5 (Black solid line)

% Add grid and title
title('Nyquist', 'FontSize', 14);
xlabel('Re', 'FontSize', 12);
ylabel('Im', 'FontSize', 12);

% Add legend for nominal and sampled systems
legend('P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'FontSize', 10, 'Location', 'Best');

hold off;


%% 
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
% bode diagram |W_1 S_0| + |W_2 T_0| s vertikální osou
figure;
semilogx(f, AW1S0_reg + AW2T0_reg, 'b', 'LineWidth', 1.5); % Součet |W_1 S_0| + |W_2 T_0|
hold on;
semilogx(f, AW1S0_reg, 'r--', 'LineWidth', 1.5); % |W_1 S_0|
semilogx(f, AW2T0_reg, 'g-.', 'LineWidth', 1.5); % |W_2 T_0|

% Najít frekvenci w, kde |W_1 S_0| + |W_2 T_0| dosahuje nejvyšší hodnoty
[max_value, max_index] = max(AW1S0_reg + AW2T0_reg); % Nejvyšší hodnota a její index
max_frequency = f(max_index); % Odpovídající frekvence

% Přidání vertikální osy
line([max_frequency max_frequency], [0 max_value], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5); % Vertikální čára
text(max_frequency, max_value, ...
    sprintf('\\omega = %.2f, max = %.4f', max_frequency, max_value), ...
    'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Nastavení grafu
title('Frekvenční charakteristika', 'FontSize', 14);
xlabel('Frekvence [rad/s]', 'FontSize', 12);
ylabel('Hodnota', 'FontSize', 12);
legend({'|W_1 S_0| + |W_2 T_0|', '|W_1 S_0|', '|W_2 T_0|', '\omega_{max}'}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;
hold off;

% Výpočet normy
normAW1S_AW2_reg = norm(AW1S0_reg + AW2T0_reg, inf);
fprintf('|| |W_1*S| + |W_2*T| ||_inf = %.4f\n\n', normAW1S_AW2_reg);


%%

% Definice T0 a S0
L = C * P0_nominal; % Otevřená smyčka
T0 = L / (1 + L);   % Funkce T0
S0 = 1 / (1 + L);   % Funkce S0

% Frekvenční rozsah
f = logspace(-2, 4, 10000); % Frekvence od 10^-2 do 10^4 rad/s

% Frekvenční odezva T0 a S0
[mag_T0, ~, omega] = bode(T0, f); % Magnituda T0 (linear), fáze a frekvence
[mag_S0, ~, ~] = bode(S0, f);     % Magnituda S0 (linear)
mag_T0 = squeeze(mag_T0);         % Převod na vektor
mag_S0 = squeeze(mag_S0);         % Převod na vektor

% Převod na decibely
mag_T0_dB = 20 * log10(mag_T0);
mag_S0_dB = 20 * log10(mag_S0);

% Najít frekvenci, kde |T0| = -3 dB
target_dB = -3; % Hodnota v dB
[~, idx] = min(abs(mag_T0_dB - target_dB)); % Index nejbližší hodnoty k -3 dB
omega_c = omega(idx); % Frekvence

% Vykreslení grafu
figure;

% T0 modře
semilogx(omega, mag_T0_dB, 'b-', 'LineWidth', 1.5); % T0 funkce
hold on;

% S0 zeleně
semilogx(omega, mag_S0_dB, 'g-', 'LineWidth', 1.5); % S0 funkce

% Horizontální červená přerušovaná čára na -3 dB
line([min(omega) max(omega)], [target_dB target_dB], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);

% Černý bod na omega_c
plot(omega_c, mag_T0_dB(idx), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(omega_c, mag_T0_dB(idx) - 5, sprintf('\\omega_c = %.4f rad/s', omega_c), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', 'k');

% Nastavení grafu
title('Bodeho diagram', 'FontSize', 14);
legend({'T_0', 'S_0', '-3 dB', '\omega_c '}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;
hold off;




%%

figure;
% otevrena smycka
L = minreal(P0_nominal*C);
FRL = squeeze(freqresp(L, f)); % frekvencni odezva otevrene smycky
iL = real(FRL);
rL = imag(FRL);
plot(rL, iL)
hold on;

[A,i] = max(AW1S + AW2T);
%i = 5500;
kruh1 = -1+AW1(i)*exp(1j*(-2*pi:0.01:2*pi));
plot(kruh1)


W2L = W2*L;
FRW2L = freqresp(W2L, f); % frekvencni odezva funkce W2*L
AW2L = (abs(squeeze(FRW2L))); % AFCH funkce L*W2
kruh2 = FRL(i) + AW2L(i)*exp(1j*(-2*pi:0.01:2*pi));
plot(kruh2)



%% Matlab regulator plot

f = logspace(-2, 4, 10000); % Frekvence od 10^-2 do 10^4 rad/s

% bode diagram |W_1S_0| + |W_2T_0|
figure;
semilogx(f, AW1S0_reg_matlab + AW2T0_reg_matlab, 'b', 'LineWidth', 1.5); % Součet |W_1 S_0| + |W_2 T_0|
hold on;
semilogx(f, AW1S0_reg_matlab, 'r--', 'LineWidth', 1.5); % |W_1 S_0|
semilogx(f, AW2T0_reg_matlab, 'g-.', 'LineWidth', 1.5); % |W_2 T_0|

% Nastavení grafu
title('Frekvenční charakteristika', 'FontSize', 14);
xlabel('Frekvence [rad/s]', 'FontSize', 12);
ylabel('Hodnota', 'FontSize', 12);
legend({'|W_1 S_0| + |W_2 T_0|', '|W_1 S_0|', '|W_2 T_0|'}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;

normAW1S_AW2_reg = norm( AW1S0_reg_matlab + AW2T0_reg_matlab, inf);
fprintf('|| |W_1*S| + |W_2*T| ||_inf = %.4f\n\n', normAW1S_AW2_reg);




L_hinf = P0_nominal * K_hinf; % Open-loop transfer function
T_hinf = feedback(L_hinf, 1); % Closed-loop transfer function
[b,a] = ss2tf(L_hinf.A, L_hinf.B, L_hinf.C, L_hinf.D);
L_Hinf = tf(b, a);



T0_reg_matlab = L_Hinf/(1+L_Hinf);
S0_reg_matlab = 1/(1+L_Hinf);

% Frekvenční odezva T0 a S0
[mag_T0, ~, omega] = bode(T0_reg_matlab, f); % Magnituda T0 (linear), fáze a frekvence
[mag_S0, ~, ~] = bode(S0_reg_matlab, f);     % Magnituda S0 (linear)
mag_T0 = squeeze(mag_T0);         % Převod na vektor
mag_S0 = squeeze(mag_S0);         % Převod na vektor

% Převod na decibely
mag_T0_dB = 20 * log10(mag_T0);
mag_S0_dB = 20 * log10(mag_S0);

% Najít frekvenci, kde |T0| = -3 dB
target_dB = -3; % Hodnota v dB
[~, idx] = min(abs(mag_T0_dB - target_dB)); % Index nejbližší hodnoty k -3 dB
omega_c = omega(idx); % Frekvence

% Vykreslení grafu
figure;

% T0 modře
semilogx(omega, mag_T0_dB, 'b-', 'LineWidth', 1.5); % T0 funkce
hold on;

% S0 zeleně
semilogx(omega, mag_S0_dB, 'g-', 'LineWidth', 1.5); % S0 funkce

% Horizontální červená přerušovaná čára na -3 dB
line([min(omega) max(omega)], [target_dB target_dB], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);

% Černý bod na omega_c
plot(omega_c, mag_T0_dB(idx), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(omega_c, mag_T0_dB(idx) - 5, sprintf('\\omega_c = %.4f rad/s', omega_c), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', 'k');

% Nastavení grafu
title('Bodeho diagram', 'FontSize', 14);
legend({'T_0', 'S_0', '-3 dB', '\omega_c '}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;
hold off;


figure
hold on;
FRL = squeeze(freqresp(L_Hinf, f)); % Frekvenční odezva otevřené smyčky
iL = imag(FRL); % Imaginární část
rL = real(FRL); % Reálná část
plot(rL, iL, 'b--', 'LineWidth', 1.5); % Nyquistova křivka pro P0

% Výpočet kruhu1 a jeho umístění na bod -1
[A, i] = max(AW1S0_reg_matlab + AW2T0_reg_matlab); % Najít index pro maximum
% i = 5033; % Vybraný index (dle zadání)
kruh1 = -1 + AW1(i) * exp(1j * (-2 * pi:0.01:2 * pi)); % Definice kruhu1
plot(real(kruh1), imag(kruh1), 'g-', 'LineWidth', 1.5); % Kruh1 zeleně
% Výpočet kruhu2 a jeho umístění na příslušnou pozici na L
W2L = W2 * L_Hinf; % Funkce W2*L
FRW2L = freqresp(W2L, f); % Frekvenční odezva W2*L
AW2L = abs(squeeze(FRW2L)); % AFCH funkce W2*L
kruh2 = FRL(i) + AW2L(i) * exp(1j * (-2 * pi:0.01:2 * pi)); % Definice kruhu2
plot(real(kruh2), imag(kruh2), 'r-', 'LineWidth', 1.5); % Kruh2 červeně
% Přidání bodů středu kruhů
plot(-1, 0, 'gx', 'MarkerSize', 8, 'LineWidth', 1.5); % Střed kruhu1 (-1, 0)
text(-1.1, 0.1, 'r=|W1|', 'FontSize', 10, 'HorizontalAlignment', 'right', 'Color', 'g');
plot(real(FRL(i)), imag(FRL(i)), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5); % Střed kruhu2
text(real(FRL(i)) - 0.1, imag(FRL(i)) - 0.1, 'r=|W2L0|', 'FontSize', 10, 'HorizontalAlignment', 'right', 'Color', 'r');
% Nastavení grafu
title('Nyquist Diagram', 'FontSize', 14);
xlabel('Real Axis', 'FontSize', 12);
ylabel('Imaginary Axis', 'FontSize', 12);
legend({'L0', 'r=|W1|', 'r=|W2L0|'}, 'FontSize', 10, 'Location', 'Best');
grid on;
axis equal;
box on;
hold off;

%%
% Frekvenční rozsah
f = logspace(-2,4,10000);
figure
semilogx(f, AW1S0_reg_matlab + AW2T0_reg_matlab, 'b', 'LineWidth', 1.5); % Součet |W_1 S_0| + |W_2 T_0|
hold on;

semilogx(f, AW1S0_reg_matlab, 'r', 'LineWidth', 1.5); % |W_1 S_0|
semilogx(f, AW2T0_reg_matlab, 'g', 'LineWidth', 1.5); % |W_2 T_0|
% Najít frekvenci w, kde |W_1 S_0| + |W_2 T_0| dosahuje nejvyšší hodnoty
[max_value, max_index] = max(AW1S0_reg_matlab + AW2T0_reg_matlab); % Nejvyšší hodnota a její index
max_frequency = f(max_index); % Odpovídající frekvence

% Přidání vertikální osy
line([max_frequency max_frequency], [0 max_value], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5); % Vertikální čára
text(max_frequency, max_value, ...
    sprintf('\\omega = %.2f, max = %.4f', max_frequency, max_value), ...
    'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Nastavení grafu
title('Frekvenční charakteristika', 'FontSize', 14);
xlabel('Frekvence [rad/s]', 'FontSize', 12);
ylabel('Hodnota', 'FontSize', 12);
legend({'|W_1 S_0| + |W_2 T_0|','|W_1 S_0|', '|W_2 T_0|', '\omega_{max}'}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;


%% porovnani regulatoru

figure;
hold on
step(T_hinf);
step(L_uzavr);
title('Porovnání jednotkových skoků regulovaných systémů');
grid on;

T0_reg_matlab = T_hinf;

% Frekvenční odezva T0 a S0
[mag_T0, ~, omega] = bode(T0_reg_matlab, f); % Magnituda T0 (linear), fáze a frekvence
mag_T0 = squeeze(mag_T0);         % Převod na vektor
       % Převod na vektor

% Převod na decibely
mag_T0_dB = 20 * log10(mag_T0);

% Najít frekvenci, kde |T0| = -3 dB
target_dB = -3; % Hodnota v dB
[~, idx] = min(abs(mag_T0_dB - target_dB)); % Index nejbližší hodnoty k -3 dB
omega_c = omega(idx); % Frekvence

% Vykreslení grafu
figure;

% T0 modře
semilogx(omega, mag_T0_dB, 'b-', 'LineWidth', 1.5); % T0 funkce
hold on;

% Horizontální červená přerušovaná čára na -3 dB
line([min(omega) max(omega)], [target_dB target_dB], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);

% Černý bod na omega_c
plot(omega_c, mag_T0_dB(idx), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(omega_c, mag_T0_dB(idx) - 5, sprintf('\\omega_c = %.4f rad/s', omega_c), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', 'k');

% Nastavení grafu
title('Bodeho diagram', 'FontSize', 14);
legend({'T_0', '-3 dB', '\omega_c '}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;
hold off;


T0_reg_matlab = L_uzavr;

% Frekvenční odezva T0 a S0
[mag_T0, ~, omega] = bode(T0_reg_matlab, f); % Magnituda T0 (linear), fáze a frekvence
mag_T0 = squeeze(mag_T0);         % Převod na vektor
       % Převod na vektor

% Převod na decibely
mag_T0_dB = 20 * log10(mag_T0);

% Najít frekvenci, kde |T0| = -3 dB
target_dB = -3; % Hodnota v dB
[~, idx] = min(abs(mag_T0_dB - target_dB)); % Index nejbližší hodnoty k -3 dB
omega_c = omega(idx); % Frekvence

% Vykreslení grafu
figure;

% T0 modře
semilogx(omega, mag_T0_dB, 'b-', 'LineWidth', 1.5); % T0 funkce
hold on;

% Horizontální červená přerušovaná čára na -3 dB
line([min(omega) max(omega)], [target_dB target_dB], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);

% Černý bod na omega_c
plot(omega_c, mag_T0_dB(idx), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(omega_c, mag_T0_dB(idx) - 5, sprintf('\\omega_c = %.4f rad/s', omega_c), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', 'k');

% Nastavení grafu
title('Bodeho diagram', 'FontSize', 14);
legend({'T_0', '-3 dB', '\omega_c '}, ...
       'FontSize', 10, 'Location', 'Best');
grid on;
hold off;



% Nyquist Plot for Stability Check
figure;
nyquist(L_hinf);
title('Nyquist Diagram of Open-Loop Transfer Function');
grid on;

figure;
nyquiststability((L.Numerator{1}),(L.Denominator{1}),0.1)
title('Nyquist Diagram of Open-Loop Transfer Function');
grid on;




%% 1) Definice přenosů systému
s = tf('s');

% a) Přenos regulovaného objektu G(s):
G = P0_nominal;  % Očekává se, že máte tuhle proměnnou v workspace (např. 2. řád apod.)

% b) 1. regulátor (z PIDlab)
C1 = C;   % Regulátor navržený v PIDlab

% c) 2. regulátor (z hinfsyn)
C2 = K_Hinf;  % Regulátor navržený pomocí hinfsyn

%% 2) Otevřená a uzavřená smyčka
L1 = series(C1, G);   % Otevřená smyčka s 1. regulátorem
L2 = series(C2, G);   % Otevřená smyčka s 2. regulátorem

T1 = feedback(L1, 1); % Uzavřená smyčka (PIDlab)
T2 = feedback(L2, 1); % Uzavřená smyčka (hinfsyn)

%% 3) Zisková a fázová bezpečnost (gain margin, phase margin)
figure('Name','Bode L1 s marginy');
margin(L1); grid on;
set(findall(gcf,'Type','Line'),'LineWidth',1.5);  % Zesílení čar
title('Bodeho diagram L_1 (PIDlab)');

figure('Name','Bode L2 s marginy');
margin(L2); grid on;
set(findall(gcf,'Type','Line'),'LineWidth',1.5);
title('Bodeho diagram L_2 (hinfsyn)');

[Gm1, Pm1, ~, ~] = margin(L1);
[Gm2, Pm2, ~, ~] = margin(L2);

disp('=====  Zisková a fázová bezpečnost  =====');
disp('--- Regulátor 1 (PIDlab) ---');
fprintf('   Gain margin (GM) = %.3f\n', Gm1);
fprintf('   Phase margin (PM) = %.3f°\n', Pm1);

disp('--- Regulátor 2 (hinfsyn) ---');
fprintf('   Gain margin (GM) = %.3f\n', Gm2);
fprintf('   Phase margin (PM) = %.3f°\n', Pm2);

%% 4) Bezpečnost ve stabilitě s_m (nejmenší vzdálenost Nyquista od -1 + j0)
figure('Name','Nyquist L1');
nyquist(L1); 
hold on;
plot(-1,0,'rx','MarkerSize',8,'LineWidth',2);
set(findall(gcf,'Type','Line'),'LineWidth',1.5);
title('Nyquist L_1 (PIDlab)');
axis equal; grid on;

figure('Name','Nyquist L2');
nyquist(L2); 
hold on;
plot(-1,0,'rx','MarkerSize',8,'LineWidth',2);
set(findall(gcf,'Type','Line'),'LineWidth',1.5);
title('Nyquist L_2 (hinfsyn)');
axis equal; grid on;

% Výpočet s_m pro L1
[re1, im1] = nyquist(L1);
re1 = squeeze(re1);
im1 = squeeze(im1);

dist1 = sqrt((re1 + 1).^2 + (im1 - 0).^2);  % vzdálenost od (-1, 0)
s_m1 = min(dist1);

% Výpočet s_m pro L2
[re2, im2] = nyquist(L2);
re2 = squeeze(re2);
im2 = squeeze(im2);

dist2 = sqrt((re2 + 1).^2 + (im2 - 0).^2);
s_m2 = min(dist2);

disp(' ');
disp('=====  Stabilitní bezpečnost s_m  =====');
disp(['   s_m (PIDlab)  = ', num2str(s_m1)]);
disp(['   s_m (hinfsyn) = ', num2str(s_m2)]);

%% 5) Porovnání kvality řízení (přechodové charakteristiky)
figure('Name','Porovnani prechodovych odezev - PIDlab vs. hinfsyn');
step(T1, T2);
grid on;
set(findall(gcf,'Type','Line'),'LineWidth',1.5);
legend('T_1 (PIDlab)','T_2 (hinfsyn)','Location','Best');
title('Porovnání přechodových charakteristik uzavřených smyček');

% Vypsání detailů stepinfo
S1 = stepinfo(T1);
S2 = stepinfo(T2);

disp(' ');
disp('=====  Kvalita řízení - stepinfo  =====');
disp('--- Regulátor 1 (PIDlab) ---');
disp(S1);
disp('--- Regulátor 2 (hinfsyn) ---');
disp(S2);

%% 6) Volitelné: Porovnání Nyquist a Bode v jednom grafu
figure('Name','Porovnani Nyquist L1 a L2');
nyquist(L1, L2);
hold on; plot(-1,0,'rx','MarkerSize',8,'LineWidth',2);
set(findall(gcf,'Type','Line'),'LineWidth',1.5);
legend('L_1 (PIDlab)','L_2 (hinfsyn)','-1 + j0','Location','Best');
axis equal; grid on;
title('Porovnání Nyquistových charakteristik');

figure('Name','Porovnani Bode L1 a L2');
bode(L1, L2);
grid on;
set(findall(gcf,'Type','Line'),'LineWidth',1.5);
legend('T_1 (PIDlab)','T_2 (hinfsyn)','Location','Best');
title('Porovnání Bodeho diagramů');

disp('Hotovo! Skript úspěšně proběhl.');




