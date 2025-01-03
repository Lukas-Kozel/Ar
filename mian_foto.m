% close all; clear all; clc;
perc = 40;
% r = ureal('r', 0.25, 'Mode', 'Range', 'Percentage', perc);
% P = tf([r], [0.0475 0.1 1.5*10*0.05]);
%C = 1;% C = tf([20 20*25], [1 300]);
% P0 = tf(P.Nominal);

K_uncertain = ureal('K', 20, 'Percentage', 15);  
T1_uncertain = ureal('T1', 0.2, 'Percentage', 15); 
T2_uncertain = ureal('T2', 0.4, 'Percentage', 15);

s = tf('s');


P0_nominal = 20/((0.2*s+1)*(0.4*s+1))
P0 = K_uncertain / ((T1_uncertain * s + 1) * (T2_uncertain * s + 1));
P0=P0_nominal

% uzavreny system s nominalni prenosem
% T = P0*C/(1 + P0*C);
T = P0 /(1+P0);

% test nominalni stability
pole(T); % stabilni system

% vahova funkce W1, nominalni kvalita - test
% W1 = tf([10], [1, 2, 10]);

% citlivostni funkce
% S = minreal(1/(1+P0*C));
S =  1 / (1 + P0);

% W1*S - budeme testovat, zda ||W1 * S || inf < 1
% W1S = minreal(W1*S);
W1S = W1*S;

f = logspace(-2,4,10000);
FRS = freqresp(S, f); % frekvencni odezva citlivostni funkce
FRW1 = freqresp(W1, f); % frekvencni odezva vahove funkce
FRW1m1 = freqresp(1/W1, f); % frekvencni odezva funkce 1/W1
FRW1S = freqresp(W1S, f); % frekvencni odezva funkce S*W1
AS = (abs(squeeze(FRS))); % AFCH citlivostni funkce
AW1 = (abs(squeeze(FRW1))); % AFCH vahove funkce
AW1m1 = (abs(squeeze(FRW1m1))); % AFCH funkce 1/W1
AW1S = (abs(squeeze(FRW1S))); % AFCH funkce S*W1

figure
semilogx(f, AS)
hold on
semilogx(f, AW1m1)
legend('S','W_1^{-1}')
grid

figure
semilogx(f, AW1S)
legend('S*W_1')
grid

figure
bodemag(S)
hold on
bodemag(1/W1)
legend('S','W_1^{-1}')
grid

norm(W1S,inf)
norm(AW1S,inf)

figure
bodemag(W1S)
legend('S*W_1')
grid

%%
% vahova funkce W2, robustni stabilita - test
%W2 = tf(perc/100);
W2T = minreal(W2*T);

FRT = freqresp(T, f); % frekvencni odezva komplementarni citlivosti
FRW2 = freqresp(W2, f); % frekvencni odezva vahove funkce
FRW2m1 = freqresp(1/W2, f); % frekvencni odezva funkce 1/W2
FRW2T = freqresp(W2T, f); % frekvencni odezva funkce T*W2
AT = (abs(squeeze(FRT))); % AFCH komplementarni citlivosti
AW2 = (abs(squeeze(FRW2))); % AFCH vahove funkce
AW2m1 = (abs(squeeze(FRW2m1))); % AFCH funkce 1/W2
AW2T = (abs(squeeze(FRW2T))); % AFCH funkce T*W2

figure
semilogx(f, AT)
hold on
semilogx(f, AW2m1)
legend('T','W_2^{-1}')
grid

figure
semilogx(f, AW2T)
legend('T*W_2')
grid

figure
bodemag(T)
hold on
bodemag(1/W2)
legend('T','W_2^{-1}')
grid

norm(AW2T, inf)
norm(W2T, inf)

figure
bodemag(W2T)
legend('T*W_2')
grid
%%

% robustni kvalita - pozor nejsou absolutni hodnoty
norm(W1S + W2T, inf)
figure
bodemag(W2T + W1S)
grid

% robustni kvalita rizeni - spravne
norm(AW1S + AW2T, inf)
figure
semilogx(f, AW1S + AW2T)
legend('|W_1S| + |W_2T|')
grid

figure
semilogx(f, 20*log10(AW1S + AW2T));
legend('|W_1S| + |W_2T|');
grid

%%
% Definice parametrů
L = minreal(P0 * C); % Přenos otevřené smyčky
FRL = squeeze(freqresp(L, f)); % Frekvenční odezva otevřené smyčky
iL = imag(FRL); % Imaginární část
rL = real(FRL); % Reálná část

% Nyquistův diagram
figure;
plot(rL, iL, 'b--', 'LineWidth', 1.5); % Nyquistova křivka pro P0
hold on;

% Výpočet kruhu1 a jeho umístění na bod -1
[A, i] = max(AW1S + AW2T); % Najít index pro maximum
i = 5033; % Vybraný index (dle zadání)
kruh1 = -1 + AW1(i) * exp(1j * (-2 * pi:0.01:2 * pi)); % Definice kruhu1
plot(real(kruh1), imag(kruh1), 'g-', 'LineWidth', 1.5); % Kruh1 zeleně

% Výpočet kruhu2 a jeho umístění na příslušnou pozici na L
W2L = W2 * L; % Funkce W2*L
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
hold off;

