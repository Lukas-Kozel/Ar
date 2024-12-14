function [] = W1_W2_test(W1, W1plot, W2, W2plot, S, T)

% W1*S - budeme testovat, zda ||W1 * S || inf < 1
W1S = minreal(W1*S);

f = logspace(-2,4,10000);
FRS = freqresp(S, f); % frekvenční odezva citlivostní funkce
FRW1 = freqresp(W1, f); % frekvenční odezva váhové funkce
FRW1m1 = freqresp(1/W1, f); % frekvenční odezva funkce 1/W1
FRW1S = freqresp(W1S, f); % frekvenční odezva funkce S*W1
AS = (abs(squeeze(FRS))); % AFCH citlivostní funkce
AW1 = (abs(squeeze(FRW1))); % AFCH váhové funkce
AW1m1 = (abs(squeeze(FRW1m1))); % AFCH funkce 1/W1
AW1S = (abs(squeeze(FRW1S))); % AFCH funkce S*W1


fprintf('-------------------------------------------------\n');
fprintf('Normy Systému:\n');


% Zobrazení normy W1*S
normW1S = norm(W1S, inf);
fprintf('||W_1 * S||_inf = %.4f\n', normW1S);

normAW1S = norm(AW1S, inf);
fprintf('ABS-> || |W_1 * S| ||_inf = %.4f\n\n', normAW1S);

if (W1plot)
    figure
    semilogx(f, AS)
    hold on
    semilogx(f, AW1m1)
    title('Frekvenční charakteristika citlivostní funkce S a W_1^{-1}')
    legend('S','W_1^{-1}')
    grid
    
    figure
    semilogx(f, AW1S)
    title('Frekvenční charakteristika funkce S*W_1')
    legend('S*W_1')
    grid
    
    figure
    bodemag(S)
    hold on
    bodemag(1/W1)
    title('Bodeho diagram funkce S a W_1^{-1}')
    legend('S','W_1^{-1}')
    grid
  
    figure
    bodemag(W1S)
    title('Bodeho diagram funkce S*W_1')
    legend('S*W_1')
    grid
end

%---------------------------------------------------------
% Váhová funkce W2, robustní stabilita - test
W2T = minreal(W2*T);

FRT = freqresp(T, f); % frekvenční odezva komplementární citlivosti
FRW2 = freqresp(W2, f); % frekvenční odezva váhové funkce
FRW2m1 = freqresp(1/W2, f); % frekvenční odezva funkce 1/W2
FRW2T = freqresp(W2T, f); % frekvenční odezva funkce T*W2
AT = (abs(squeeze(FRT))); % AFCH komplementární citlivosti
AW2 = (abs(squeeze(FRW2))); % AFCH váhové funkce
AW2m1 = (abs(squeeze(FRW2m1))); % AFCH funkce 1/W2
AW2T = (abs(squeeze(FRW2T))); % AFCH funkce T*W2

% Zobrazení normy T*W2
normW2T = norm(W2T, inf);
fprintf('||W_2 * T||_inf = %.4f\n', normW2T);

% Zobrazení normy T*W2
normAW2T = norm(AW2T, inf);
fprintf('ABS-> || |W_2 * T| ||_inf = %.4f\n\n', normAW2T);



if (W2plot)
    figure
    semilogx(f, AT)
    hold on
    semilogx(f, AW2m1)
    title('Frekvenční charakteristika komplementární citlivosti T a W_2^{-1}')
    legend('T','W_2^{-1}')
    grid
    
    figure
    semilogx(f, AW2T)
    title('Frekvenční charakteristika funkce T*W_2')
    legend('T*W_2')
    grid
    
    figure
    bodemag(T)
    hold on
    bodemag(1/W2)
    title('Bodeho diagram funkce T a W_2^{-1}')
    legend('T','W_2^{-1}')
    grid
    
    figure
    bodemag(W2T)
    title('Bodeho diagram funkce T*W_2')
    legend('T*W_2')
    grid
end

% % Robustní kvalita - pozor nejsou absolutní hodnoty
% normW1S_W2T = norm(W1S + W2T, inf);
% fprintf('||W_1*S + W_2*T||_inf = %.4f\n', normW1S_W2T);

% Robustní kvalita řízení - správně
normAW1S_AW2T = norm(AW1S + AW2T, inf);
fprintf('|| |W_1*S| + |W_2*T| ||_inf = %.4f\n\n', normAW1S_AW2T);


figure
bodemag(W2T + W1S)
title('Bodeho diagram funkce W_2*T + W_1*S')
grid


figure
semilogx(f, AW1S + AW2T)
title('Frekvenční charakteristika |W_1S| + |W_2T|')
legend('|W_1S| + |W_2T|')
grid

figure
semilogx(f, 20*log10(AW1S + AW2T));
title('Logaritmická frekvenční charakteristika |W_1S| + |W_2T| v dB')
legend('|W_1S| + |W_2T|')
grid

end
