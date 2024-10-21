K_uncertain = ureal('K', 20, 'Percentage', 15);  
T1_uncertain = ureal('T1', 0.2, 'Percentage', 15); 
T2_uncertain = ureal('T2', 0.4, 'Percentage', 15);

s = tf('s');

P0_uncertain = K_uncertain / ((T1_uncertain * s + 1) * (T2_uncertain * s + 1));

S_uncertain = 1 / (1 + P0_uncertain);


W1 = (s+0.1)/(s+15) %allright i guess
bode(W1*S_uncertain);
grid on;
