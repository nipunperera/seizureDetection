close all;

t1 = 0:0.01:50;
t2 = 50.01:0.01:100;
t = [t1 t2];
x = [a b];
a = sin(2*pi*0.2*t1);
b = sin(2*pi*0.4*t2);
plot([t1 t2], [a b])

dirDec = 'r';         % Direction of decomposition
level  = 1;           % Level of decomposition
wname  = 'sym8';      % Debauchies wavelet

[wpta,~,Falign] = modwpt([a b], 'timealign', true);

contour([t1 t2],Falign.*(1/0.01),abs(wpta).^2);
grid on;
xlabel('Time');
ylabel('Hz');
title('Time-Frequency Plot (Aligned)');

% decROW = mdwtdec(dirDec,[a b],level,wname);
% cfs = cat(2,[decROW.cd{:},decROW.ca]);

% A1_ROW  = mdwtrec(decROW,'a',1);
% D1_ROW  = mdwtrec(decROW,'d',1);
