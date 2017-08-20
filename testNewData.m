parfor i = 10:16
    fileName = sprintf('chb01/chb01_%d.edf', i);
    [~,recordsn(:,:,i-9)] = edfread(fileName);
end

% Each channel is represented by a column for filtering 
datan = permute(recordsn, [2, 1, 3]);

% Create low pass filter 
Fs = 256;  % Sampling Frequency

N      = 50;      % Order
Fc     = 25;      % Cutoff Frequency
DpassU = 0.01;    % Upper Passband Ripple
DpassL = 0.01;    % Lower Passband Ripple
DstopU = 0.0001;  % Upper Stopband Attenuation
DstopL = 0.0001;  % Lower Stopband Attenuation

% Calculate the coefficients using the FIRCLS function.
b  = fircls(N, [0 Fc Fs/2]/(Fs/2), [1 0], [1+DpassU DstopU], [1-DpassL ...
            -DstopL]);
Hd = dfilt.dffir(b);

LPFilteredn = filter(Hd, datan);
recordsn = permute(LPFilteredn,[2 1 3]);

rec15 = recordsn(:,:,6);
rec16 = recordsn(:,:,7);

matrix15 = ones(32, 900, 23);
matrix16 = ones(32, 900, 23);

wdSize = 1024;
level = 5;
wName = 'db4';

%% For third record
for j = 1:size(rec15, 2) / wdSize
    for k = 1:size(rec15, 1)
        signal = rec15(k,1 + (wdSize * (j - 1)):wdSize * j);
        wpt = wpdec(signal,level,wName);    % Wavelet packet tree
        for m = 1:2^level
            wpCoef = wpcoef(wpt, [level, m - 1]);
            matrix15(m, j, k) = log(sqrt(mean(wpCoef.*wpCoef)));
        end
    end
end

%% For third record
for j = 1:size(rec16, 2) / wdSize
    for k = 1:size(rec16, 1)
        signal = rec16(k,1 + (wdSize * (j - 1)):wdSize * j);
        wpt = wpdec(signal,level,wName);    % Wavelet packet tree
        for m = 1:2^level
            wpCoef = wpcoef(wpt, [level, m - 1]);
            matrix16(m, j, k) = log(sqrt(mean(wpCoef.*wpCoef)));
        end
    end
end
%%
featureVector15 = [];
featureVec16 = [];
for i = 1:length(promChan)
    featureVector15 = [featureVector15; matrix15(promWvPkt(i),:,promChan(i))];
    featureVec16 = [featureVec16; matrix16(promWvPkt(i),:,promChan(i))];
end