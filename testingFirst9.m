%% Read EEG Data for patient 01
clear
parfor i =28:42
    fileName = sprintf('chb01/chb01_%d.edf', i);
    [~,records(:,:,i-27)] = edfread(fileName);
end

% Channels being used for detection
chansOfInterest = [2 6 22 9];
waveletPkts = [13 9 13 9];

% Each channel is represented by a column for filtering 
data = permute(records(chansOfInterest,:,:), [2, 1, 3]);

% Indices corresponding to seizure events 
seizureInd = [900*2+(749:759) 900*3+(367:373)];
% 900*14+(433:443) 900*15+(254:266) 900*17+(430:452) 900*20+(82:105) 900*25+(465:491)];

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

LPFiltered = filter(Hd, data);
records = permute(LPFiltered,[2 1 3]);
%%
numFeat = 2;
featMatrix = ones(numFeat*size(records, 1), 900, size(records, 3));       % Matrix to store the features 
                                    % First dimension - WP Coeffs of each
                                    % window
                                    % Last dimension - Record Number

wdSize = 1024;
level = 5;
wName = 'db4';

%% Extraction of features (Wavelet Packet Coefficients)
for i = 1:size(records, 3)
    for j = 1:size(records, 2) / wdSize
        for k = 1:size(records, 1)
            signal = records(k,1 + (wdSize * (j - 1)):wdSize * j, i);
            wpt = wpdec(signal,level,wName);    % Wavelet packet tree
            wpCoef = wpcoef(wpt, [level, waveletPkts(k) - 1]);
            featMatrix((2*k) - 1, j, i) = log(sqrt(mean(wpCoef.*wpCoef)));
            featMatrix((2*k), j, i) = var(wpCoef);
        end
    end
end
%%
featMatrix2D = reshape(featMatrix, size(featMatrix, 1), size(featMatrix, 2)*size(featMatrix, 3)); 