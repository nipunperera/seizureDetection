close all
clear

% [~, records(:,:,1)] = edfread('chb01/chb01_01.edf');
% [~, records(:,:,2)] = edfread('chb01/chb01_02.edf');
% [~, records(:,:,3)] = edfread('chb01/chb01_03.edf');
% [~, records(:,:,4)] = edfread('chb01/chb01_04.edf');
% [~, records(:,:,5)] = edfread('chb01/chb01_05.edf');
% [~, records(:,:,6)] = edfread('chb01/chb01_06.edf');
% [~, records(:,:,7)] = edfread('chb01/chb01_07.edf');
% [~, records(:,:,8)] = edfread('chb01/chb01_08.edf');
% [~, records(:,:,9)] = edfread('chb01/chb01_09.edf');
% [~, records(:,:,10)] = edfread('chb01/chb01_10.edf');
% [~, records(:,:,11)] = edfread('chb01/chb01_11.edf');
% [~, records(:,:,12)] = edfread('chb01/chb01_12.edf');
% [~, records(:,:,13)] = edfread('chb01/chb01_13.edf');
% [~, records(:,:,14)] = edfread('chb01/chb01_14.edf');
% [~, records(:,:,15)] = edfread('chb01/chb01_15.edf');
% [~, records(:,:,16)] = edfread('chb01/chb01_16.edf');
% [~, records(:,:,17)] = edfread('chb01/chb01_17.edf');
% [~, records(:,:,18)] = edfread('chb01/chb01_18.edf');
% [~, records(:,:,19)] = edfread('chb01/chb01_19.edf');
% [~, records(:,:,20)] = edfread('chb01/chb01_20.edf');
% [~, records(:,:,21)] = edfread('chb01/chb01_21.edf');
% [~, records(:,:,22)] = edfread('chb01/chb01_22.edf');
% [~, records(:,:,23)] = edfread('chb01/chb01_23.edf');
% [~, records(:,:,24)] = edfread('chb01/chb01_24.edf');
% [~, records(:,:,25)] = edfread('chb01/chb01_25.edf');
% [~, records(:,:,26)] = edfread('chb01/chb01_26.edf');
% [~, records(:,:,27)] = edfread('chb01/chb01_27.edf');
% [~, records(:,:,28)] = edfread('chb01/chb01_29.edf');
% [~, records(:,:,29)] = edfread('chb01/chb01_30.edf');
% [~, records(:,:,30)] = edfread('chb01/chb01_31.edf');
% [~, records(:,:,31)] = edfread('chb01/chb01_32.edf');
% [~, records(:,:,32)] = edfread('chb01/chb01_33.edf');
% [~, records(:,:,33)] = edfread('chb01/chb01_34.edf');
% [~, records(:,:,34)] = edfread('chb01/chb01_36.edf');
% [~, records(:,:,35)] = edfread('chb01/chb01_37.edf');
% [~, records(:,:,36)] = edfread('chb01/chb01_38.edf');
% [~, records(:,:,37)] = edfread('chb01/chb01_39.edf');
% [~, records(:,:,38)] = edfread('chb01/chb01_40.edf');
% [~, records(:,:,39)] = edfread('chb01/chb01_41.edf');
% [~, records(:,:,40)] = edfread('chb01/chb01_42.edf');
% [~, records(:,:,41)] = edfread('chb01/chb01_43.edf');
% [~, records(:,:,42)] = edfread('chb01/chb01_46.edf');

% ----------------------------Seizure occurences---------------------------

% File Name: chb01_03.edf           (Record 3)
% Seizure Start Time: 2996 seconds  - 749
% Seizure End Time: 3036 seconds    - 759

% File Name: chb01_04.edf           (Record 4)
% Seizure Start Time: 1467 seconds  - 367
% Seizure End Time: 1494 seconds    - 373

% File Name: chb01_15.edf           (Record 15)
% Seizure Start Time: 1732 seconds  - 433
% Seizure End Time: 1772 seconds    - 443
% 
% File Name: chb01_16.edf           (Record 16)
% Seizure Start Time: 1015 seconds  - 254
% Seizure End Time: 1066 seconds    - 266

% File Name: chb01_18.edf           (Record 18)
% Seizure Start Time: 1720 seconds  - 430
% Seizure End Time: 1810 seconds    - 452

% File Name: chb01_21.edf           (Record 21)
% Seizure Start Time: 327 seconds   - 82
% Seizure End Time: 420 seconds     - 105

% File Name: chb01_26.edf           (Record 26)
% Seizure Start Time: 1862 seconds  - 465
% Seizure End Time: 1963 seconds    - 491



% Extraction of spatial features 
% Prominent channels are selected
% Fp1-F7, Fp1-Fp3, F3-C3, Fp2-F4, F4-C4, C4-P4, P4-O2, Fp2-F8, F8-T8,
% T8-P8, P8-O2, Fz-Cz, Cz-Pz, FT9-FT10, FT10-T8, T8-P8 
% 1, 5, 6, 9, 10, 1, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23

% Prominent channels for patient 01
promChannels = [1, 5, 6, 9, 10, 1, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23];

%--------------------Read EEG Data for patient 01--------------------------
% Read first batch of data
parfor i = 1:9
    fileName = sprintf('chb01/chb01_0%d.edf', i);
    [~,records(:,:,i)] = edfread(fileName);
end


data = permute(records, [2, 1, 3]);
data = data(:,promChannels,:);

% Indices corresponding to seizure events 
seizureInd = [900*3+(749:759) 900*4+(367:373) 900*15+(433:443)...
    900*16+(254:266) 900*18+(430:452) 900*21+(82:105) 900*26+(465:491)];

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

%%
featureVector = [];

% Perform Short Time Fourier Transform for time intervals of 2s
for i = 1:size(LPFiltered, 3)
    for j = 1:size(LPFiltered, 1)/1024
        [pxx, f] = periodogram(LPFiltered(1 + (1024*(j - 1)):1024*j,:,i),...
            [], [], 256);
        power0_16 = bandpower(pxx(1:65, :), f(1:65), 'psd');
        power16_25 = bandpower(pxx(65:101, :), f(65:101), 'psd');
        powerVector = [power0_16 power16_25];
        featureVector = [featureVector powerVector'];
    end
end

% Preparing the training dataset and target vectors
targetVectorNonSeiz = zeros(1, 8100 - length(seizureInd(1:18)))';
targetVectorSeiz = ones(1, length(seizureInd(1:18)))';

featureVectorSeiz = featureVector(:,seizureInd(1:18));
featureVectorSeiz = featureVectorSeiz';
featureVectorNonSeiz = featureVector;
featureVectorNonSeiz(:,seizureInd(1:18)) = [];
featureVectorNonSeiz = featureVectorNonSeiz';

% Create and train SVM model
% Randomly select 200 samples from non-seizure data
randSampleInd = randi(size(featureVectorNonSeiz, 1), 1, 200);
finalFeatureVector = [featureVectorSeiz;featureVectorNonSeiz(randSampleInd,:)];
finalTargetVector = [targetVectorSeiz;targetVectorNonSeiz(randSampleInd,:)];

SVMmodel = fitcsvm(finalFeatureVector, finalTargetVector, 'Standardize',...
    true,'KernelFunction','RBF','KernelScale','auto');

%%
% Testing new data
newX = featureVector(seizureInd(1:18), :);
[label,score] = predict(SVMmodel,newX);

% Visualization of data points
figure;
plot(featureVector(seizureInd(1:18), 1),featureVector(seizureInd(1:18), 17),'r.','MarkerSize',10)
hold on
plot(featureVector(101:200, 1),featureVector(101:200, 17),'b.','MarkerSize',10)
xlim([0 10000])
ylim([0 500])

xlabel('Energy of 0 - 16Hz Band')
ylabel('Energy of 16 - 25Hz Band')
legend('Seizure', 'Non-Seizure')
