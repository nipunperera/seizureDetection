close all
clear

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

% Indices corresponding to seizure events
seizureInd = [900*2+(749:759) 900*3+(367:373) 900*17+(430:452)];
% 900*3+(367:373) 900*15+(254:266) 900*17+(430:452)

% Reading data and generating feature matrices
featureVector1 = generateSpectralFeat(0);
featureVector2 = generateSpectralFeat(10);

featureVector = [featureVector1 featureVector2];

%% Preparing the training dataset and target vectors
targetVectorNonSeiz = zeros(1, size(featureVector, 2) - length(seizureInd))';
targetVectorSeiz = ones(1, length(seizureInd))';

featureVectorSeiz = featureVector(:,seizureInd);
featureVectorSeiz = featureVectorSeiz';
featureVectorNonSeiz = featureVector;
featureVectorNonSeiz(:,seizureInd) = [];
featureVectorNonSeiz = featureVectorNonSeiz';

%% Create and train SVM model
% Randomly select 1000 samples from non-seizure data
randSampleInd = randi(size(featureVectorNonSeiz, 1), 1, 1000);
finalFeatureVector = [featureVectorSeiz;featureVectorNonSeiz(randSampleInd,:)];
finalTargetVector = [targetVectorSeiz;targetVectorNonSeiz(randSampleInd,:)];

SVMmodel = fitcsvm(finalFeatureVector, finalTargetVector, 'Standardize',...
    true,'KernelFunction','RBF','KernelScale','auto');

%%
% Testing new data
newX = featureVector(:,1:900)';
[label,score] = predict(SVMmodel,newX);


% %% Visualization of data points
% figure;
% plot(featureVector(seizureInd(1:18), 1),featureVector(seizureInd(1:18), 17),'r.','MarkerSize',10)
% hold on
% plot(featureVector(101:200, 1),featureVector(101:200, 17),'b.','MarkerSize',10)
% xlim([0 10000])
% ylim([0 500])
% 
% xlabel('Energy of 0 - 16Hz Band')
% ylabel('Energy of 16 - 25Hz Band')
% legend('Seizure', 'Non-Seizure')


%% Read EEG Data for patient 01 (In blocks of 10 records)
% This functions reads blocks of 10 records and a low pass filtering is 
% applied and only channels of interest specified by the user is considered

function featureVector = generateSpectralFeat(index)

parfor i = 1:10
    fileName = sprintf('chb01/chb01_%d.edf', i+index);
    [~,records(:,:,i)] = edfread(fileName);
end

% Prominent channels for patient 01
promChannels = [1, 5, 6, 9, 10, 1, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23];

% Each channel is represented by a column for filtering
data = permute(records, [2, 1, 3]);
data = data(:,promChannels,:);

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

featureVector = [];

% Perform Short Time Fourier Transform for time intervals of 4s
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
end