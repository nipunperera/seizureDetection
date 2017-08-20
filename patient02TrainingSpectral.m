close all
clear

% ----------------------------Seizure occurences---------------------------

% File Name: chb03_01.edf
% File Start Time: 13:23:36
% File End Time: 14:23:36
% Number of Seizures in File: 1
% Seizure Start Time: 362 seconds  - 91
% Seizure End Time: 414 seconds - 104

% File Name: chb03_02.edf
% File Start Time: 14:23:39
% File End Time: 15:23:39
% Number of Seizures in File: 1
% Seizure Start Time: 731 seconds - 183
% Seizure End Time: 796 seconds - 199

% File Name: chb03_03.edf
% File Start Time: 15:23:47
% File End Time: 16:23:47
% Number of Seizures in File: 1
% Seizure Start Time: 432 seconds - 108
% Seizure End Time: 501 seconds - 125

% File Name: chb03_04.edf
% File Start Time: 16:23:54
% File End Time: 17:23:54
% Number of Seizures in File: 1
% Seizure Start Time: 2162 seconds - 541
% Seizure End Time: 2214 seconds - 554

% File Name: chb03_05.edf
% File Start Time: 01:51:23
% File End Time: 2:51:23
% Number of Seizures in File: 1
% Seizure Start Time: 1982 seconds - 496
% Seizure End Time: 2029 seconds - 507

% File Name: chb03_06.edf
% File Start Time: 04:51:45
% File End Time: 5:51:45
% Number of Seizures in File: 1
% Seizure Start Time: 1725 seconds - 437
% Seizure End Time: 1778 seconds - 445

% Extraction of spatial features
% Prominent channels are selected
% 1, 2, 3, 4, 5, 6, 9, 13, 19, 20, 21, 22, 23

% Indices corresponding to seizure events
seizureInd = [(91:104) 900*1+(183:199) 900*2+(108:125) 900*3+(541:554)...
    900*4+(496:507) 900*5+(437:445)];

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
% newX = featureVector(:,1:900)';
% [label,score] = predict(SVMmodel,newX);


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
    fileName = sprintf('chb03/chb03_%d.edf', i+index);
    [~,records(:,:,i)] = edfread(fileName);
end

% Prominent channels for patient 01
promChannels = [1, 2, 3, 4, 5, 6, 9, 13, 19, 20, 21, 22, 23];

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