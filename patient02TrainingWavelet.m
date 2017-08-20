%%---------------------Seizure detection in continuous EEG-----------------
% In this method, seizure events are detected and classified using discrete
% wavelet packet transform. 
% Seizure events of the ECG are identified by considering 4 second
% intervals in the EEG signal. Prominent channels and the wavelet packet
% are extracted by significance.m based on a one way ANOVA test

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

% Prominent channels for the first patient
% Channels 2, 19, 20, 21
% Corresponding wavelet packets 2, 2, 2, 2

clear

% Indices corresponding to seizure events
seizureInd = [(91:104) 900*1+(183:199) 900*2+(108:125) 900*3+(541:554)...
    900*4+(496:507) 900*5+(437:445)];

%% Preparing feature vector and targets
featMatrix1 = generateFeatMat(0);
featMatrix2 = generateFeatMat(10);

featMatrix = cat(3, featMatrix1, featMatrix2);
featMatrix2D = reshape(featMatrix, size(featMatrix, 1), ...
    size(featMatrix, 2)*size(featMatrix, 3));

% Preparing the training dataset and target vectors
targetVectorNonSeiz = zeros(1, size(featMatrix2D, 2) - length(seizureInd))';
targetVectorSeiz = ones(1, length(seizureInd))';

featureVectorSeiz = featMatrix2D(:,seizureInd);
featureVectorSeiz = featureVectorSeiz';

featureVectorNonSeiz = featMatrix2D;
featureVectorNonSeiz(:,seizureInd) = [];
featureVectorNonSeiz = featureVectorNonSeiz';

% Create and train SVM model
% Randomly select 1000 samples from non-seizure data
randSampleInd = randi(size(featureVectorNonSeiz, 1), 1, 1000);
finalFeatureVector = [featureVectorSeiz;featureVectorNonSeiz(randSampleInd,:)];
finalTargetVector = [targetVectorSeiz;targetVectorNonSeiz(randSampleInd,:)];

SVMmodel = fitcsvm(finalFeatureVector, finalTargetVector, 'Standardize',...
    true,'KernelFunction','RBF','KernelScale','auto');

%% Read EEG Data for patient 01 (In blocks of 10 records)
function featMatrix = generateFeatMat(index)

parfor i = 1:10
    fileName = sprintf('chb03/chb03_%d.edf', i+index);
    [~,records(:,:,i)] = edfread(fileName);
end

% Channels being used for detection
chansOfInterest = [2 19 20 21];
waveletPkts = [2 2 2 2];

% Each channel is represented by a column for filtering
data = permute(records(chansOfInterest,:,:), [2, 1, 3]);

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

wdSize = 1024;
level = 5;
wName = 'db4';

% featMatrix = ones(numFeat*size(records, 1), size(records, 2)/wdSize,...
%     size(records, 1));      % Matrix to store the features
%                             % First dimension - WP Coeffs of each window
%                             % Last dimension - Record Number


% Extraction of features (Wavelet Packet Coefficients)
for i = 1:size(records, 3)
    for j = 1:size(records, 2) / wdSize
        for k = 1:size(records, 1)
            signal = records(k,1 + (wdSize * (j - 1)):wdSize * j, i);
            wpt = wpdec(signal,level,wName);
            wpCoef = wpcoef(wpt, [level, waveletPkts(k) - 1]);
            featMatrix(k, j, i) = log(sqrt(mean(wpCoef.*wpCoef)));
            % featMatrix((2*k), j, i) = var(wpCoef);
        end
    end
end
end