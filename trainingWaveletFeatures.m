%%---------------------Seizure detection in continuous EEG-----------------
% In this method, seizure events are detected and classified using discrete
% wavelet packet transform. 
% Seizure events of the ECG are identified by considering 4 second
% intervals in the EEG signal. Prominent channels and the wavelet packet
% are extracted by significance.m based on a one way ANOVA test

% Prominent channels for the first patient
% Channels 14, 15, 21, 22
% Corresponding wavelet packets 2, 2, 2, 2

clear
%% Indices corresponding to seizure events
seizureInd = [900*2+(749:759) 900*3+(367:373) 900*14+(433:443)...
    900*15+(254:266) 900*17+(430:452)];

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
    fileName = sprintf('chb01/chb01_%d.edf', i+index);
    [~,records(:,:,i)] = edfread(fileName);
end

% Channels being used for detection
chansOfInterest = [14 15 21 22];
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
