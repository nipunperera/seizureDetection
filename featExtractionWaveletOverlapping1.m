%%---------------------Seizure detection in continuous EEG-----------------
% In this method, seizure events are detected and classified using discrete
% wavelet packet transform. 
% Seizure events of the ECG are identified by considering 4 second
% intervals in the EEG signal. Prominent channels and the wavelet packet
% are extracted by significance.m based on a one way ANOVA test

% Prominent channels for the first patient
% Channels 14, 15, 21, 22
% Corresponding wavelet packets 2, 2, 2, 2

% ----------------------------Seizure occurences---------------------------

% File Name: chb01_03.edf           (Record 3)
% Seizure Start Time: 2996 seconds      - 2993
% Seizure End Time: 3036 seconds        - 3033

% File Name: chb01_04.edf           (Record 4)
% Seizure Start Time: 1467 seconds      - 1464
% Seizure End Time: 1494 seconds        - 1491

% File Name: chb01_15.edf           (Record 15)
% Seizure Start Time: 1732 seconds      - 1729
% Seizure End Time: 1772 seconds        - 1769
%
% File Name: chb01_16.edf           (Record 16)
% Seizure Start Time: 1015 seconds      - 1012
% Seizure End Time: 1066 seconds        - 1063

% File Name: chb01_18.edf           (Record 18)
% Seizure Start Time: 1720 seconds      - 1717
% Seizure End Time: 1810 seconds        - 1807

% File Name: chb01_21.edf           (Record 21)
% Seizure Start Time: 327 seconds       - 324
% Seizure End Time: 420 seconds         - 417

% File Name: chb01_26.edf           (Record 26)
% Seizure Start Time: 1862 seconds
% Seizure End Time: 1963 seconds

clear
%% Indices corresponding to seizure events
seizureInd = [3597*2+(2993:3033) 3597*3+(1464:1491) 3597*14+(1729:1769)...
    3597*15+(1012:1063) 3597*17+(1717:1807)];
seizureInd = [];
%% Preparing feature vector and targets
featMatrix1 = generateFeatMat(19);
%featMatrix2 = generateFeatMat(5);
%featMatrix3 = generateFeatMat(10);
%featMatrix4 = generateFeatMat(15);

featMatrix = featMatrix1;
%featMatrix = cat(3, featMatrix1, featMatrix2, featMatrix3, featMatrix4);
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
% Randomly select 3000 samples from non-seizure data
randSampleInd = randi(size(featureVectorNonSeiz, 1), 1, 3000);
finalFeatureVector = [featureVectorSeiz;featureVectorNonSeiz(randSampleInd,:)];
finalTargetVector = [targetVectorSeiz;targetVectorNonSeiz(randSampleInd,:)];

SVMmodel = fitcsvm(finalFeatureVector, finalTargetVector, 'Standardize',...
    true,'KernelFunction','RBF','KernelScale','auto');

%% Read EEG Data for patient 01 (In blocks of 10 records)
function featMatrix = generateFeatMat(index)

parfor i = 1:1
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
disp('Computing wavelet features...')
for i = 1:size(records, 3)
    for j = 4:size(records, 2) / Fs
        if mod(j, 100) == 0
            fprintf('%d seconds processed....\n', j);
        end
        for k = 1:size(records, 1)
            signal = records(k,1 + (Fs * (j - 4)):Fs * j, i);
            wpt = wpdec(signal,level,wName);
            wpCoef = wpcoef(wpt, [level, waveletPkts(k) - 1]);
            featMatrix(k, j-3, i) = log(sqrt(mean(wpCoef.*wpCoef)));
            % featMatrix((2*k), j, i) = var(wpCoef);
        end
    end
end
end
