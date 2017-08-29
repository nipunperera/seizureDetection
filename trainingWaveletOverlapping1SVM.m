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
%% Loading feature vectors
fileName = sprintf('WaveletFeatures/Patient01/featMatrix');
load(fileName);

%% Preparing feature vector and targets
featMatrix2D = reshape(featMatrix, size(featMatrix, 1), ...
    size(featMatrix, 2)*size(featMatrix, 3));

%% Indices corresponding to seizure events
seizureInd = [3597*2+(2993:3033) 3597*15+(1012:1063) 3597*17+(1717:1807)];

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

%% Predicting and calculating accuracy

seizureInd = [3597*3+(1464:1491) 3597*14+(1729:1769)];

[label,score] = predict(SVMmodel,featMatrix2D(:,:)');

trueSeizEvents = length(seizureInd);
detectedSeizEvents = sum(label(seizureInd));

label(seizureInd) = 0;
falseDetections = find(label);

accuracy = detectedSeizEvents * 100 / trueSeizEvents;

display(accuracy);

