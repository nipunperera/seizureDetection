%%---------------Seizure Detection in Continuous EEG-----------------------
% This script uses the calcualted feature vector and train the system using
% Neural Networks

clear
%% Loading feature vectors
fileName = sprintf('WaveletFeatures/Patient01/featMatrix');
load(fileName);

%% Preparing feature vector and targets
featMatrix2D = reshape(featMatrix, size(featMatrix, 1), ...
    size(featMatrix, 2)*size(featMatrix, 3));

%% Indices corresponding to seizure events
seizureInd = [3597*2+(2993:3033) 3597*3+(1464:1491) 3597*14+(1729:1769)...
    3597*15+(1012:1063) 3597*17+(1717:1807)];

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
randSampleInd = randi(size(featureVectorNonSeiz, 1), 1, 10000);
finalFeatureVector = [featureVectorSeiz;featureVectorNonSeiz(randSampleInd,:)];
finalTargetVector = [targetVectorSeiz;targetVectorNonSeiz(randSampleInd,:)];

finalFeatureVector = finalFeatureVector';
finalTargetVector = finalTargetVector';

% Create a Fitting Network
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);

% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
 
% Train the Network
[net,tr] = train(net, finalFeatureVector, finalTargetVector);
 
% Test the Network
outputs = net(finalFeatureVector);
errors = gsubtract(outputs, finalTargetVector);
performance = perform(net, finalTargetVector, outputs);
 
% View the Network
% view(net)

%% Predicting and calculating accuracy

seizureInd = [3597*2+(2993:3033) 3597*3+(1464:1491) 3597*14+(1729:1769)...
    3597*15+(1012:1063) 3597*17+(1717:1807)];

label = sim(net,featMatrix2D);
label = hardlim(label - 0.5);
trueSeizEvents = length(seizureInd);
detectedSeizEvents = sum(label(seizureInd));

% label(seizureInd) = 0;
% falseDetections = find(label);

accuracy = detectedSeizEvents * 100 / trueSeizEvents;

display(accuracy);