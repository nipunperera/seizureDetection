% This script trains a classifier based on spectral features using the 
% saved feature vector. 
% This can be used to change the amount of training data and evaluate the
% accuracy of the classifier

clear 
close all

load spectralFeatureVectorPatient3
 
% Indices corresponding to seizure events
% Patient 01
% seizureInd = [900*2+(749:759) 900*3+(367:373) 900*14+(433:443)...
%     900*15+(254:266) 900*17+(430:452)];

% Patient 03
seizureInd = [900*1+(183:199) 900*4+(496:507) 900*2+(108:125) 900*3+(541:554)];
% seizureInd = [(91:104) 900*5+(437:445)];

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


%% Predicting and calculating accuracy

% Patient 01
% seizureInd = [900*2+(749:759) 900*3+(367:373) 900*14+(433:443)...
%     900*15+(254:266) 900*17+(430:452)];

% Patient 03
% seizureInd = [ 900*1+(183:199) 900*2+(108:125) 900*3+(541:554) 900*4+(496:507)];
seizureInd = [(91:104) 900*5+(437:445)];

[label,score] = predict(SVMmodel,featureVector(:,:)');

trueSeizEvents = length(seizureInd);
detectedSeizEvents = sum(label(seizureInd));

% label(seizureInd) = 0;
% falseDetections = find(label);

accuracy = detectedSeizEvents * 100 / trueSeizEvents;

display(accuracy);
