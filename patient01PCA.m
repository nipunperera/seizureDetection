close all
clear
%% Wavelet Transform based Seizure Detection

% The EEG data is acquired from CHB-MIT Scalp EEG Database
% Each recording is one hour long
% 256Hz sampling rate has been used 

% --------------------Seizure occurences for patient 01--------------------

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

% % Each channel is represented by a column for filtering 
% data = permute(records, [2, 1, 3]);
% data = data(:,promChannels,:);

% Indices corresponding to seizure events 
seizureInd = [900*3+(749:759) 900*4+(367:373) 900*15+(433:443)...
    900*16+(254:266) 900*18+(430:452) 900*21+(82:105) 900*26+(465:491)];

%%
% 

dirDec = 'r';         % Direction of decomposition
level  = 10;           % Level of decomposition
wname  = 'db4';      % Debauchies wavelet


decROW = mdwtdec(dirDec,records(:,:,1),level,wname);

cfs = cat(2,[decROW.cd{:},decROW.ca]);
