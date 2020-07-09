%the struct runData holds data about patients and where the different event
%types are stored
data_p_path = 'E:\Data_p\';

patients = {'p488','p487','p489','p486','p499','p485','p490','p496','p498','p510','p505','p515','p510','p520','p497'};
expNames = {'EXP4','EXP3','EXP3','EXP8','EXP8','EXP8','EXP3','EXP8','EXP3','EXP7','EXP1','EXP9','EXP1','EXP3','EXP2'};
sleepScoreFileName = {'sleepScore_manualValidated_p488_4_LPHG2',...
    'sleepScore_manualValidated_p487_3_ROF2',...
    'sleepScore_manualValidated_p489_3_RPHG3',...
    'sleepScore_manualValidated_p486_8_RA1',...
    'sleepScore_manualValidated_p499_8_RPT2',...
    'sleepScore_manualValidated_p485_8_ROF7_0016',...
    'sleepScore_manualValidated_p490_3_RSTG1',...
    'sleepScore_manualValidated_p496_8_RSO4',...
    'sleepScore_manualValidated_p498_3_ROF1',...
    'sleepScore_manualValidated_p510_7_RPSM2',...
    'sleepScore_manualValidated_p505_1_RSTG1',...
    'sleepScore_manualValidated_p515_9_LOF3',...
    'sleepScore_manualValidated_p510_1_RPSM2',...
    'sleepScore_manualValidated_p520_3_CSC26',...
    'sleepScore_manualValidated_p497_2_RPHG4'};


%note macro montage is not necessary to run AnalyzeCoupling
% macroMontageNames = {'MM488','MM487','MM489','MM486','MM499','MM485','MM490','MM496','MM498','MM510','MM505','MM515','MM510_1'};

%channels on which detections will be performed (just an example)
channelsPerPatient = {[78  1  23  ],... % p488
    [16  10  39 31 61 ],... % p487
    [1  8  36  29  44 ],... % p489
    [8 46 ],... % p486
    [46 38 24],... % p499
    [52 ],... % p485
    [],... % 490 - RANH - excluded from analysis
    [1],...% 496
    [32 40 48],... % 498
    [9 47],... % p510
    [2],... % p505
    [9 39 65],... %p515
    [9 47],...% p510
    [17 9 ],...% 520 
    [9 16 1 40 48]};  % 497


%for bipolar ripple detection - in every row the first index is the channel in which ripple
%detection is required and the second is the reference channel
biPolarCouplesPerPatient = {[],[],[],[],[],[],[],[],[32 35; 33 35; 39 42; 40 42; 47 50; 48 50; ],[],[2 4],[9 12; 39 42; 40 42; 65 67; 66 67],[],[],[]};

%building the run data, not that for all file names of detections the
%methods assume the name is <provided name according to runData><channel
%num>, e.g. if runData(iPatient).SWStaresinaFileName='c:\slow_wave', then
%the slow waves file for channel 1 is 'c:\slow_wave1.mat'.
runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients 
    runData(iPatient).patientName = patients{iPatient};
    %The folder where the raw data is stored - you will need to change it
    runData(iPatient).DataFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO'];
       
    %The folder+filename from which spindles are going to be loaded
    %The folder+filename from which spindles are going to be loaded (should be the same
    %as the line above if the detections are first run and saved into SpindlesStaresinaFileNames
    spindleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleResults'];
    runData(iPatient).SpindlesFileNames = fullfile(spindleFolder,'spindleTimes');
    if isempty(dir(spindleFolder))
        mkdir(spindleFolder)
    end
    
    %The folder+filename into which the bipolar ripples detections results is going
    %to be stored
    rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\rippleBipolarResults'];
    runData(iPatient).RipplesBipolarFileNames = fullfile(rippleFolder,'rippleTimes');
    if isempty(dir(rippleFolder))
        mkdir(rippleFolder)
    end
    
    %The folder+filename from which ripples are going to be loaded (should be the same
    %as the line above if the bipolar detections are first run and saved into
    %RipplesBipolarFileNames)
    rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\rippleResults'];
    runData(iPatient).RipplesFileNames = fullfile(rippleFolder,'rippleTimes');
    if isempty(dir(rippleFolder))
        mkdir(rippleFolder)
    end
    
    %list of couples for bipolar ripple detection - where each row has the channel
    %in which the detection is performed in the first index, and the
    %reference channel in the second
    runData(iPatient).biPolarCouples = biPolarCouplesPerPatient{iPatient};
    %The folder+filename into which the spikes results is going to be stored or is
    %already stored if the spikes detection was already run (the folder should
    %exist)
    runData(iPatient).SpikesFileNames = fullfile(runData(iPatient).DataFolder, ...
        sprintf('MacroInterictalSpikeTimesFor_%s_%s_',patients{iPatient},expNames{iPatient}));
   
    %name of the EXP data for the patient
    runData(iPatient).ExpDataFileName = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
    % Extract stimulation times
    mm = matfile(runData(iPatient).ExpDataFileName); EXP_DATA = mm.EXP_DATA;
    runData(iPatient).stimulation_times = EXP_DATA.stimTiming.validatedTTL_NLX;
    clear mm EXP_DATA
    
    %name of the sleep scoring mat file for the patient
    runData(iPatient).sleepScoringFileName = [runData(iPatient).DataFolder,'\',sleepScoreFileName{iPatient},'.mat'];
    %channels that the ripples analyses will be performed on
    runData(iPatient).channelsToRunOn = channelsPerPatient{iPatient};   
    %macromontage file name per patient
    runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
    %The micro montage filename is required for micro channels analysis
    %which is not performed in this example
    %runData(iPatient).microMontageFileName = ['C:\Users\user\google drive\MayaProject\montages\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
    
    %name of file with single units info
    runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
end

%% an example for saving ripples using the wrapper AnalyzeCoupling.saveDetectionResults
ac = AnalyzeCoupling;

%setting which detections to run - in this example, bipolar ripples will be
%detected for channelsToRunOn per patient
whatToRun.runSpikes = false;
whatToRun.runRipples = true;
whatToRun.runRipplesBiPolar = false;
whatToRun.runSpindles = false;
whatToRun.runSpindlesStaresina = false;
whatToRun.runSWStaresina = false;
whatToRun.runSWMaingret = false;

%saving detections (in this example, bipolar ripples detection)
ac.saveDetectionResults(runData, whatToRun);

%% an example for detecting ripples directly using RippleDetector (it's the same thing the wrapper does inside)
rd = RippleDetector;

%an example of using the ripple detection directly and not with the wrapper
%(on the first channel of the first patient for this example)
currChan = runData(iPatient).channelsToRunOn(1);

%loading - sleep scoring, IIS, data
sleepScoring = load(runData(1).sleepScoringFileName);
sleepScoring = sleepScoring.sleep_score_vec;
peakTimes = load([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat']);
peakTimes = peakTimes.peakTimes;
currData = load([runData(iPatient).DataFolder,'\CSC',num2str(currChan),'.mat']);
currData = currData.data;
%detecting the ripples
[ripplesTimes, ripplesStartEnd] = rd.detectRipple(currData, sleepScoring, peakTimes);
%plotting the single ripples and saving the figures
rd.plotRipples(currData,ripplesTimes,'D:\singleRipples');

%% ripple related analyses

% Plot ripple summary figures per patient
for ii = 1:length(runData)
    resultsFile = fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\',...
        sprintf('resultsRipplesMACRO_%s.mat',runData(ii).patientName));
    resultsData = rd.runRippleData(runData(ii),resultsFile);
    figureFolder = 'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACROrippleFigures';
    rd.plotResultsRipplesData(resultsData,figureFolder);   
end

% Alternatively - run for all patients 
%runs the "ripple overview" analysis for the channels specified in
%channelsToRunOn
resultsData = rd.runRippleData(runData,...
    'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\resultsRipplesMACRO.mat');
%plots the results
rd.plotResultsRipplesData(resultsData,'E:\Data_p\ClosedLoopDataset\rippleDetResults\\ripplesFigures');

%runs the ripples-spikes correlation analysis for the channels specified in
%channelsToRunOn
resultsRS = rd.runRipSpikes(runData,'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\resultsRipplesRS.mat');
%plots the results
rd.plotResultsSpikes(resultsRS,'E:\Data_p\ClosedLoopDataset\rippleDetResults\ripplesSpikesFigures');

