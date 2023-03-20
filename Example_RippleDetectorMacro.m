%the struct runData holds data about patients and where the different event
%types are stored
data_p_path = 'E:\Data_p\';

patients = {'p1'};
expNames = {'EXP1'};
sleepScoreFileName = {'sleepScore_manualValidated_p001_1_LPHG2'};

%channels on which detections will be performed (just an example)
channelsPerPatient = {[78  1  23  ],... % p1
                      }; 

%for bipolar ripple detection - in every row the first index is the channel in which ripple
%detection is required and the second is the reference channel
biPolarCouplesPerPatient = {[78 81;1 4;23 26],... % p1
}; 
bipolarDet = 0; % determines if uni-polar or bi-polar detection is used

%building the run data, not that for all file names of detections the
%methods assume the name is <provided name according to runData><channel
%num>, e.g. if runData(iPatient).SWStaresinaFileName='c:\slow_wave', then
%the slow waves file for channel 1 is 'c:\slow_wave1.mat'.
runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients 
    runData(iPatient).patientName = patients{iPatient};
    %The folder where the raw data is stored - you will need to change it
    runData(iPatient).DataFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\MACRO'];
       
    %The folder+filename from which spindles are going to be loaded
    %The folder+filename from which spindles are going to be loaded (should be the same
    %as the line above if the detections are first run and saved into SpindlesStaresinaFileNames
    spindleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\MACRO\spindleResults'];
    runData(iPatient).SpindlesFileNames = fullfile(spindleFolder,'spindleTimes');
    if isempty(dir(spindleFolder))
        mkdir(spindleFolder)
    end
    
    %The folder+filename into which the bipolar ripples detections results is going
    %to be stored
    rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\MACRO\rippleBipolarResults'];
    runData(iPatient).RipplesBipolarFileNames = fullfile(rippleFolder,'rippleTimes');
    if isempty(dir(rippleFolder))
        mkdir(rippleFolder)
    end
    
    %The folder+filename from which ripples are going to be loaded (should be the same
    %as the line above if the bipolar detections are first run and saved into
    %RipplesBipolarFileNames)
    rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\MACRO\rippleBipolarResults'];
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
    
    %name of the sleep scoring mat file for the patient
    runData(iPatient).sleepScoringFileName = [runData(iPatient).DataFolder,'\',sleepScoreFileName{iPatient},'.mat'];
    %channels that the ripples analyses will be performed on

    % For ripple detection
    if ~bipolarDet
        runData(iPatient).channelsToRunOn = biPolarCouplesPerPatient{iPatient}(:,1);   
    else
    % For analysis
        runData(iPatient).channelsToRunOn = biPolarCouplesPerPatient{iPatient};   
    end
    
    %macromontage file name per patient
    runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
    
    %name of file with single units info
    runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
end


%% an example for saving ripples using the wrapper AnalyzeCoupling.saveDetectionResults
as = AnalyzeSleepOsc;

%setting which detections to run - 
whatToRun.runSpikes = false;
whatToRun.runRipples = true;
whatToRun.runRipplesBiPolar = false;
whatToRun.runSpindles = false;
whatToRun.runSpindlesStaresina = false;
whatToRun.runSWStaresina = false;
whatToRun.runSWMaingret = false;
whatToRun.HighFreqSpindles = false;

%saving detections (in this example, bipolar ripples detection)
parfor ii = 1:length(runData)
    as.saveDetectionResults(runData(ii), whatToRun);
end

%% an example for detecting ripples directly using RippleDetector (it's the same thing the wrapper does inside)
rd = RippleDetector_class;

%an example of using the ripple detection directly and not with the wrapper
%(on the first channel of the first patient for this example)
currChan = runData(iPatient).channelsToRunOn(1);

%loading - sleep scoring, IED detection times, data
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


