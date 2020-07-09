%this struct can be longer than 1, where each element represents a patient
runDataForMicroSaving(1).microChannelsFolderToLoad = 'D:\microFolder\p488\fullData';
runDataForMicroSaving(1).microChannelsFolderToSave = 'D:\microFolder\p488\newData';
% runDataForMicroSaving(1).microChannelsToRunOn = [14]; Optional - if not
% provided will run on all the channels (files in the form: CSC#_#) in the
% folder microChannelsFolderToLoad.

%Note the code assumes the format of the filename of the parts is
%CSC<chan_number>_<part number in 3 digits exactly>. if it's different the
%properties microChannels.fileNamePrefix (currently 'CSC') and
%microChannels.fileNumberLength (currently 3) should be changed accordingly

mc = microChannels;
%These are the defaults, change them if it's not true
%mc.samplingRate = 1000; required sampling rate
%mc.oldSamplingRate = 40000;
mc.mergeAndDownsample(runDataForMicroSaving);

%% preparing the input struct to the methods
%building runData - the struct required by the methods for running the
%analyses

patients = {'p487'};
expNames = {'EXP3'};
sleepScoreFileName = {'sleepScore_manualValidated_p487_3_ROF2'};
%channels (micro) for which you want to detect ripples (this is just an
%example)
channelsPerPatient = {[1:16 35:48]};

%this is optional for the micro analysis - if you provide a cell of area
%names it will only run the analysis on these areas. If left empty or the
%field doesn't exist it will run on all areas (which will require having
%first saved ripples files for all the channels)
areasPerPatient = {{'REC', 'RAH', 'LPHG', 'LAH', 'LEC'}};
areasPerPatient = {{'RAH'}};

%this field is not relevant to micro analysis
% channelCouplesPerPatient = {[78 55; 1 9]};

% macroMontageNames = {'MM488'};

nPatients = length(patients);
for iPatient = 1:nPatients
    runData(iPatient).patientName = patients{iPatient};
    %The folder where the raw data is stored - you will need to change it
    %to the folder name of micro raw data
    runData(iPatient).DataFolder = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\MICRO'];
    
    %The folder where the ripples detections results is going to be stored
    %(it should exist, the code doesn't create it)
    runData(iPatient).RipplesFileNames = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\MICRO\rippleResults\rippleTimes'];
    
    %The folder from which the ripples will be read in the
    %runRipSpikesMicro method
    runData(iPatient).microRipplesFileNames = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\MICRO\rippleResults\rippleTimes'];

    %folder+filename where spikes are stored - i don't know if you want to remove
    %spikes from micro channels (if it's relevant). There are three
    %options: 1. enter here a folder that already stores spikes detections
    %results (mat files with the parameter peakTimes, the name for each file should be runData(iPatient).SpikesFileNames<#channel number>), 
    %2. enter here a folder into which you want to save newly generated spike detection results (which will be detected along with the ripples),
    %3. leave empty if it's not relevant for the analysis
    runData(iPatient).SpikesFileNames = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spikesResultsShdema\SpikesResults'];
    runData(iPatient).SpikesFileNames = [];
    %file name of ExpData - required
    runData(iPatient).ExpDataFileName = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
    %file name of sleep scoring - required (if left empty the entire night
    %will be used)
    runData(iPatient).sleepScoringFileName = ['D:\data_p\SleepScore_v1\',sleepScoreFileName{iPatient},'.mat'];
    %micro montage file name - required
    runData(iPatient).microMontageFileName = ['C:\Users\user\google drive\MayaProject\montages\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
    %spike data file name - required
    runData(iPatient).spikeData = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
    
    runData(iPatient).channelsToRunOn = channelsPerPatient{iPatient};
    runData(iPatient).areasToRunOn = areasPerPatient{iPatient};
    
    %The next fields are not used in ripples-spikes analyses
%     runData(iPatient).SWStaresinaFileName = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\SWStaresinaResults\SWStaresina'];
%     runData(iPatient).SWMaingretFileName = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\SWMaingretResults\SWMaingret'];
%     runData(iPatient).SpindlesFileNames = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleResults\spindleTimes'];
%     runData(iPatient).channelCouples = channelCouplesPerPatient{iPatient};
%     runData(iPatient).macroMontageFileName = ['C:\Users\user\google drive\MayaProject\badChans\newMontages\',macroMontageNames{iPatient},'.mat'];
%     
end

%% running

%first step - detecting ripples
whatToRun.runSpikes = false; %OR: whatToRun.runSpikes = true if you want to detect spikes and save the results into the folder written into SpikesFileNames field
whatToRun.runRipples = true;
whatToRun.runSpindles = false;
whatToRun.runSWStaresina = false;
whatToRun.runSWMaingret = false;

ac = AnalyzeCoupling;
ac.dataFilePrefix = 'CSC'; %Enter here the file name for the raw micro data, the default is CSC which is the name of Macro raw data, if it's the same for micro then delete this line
%saves ripples detection (and spikes if specified as required)
ac.saveDetectionResults(runData, whatToRun);

%second step - running analysis, if areasPerPatient is empty it will run on
%all areas
rd = RippleDetector;
%second parameter is a filename for saving results, can be left empty (and
%it will not save it)
runData(1).areasToRunOn = {'REC'};
results = rd.runRipSpikesMicro(runData(1), 'D:\results\resultsMicro487.mat');

%third step - producing figures, the second parameter is the folder for
%saving the figures, can be left empty and the figures will be presented on
%screen (note this folder should include subfolders with the relevant
%patient names into which it tries to save the figures, I should probably
%add a creation of these folders if they don't exist to the code but it's not there now...)
rd.plotResultsSpikes(results);