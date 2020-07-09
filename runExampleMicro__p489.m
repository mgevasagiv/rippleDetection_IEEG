%this struct can be longer than 1, where each element represents a patient
runDataForMicroSaving(1).microChannelsFolderToLoad = 'E:\Data_p\p489\EXP3\spikeSorting\MICRO_FILES';
runDataForMicroSaving(1).microChannelsFolderToSave = 'E:\Data_p\p489\EXP3\Denoised_Downsampled_InMicroVolt\MICRO';
runDataForMicroSaving(1).microChannelsToRunOn = [17:24]; % Optional - if not
% provided will run on all the channels (files in the form: CSC#_#) in the
% folder microChannelsFolderToLoad.
runDataForMicroSaving(1).microChannelsToRunOn = [1:16,49:56,25:32,57:64];% LEC, LMH, LPHG, RPHG

%Note the code assumes the format of the filename of the parts is
%CSC<chan_number>_<part number in 3 digits exactly>. if it's different the
%properties microChannels.fileNamePrefix (currently 'CSC') and
%microChannels.fileNumberLength (currently 3) should be changed accordingly

mc = microChannels;
IIS_det = SpikeWaveDetector;

%These are the defaults, change them if it's not true
mc.samplingRate = 1000; % required sampling rate
mc.oldSamplingRate = 40000;
mc.mergeAndDownsample(runDataForMicroSaving);

%% preparing the input struct to the methods
%building runData - the struct required by the methods for running the
%analyses

patients = {'p489'};
expNames = {'EXP3'};
sleepScoreFileName = {'sleepScore_manualValidated_p489_3_RPHG3'};
MACRO_ch = {[1]};
%channels (micro) for which you want to detect ripples (this is just an
%example)
MICROchannelsPerPatient = {{[1:8]},{[9:16]},{[17:24]},{[25:32]},{[49:56]},{[57:64]}};

%this is optional for the micro analysis - if you provide a cell of area
%names it will only run the analysis on these areas. If left empty or the
%field doesn't exist it will run on all areas (which will require having
%first saved ripples files for all the channels)
areasPerPatient = {{'RAH'}};

%this field is not relevant to micro analysis
channelCouplesPerPatient = {[1 19; 8 19; 36 43]};

data_p_path = 'E:\Data_p\';

nPatients = length(patients);

for iPatient = 1:nPatients
    
    nChan = length(MICROchannelsPerPatient{iPatient});
    
    for iChan = [3]
        
        runData(iPatient).patientName = patients{iPatient};
        %The folder where the raw data is stored - you will need to change it
        %to the folder name of micro raw data
        runData(iPatient).DataFolder = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO'];
        runData(iPatient).microChannelsFolderToLoad = runData(iPatient).DataFolder;
        
        %The folder where the ripples detections results is going to be stored
        %(it should exist, the code doesn't create it)
        runData(iPatient).RipplesFileNames = [data_p_path ,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO\rippleResults\rippleTimes'];
        runData(iPatient).microRipplesFileNames = runData(iPatient).RipplesFileNames;
        a = dir(runData(iPatient).RipplesFileNames);
        if isempty(a)
            mkdir(runData(iPatient).RipplesFileNames);
        end
        
        micro_ch = cell2mat(MICROchannelsPerPatient{iPatient, iChan});

        %folder+filename where spikes are stored - i don't know if you want to remove
        %spikes from micro channels (if it's relevant). There are three
        %options: 1. enter here a folder that already stores spikes detections
        %results (mat files with the parameter peakTimes, the name for each file should be runData(iPatient).SpikesFileNames<#channel number>),
        %2. enter here a folder into which you want to save newly generated spike detection results (which will be detected along with the ripples),
        %3. leave empty if it's not relevant for the analysis
        
        %micro montage file name - required
        runData(iPatient).microMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
        mfile = matfile(runData(iPatient).microMontageFileName);
        MicroMontage = mfile.Montage;

        %miacro montage file name - required
        runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
        mfile = matfile(runData(iPatient).macroMontageFileName);
        MacroMontage = mfile.MacroMontage;
        MACRO_ch(iPatient,iChan) = getDeepMacros(MacroMontage, MicroMontage(micro_ch(1)).Area);
        
        runData(iPatient).SpikesFileNames = [];
        if ~isempty(MACRO_ch(iPatient,iChan))
            spikeFilename = fullfile([data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\'],...
                sprintf('InterictalSpikeTimesFor_%s_%s_C%03d.mat',patients{iPatient},expNames{iPatient},MACRO_ch(iPatient,iChan)));
            a = dir(spikeFilename);
            if isempty(a)
                filename = fullfile([data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\'], ...
                    sprintf('CSC%d.mat',MACRO_ch(iPatient,iChan)) );
                mlink = matfile(filename);
                [peakTimes, passedConditions]= detectTimes(IIS_det, mlink.data,true);
                save(spikeFilename,'peakTimes','passedConditions','IIS_det')
            end
            
            runData(iPatient).SpikesFileNames = fullfile(runData(iPatient).DataFolder, ...
                sprintf('MacroInterictalSpikeTimesFor_%s_%s_C',patients{iPatient},expNames{iPatient}));
            
            
            for ii = micro_ch
                copyfile(spikeFilename,fullfile(runData(iPatient).DataFolder, ...
                    sprintf('MacroInterictalSpikeTimesFor_%s_%s_C%d.mat',patients{iPatient},expNames{iPatient},ii)));
            end
        end
        
        %file name of ExpData - required
        runData(iPatient).ExpDataFileName = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
        %file name of sleep scoring - required (if left empty the entire night
        %will be used)
        runData(iPatient).sleepScoringFileName = [data_p_path,'SleepScore_v1\',sleepScoreFileName{iPatient},'.mat'];
        %micro montage file name - required
        runData(iPatient).microMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
        %spike data file name - required
        runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
        
        runData(iPatient).channelsToRunOn = cell2mat(MICROchannelsPerPatient{iPatient,iChan});
        runData(iPatient).areasToRunOn = areasPerPatient{iPatient};
        
        %The next fields are not used in ripples-spikes analyses
        runData(iPatient).SWStaresinaFileName = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\SWStaresinaResults\SWStaresina'];
        runData(iPatient).SWMaingretFileName = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\SWMaingretResults\SWMaingret'];
        runData(iPatient).SpindlesFileNames = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleResults\spindleTimes'];
        runData(iPatient).channelCouples = channelCouplesPerPatient{iPatient};
        
        
        
        %first step - detecting ripples
        whatToRun.runSpikes = false; %OR: whatToRun.runSpikes = true if you want to detect spikes and save the results into the folder written into SpikesFileNames field
        whatToRun.runRipples = true;
        whatToRun.runSpindles = false;
        whatToRun.runSWStaresina = false;
        whatToRun.runSWMaingret = false;
        whatToRun.runSpindlesStaresina = false;
        whatToRun.runRipplesBiPolar = false;
        
        ac = AnalyzeCoupling;
        ac.dataFilePrefix = 'CSC'; %Enter here the file name for the raw micro data, the default is CSC which is the name of Macro raw data, if it's the same for micro then delete this line
        %saves ripples detection (and spikes if specified as required)
        ac.saveDetectionResults(runData, whatToRun);
        
    end
end



%% running

%first step - detecting ripples
whatToRun.runSpikes = false; %OR: whatToRun.runSpikes = true if you want to detect spikes and save the results into the folder written into SpikesFileNames field
whatToRun.runRipples = true;
whatToRun.runSpindles = false;
whatToRun.runSWStaresina = false;
whatToRun.runSWMaingret = false;
whatToRun.runSpindlesStaresina = false;
whatToRun.runRipplesBiPolar = false;

ac = AnalyzeCoupling;
ac.dataFilePrefix = 'CSC'; %Enter here the file name for the raw micro data, the default is CSC which is the name of Macro raw data, if it's the same for micro then delete this line
%saves ripples detection (and spikes if specified as required)
ac.saveDetectionResults(runData, whatToRun);

%second step - running analysis, if areasPerPatient is empty it will run on
%all areas
rd = RippleDetector;
%second parameter is a filename for saving results, can be left empty (and
%it will not save it)
results = rd.runRipSpikesMicro(runData, 'E:\Data_p\ClosedLoopDataset\rippleDetResults\results.mat');

%third step - producing figures, the second parameter is the folder for
%saving the figures, can be left empty and the figures will be presented on
%screen (note this folder should include subfolders with the relevant
%patient names into which it tries to save the figures, I should probably
%add a creation of these folders if they don't exist to the code but it's not there now...)
rd.plotResultsSpikes(results,'E:\Data_p\ClosedLoopDataset\rippleDetResults\figures');


% This plots MACRO-style ripples - one per panel -
for iChannel = 1:length(runData(iPatient).channelsToRunOn)
    disp(['channel ',num2str(runData(iPatient).channelsToRunOn(iChannel))]);
    currChan = runData(iPatient).channelsToRunOn(iChannel);
    try
        currData = load([runData(iPatient).DataFolder,'\',rd.dataFilePrefix ,num2str(currChan),'.mat']);
        currData = currData.data;
    catch
        disp([runData(iPatient).DataFolder,'\',obj.dataFilePrefix,num2str(currChan),'.mat doesnt exist']);
        continue;
    end
    
    try
        ripplesTimes = load([runData(iPatient).RipplesFileNames ,num2str(currChan),'.mat']); ripplesTimes = ripplesTimes.ripplesTimes;
    catch
        disp([runData(iPatient).RipplesFileNames ,num2str(currChan),'.mat doesnt exist']);
        continue;
    end
    rd.plotRipples(currData, ripplesTimes, 'E:\Data_p\ClosedLoopDataset\rippleDetResults\figures');
end

AreaStr = 'LMH';
folderToSave = fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\figures',runData(iPatient).patientName,AreaStr);
mkdir(folderToSave);
rd.plotRipplesMicro(runData(iPatient), AreaStr, folderToSave)

rd.plotResultsSpikes(results, folderToSave)
