%the struct runData holds data about patients and where the different event
%types are stored

patients = {'p487'};
expNames = {'EXP3'};
sleepScoreFileName = {'sleepScore_manualValidated_p487_3_ROF2'};
%note macro montage is not necessary to run AnalyzeCoupling
macroMontageNames = {'MM487'};

%ripple detection is per area - the field areasToRunOn specifies which
%areas are run on per patient
areasPerPatient = {{'RAH','REC'}};
%ripples that area detected in a reference area at the same time will be
%discarded, for each patient this list should be at the same length as
%areasPerPatient to correspond to each area. If there is only one element
%in the list all areas will be referenced to it, if the list is empty no
%reference will be used
referenceAreasPerPatient = {{'LAH','LEC'}};

%this field will usually not be used, it's only relevant if for some reason
%there is a need to do detection on specific micro channels - i.e. by
%channel and not by area
channelsPerPatient = {[1 2 3]};
%this field can include micro channels that should be disregarded when
%finding a ripple per area (i.e. the merging process will not include them, 
%the saving method will also not save single ripples for them) and the
%plotting method will not plot them
noisyChannelsPerPatient = {[9]};

runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients 
    runData(iPatient).patientName = patients{iPatient};
    runData(iPatient).ExpDataFileName = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
    %sets for which areas the detection and analysis is performed
    runData(iPatient).areasToRunOn = areasPerPatient{iPatient};
    %reference areas per patient (for the detection phase)
    runData(iPatient).referenceAreas = referenceAreasPerPatient{iPatient};
    %The folder where the raw data is stored - you will need to change it
    runData(iPatient).MicroDataFolder = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\MICRO'];       
    %The folder+filename from which micro ripples are going to be loaded
    runData(iPatient).MicroRipplesFileNames = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\MICRO\rippleResults\rippleTimes'];
    %The folder+filename into which the spikes results is going to be stored or is
    %already stored if the spikes detection was already run (spikes for
    %micro channels are per area - so they should be saved as
    %<runData(iPatient).MicroSpikesFileNames><area name>.mat
    runData(iPatient).MicroSpikesFileNames = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\MICRO\SpikesResults\SpikesResults'];
    %name of the sleep scoring mat file for the patient
    runData(iPatient).sleepScoringFileName = ['D:\data_p\SleepScore_v1\',sleepScoreFileName{iPatient},'.mat'];    
    %this is usually not relevant - only if we want the detection to be for
    %specific channels
    runData(iPatient).channelsToRunOn = channelsPerPatient{iPatient};   
    %micro montage file name - required
    runData(iPatient).microMontageFileName = ['C:\Users\user\google drive\MayaProject\montages\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
    %name of file with single units info
    runData(iPatient).spikeData = ['D:\data_p\',patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
    runData(iPatient).noisyChannels = noisyChannelsPerPatient{iPatient};
end

%% an example for saving ripples using the wrapper RippleDetector.saveRipplesDetectionsMicro

rd = RippleDetector;

%saving detections for the specified areas
rd.saveRipplesDetectionsMicro(runData);

%% plotting micro ripples, second input parameter is the area name, third is a folder to save the figures to
rd.plotRipplesMicro(runData(1), 'RAH', 'LAH', 'D:\SingleRipplesFolder');

%% ripple related analyses - ripple spike correlation

rd = RippleDetector;
%if areasToRunOn is empty it will run on
%all areas, in this example we set areasToRunOn to 'RAH' so the analysis
%will be performed only on RAH
runData(1).areasToRunOn = {'RAH'};
%second parameter is a filename for saving results, can be left empty (and
%it will not save it)
results = rd.runRipSpikesMicro(runData(1), 'D:\results\resultsMicro487.mat');
results = rd.runRipSpikesMicro(runData(1));


%plotting the figures, the second parameter is the folder for
%saving the figures, can be left empty and the figures will be presented on
%screen
rd.plotResultsSpikes(results,'D:\results\figuresFolder');
