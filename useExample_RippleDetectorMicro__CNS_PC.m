% source code for interictal spike detection
IIS_det = SpikeWaveDetector;

% database to run on - 
patients = {'p486','p487','p488','p489','p498'};
expNames = {'EXP8','EXP3','EXP4','EXP3','EXP3'};
sleepScoreFileName = {    'sleepScore_manualValidated_p486_8_RA1',...    
    'sleepScore_manualValidated_p487_3_ROF2',...
    'sleepScore_manualValidated_p488_4_LPHG2',...
    'sleepScore_manualValidated_p489_3_RPHG3'...
    'sleepScore_manualValidated_p498_3_ROF1'};

%channels (micro) for which you want to detect ripples (this is just an
%example)
MICROchannelsPerPatient{1} = {{[1:8]},{[9:16]},{[17:24]},{[25:32]},{[33:40]},{[41:48]},{[49:56]},{[57:64]},...
    {[65:72]},{[73:80]},{[81:88]},{[89:96]}}; % p486
MICROchannelsPerPatient{2} = {{[1:8]},{[9:16]},{[17:24]},{[25:32]},{[33:40]},{[41:48]},{[49:56]},{[57:64]},...
    {[65:72]},{[73:80]}}; % p487
MICROchannelsPerPatient{3} = {{[1:8]},{[9:16]},{[17:24]},{[25:32]},{[33:40]},{[41:48]},{[49:56]},{[57:64]},...
    {[65:72]},{[73:80]},{[81:88]}};
MICROchannelsPerPatient{4} = {{[1:8]},{[9:16]},{[17:24]},{[25:32]},{[49:56]},{[57:64]}};% LEC, LMH, LPHG, RPHG
MICROchannelsPerPatient{5} = {{[1:8]},{[17:24]},{[25:32]},{[33:40]},{[49:56]},{[57:64]}};
%this is optional for the micro analysis - if you provide a cell of area
%names it will only run the analysis on these areas. If left empty or the
%field doesn't exist it will run on all areas (which will require having
%first saved ripples files for all the channels)
% areasPerPatient = {{'REC','RPHG','LPHG'},{'RAH'},{'RMH','RPHG','LAH','LPHG'}};

%this field is not relevant to micro analysis
channelCouplesPerPatient = {[1 19; 8 19; 36 43],[],[],[],[]};

data_p_path = 'E:\Data_p\';
PLOT_SINGLE_RIPPLES = 1;
nPatients = length(patients);

for iPatient = 1:nPatients
    %the struct runData holds data about patients and where the different event
    %types are stored

    nChan = length(MICROchannelsPerPatient{iPatient});
    runData(iPatient).patientName = patients{iPatient};
    %The folder where the raw data is stored - you will need to change it
    %to the folder name of micro raw data
    runData(iPatient).MicroDataFolder = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO'];
    runData(iPatient).MacroDataFolder = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO'];
    
    runData(iPatient).microChannelsFolderToLoad = runData(iPatient).MicroDataFolder;
    runData(iPatient).DataFolder = runData(iPatient).MicroDataFolder;
    
    runData(iPatient).MicroSpikesFileNames = fullfile(runData(iPatient).MicroDataFolder, ...
        sprintf('MacroInterictalSpikeTimesFor_%s_%s_',patients{iPatient},expNames{iPatient}));
    
 
    %file name of ExpData - required
    runData(iPatient).ExpDataFileName = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\',sprintf('%s_%s_dataset.mat', patients{iPatient},expNames{iPatient})];
    %file name of sleep scoring - required (if left empty the entire night
    %will be used)
    runData(iPatient).sleepScoringFileName = [runData(iPatient).MacroDataFolder,'\',sleepScoreFileName{iPatient},'.mat'];
    %micro montage file name - required
    runData(iPatient).microMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
    %spike data file name - required
    runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
 
    %The next fields are not used in ripples-spikes analyses
    % but are used for calculating the oscilations in the MACRO channels
    runData(iPatient).SWStaresinaFileName = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\SWStaresinaResults\SWStaresina'];
    runData(iPatient).SWMaingretFileName = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\SWMaingretResults\SWMaingret'];
    runData(iPatient).SpindlesFileNames = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleResults\spindleTimes'];
    runData(iPatient).channelCouples = channelCouplesPerPatient{iPatient};
    
    %The folder where the ripples detections results is going to be stored
    %(it should exist, the code doesn't create it)
    runData(iPatient).RipplesFolder = [data_p_path ,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO\rippleResults\'];
    runData(iPatient).RipplesFileNames = fullfile(runData(iPatient).RipplesFolder,'rippleTimes');
    runData(iPatient).MicroRipplesFileNames = runData(iPatient).RipplesFileNames;
    a = dir(runData(iPatient).RipplesFolder);
    if isempty(a); mkdir(runData(iPatient).RipplesFolder); end
    
    results = [];

    for iChan = 1:nChan
            
        micro_ch = cell2mat(MICROchannelsPerPatient{iPatient}{iChan});
        
        %folder+filename where spikes are stored - i don't know if you want to remove
        %spikes from micro channels (if it's relevant). There are three
        %options: 1. enter here a folder that already stores spikes detections
        %results (mat files with the parameter peakTimes, the name for each file should be runData(iPatient).SpikesFileNames<#channel number>),
        %2. enter here a folder into which you want to save newly generated spike detection results (which will be detected along with the ripples),
        %3. leave empty if it's not relevant for the analysis
        
        % micro montage file name - required
        runData(iPatient).microMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
        % macro montage file name - required
        runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
       
        mfile = matfile(runData(iPatient).microMontageFileName);
        MicroMontage = mfile.Montage;
        ch_area = MicroMontage(micro_ch(1)).Area

        MicroSpikesFileNames_ch = [runData(iPatient).MicroSpikesFileNames,ch_area,'.mat'];
        
        mfile = matfile(runData(iPatient).macroMontageFileName);
        MacroMontage = mfile.MacroMontage;
        MACRO_ch(iPatient,iChan) = getDeepMacros(MacroMontage, ch_area);
        
        mfile = matfile(runData(iPatient).ExpDataFileName);
        EXP_DATA = mfile.EXP_DATA;
        stimTimestamps_ms = EXP_DATA.stimTiming.validatedTTL_NLX;
        
        runData(iPatient).areasToRunOn{iChan} = ch_area;
        runData(iPatient).SpikesFileNames = runData(iPatient).MicroSpikesFileNames;
        
        if ~isempty(MACRO_ch(iPatient,iChan))
            MicroSpikesFileNames_ch = [runData(iPatient).MicroSpikesFileNames,ch_area,'.mat'];
            
            a = dir(MicroSpikesFileNames_ch);
            if isempty(a)
                filename = fullfile([data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\'], ...
                    sprintf('CSC%d.mat',MACRO_ch(iPatient,iChan)) );
                mlink = matfile(filename);
                [peakTimes, passedConditions]= detectTimes(IIS_det, mlink.data,true);
                %% Add stimulation times into peak time struct
                peakTimes = [peakTimes, stimTimestamps_ms];
                
                save(MicroSpikesFileNames_ch,'peakTimes','passedConditions','IIS_det')
                for ii_c = 1:length(micro_ch)
                    copyfile(MicroSpikesFileNames_ch, [runData(iPatient).MicroSpikesFileNames,num2str(micro_ch(ii_c)),'.mat']);
                end
            end
            
            % To analyze the realtionship with spindles and SWs, the detections
            % need to exist for the deep channels
            whatToRun.runSpikes = true;
            whatToRun.runRipples = false;
            whatToRun.runSpindles = true;
            whatToRun.runSWStaresina = true;
            whatToRun.runSWMaingret = false;
            whatToRun.runSpindlesStaresina = false;
            whatToRun.runRipplesBiPolar = false;
            
            ac = AnalyzeCoupling;
            ac.dataFilePrefix = 'CSC'; %Enter here the file name for the raw micro data, the default is CSC which is the name of Macro raw data, if it's the same for micro then delete this line
            % saves ripples detection (and spikes if specified as required)
            % ripples are detected obased on channels numbers for both MACRO and MICRO ripples (runData(iPatient).channelsToRunOn)
            runDataMACRO = runData(iPatient);
            runDataMACRO.DataFolder = runData(iPatient).MacroDataFolder;
            runDataMACRO.SpikesFileNames = fullfile(runDataMACRO.DataFolder,...
                sprintf('MacroInterictalSpikeTimesFor_%s_%s_',patients{iPatient},expNames{iPatient}));
            runDataMACRO.channelsToRunOn = MACRO_ch(iPatient,iChan);
            ac.saveDetectionResults(runDataMACRO, whatToRun);
            
        end
        
    end % run spike/sw/spindle detection for deep macro-channels
    
    
    
    
    
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
    % saves ripples detection (and spikes if specified as required)
    % ripples are detected obased on channels numbers for both MACRO and MICRO ripples (runData(iPatient).channelsToRunOn) 
    runData(iPatient).channelsToRunOn = MICROchannelsPerPatient{iPatient}{1}{1}(1):MICROchannelsPerPatient{iPatient}{end}{1}(end);
    ac.saveDetectionResults(runData(iPatient), whatToRun);
    
    
    rd = RippleDetector; 
    % saving detections for the specified areas (this is diffferent for
    % MICRO ripples )
    runByChannel = false;
    useExistingRipples = true;
    % runData(iPatient).noisyChannels =
    rd.saveRipplesDetectionsMicro(runData(iPatient), runByChannel, useExistingRipples);
    
    if PLOT_SINGLE_RIPPLES
        %% plotting micro ripples, second input parameter is the area name, third is a folder to save the figures to
        ch_area = 'LEC';
        mkdir(fullfile(runData(iPatient).RipplesFolder,'rippleFigures',ch_area));
        refArea = [];
        rd.plotRipplesMicro(runData(iPatient), ch_area, refArea,...
            fullfile(runData(iPatient).RipplesFolder,'rippleFigures',ch_area));
    end
    
    summary_plots_folder = fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\figures\');
    pt_figure_folder = fullfile(summary_plots_folder,sprintf('%s',runData(iPatient).patientName));
    if isempty(dir(pt_figure_folder)); mkdir(pt_figure_folder); end
   
    
    %% Ripple - unit spiking correlation
    fileNameResults = fullfile(runData(iPatient).RipplesFolder,sprintf('MICRO_ripplesSpikes_%s_%s.mat',...
        runData(iPatient).patientName,expNames{iPatient}));
    results = rd.runRipSpikesMicro(runData(iPatient), fileNameResults);       
    %plotting the figures, the second parameter is the folder for
    %saving the figures, can be left empty and the figures will be presented on
    %screen
    rd.plotResultsSpikes(results,summary_plots_folder);

    
    
    %% Ripple - summary of spectral content 
    
    fileNameResults = fullfile(runData(iPatient).RipplesFolder,sprintf('MICRO_ripplesSummary_%s.mat',runData(iPatient).patientName));
    results = [];
    runData(iPatient).channelsToRunOn = [];
    results = rd.runRippleDataMicro(runData(iPatient), fileNameResults);
    rd.plotResultsRipplesDataMicro(results,summary_plots_folder);
    
    
    
end
