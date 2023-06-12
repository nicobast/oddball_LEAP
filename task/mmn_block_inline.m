function [ logOut ]=mmn_block(pres,track,log,vars,varNames,curSample)

% fold variable names and values into a struct for ease of
% access/readability
design=cell2struct(vars,varNames,2);
design.BlockNo=curSample;

% DESIGN
trials=ECKList;
trials.Name='auditory mmn trial list';
trials.Presenter=pres;
trials.Tracker=track;
trials.Log=log;

% load mmn stimuli data
load([pres.Paths.mmn_stim, filesep, 'stim_table.mat']);

% load sounds
[sndStd, audioFs]       =...
    wavread([pres.Paths.mmn_stim, filesep, '1000Hz_50ms.wav']);
[sndDevFreq, audioFs]   =...
    wavread([pres.Paths.mmn_stim, filesep, '1500Hz_50ms.wav']);
[sndDevDur, audioFs]    =...
    wavread([pres.Paths.mmn_stim, filesep, '1000Hz_100ms.wav']);
[sndDevComb, audioFs]   =...
    wavread([pres.Paths.mmn_stim, filesep, '1500Hz_100ms.wav']);

% create buffers
% global sndPtr bufStd bufFreqDev bufDurDev bufCombDev mov

% load movie
mov = ECKLookupObject(pres.Movies, 'mmnpingu_silent.m4v');

% open sound device
sndPtr                  =   PsychPortAudio('Open',[],[],1,audioFs,1);

bufStd                  =   PsychPortAudio('CreateBuffer', sndPtr, sndStd');
bufFreqDev              =   PsychPortAudio('CreateBuffer', sndPtr, sndDevFreq');
bufDurDev               =   PsychPortAudio('CreateBuffer', sndPtr, sndDevDur');
bufCombDev              =   PsychPortAudio('CreateBuffer', sndPtr, sndDevComb');

stim                    =   stim_table(design.StartTrial:design.StartTrial +...
                            design.NumTrials - 1, :);
freq                    =   stim(:, 1);
dur                     =   stim(:, 2);
eeg_code                =   stim(:, 3);

trials.ImportVariables(varNames, vars, 'ms_block');
trials.NumSamples=1;
trials.StartSample=1;
trials.Order='SEQUENTIAL';
trials.DisableLogSaveOnTrial = true;

mov.Play(pres.WindowPtr);

profile on

curTrial = 1;
while curTrial <= design.NumTrials && ~pres.QuitRequested

    % sync EEG
    pres.EEGSync
    
    % random ISI
    ISI = 500 + (rand * 100);

    % determine stimulus type, and fill sound buffer accordingly
    if freq(curTrial) == 1000 && dur(curTrial) == 50
        % std
        PsychPortAudio('FillBuffer', sndPtr, bufStd);
    elseif freq(curTrial) == 1500 && dur(curTrial) == 50
        % freq dev
        PsychPortAudio('FillBuffer', sndPtr, bufFreqDev);
    elseif freq(curTrial) == 1000 && dur(curTrial) == 100
        % freq dev
        PsychPortAudio('FillBuffer', sndPtr, bufDurDev);
    elseif freq(curTrial) == 1500 && dur(curTrial) == 100
        % comb dev
        PsychPortAudio('FillBuffer', sndPtr, bufCombDev);
    end

    pres.EEGSendEvent(eeg_code(curTrial));

    % play sound
    sndOnsetTime = GetSecs;
    PsychPortAudio('Start', sndPtr, 1, [], 0);
    fprintf('Start audio: %.3f', GetSecs - sndOnsetTime);
    
    % loop for duration of sound
    while (GetSecs - sndOnsetTime) < (dur(curTrial) / 1000)
        pres.DrawImageFullscreen(mov.Frame);
        pres.RefreshDisplay;
    end    
    
    fprintf('Scheduled: %.3f, achieved: %.3f\n', dur(curTrial) / 1000, (GetSecs - sndOnsetTime) * 1000);
    
    PsychPortAudio('Stop', sndPtr);

    % loop for duration and ISI
    isiOnset = GetSecs;
    while GetSecs - isiOnset < ISI / 1000
        pres.DrawImageFullscreen(mov.Frame);
        pres.RefreshDisplay;
    end

    trialOffsetTime = GetSecs;
% 
%     logOut.Data=horzcat({...
%         num2str(trialOnsetTime),...
%         num2str(trialOffsetTime),...
%         num2str(sndOnsetTime),...
%         num2str(ISI),...
%         },vars);
% 
%     logOut.Headings=horzcat({...
%         'TrialOnsetTime',...
%         'TrialOffsetTime',...
%         'SoundOnsetTime',...
%         'ISI',...
%         }, varNames);

    % check for pause
    if pres.PauseRequested
        mov.Stop;
        fprintf('\nPAUSED - press "p" to continue...\n');
        while pres.PauseRequested && ~pres.QuitRequested
            pres.RefreshDisplay;
        end
        mov.Play(pres.WindowPtr);
    end
    
end

profile off

mov.Stop;

% Close the audio device:
PsychPortAudio('Close', sndPtr);

logOut=[];

end