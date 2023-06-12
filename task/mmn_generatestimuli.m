function [stim_table] = mmn_generatestimuli(numTrials, propFreqDev, propDurDev...
    , propCombDev, stdFreq, devFreq, stdDur, devDur)

    rng('shuffle');
    
    % eeg codes
    eegStd              =   201;
    eegDevFreq          =   202;
    eegDevDur           =   203;
    eegDevComb          =   204;
    
    % min safe gap to insert repeats
    minGap = 5;
    
    % calculate proportion of standards (freq, dev)
    propStd             =   1 - propFreqDev - propDurDev - propCombDev;

    % fill array with standard freqs and durations
    arr                 =   repmat([stdFreq, stdDur], [numTrials, 1]);

    % generate prob. table   
    probDevFreq         =   zeros(numTrials, 1);
    probDurDev          =   zeros(numTrials, 1);
    probCombDev         =   zeros(numTrials, 1);
    
    while abs((sum(probDevFreq) / numTrials) - propFreqDev) > .002
        probSeed(:,1)   =   rand(numTrials, 1); 
        probDevFreq     =   probSeed(:, 1) < propFreqDev;
    end
    
    while abs((sum(probDurDev) /numTrials) - propDurDev) > .002
        probSeed(:,2)   =   rand(numTrials, 1); 
        probDurDev      =   probSeed(:, 1) < propDurDev;
    end
    
    while abs((sum(probCombDev) / numTrials) - propCombDev) > .002
        probSeed(:,3)   =   rand(numTrials, 1); 
        probCombDev     =   probSeed(:, 1) < propCombDev;
    end
    
    probDevFreq         =   probSeed(:, 1) < propFreqDev;
    probDurDev          =   probSeed(:, 2) < propDurDev;
    probCombDev         =   probSeed(:, 3) < propCombDev;

    % replace stds with devs according to prob. table
    arr(probDevFreq | probCombDev, 1) =  devFreq;
    arr(probDurDev | probCombDev, 2) = devDur;
      
    % find repeats
    freqRep             =   findcontig(arr(:,1), devFreq);
    durRep              =   findcontig(arr(:,2), devDur);

    % look for repeats of two or more
    freqRepFound        =   freqRep(freqRep(:,3) >= 2, :);
    durRepFound         =   durRep(durRep(:,3) >= 2, :);

    % find clear chunks of stds as targets for replacement
    freqStdRep          =   findcontig(arr(:, 1), stdFreq);
    freqStdFound        =   freqStdRep(freqStdRep(:,3) >= 3, :);
    durStdRep           =   findcontig(arr(:, 2), stdDur);
    durStdFound         =   durStdRep(durStdRep(:,3) >= 3, :);

    %  loop until no repeats are found
    while ~isempty(freqRepFound) || ~isempty(durRepFound)

        % loop through dev repeats and move to safe chunks of std repeats, do this
        % separately for freq and dev
        for curRep = 1:size(freqRepFound, 1)

            % loop through each repeat (as may be more than two in a row)
            for curRepEntry = 1:freqRepFound(curRep, 3) - 1

                % store current repeat, replace with std
                curIdx          =   freqRepFound(curRep, 1) + curRepEntry;
                curRepItem      =   arr(curIdx, :);
                arr(curIdx, :)  =   [stdFreq, stdDur];

                % choose location for replacement
                repIdx          =   freqStdFound(randi(size(freqStdFound, 1)), 1) + 1;
                arr(repIdx, :)  =   curRepItem;

                % update std free chunks
                freqStdRep      =   findcontig(arr(:, 1), stdFreq);
                freqStdFound    =   freqStdRep(freqStdRep(:,3) >= 3, :);
                durStdRep       =   findcontig(arr(:, 2), stdDur);
                durStdFound     =   durStdRep(durStdRep(:,3) >= 3, :);

            end

        end

        for curRep = 1:size(durRepFound, 1)

            % loop through each repeat (as may be more than two in a row)
            for curRepEntry = 1:durRepFound(curRep, 3) - 1

                % store current repeat, replace with std
                curIdx          =   durRepFound(curRep, 1) + curRepEntry;
                curRepItem      =   arr(curIdx, :);
                arr(curIdx, :)  =   [stdFreq, stdDur];

                % choose location for replacement
                repIdx          =   freqStdFound(randi(size(freqStdFound, 1)), 1) + 1;
                arr(repIdx, :)  =   curRepItem;

                % update std free chunks
                freqStdRep      =   findcontig(arr(:, 1), stdFreq);
                freqStdFound    =   freqStdRep(freqStdRep(:,3) >= minGap, :);
                durStdRep       =   findcontig(arr(:, 2), stdDur);
                durStdFound     =   durStdRep(durStdRep(:,3) >= minGap, :);
            end

        end

        % find repeats
        freqRep             =   findcontig(arr(:,1), devFreq);
        durRep              =   findcontig(arr(:,2), devDur);

        % look for repeats of two or more
        freqRepFound        =   freqRep(freqRep(:,3) >= 2, :);
        durRepFound         =   durRep(durRep(:,3) >= 2, :);

    end
    
    % insert EEG codes
    arr(:,3) = eegStd;
    arr(arr(:, 1) == devFreq & arr(:,2) == stdDur, 3) = eegDevFreq;
    arr(arr(:, 1) == stdFreq & arr(:,2) == devDur, 3) = eegDevDur;
    arr(arr(:, 1) == devFreq & arr(:,2) == devDur, 3) = eegDevComb;
    
    stim_table = arr;
    paths = ECKPaths;
    save([paths.mmn_stim, filesep, 'stim_table.mat']);

end