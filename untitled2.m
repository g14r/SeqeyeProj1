% cnt = zeros(size(ANA.TN));
% for tn = 1:length(ANA.TN)
%     if max(ANA.AllPressIdx(tn , :)) > length(ANA.xEye{tn})
%         %if sum(isnan(ANA.EyePressTimePos)) == 14
%         cnt(tn) = 1;
%     end
% end
% sum(cnt)
%
%
% cnt = zeros(size(ANA.TN));
% for tn  =1:length(ANA.TN)
%     if length(ANA.xEyePosAng{tn,1}) ~= length(ANA.PressTimeSeries{tn})
%         cnt(tn) = 1;
%     end
% end
% sum(cnt)

for tn = 1:length(ANA.TN)
    % detrend the IPIs
    ANA.IPI(tn , :) = detrend(ANA.IPI(tn , :),'linear',7) + mean(ANA.IPI(tn , :));
    
    ANA.xEye{tn}(ANA.Pupil{tn}<= 5500) = NaN;
    ANA.yEye{tn}(ANA.Pupil{tn}<= 5500) = NaN;
    
    % Calculate the eye start position during state = 3
    idx = (ANA.state{tn} == 3);
    % the first 1.5 seconds is always start fixation. We're only considerin 300ms before the seqence is presented
    ANA.EyeStartPos(tn,1)  = nanmedian(ANA.xEye{tn}(idx,:));
    ANA.EyeStartPos(tn,2)  = nanmedian(ANA.yEye{tn}(idx,:));
    
    
    % Calculate the eye end position during state = 7
    idx = (ANA.state{tn} == 7);
    % the first 1.5 seconds is always start fixation. We're only considerin 300ms before the seqence is presented
    ANA.EyeEndPos(tn,1)  = nanmedian(ANA.xEye{tn}(idx,:));
    ANA.EyeEndPos(tn,2)  = nanmedian(ANA.yEye{tn}(idx,:));
    
    
    % convert eye movements to cms knowing that the distance
    % between the two pluses is 29.9 cm (2.3 cm between digits)
    
    ANA.pix2cm(tn,1) = abs(ANA.EyeStartPos(tn,1) - ANA.EyeEndPos(tn,1))/((sum(~isnan(ANA.AllPress(tn,:)))-1) * 2.3); % pixels per 1cm calculates for every trial
    
    ANA.xEyePosCM{tn , 1}(1,:) = [(ANA.xEye{tn} - repmat(ANA.EyeStartPos(tn,1),length(ANA.xEye{tn}) , 1))...
        ./repmat(ANA.pix2cm(tn,1),length(ANA.xEye{tn}) , 1)];
    ANA.yEyePosCM{tn , 1}(1,:) = [-(ANA.yEye{tn} - repmat(ANA.EyeStartPos(tn,2),length(ANA.yEye{tn}) , 1))...
        ./repmat(ANA.pix2cm(tn,1),length(ANA.yEye{tn}) , 1)];
    
    ANA.xEyePosAng{tn ,1}(1,:) = rad2deg(atan((ANA.xEyePosCM{tn , 1}-15) /46));
    
    
    
    % whare eyes are at the time of press
    for fing = 1:size(ANA.AllPressIdx , 2)
        if ~isnan(ANA.AllPressIdx(tn,fing)) & ANA.AllPressIdx(tn,fing)<length(ANA.xEyePosCM{tn})
            ANA.EyePressTimePos(tn,fing) = nanmedian(ANA.xEyePosCM{tn}(1, ANA.AllPressIdx(tn,fing)-10 : ANA.AllPressIdx(tn,fing)));
        else
            ANA.EyePressTimePos(tn,fing) = NaN;
        end
    end
    
    window  = 50;
    for w = window + 1 : length(ANA.xEyePosAng{tn , 1}) - window
        id = [w - window  w + window];
        if id (2) < length(ANA.xEyePosAng{tn,1})
            ANA.xEyeAngVelocity{tn ,1}(w-window , 1) = abs(ANA.xEyePosAng{tn}(id(1)) - ANA.xEyePosAng{tn}(id(2)))/300;
        end
    end
    
    
    % find eye angular velocity over the 200 ms (100 samples) vicinity of presses
    ANA.xEyePressAngVelocity(tn , :) = zeros(1,14);
    ANA.xEyePressAngVelocity(tn , isnan(ANA.AllPressIdx(tn , :))) = NaN;
    ANA.xEyeVelEstCnkPlcmnt(tn,:) = NaN * ones(1,14);
    for p = 1:sum(~isnan(ANA.AllPressIdx(tn , :)))
        id = [ANA.AllPressIdx(tn, p) - 50  ANA.AllPressIdx(tn, p) + 50 ];
        if id (2) < length(ANA.xEyePosAng{tn,1})
            ANA.xEyePressAngVelocity(tn , p) = abs(ANA.xEyePosAng{tn}(id(1)) - ANA.xEyePosAng{tn}(id(2)))/300;
        end
    end
    
    
    
    % create a smooth press time series
    presses = [0 : length(find(ANA.AllPress(tn , :)))+1];
    chunks  = [1 , ANA.ChnkPlcmnt(tn , :) , 1];
    ANA.PressTimeSeries{tn,1} = zeros(length(ANA.xEyePosCM{tn}) , 1);
    IND = [0 , ANA.AllPressIdx(tn , 1:sum(~isnan(ANA.AllPressIdx(tn,:)))) , length(ANA.PressTimeSeries{tn,1})];
    for ind = 1:length(IND) - 1
        ANA.PressTimeSeries{tn,1}(IND(ind)+1 : IND(ind + 1)) = linspace(presses(ind) , presses(ind+1) , IND(ind + 1) - IND(ind));
    end
    
    ANA.ChnkPlcEyeVelCorr(tn , 1) = NaN;
    % Create a smooth imposed-chunk timeseries
    ANA.IPIChnkPlcmnt(tn , :) = [NaN NaN NaN NaN];
    ANA.ChunkTimeSeries{tn,1} = zeros(size(ANA.xEyePosCM{tn}));
    if ismember(ANA.seqNumb(tn), [1:6])   % Intermixed or CLAT
        
        ANA.IPIwithin(tn,1)  = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 2));
        ANA.IPIbetween(tn,1) = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 1));
        
        ANA.IPIChnkPlcmnt(tn , 1)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 0)); % basically between
        ANA.IPIChnkPlcmnt(tn , 2)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 1));
        ANA.IPIChnkPlcmnt(tn , 3)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 2));
        ANA.IPIChnkPlcmnt(tn , 4)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 3));
        
        
        for ind = 1:length(IND) - 1
            ANA.ChunkTimeSeries{tn,1}(IND(ind)+1 : IND(ind + 1)) = linspace(chunks(ind) , chunks(ind+1) , IND(ind + 1) - IND(ind));
        end
      
        % correlate the chunk placement time series with the eye velocity
        
        id = ~isnan(ANA.xEyeAngVelocity{tn});
        A = ANA.ChunkTimeSeries{tn,1}(window+1 : end-window);
        A = A(id);
        B = ANA.xEyeAngVelocity{tn,1}(id);
        pressid = ANA.AllPressIdx(tn,1)- window : ANA.AllPressIdx(tn,sum(~isnan(ANA.AllPressIdx(tn  , :)))) - window;
        if pressid(end)<length(A)
            temp = corrcoef(A(pressid) , B(pressid));
            ANA.ChnkPlcEyeVelCorr(tn  ,1) = temp(2);
        end
        
        
    else
        ANA.IPIwithin(tn,1)  = NaN;
        ANA.IPIbetween(tn,1) = NaN;
    end
    
    % estimate chunk bounries as the ones that are 75% of a std above the rest of the IPIs
    
    thresh = .2 * std(ANA.IPI(tn , :));
    [dum , estChnkBndry] = findpeaks(ANA.IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
    
    if ~isempty(estChnkBndry)
        goodpeak = ones(1, length(estChnkBndry));
        for cb = 1:length(estChnkBndry)
            if ANA.IPI(tn , estChnkBndry(cb)) < mean(ANA.IPI(tn  ,:)) + thresh
                goodpeak(cb) = 0;
            end
        end
        if sum(goodpeak)
            estChnkBndry = estChnkBndry(logical(goodpeak));
        else
            estChnkBndry = [];
        end
    end
    
    
    ANA.estChnkBndry(tn , :) = zeros(1, 14);  % first presses of chunks will be 1
    ANA.estChnkPlcmnt(tn , :) = zeros(1, 14);  % chunk placements of digits
    ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
    ANA.estChnkBndry(tn , 1) = 1;
    dum = find(~ANA.estChnkBndry(tn,:));
    
    ANA.ChnkTimeSeriesCorr(tn,1) = NaN;
    ANA.ChnkPlcmntDist(tn,:)     = [NaN NaN NaN NaN NaN NaN];
    ANA.ChnkPlcmntCorr(tn,:)     = [NaN NaN NaN NaN NaN NaN];
    ANA.estChnkPlcEyeVelCorr(tn , 1) = NaN;
    if sum(ANA.estChnkBndry(tn,:))>0
        ANA.estChnkPlcmnt(tn,logical(ANA.estChnkBndry(tn,:))) = 1;
        ANA.estChnkPlcmnt(tn,1) = 1;
        dum = find(~ANA.estChnkPlcmnt(tn,:));
        for h = 1:length(dum)
            ANA.estChnkPlcmnt(tn,dum(h)) = ANA.estChnkPlcmnt(tn,dum(h)-1) + 1;
        end
        
        ANA.estChunkTimeSeries{tn,1} = zeros(size(ANA.xEyePosCM{tn}));
        ANA.estChunkTimeSeries{tn,1}(ANA.AllPressIdx(tn,1:sum(~isnan(ANA.AllPressIdx(tn,:))))) = ANA.estChnkPlcmnt(tn , 1:sum(~isnan(ANA.AllPressIdx(tn,:))));
        chunks  = [1 , ANA.estChnkPlcmnt(tn , :) , 1];
        % generate the timeseries for the estimated chunking patterns
        for ind = 1:length(IND) - 1
            ANA.estChunkTimeSeries{tn,1}(IND(ind)+1 : IND(ind + 1)) = linspace(chunks(ind) , chunks(ind+1) , IND(ind + 1) - IND(ind));
        end
        
        
        if ismember(ANA.seqNumb(tn), [1:6])
            for st = 1:6
                temp = corrcoef(ANA.estChunkTimeSeries{tn , 1} , ANA.ChunkTimeSeries{tn,1});
                ANA.ChnkTimeSeriesCorr(tn,1) = temp(2);
                
                temp = corrcoef(ANA.estChnkPlcmnt(tn , :) , ChnkPlcmnt(st,:));
                ANA.ChnkPlcmntCorr(tn,st) = temp(2);
                
                ANA.ChnkPlcmntDist(tn,st) = sqrt(sum((ANA.estChnkPlcmnt(tn , :) - ChnkPlcmnt(st,:)).^2));
            end
        end
    else
        ANA.estChunkTimeSeries{tn , 1} =  zeros(length(ANA.xEyePosCM{tn}) , 1);
    end
    % find mean eye angular velocity in terms of imposed chunk placement
    for cp  =  1:4
        ANA.xEyeVelCnkPlcmnt(tn,cp)  = nanmean(ANA.xEyePressAngVelocity(tn,ANA.ChnkPlcmnt(tn,:) == cp));
    end
    
    % find mean eye angular velocity in terms of estimated chunk placement
    for cp = 1:length(unique(ANA.estChnkPlcmnt(tn , :)))
        ANA.xEyeVelEstCnkPlcmnt(tn,cp)  = nanmean(ANA.xEyePressAngVelocity(tn,ANA.estChnkPlcmnt(tn,:) == cp));
    end
    % correlate the chunk placement time series with the eye velocity
    id = ~isnan(ANA.xEyeAngVelocity{tn});
    A1 = ANA.estChunkTimeSeries{tn,1}(window+1 : end-window);
    A1 = A1(id);
    B1 = ANA.xEyeAngVelocity{tn,1}(id);
    pressid = ANA.AllPressIdx(tn,1)- window : ANA.AllPressIdx(tn,sum(~isnan(ANA.AllPressIdx(tn  , :)))) - window;
    if pressid(end)<length(A1)
        temp = corrcoef(A1(pressid) , B1(pressid));
        ANA.estChnkPlcEyeVelCorr(tn  ,1) = temp(2);
    end
    
    
end