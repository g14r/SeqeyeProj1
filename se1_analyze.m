function Dout=se1_analyze (what , getdat , SubjCodes , Dall)

prefix = 'se1_';

%baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye1/analyze';
%baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/se1_data/analyze';
baseDir = '/Volumes/motorcontrol/data/SeqEye1/SEp/data/analyze';  % server
subj_name = {'SZ1','JN2' ,'SP1','AT1','DW1','KL1','JG1','GP1','SK1' ,'NM1', 'VR1', 'PB1'};
load([baseDir , '/CMB.mat'])
%load([baseDir , '/se1_all.mat'])
% D   = load('/Users/nedakordjazi/Documents/SeqEye/SequenceHierarchical/Analysis/sh3_avrgPattern.mat');
% MTW = 2 - D.MT(126:end)/max(D.MT(126:end));
% possibleDuo = D.Sequence(126:end , 1:2);
% create an emty structure with the same fields as Dall (empty b/c isError can never be 2)

switch what
    case 'all_subj'
        for i=1:length(subj_name)
            clear ANA
            if getdat
                ANA = se1_subj(subj_name{i} , 0);
            else
                ANA = getrow(Dall , Dall.SN == i);
            end
            BlPerDay = {[1 :14] [15 :29] [30 :44] [45 :59]};
            ANA.Day = zeros(size(ANA.BN));
            for d = 1:length(BlPerDay)
                ANA.Day(ismember(ANA.BN , BlPerDay{d})) = d;
            end
            GroupCode = str2num(subj_name{i}(3));
            ANA.SN(1:length(ANA.BN) , :) = i;
            ANA.Group(1:length(ANA.BN),:) = GroupCode;
            ChnkArrng = zeros(6,14);
            Chnkplcmnt = zeros(6,14);
            if GroupCode == 1
                ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group1 , :);
                for j = 1:size(ChnkArrang , 1)
                    temp = [];
                    temp1 = [];
                    for k = 1:length(find(ChnkArrang(j , :)))
                        temp = [temp k*ones(1,ChnkArrang(j,k))];
                        temp1 = [temp1 1:ChnkArrang(j,k)];
                    end
                    ChnkArrng(j , :) = temp;
                    ChnkPlcmnt(j,:)  = temp1;
                    
                    ICP(j , :) = diff(ChnkPlcmnt(j,:),1,2);
                    ICP(j , (ICP(j , :)<0)) = 0;
                    for k = 1:length(ICP(j , :))-1
                        if ICP(j,k) & ICP(j,k+1)
                            ICP(j,k+1) = ICP(j,k) + ICP(j,k+1);
                        end
                    end
                end
                IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
                IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
            elseif GroupCode == 2
                ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group2 , :);
                for j = 1:size(ChnkArrang , 1)
                    temp = [];
                    temp1 = [];
                    for k = 1:length(find(ChnkArrang(j , :)))
                        temp = [temp k*ones(1,ChnkArrang(j,k))];
                        temp1 = [temp1 1:ChnkArrang(j,k)];
                    end
                    ChnkArrng(j , :) = temp;
                    ChnkPlcmnt(j,:)  = temp1;
                    
                    ICP(j , :) = diff(ChnkPlcmnt(j,:),1,2);
                    ICP(j , (ICP(j , :)<0)) = 0;
                    for k = 1:length(ICP(j , :))-1
                        if ICP(j,k) & ICP(j,k+1)
                            ICP(j,k+1) = ICP(j,k) + ICP(j,k+1);
                        end
                    end
                end
                IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
                IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
            end
            ANA.IPI = diff(ANA.AllPressTimes,1,2);
            ANA.badPress = ANA.AllPress ~= ANA.AllResponse;
            
            
            ANA.ChnkArrang     = zeros(size(ANA.AllPressTimes));
            ANA.ChnkPlcmnt     = zeros(size(ANA.AllPressTimes));
            ANA.IPIarrangement = zeros(size(ANA.IPI));
            ANA.IPIChnkPlcmntArr  = zeros(size(ANA.IPI));
            for cl = 1:6
                ANA.ChnkArrang(ANA.seqNumb == cl , :) = repmat(ChnkArrng(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
                ANA.ChnkPlcmnt(ANA.seqNumb == cl , :) = repmat(ChnkPlcmnt(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
                ANA.IPIarrangement(ANA.seqNumb == cl , :) = repmat(IPIarrangement(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
                ANA.IPIChnkPlcmntArr(ANA.seqNumb == cl , :) = repmat(ICP(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
            end
            
            uBN = unique(ANA.BN);
            ANA.Rep = ones(size(ANA.BN));
            
            for bn = 1:length(uBN)
                
                switch bn
                    case {1 ,2 ,3} % Random
                        ANA.BT(ANA.BN == uBN(bn) , :) = zeros(length(find(ANA.BN == uBN(bn))),1);
                    case {4,5,6,7,8,9,10,11,12,13,14,15,16,30,31,45,46} % Chunk
                        ANA.BT(ANA.BN == uBN(bn) , :) = ones(length(find(ANA.BN == uBN(bn))),1);
                    case {17,18,19,27,28,29,32,33,34,42,43,44,47,48,49,57,58,59} % Intermixed
                        ANA.BT(ANA.BN == uBN(bn) , :) = 2*ones(length(find(ANA.BN == uBN(bn))),1);
                        ANA.Rep(ANA.BN == uBN(bn) , :) = repmat([1;2] , 0.5* length(find(ANA.BN == uBN(bn))), 1);
                    case {20,21,22,23,24,25,26,35,36,37,38,39,40,41,50,51,52,53,54,55,56} % CLAT
                        ANA.BT(ANA.BN == uBN(bn) , :) = 3*ones(length(find(ANA.BN == uBN(bn))),1);
                        ANA.Rep(ANA.BN == uBN(bn) , :) = repmat([1;2] , 0.5* length(find(ANA.BN == uBN(bn))), 1);
                end
            end
            
            for tn = 1:length(ANA.TN)
                % detrend the IPIs
                %ANA.IPI(tn , :) = detrend(ANA.IPI(tn , :),'linear',7) + mean(ANA.IPI(tn , :));
                ANA.xEye{tn} = nanfastsmooth(ANA.xEye{tn} , 9);
                ANA.xEye{tn}(ANA.Pupil{tn}<= .95 *nanmedian(ANA.Pupil{tn})) = NaN;
                ANA.yEye{tn}(ANA.Pupil{tn}<= .95 *nanmedian(ANA.Pupil{tn})) = NaN;
                
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
                ANA.pix2cm(tn,1) = abs(ANA.EyeStartPos(tn,1) - ANA.EyeEndPos(tn,1))/((sum(~isnan(ANA.AllPress(tn,:)))-1) * 2.3); % pixels per 1cm calculates for every trial
            end
            % there is high discrepency in the statrting position and
            % ending position betwen trails so it make sense to use the
            % median of the starting position and also the median of pix2cm
            % for each block
            uBN = unique(ANA.BN);
            for bn = 1:length(uBN)
                id = ANA.BN == uBN(bn) & ~ANA.isError;
                l  =length(find(ANA.BN == uBN(bn)));
                ANA.EyeStartPos(ANA.BN == uBN(bn) , :) = repmat(nanmedian(ANA.EyeStartPos(id , :)), l , 1);
                ANA.EyeEndPos(ANA.BN == uBN(bn) , :)   = repmat(nanmedian(ANA.EyeEndPos(id , :)), l , 1);
                ANA.pix2cm(ANA.BN == uBN(bn) , :)   = repmat(nanmedian(ANA.pix2cm(id , :)), l , 1);
            end
            if ~isfield (ANA , 'PPDx')
                for tn = 1:length(ANA.TN)
                    ppdx = (ANA.EyeEndPos(tn,1) - ANA.EyeStartPos(tn,1))/36;
                    ANA.PPDx{tn,1}(:,1) = ppdx*ones(length(ANA.xEye{tn}) , 1);
                end
            end
            

            ANA.isgood = ones(length(ANA.TN) , 1);
            window  = 12;
            for tn = 1:length(ANA.TN)
                switch ANA.seqNumb(tn,1)
                    case {0,1,2,3,4,5,6}
                        ANA.seqlength(tn,1) = 14;
                    case {102 , 202}
                        ANA.seqlength(tn,1) = 2;
                    case {103 , 203}
                        ANA.seqlength(tn,1) = 3;
                    case {104 , 204}
                        ANA.seqlength(tn,1) = 4;
                end
                
                % convert eye movements to cms knowing that the distance between the two pluses is 29.9 cm (2.3 cm between digits)
                ANA.xEyePosDigit{tn , 1}(:,1) = [(ANA.xEye{tn} - repmat(ANA.EyeStartPos(tn,1),length(ANA.xEye{tn}) , 1))...
                    ./repmat(ANA.pix2cm(tn,1),length(ANA.xEye{tn}) , 1)];
                ANA.xEyePosDigit{tn , 1} = 1 + ANA.xEyePosDigit{tn , 1}/2.3;
                
                if ANA.AllPressIdx(tn,ANA.seqlength(tn)) + window + 5 > length(ANA.xEyePosDigit{tn}) | ANA.seqlength(tn) ~= sum(~isnan(ANA.AllPressIdx(tn , :)))
                    ANA.isgood(tn) = 0;
                end
                
                
                
                
                % whare eyes are at the time of press
                for fing = 1:size(ANA.AllPressIdx , 2)
                    if ~isnan(ANA.AllPressIdx(tn,fing)) & ANA.isgood(tn) & ~ANA.isError(tn)
                        ANA.EyePressTimePos(tn,fing) = ceil(nanmedian(ANA.xEyePosDigit{tn}(ANA.AllPressIdx(tn,fing)-window : ANA.AllPressIdx(tn,fing) + window)));
                        if ANA.EyePressTimePos(tn,fing) <= 1
                            ANA.EyePressTimePos(tn,fing) = 1;
                        elseif ANA.EyePressTimePos(tn,fing) >= 14
                            ANA.EyePressTimePos(tn,fing) = 14;
                        end
                    else
                        ANA.EyePressTimePos(tn,fing) = NaN;
                    end
                end
                
                
                % My way of eye velocity calculation
                %                 for w = window + 1 : length(ANA.xEyePosAng{tn , 1}) - window
                %                     id = [w - window  w + window];
                %                     if id (2) < length(ANA.xEyePosAng{tn,1})
                %                         ANA.xEyeAngVelocity{tn ,1}(w-window , 1) = abs(ANA.xEyePosAng{tn}(id(1)) - ANA.xEyePosAng{tn}(id(2)))/(window*2*2);
                %                     end
                %                 end
                
                
                % find eye angular velocity over the 200 ms (100 samples) vicinity of presses
                ANA.xEyePressAngVelocity(tn , :) = NaN * ones(1,14);
                ANA.xEyeVelEstCnkPlcmnt(tn,:) = NaN * ones(1,14);
                
                % create a smooth press time series
                presses = [0 : length(find(ANA.AllPress(tn , :)))+1];
                chunks  = [1 , ANA.ChnkPlcmnt(tn , :) , 1];
                ANA.PressPressVelocity(tn , :) = NaN * ones(1,14);
                ANA.SaccPerSec(tn,1) = NaN;
                if ANA.isgood(tn)
                    ANA.PressTimeSeries{tn,1} = zeros(size(ANA.xEye{tn}));
                    IND = [0 , ANA.AllPressIdx(tn , 1:ANA.seqlength(tn)) , length(ANA.PressTimeSeries{tn,1})];
                    ANA.PressTimeSeries{tn,1}(1 : IND(2)) = NaN*ANA.PressTimeSeries{tn,1}(1 : IND(2));
                    ANA.PressTimeSeries{tn,1}(IND(end-1) : IND(end)) = NaN*ANA.PressTimeSeries{tn,1}(IND(end-1) : IND(end));
                    
                    for ind = 2:length(IND) - 2
                        temp = linspace(presses(ind) , presses(ind+1) , IND(ind + 1) - IND(ind)+1);
                        ANA.PressTimeSeries{tn,1}(IND(ind)+1 : IND(ind + 1)) = temp(1:end-1);
                    end
                    temp = [NaN; diff(ANA.PressTimeSeries{tn,1})];
                    for w = 3 : length(ANA.PressTimeSeries{tn , 1}) - 1
                        id = [w - 1  w + 1];
                        if id (2) < length(ANA.PressTimeSeries{tn,1})
                            ANA.pressVelocity{tn ,1}(w-2 , 1) = 1000*.25*(ANA.PressTimeSeries{tn,1}(id(2)) - ANA.PressTimeSeries{tn,1}(id(1)));
                        end
                    end

                    for p = 1:ANA.seqlength(tn)
                        id = [ANA.AllPressIdx(tn, p) - window  ANA.AllPressIdx(tn, p) + window ];
                        ANA.PressPressVelocity(tn , p) = nanmean(ANA.pressVelocity{tn,1}(id(1) : id(2)));
                    end
                else
                    ANA.pressVelocity{tn ,1} = [];
                    ANA.PressTimeSeries{tn,1} = [];
                end
                
                
                %         ANA.ChnkPlcEyeVelCorr(tn , 1) = NaN;
                % Create a smooth imposed-chunk timeseries
                ANA.IPIChnkPlcmnt(tn , :) = [NaN NaN NaN NaN];
                %                 ANA.ChunkTimeSeries{tn,1} = zeros(size(ANA.xEyePosDigit{tn}));
                if ismember(ANA.seqNumb(tn), [1:6])   % Intermixed or CLAT
                    
                    ANA.IPIwithin(tn,1)  = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 2));
                    ANA.IPIbetween(tn,1) = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 1));
                    
                    ANA.IPIChnkPlcmnt(tn , 1)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 0)); % basically between
                    ANA.IPIChnkPlcmnt(tn , 2)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 1));
                    ANA.IPIChnkPlcmnt(tn , 3)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 2));
                    ANA.IPIChnkPlcmnt(tn , 4)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 3));
                    
                else
                    ANA.IPIwithin(tn,1)  = NaN;
                    ANA.IPIbetween(tn,1) = NaN;
                end
                
                
                % estimate chunk bounries as the ones that are 20% of a std above the rest of the IPIs
                
                thresh = .3 * std(ANA.IPI(tn , :));
                [dum , estChnkBndry] = findpeaks(ANA.IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
                
                if ~isempty(estChnkBndry)
                    goodpeak = ones(1, length(estChnkBndry));
                    for cb = 1:length(estChnkBndry)
                        if ANA.IPI(tn , estChnkBndry(cb)) < nanmean(ANA.IPI(tn  ,:)) + thresh
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
                
                %                 ANA.ChnkPlcmntDist(tn,:)     = [NaN NaN NaN NaN NaN NaN];
                %                 ANA.ChnkPlcmntCorr(tn,:)     = [NaN NaN NaN NaN NaN NaN];
                %                 ANA.estChnkPlcEyeVelCorr(tn , 1) = NaN;
                if sum(ANA.estChnkBndry(tn,:))>0
                    ANA.estChnkPlcmnt(tn,logical(ANA.estChnkBndry(tn,:))) = 1;
                    ANA.estChnkPlcmnt(tn,1) = 1;
                    dum = find(~ANA.estChnkPlcmnt(tn,:));
                    for h = 1:length(dum)
                        ANA.estChnkPlcmnt(tn,dum(h)) = ANA.estChnkPlcmnt(tn,dum(h)-1) + 1;
                    end
                end
                
                
                %% Eye stuff
                
                ANA.xEyePosAng{tn ,1}(:,1) = ANA.xEye{tn}./ANA.PPDx{tn};   % rad2deg(atan((ANA.xEyePosDigit{tn , 1}-15) /46));
                for samp = 3:length(ANA.xEye{tn})-2
                    ANA.xEyeAngVelocity{tn,1}(samp-1 , :) = -500*(ANA.xEye{tn}(samp - 1) - ANA.xEye{tn}(samp + 1))/(2*ANA.PPDx{tn}(samp));
                end
                for samp = 3:length(ANA.xEye{tn})-2
                    ANA.xEyeAngAccel{tn,1}(samp-1 , :) = 500^2*(ANA.xEye{tn}(samp - 2) - 2*ANA.xEye{tn}(samp) + ANA.xEye{tn}(samp +2))/(4*ANA.PPDx{tn}(samp));
                end
                for p = 1:sum(~isnan(ANA.AllPressIdx(tn , :)))
                    id = [ANA.AllPressIdx(tn, p) - window-2  ANA.AllPressIdx(tn, p) + window-2 ];
                    if id (2) < length(ANA.xEyeAngVelocity{tn,1})
                        ANA.xEyePressAngVelocity(tn , p) = nanmean(ANA.xEyeAngVelocity{tn}(id(1) : id(2)));
                    end
                end
                
                % find mean eye angular velocity in terms of imposed chunk placement
                for cp  =  1:4
                    ANA.xEyeVelCnkPlcmnt(tn,cp)  = nanmean(ANA.xEyePressAngVelocity(tn,ANA.ChnkPlcmnt(tn,:) == cp));
                end
                
                % find mean eye angular velocity in terms of estimated chunk placement
                for cp = 1:length(unique(ANA.estChnkPlcmnt(tn , :)))
                    ANA.xEyeVelEstCnkPlcmnt(tn,cp)  = nanmean(ANA.xEyePressAngVelocity(tn,ANA.estChnkPlcmnt(tn,:) == cp));
                end
                
                
                signal = abs(ANA.xEyeAngVelocity{tn});
                if ANA.isgood(tn) & ismember(ANA.seqNumb(tn) , [0:6])
                    tempVel = ANA.xEyeAngVelocity{tn}((ANA.AllPressIdx(tn , 1)+2)-window : (ANA.AllPressIdx(tn , ANA.seqlength(tn))+2) + window);
                    tempPos = ANA.xEyePosAng{tn}((ANA.AllPressIdx(tn , 1)+2)-window : (ANA.AllPressIdx(tn , ANA.seqlength(tn))+2)+window);
                    ANA.SaccFlag{tn,1} = ones(size(tempVel));
                    ANA.isSaccWhilePress(tn , :) = NaN * ones(1,14);
                    [pks , locs] = findpeaks(abs(tempVel) , 'MinPeakHeight',50,'MinPeakDistance',50);
                    if length(locs) == 0
                        ANA.SaccDuration{tn,1} = [];
                        ANA.SaccPeakVel{tn,1} = [];
                        ANA.SaccAmplitude{tn,1} = [];
                        ANA.NumSaccs(tn,1) = NaN;
                    else
                        for l = 1:length(locs)
                            cntr = l;
                            peakHeight = tempVel(locs(l));
                            flrs = flipud(tempVel);
                            x1{l} = find(flrs(end-locs(l)+1:end) <= 0.05*peakHeight, 1, 'first');
                            x2{l} = find(tempVel(locs(l):end) <= 0.05*peakHeight, 1, 'first');
                            
                            n = 1;
                            while isempty(x1{l})
                                fac = .1*n;
                                x1{l} = find(flrs(end-locs(l)+1:end) <= fac*peakHeight, 1, 'first');
                                n = n+1;
                            end
                            n = 1;
                            while isempty(x2{l})
                                fac = .1*n;
                                x2{l} = find(tempVel(locs(l):end) <= fac*peakHeight, 1, 'first');
                                n = n+1;
                            end
                            
                            ANA.SaccFlag{tn,1}(locs(l) - x1{l}+1 : locs(l) + x2{l}-1) = 0;
                            ANA.SaccDuration{tn,1}(cntr , 1) = (x1{l}+x2{l})*2;
                            
                            
                            
                            ANA.SaccPeakVel{tn,1}(cntr , 1) = pks(l);
                            ANA.SaccAmplitude{tn,1}(cntr , 1) = abs(tempPos(locs(l) - x1{l} +1) - tempPos(locs(l) + x2{l}-1));
                        end
                        ANA.NumSaccs(tn,1) = l;
                        
                        if ismember(ANA.seqNumb(tn) , [0:6]) & ANA.isgood(tn) & ~ANA.isError(tn)
                            for fing = 1:ANA.seqlength(tn)
                                %                                 if fing == 1
                                %                                     fid = ANA.AllPressIdx(tn , fing) - ANA.AllPressIdx(tn , 1)+1;
                                %                                     ANA.isSaccWhilePress(tn , fing) = nanmedian(ANA.SaccFlag{tn}(fid : fid+window));
                                %                                 else
                                fid = ANA.AllPressIdx(tn , fing) - ANA.AllPressIdx(tn , 1)+1 + window;
                                ANA.isSaccWhilePress(tn , fing) = nanmedian(ANA.SaccFlag{tn}(fid-window : fid+window));
                                %                                 end
                            end
                        end
                        
                        clear x1 x2
                    end
                    ANA.SaccFlag{tn,1}(isnan(tempVel)) = NaN;
                    tempdiff = diff(ANA.SaccFlag{tn,1});
                    a = find(tempdiff == 1  );
                    b = find(tempdiff == -1);
                    if length(a) > length(b)
                        a = a(1+length(a) - length(b) : end);
                    elseif length(b)>length(a)
                        b = b(1:end -(length(b) - length(a)));
                    elseif size(a) == size(b)
                         ANA.EyeFixDuration{tn  ,1} = 2*(b(2:end) - a(1:end-1));
                    else
                        a = a';
                        ANA.EyeFixDuration{tn  ,1} = NaN;
                    end
                    
                    ANA.EyeFixDigit{tn , 1} = ANA.xEyePosDigit{tn}((ANA.AllPressIdx(tn , 1)+2)-window :ANA.AllPressIdx(tn , ANA.seqlength(tn)) + 2 + window) .* ANA.SaccFlag{tn};
                    for p = 1:14
                        id = ANA.EyeFixDigit{tn , 1}<=p & ANA.EyeFixDigit{tn , 1}>p-1;
                        ANA.EyeFixDigit{tn , 1}(id) = p;
                    end
                    %                     if tn == 49
                    %                         keyboard
                    %                     end
                    ANA.EyeFixDigit{tn}(ANA.EyeFixDigit{tn} > 14) = 14;
                    ANA.EyeFixDigit{tn}(ANA.EyeFixDigit{tn} <= 0) = 0;
                    ANA.EyeFixDigit{tn}(ANA.EyeFixDigit{tn} == 0) = NaN;
                    ANA.SaccPerSec(tn,1) = 1000 *ANA.NumSaccs(tn,1)/(ANA.AllPressTimes(tn , ANA.seqlength(tn,1)) - ANA.AllPressTimes(tn , 1));
                    
                else
                    ANA.SaccFlag{tn,1} = [];
                    ANA.EyeFixDuration{tn  ,1} = [];
                    ANA.SaccPeakVel{tn,1} = [];
                    ANA.SaccAmplitude{tn,1} = [];
                    ANA.SaccDuration{tn ,1} = [];
                end
                
            end
            if i ==1
                Dout = ANA;
            else
                Dout = addstruct(Dout , ANA);
            end
            
        end
%%
%%
%%
%%
    case 'single_subj'
        subjnums = find(ismember(subj_name , SubjCodes));
         for i = subjnums
            clear ANA
            if getdat
                ANA = se1_subj(subj_name{i} , 0);
            else
                ANA = getrow(Dall , Dall.SN == i);
            end
            BlPerDay = {[1 :14] [15 :29] [30 :44] [45 :59]};
            ANA.Day = zeros(size(ANA.BN));
            for d = 1:length(BlPerDay)
                ANA.Day(ismember(ANA.BN , BlPerDay{d})) = d;
            end
            GroupCode = str2num(subj_name{i}(3));
            ANA.SN(1:length(ANA.BN) , :) = i;
            ANA.Group(1:length(ANA.BN),:) = GroupCode;
            ChnkArrng = zeros(6,14);
            Chnkplcmnt = zeros(6,14);
            if GroupCode == 1
                ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group1 , :);
                for j = 1:size(ChnkArrang , 1)
                    temp = [];
                    temp1 = [];
                    for k = 1:length(find(ChnkArrang(j , :)))
                        temp = [temp k*ones(1,ChnkArrang(j,k))];
                        temp1 = [temp1 1:ChnkArrang(j,k)];
                    end
                    ChnkArrng(j , :) = temp;
                    ChnkPlcmnt(j,:)  = temp1;
                    
                    ICP(j , :) = diff(ChnkPlcmnt(j,:),1,2);
                    ICP(j , (ICP(j , :)<0)) = 0;
                    for k = 1:length(ICP(j , :))-1
                        if ICP(j,k) & ICP(j,k+1)
                            ICP(j,k+1) = ICP(j,k) + ICP(j,k+1);
                        end
                    end
                end
                IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
                IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
            elseif GroupCode == 2
                ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group2 , :);
                for j = 1:size(ChnkArrang , 1)
                    temp = [];
                    temp1 = [];
                    for k = 1:length(find(ChnkArrang(j , :)))
                        temp = [temp k*ones(1,ChnkArrang(j,k))];
                        temp1 = [temp1 1:ChnkArrang(j,k)];
                    end
                    ChnkArrng(j , :) = temp;
                    ChnkPlcmnt(j,:)  = temp1;
                    
                    ICP(j , :) = diff(ChnkPlcmnt(j,:),1,2);
                    ICP(j , (ICP(j , :)<0)) = 0;
                    for k = 1:length(ICP(j , :))-1
                        if ICP(j,k) & ICP(j,k+1)
                            ICP(j,k+1) = ICP(j,k) + ICP(j,k+1);
                        end
                    end
                end
                IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
                IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
            end
            ANA.IPI = diff(ANA.AllPressTimes,1,2);
            ANA.badPress = ANA.AllPress ~= ANA.AllResponse;
            
            
            ANA.ChnkArrang     = zeros(size(ANA.AllPressTimes));
            ANA.ChnkPlcmnt     = zeros(size(ANA.AllPressTimes));
            ANA.IPIarrangement = zeros(size(ANA.IPI));
            ANA.IPIChnkPlcmntArr  = zeros(size(ANA.IPI));
            for cl = 1:6
                ANA.ChnkArrang(ANA.seqNumb == cl , :) = repmat(ChnkArrng(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
                ANA.ChnkPlcmnt(ANA.seqNumb == cl , :) = repmat(ChnkPlcmnt(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
                ANA.IPIarrangement(ANA.seqNumb == cl , :) = repmat(IPIarrangement(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
                ANA.IPIChnkPlcmntArr(ANA.seqNumb == cl , :) = repmat(ICP(cl , :) , length(find(ANA.seqNumb == cl)) , 1);
            end
            
            uBN = unique(ANA.BN);
            ANA.Rep = ones(size(ANA.BN));
            
            for bn = 1:length(uBN)
                
                switch bn
                    case {1 ,2 ,3} % Random
                        ANA.BT(ANA.BN == uBN(bn) , :) = zeros(length(find(ANA.BN == uBN(bn))),1);
                    case {4,5,6,7,8,9,10,11,12,13,14,15,16,30,31,45,46} % Chunk
                        ANA.BT(ANA.BN == uBN(bn) , :) = ones(length(find(ANA.BN == uBN(bn))),1);
                    case {17,18,19,27,28,29,32,33,34,42,43,44,47,48,49,57,58,59} % Intermixed
                        ANA.BT(ANA.BN == uBN(bn) , :) = 2*ones(length(find(ANA.BN == uBN(bn))),1);
                        ANA.Rep(ANA.BN == uBN(bn) , :) = repmat([1;2] , 0.5* length(find(ANA.BN == uBN(bn))), 1);
                    case {20,21,22,23,24,25,26,35,36,37,38,39,40,41,50,51,52,53,54,55,56} % CLAT
                        ANA.BT(ANA.BN == uBN(bn) , :) = 3*ones(length(find(ANA.BN == uBN(bn))),1);
                        ANA.Rep(ANA.BN == uBN(bn) , :) = repmat([1;2] , 0.5* length(find(ANA.BN == uBN(bn))), 1);
                end
            end
            
            for tn = 1:length(ANA.TN)
                % detrend the IPIs
                %ANA.IPI(tn , :) = detrend(ANA.IPI(tn , :),'linear',7) + mean(ANA.IPI(tn , :));
                ANA.xEye{tn} = nanfastsmooth(ANA.xEye{tn} , 9);
                ANA.xEye{tn}(ANA.Pupil{tn}<= .95 *nanmedian(ANA.Pupil{tn})) = NaN;
                ANA.yEye{tn}(ANA.Pupil{tn}<= .95 *nanmedian(ANA.Pupil{tn})) = NaN;
                
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
                ANA.pix2cm(tn,1) = abs(ANA.EyeStartPos(tn,1) - ANA.EyeEndPos(tn,1))/((sum(~isnan(ANA.AllPress(tn,:)))-1) * 2.3); % pixels per 1cm calculates for every trial
            end
            % there is high discrepency in the statrting position and
            % ending position betwen trails so it make sense to use the
            % median of the starting position and also the median of pix2cm
            % for each block
            uBN = unique(ANA.BN);
            for bn = 1:length(uBN)
                id = ANA.BN == uBN(bn) & ~ANA.isError;
                l  =length(find(ANA.BN == uBN(bn)));
                ANA.EyeStartPos(ANA.BN == uBN(bn) , :) = repmat(nanmedian(ANA.EyeStartPos(id , :)), l , 1);
                ANA.EyeEndPos(ANA.BN == uBN(bn) , :)   = repmat(nanmedian(ANA.EyeEndPos(id , :)), l , 1);
                ANA.pix2cm(ANA.BN == uBN(bn) , :)   = repmat(nanmedian(ANA.pix2cm(id , :)), l , 1);
            end
            if ~isfield (ANA , 'PPDx')
                for tn = 1:length(ANA.TN)
                    ppdx = (ANA.EyeEndPos(tn,1) - ANA.EyeStartPos(tn,1))/36;
                    ANA.PPDx{tn,1}(:,1) = ppdx*ones(length(ANA.xEye{tn}) , 1);
                end
            end
            

            ANA.isgood = ones(length(ANA.TN) , 1);
            window  = 12;
            for tn = 1:length(ANA.TN)
                switch ANA.seqNumb(tn,1)
                    case {0,1,2,3,4,5,6}
                        ANA.seqlength(tn,1) = 14;
                    case {102 , 202}
                        ANA.seqlength(tn,1) = 2;
                    case {103 , 203}
                        ANA.seqlength(tn,1) = 3;
                    case {104 , 204}
                        ANA.seqlength(tn,1) = 4;
                end
                
                % convert eye movements to cms knowing that the distance between the two pluses is 29.9 cm (2.3 cm between digits)
                ANA.xEyePosDigit{tn , 1}(:,1) = [(ANA.xEye{tn} - repmat(ANA.EyeStartPos(tn,1),length(ANA.xEye{tn}) , 1))...
                    ./repmat(ANA.pix2cm(tn,1),length(ANA.xEye{tn}) , 1)];
                ANA.xEyePosDigit{tn , 1} = 1 + ANA.xEyePosDigit{tn , 1}/2.3;
                
                if ANA.AllPressIdx(tn,ANA.seqlength(tn)) + window + 5 > length(ANA.xEyePosDigit{tn}) | ANA.seqlength(tn) ~= sum(~isnan(ANA.AllPressIdx(tn , :)))
                    ANA.isgood(tn) = 0;
                end
                
                
                
                
                % whare eyes are at the time of press
                for fing = 1:size(ANA.AllPressIdx , 2)
                    if ~isnan(ANA.AllPressIdx(tn,fing)) & ANA.isgood(tn) & ~ANA.isError(tn)
                        ANA.EyePressTimePos(tn,fing) = ceil(nanmedian(ANA.xEyePosDigit{tn}(ANA.AllPressIdx(tn,fing)-window : ANA.AllPressIdx(tn,fing) + window)));
                        if ANA.EyePressTimePos(tn,fing) <= 1
                            ANA.EyePressTimePos(tn,fing) = 1;
                        elseif ANA.EyePressTimePos(tn,fing) >= 14
                            ANA.EyePressTimePos(tn,fing) = 14;
                        end
                    else
                        ANA.EyePressTimePos(tn,fing) = NaN;
                    end
                end
                
                
                % My way of eye velocity calculation
                %                 for w = window + 1 : length(ANA.xEyePosAng{tn , 1}) - window
                %                     id = [w - window  w + window];
                %                     if id (2) < length(ANA.xEyePosAng{tn,1})
                %                         ANA.xEyeAngVelocity{tn ,1}(w-window , 1) = abs(ANA.xEyePosAng{tn}(id(1)) - ANA.xEyePosAng{tn}(id(2)))/(window*2*2);
                %                     end
                %                 end
                
                
                % find eye angular velocity over the 200 ms (100 samples) vicinity of presses
                ANA.xEyePressAngVelocity(tn , :) = NaN * ones(1,14);
                ANA.xEyeVelEstCnkPlcmnt(tn,:) = NaN * ones(1,14);
                
                % create a smooth press time series
                presses = [0 : length(find(ANA.AllPress(tn , :)))+1];
                chunks  = [1 , ANA.ChnkPlcmnt(tn , :) , 1];
                ANA.PressPressVelocity(tn , :) = NaN * ones(1,14);
                ANA.SaccPerSec(tn,1) = NaN;
                if ANA.isgood(tn)
                    ANA.PressTimeSeries{tn,1} = zeros(size(ANA.xEye{tn}));
                    IND = [0 , ANA.AllPressIdx(tn , 1:ANA.seqlength(tn)) , length(ANA.PressTimeSeries{tn,1})];
                    ANA.PressTimeSeries{tn,1}(1 : IND(2)) = NaN*ANA.PressTimeSeries{tn,1}(1 : IND(2));
                    ANA.PressTimeSeries{tn,1}(IND(end-1) : IND(end)) = NaN*ANA.PressTimeSeries{tn,1}(IND(end-1) : IND(end));
                    
                    for ind = 2:length(IND) - 2
                        temp = linspace(presses(ind) , presses(ind+1) , IND(ind + 1) - IND(ind)+1);
                        ANA.PressTimeSeries{tn,1}(IND(ind)+1 : IND(ind + 1)) = temp(1:end-1);
                    end
                    temp = [NaN; diff(ANA.PressTimeSeries{tn,1})];
                    for w = 3 : length(ANA.PressTimeSeries{tn , 1}) - 1
                        id = [w - 1  w + 1];
                        if id (2) < length(ANA.PressTimeSeries{tn,1})
                            ANA.pressVelocity{tn ,1}(w-2 , 1) = 1000*.25*(ANA.PressTimeSeries{tn,1}(id(2)) - ANA.PressTimeSeries{tn,1}(id(1)));
                        end
                    end

                    for p = 1:ANA.seqlength(tn)
                        id = [ANA.AllPressIdx(tn, p) - window  ANA.AllPressIdx(tn, p) + window ];
                        ANA.PressPressVelocity(tn , p) = nanmean(ANA.pressVelocity{tn,1}(id(1) : id(2)));
                    end
                else
                    ANA.pressVelocity{tn ,1} = [];
                    ANA.PressTimeSeries{tn,1} = [];
                end
                
                
                %         ANA.ChnkPlcEyeVelCorr(tn , 1) = NaN;
                % Create a smooth imposed-chunk timeseries
                ANA.IPIChnkPlcmnt(tn , :) = [NaN NaN NaN NaN];
                %                 ANA.ChunkTimeSeries{tn,1} = zeros(size(ANA.xEyePosDigit{tn}));
                if ismember(ANA.seqNumb(tn), [1:6])   % Intermixed or CLAT
                    
                    ANA.IPIwithin(tn,1)  = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 2));
                    ANA.IPIbetween(tn,1) = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 1));
                    
                    ANA.IPIChnkPlcmnt(tn , 1)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 0)); % basically between
                    ANA.IPIChnkPlcmnt(tn , 2)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 1));
                    ANA.IPIChnkPlcmnt(tn , 3)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 2));
                    ANA.IPIChnkPlcmnt(tn , 4)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 3));
                    
                else
                    ANA.IPIwithin(tn,1)  = NaN;
                    ANA.IPIbetween(tn,1) = NaN;
                end
                
                
                % estimate chunk bounries as the ones that are 20% of a std above the rest of the IPIs
                
                thresh = .3 * std(ANA.IPI(tn , :));
                [dum , estChnkBndry] = findpeaks(ANA.IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
                
                if ~isempty(estChnkBndry)
                    goodpeak = ones(1, length(estChnkBndry));
                    for cb = 1:length(estChnkBndry)
                        if ANA.IPI(tn , estChnkBndry(cb)) < nanmean(ANA.IPI(tn  ,:)) + thresh
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
                
                %                 ANA.ChnkPlcmntDist(tn,:)     = [NaN NaN NaN NaN NaN NaN];
                %                 ANA.ChnkPlcmntCorr(tn,:)     = [NaN NaN NaN NaN NaN NaN];
                %                 ANA.estChnkPlcEyeVelCorr(tn , 1) = NaN;
                if sum(ANA.estChnkBndry(tn,:))>0
                    ANA.estChnkPlcmnt(tn,logical(ANA.estChnkBndry(tn,:))) = 1;
                    ANA.estChnkPlcmnt(tn,1) = 1;
                    dum = find(~ANA.estChnkPlcmnt(tn,:));
                    for h = 1:length(dum)
                        ANA.estChnkPlcmnt(tn,dum(h)) = ANA.estChnkPlcmnt(tn,dum(h)-1) + 1;
                    end
                end
                
                
                %% Eye stuff
                
                ANA.xEyePosAng{tn ,1}(:,1) = ANA.xEye{tn}./ANA.PPDx{tn};   % rad2deg(atan((ANA.xEyePosDigit{tn , 1}-15) /46));
                for samp = 3:length(ANA.xEye{tn})-2
                    ANA.xEyeAngVelocity{tn,1}(samp-1 , :) = -500*(ANA.xEye{tn}(samp - 1) - ANA.xEye{tn}(samp + 1))/(2*ANA.PPDx{tn}(samp));
                end
                for samp = 3:length(ANA.xEye{tn})-2
                    ANA.xEyeAngAccel{tn,1}(samp-1 , :) = 500^2*(ANA.xEye{tn}(samp - 2) - 2*ANA.xEye{tn}(samp) + ANA.xEye{tn}(samp +2))/(4*ANA.PPDx{tn}(samp));
                end
                for p = 1:sum(~isnan(ANA.AllPressIdx(tn , :)))
                    id = [ANA.AllPressIdx(tn, p) - window-2  ANA.AllPressIdx(tn, p) + window-2 ];
                    if id (2) < length(ANA.xEyeAngVelocity{tn,1})
                        ANA.xEyePressAngVelocity(tn , p) = nanmean(ANA.xEyeAngVelocity{tn}(id(1) : id(2)));
                    end
                end
                
                % find mean eye angular velocity in terms of imposed chunk placement
                for cp  =  1:4
                    ANA.xEyeVelCnkPlcmnt(tn,cp)  = nanmean(ANA.xEyePressAngVelocity(tn,ANA.ChnkPlcmnt(tn,:) == cp));
                end
                
                % find mean eye angular velocity in terms of estimated chunk placement
                for cp = 1:length(unique(ANA.estChnkPlcmnt(tn , :)))
                    ANA.xEyeVelEstCnkPlcmnt(tn,cp)  = nanmean(ANA.xEyePressAngVelocity(tn,ANA.estChnkPlcmnt(tn,:) == cp));
                end
                
                
                signal = abs(ANA.xEyeAngVelocity{tn});
                if ANA.isgood(tn) & ismember(ANA.seqNumb(tn) , [0:6])
                    tempVel = ANA.xEyeAngVelocity{tn}((ANA.AllPressIdx(tn , 1)+2)-window : (ANA.AllPressIdx(tn , ANA.seqlength(tn))+2) + window);
                    tempPos = ANA.xEyePosAng{tn}((ANA.AllPressIdx(tn , 1)+2)-window : (ANA.AllPressIdx(tn , ANA.seqlength(tn))+2)+window);
                    ANA.SaccFlag{tn,1} = ones(size(tempVel));
                    ANA.isSaccWhilePress(tn , :) = NaN * ones(1,14);
                    [pks , locs] = findpeaks(abs(tempVel) , 'MinPeakHeight',50,'MinPeakDistance',50);
                    if length(locs) == 0
                        ANA.SaccDuration{tn,1} = [];
                        ANA.SaccPeakVel{tn,1} = [];
                        ANA.SaccAmplitude{tn,1} = [];
                        ANA.NumSaccs(tn,1) = NaN;
                    else
                        for l = 1:length(locs)
                            cntr = l;
                            peakHeight = tempVel(locs(l));
                            flrs = flipud(tempVel);
                            x1{l} = find(flrs(end-locs(l)+1:end) <= 0.05*peakHeight, 1, 'first');
                            x2{l} = find(tempVel(locs(l):end) <= 0.05*peakHeight, 1, 'first');
                            
                            n = 1;
                            while isempty(x1{l})
                                fac = .1*n;
                                x1{l} = find(flrs(end-locs(l)+1:end) <= fac*peakHeight, 1, 'first');
                                n = n+1;
                            end
                            n = 1;
                            while isempty(x2{l})
                                fac = .1*n;
                                x2{l} = find(tempVel(locs(l):end) <= fac*peakHeight, 1, 'first');
                                n = n+1;
                            end
                            
                            ANA.SaccFlag{tn,1}(locs(l) - x1{l}+1 : locs(l) + x2{l}-1) = 0;
                            ANA.SaccDuration{tn,1}(cntr , 1) = (x1{l}+x2{l})*2;
                            
                            
                            
                            ANA.SaccPeakVel{tn,1}(cntr , 1) = pks(l);
                            ANA.SaccAmplitude{tn,1}(cntr , 1) = abs(tempPos(locs(l) - x1{l} +1) - tempPos(locs(l) + x2{l}-1));
                        end
                        ANA.NumSaccs(tn,1) = l;
                        
                        if ismember(ANA.seqNumb(tn) , [0:6]) & ANA.isgood(tn) & ~ANA.isError(tn)
                            for fing = 1:ANA.seqlength(tn)
%                                 if fing == 1
%                                     fid = ANA.AllPressIdx(tn , fing) - ANA.AllPressIdx(tn , 1)+1;
%                                     ANA.isSaccWhilePress(tn , fing) = nanmedian(ANA.SaccFlag{tn}(fid : fid+window));
%                                 else
                                    fid = ANA.AllPressIdx(tn , fing) - ANA.AllPressIdx(tn , 1)+1 + window;
                                    ANA.isSaccWhilePress(tn , fing) = nanmedian(ANA.SaccFlag{tn}(fid-window : fid+window));
%                                 end
                            end
                        end
                        
                        clear x1 x2
                    end
                    ANA.SaccFlag{tn,1}(isnan(tempVel)) = NaN;
                    tempdiff = diff(ANA.SaccFlag{tn,1});
                    a = find(tempdiff == 1);
                    b = find(tempdiff == -1);
                    if length(a) > length(b)
                        a = a(1+length(a) - length(b) : end);
                    elseif length(b)>length(a)
                        b = b(1:end -(length(b) - length(a)));
                    end
                    if size(a) == size(b)
                         ANA.EyeFixDuration{tn  ,1} = 2*(b(2:end) - a(1:end-1));
                    else
                        a = a';
                        ANA.EyeFixDuration{tn  ,1} = NaN;
                    end
                    
                    ANA.EyeFixDigit{tn , 1} = ANA.xEyePosDigit{tn}((ANA.AllPressIdx(tn , 1)+2)-window :ANA.AllPressIdx(tn , ANA.seqlength(tn)) + 2 + window) .* ANA.SaccFlag{tn};
                    for p = 1:14
                        id = ANA.EyeFixDigit{tn , 1}<=p & ANA.EyeFixDigit{tn , 1}>p-1;
                        ANA.EyeFixDigit{tn , 1}(id) = p;
                    end
                    %                     if tn == 49
                    %                         keyboard
                    %                     end
                    ANA.EyeFixDigit{tn}(ANA.EyeFixDigit{tn} > 14) = 14;
                    ANA.EyeFixDigit{tn}(ANA.EyeFixDigit{tn} <= 0) = 0;
                    ANA.EyeFixDigit{tn}(ANA.EyeFixDigit{tn} == 0) = NaN;
                    ANA.SaccPerSec(tn,1) = 1000 *ANA.NumSaccs(tn,1)/(ANA.AllPressTimes(tn , ANA.seqlength(tn,1)) - ANA.AllPressTimes(tn , 1));
                    
                else
                    ANA.SaccFlag{tn,1} = [];
                    ANA.EyeFixDuration{tn  ,1} = [];
                    ANA.SaccPeakVel{tn,1} = [];
                    ANA.SaccAmplitude{tn,1} = [];
                    ANA.SaccDuration{tn ,1} = [];
                end
                
            end
            if i ==min(subjnums)
                Dout = ANA;
            else
                Dout = addstruct(Dout , ANA);
            end
            
        end
end
