ANA = Dall;
for tn = 1:length(ANA.TN)
    % detrend the IPIs
    %ANA.IPI(tn , :) = detrend(ANA.IPI(tn , :),'linear',7) + mean(ANA.IPI(tn , :));
    ANA.xEye{tn} = nanfastsmooth(ANA.xEye{tn} , 9);
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
    ANA.pix2cm(tn,1) = abs(ANA.EyeStartPos(tn,1) - ANA.EyeEndPos(tn,1))/((sum(~isnan(ANA.AllPress(tn,:)))-1) * 2.3); % pixels per 1cm calculates for every trial
end
if ~isfield (ANA , 'PPDx')
    for tn = 1:length(ANA.TN)
    ppdx = (Dall.EyeEndPos(tn) - Dall.EyeStartPos(tn))/36;
    ANA.PPDx{tn,1}(:,1) = ppdx*ones(length(ANA.xEye{tn}) , 1);
    end
end
for tn = 1:length(ANA.TN)
    
    ANA.xEyePosAng{tn ,1}(1,:) = ANA.xEye{tn}./ANA.PPDx{tn};   % rad2deg(atan((ANA.xEyePosCM{tn , 1}-15) /46));
    for samp = 2:length(ANA.xEye{tn})-2
        ANA.xEyeAngVelocity{tn}(samp-1 , :) = -500*(ANA.xEye{tn}(samp - 1) - ANA.xEye{tn}(samp + 1))/(2*ANA.PPDx{tn}(samp));
    end
    for samp = 3:length(ANA.xEye{tn})-2
        ANA.xEyeAngAccel{tn}(samp-1 , :) = 500^2*(ANA.xEye{tn}(samp - 2) - 2*ANA.xEye{tn}(samp) + ANA.xEye{tn}(samp +2))/(4*ANA.PPDx{tn}(samp));
    end
    [pks , locs] = findpeaks(abs(ANA.xEyeAngVelocity{tn}) , 'MinPeakHeight',30,'MinPeakDistance',50);
    locs = locs(pks<400);
    locs = locs(2:end - 1);
    pks  = pks(pks<400);
    pks = pks(2:end - 1);
    
    signal = abs(ANA.xEyeAngVelocity{tn});
    ANA.SaccFlag{tn} = zeros(size(ANA.xEyeAngVelocity{tn}));
    cntr = 1;
    for l = 1:length(locs)
        
        peakHeight = signal(locs(l));
        flrs = flipud(signal);
        x1 = find(flrs(end-locs(l)+1:end) <= 0.05*peakHeight, 1, 'first');
        x2 = find(signal(locs(l):end) <= 0.05*peakHeight, 1, 'first');
        if max(ANA.xEyeAngAccel{tn}(locs(l) - x1 : locs(l) + x2-1))>8000
            n = 1;
            while isempty(x1)
                fac = .1*n;
                x1 = find(flrs(end-locs(l)+1:end) <= fac*peakHeight, 1, 'first');
                n = n+1;
            end
            n = 1;
            while isempty(x2)
                fac = .1*n;
                x2 = find(signal(locs(l):end) <= fac*peakHeight, 1, 'first');
                n = n+1;
            end
            ANA.SaccFlag{tn}(locs(l) - x1 : locs(l) + x2-1) = 1;
            ANA.SaccWidth{tn}(cntr , 1) = x1+x2;
            ANA.SaccPeakVel{tn}(cntr , 1) = pks(l);
            ANA.SaccAmplitude{tn}(cntr , 1) = abs(ANA.xEyePosAng{tn}(locs(l) - x1) - ANA.xEyePosAng{tn}(locs(l) + x2-1));
            cntr = cntr + 1;
        end
    end
%     close all
%     plot(ANA.SaccAmplitude{tn} , ANA.SaccPeakVel{tn} , '*')
    tn
end
ANA.SaccFlag{tn}(isnan(ANA.xEyeAngVelocity{tn}) | abs(ANA.xEyeAngVelocity{tn})>400) = NaN;
ANA.SaccPerSec{tn} = 1000 *(cntr-1)/(ANA.AllPressTimes(tn , length(find(ANA.AllPressTimes(tn , :)))) - ANA.AllPressTimes(tn , 1));
