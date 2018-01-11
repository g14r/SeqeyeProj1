function out  = se1_visualize(Dall , subjnum, what, distance, calc , day , rep, GroupCode)
%%  distances:
% 'euclidean'	      Euclidean distance (default).
% 'squaredeuclidean'  Squared Euclidean distance. (This option is provided for efficiency only. It does not satisfy the triangle inequality.)
% 'seuclidean'	      Standardized Euclidean distance. Each coordinate difference between rows in X and Y is scaled by dividing by the corresponding element of the standard deviation computed from X, S=nanstd(X). To specify another value for S, use D = PDIST2(X,Y,'seuclidean',S).
% 'cityblock'	      City block metric.
% 'minkowski'         Minkowski distance. The default exponent is 2. To compute the distance with a different exponent, use D = pdist2(X,Y,'minkowski',P), where the exponent P is a scalar positive value.
% 'Chebychev'	      Chebychev distance (maximum coordinate difference).
% 'mahalanobis'	      Mahalanobis distance, using the sample covariance of X as computed by nancov. To compute the distance with a different covariance, use D = pdist2(X,Y,'mahalanobis',C) where the matrix C is symmetric and positive definite.
% 'cosine'	          One minus the cosine of the included angle between points (treated as vectors).
% 'correlation'	      One minus the sample correlation between points (treated as sequences of values).
% 'spearman'	      One minus the sample Spearman's rank correlation between observations, treated as sequences of values.
% 'hamming'           Hamming distance, the percentage of coordinates that differ.
% 'jaccard'           One minus the Jaccard coefficient, the percentage of nonzero coordinates that differ.%     case 'chunk_est_instance'
%%  Cases
%     case 'chunk_est_instance'
%     case 'eye_vel_instance'
%     case 'eye_mainSequence'
%     case 'eye_vs_press_instance'
%     case 'eyePos_eyeVel'
%     case 'IPI_dist'
%     case 'IPI_ttest_rand'
%     case 'chunk_dist'
%     case 'Mychunk_dist'
%     case 'Avg_pattern_sh3'
%     case 'eye_vel_chunkplace'
%     case 'eye_vel_seqplace'
%     case 'eye_pos_seqplace'


%     case 'eyepress_pos_traces'
%     case 'eyepress_pos_distances'
%     case 'eyepress_pos_avg_distances'


%     case 'eyepress_vel_traces'
%     case 'eyepress_vel_distances'
%     case 'eyepress_vel_avg_distances'


%     case 'run_f-tests_mt'
%     case 'dtw'
%     case 'eyeEndPos_ChunkLength'
%     case 'eyeFrstPos_ChunkLength'


%     case 'crossvaldist_pos'
%     case 'crossvaldist_vel'


%     case 'crossvaldist_pos_presstim'
%     case 'crossvaldist_vel_presstime'
%     case 'crossval_IPI_dist'
%     case 'crossvaldist_chunk'
%     case 'perveiw_benefit'
%     case 'saccades'

%%

prefix = 'se1_';
% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye1/analyze';
baseDir = '/Users/nkordjazi/Documents/SeqEye/SeqEye1/analyze';
subj_name = {'SZ1','JN2' ,'SP1','AT1','DW1','KL1','JG1','GP1','SK1' ,'NM1','VR1','PB1' , 'XX1'};
if subjnum == length(subj_name)
    %     subjnum = [1 3:length(subj_name)-1];
    subjnum = 1:length(subj_name)-1;
end


load([baseDir , '/CMB.mat'])
%load([baseDir , '/se1_all.mat'])
days = {1 ,2 ,3 ,4 ,[2:4] ,[3:4] [1:4]};
% subjnum = 2;

tid = {[1:250] [251 :750] [751:1000] [1:1000]};
switch what
    case 'chunk_est_instance'
        
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        prompt = 'Which trial?';
        tnoi = input(prompt);
        a = ANA.IPI(tnoi , :);
        
        % ad = detrend(a,'linear',7) + mean(a);
        figure('color' , [1 1 1])
        plot(a , 'LineWidth' , 6)
        hold  on
        % plot(ad)
        peakid = find(ANA.estChnkBndry(tnoi , 2:end));
        plot(peakid , a(peakid) , 'o' , 'MarkerSize' , 4, 'LineWidth' , 10)
        line([0 13] , [mean(a)+.3*std(a) mean(a)+.3*std(a)] , 'LineWidth' , 3 , 'color' , 'r')
        line([0 13] , [mean(a)-.3*std(a) mean(a)-.3*std(a)], 'LineWidth' , 3,'color' , 'r')
        line([0 13] , [mean(a) mean(a)] , 'LineStyle' , ':', 'LineWidth' , 3,'color' , 'r')
        xlabel('Inter Press Interval' , 'FontSize' , 20)
        ylabel('msec', 'FontSize' , 20)
        title('Chunk boundary detection' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:13];
        ax.XLim = [1 13];
        ax.Box = 'off';
        grid on
        out = [];
    case 'eye_vel_instance'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        prompt = 'Which trial?';
        tnoi = input(prompt);
        a = ANA.IPI(tnoi , :);
        b = ANA.xEyeAngVelocity{tnoi};
        figure('color' , [1 1 1])
        subplot(2,1,1)
        plot(a, 'LineWidth' , 3)
        hold  on
        % plot(ad)
        peakid = find(ANA.estChnkBndry(tnoi , 2:end));
        title('IPI')
        plot(peakid , a(peakid) , 'o')
        line([0 13] , [mean(a)+.2*std(a) mean(a)+.2*std(a)])
        line([0 13] , [mean(a)-.2*std(a) mean(a)-.2*std(a)])
        line([0 13] , [mean(a) mean(a)] , 'LineStyle' , ':')
        xlabel('IPI')
        ylabel('msec')
        
        subplot(2,1,2)
        plot(b, 'LineWidth' , 3)
        hold  on
        for m = 1:length(find(ANA.AllPressIdx(tnoi , :)))
            x = ANA.AllPressIdx(tnoi , m) - 50;
            colorid  = ANA.ChnkPlcmnt(tnoi , m);
            color = {'r' 'g' , 'c' 'black'};
            line([x x] , [0 max(b)], 'LineStyle' , ':' , 'LineWidth' , 3 ,'color' , color{colorid});
            ylabel('deg/sec')
        end
        title('Eye angular velocity')
        out = [];
    case 'eye_mainSequence'
        figure('color' , 'white')
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        for tn  = 1:length(ANA.TN)
            plot(ANA.SaccAmplitude{tn} , ANA.SaccPeakVel{tn} , '*')
            xlabel('Saccade Amiplitude (deg)')
            ylabel('Saccade Peak Velocity (deg/sec)')
            title('Saccade Main Sequence')
            drawnow
            pause()
        end
        out = [];
    case 'eye_vs_press_instance'
        
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        prompt = 'Which trial?';
        tn = input(prompt);
        window  = 50;
        ANA = getrow(ANA , ANA.TN ==  ANA.TN(tn) & ANA.BN ==  ANA.BN(tn));
        if ANA.AllPressIdx(14) < length(ANA.xEyePosDigit{1})
            idd   = linspace(1 , ANA.AllPressIdx(1,14)-ANA.AllPressIdx(1,1)+1 , 1000);
            press = inpaint_nans(ANA.PressTimeSeries{1}(ANA.AllPressIdx(1,1):ANA.AllPressIdx(1,14)));
            press = interp1([1 :ANA.AllPressIdx(1,14)-ANA.AllPressIdx(1,1)+1] , press , idd);
            
            idd = floor(linspace(ANA.AllPressIdx(1) , ANA.AllPressIdx(14) , 1001));
            for j = 1:length(idd)-1
                eye(j) = nanmean(ANA.xEyePosDigit{1}(idd(j) : idd(j+1)));
                %                 eyevel(j) = nanmean(ANA.xEyeAngVelocity{1}(iddv(j)-window : iddv(j+1)-window));
                %                 pressvel(j) = nanmean(ANA.pressVelocity{1}(iddv(j)-window : iddv(j+1)-window));
            end
            eyevel = eye(10 : end) - eye(1:end-10+1);
            pressvel = diff(press);
        end
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        plot(inpaint_nans(eye)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title('The Time Normalized Eye Position Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        hold on
        for i = 1:14
            [~, pt(i)] = min(abs((press - i)));
            line([pt(i) , pt(i)] , [1 14] , 'LineStyle' , ':' , 'LineWidth' , 3 , 'color' , 'r')
        end
        
        
        subplot(2,1,2)
        plot(inpaint_nans(eyevel)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Deg/Sec' , 'FontSize' , 20)
        title('The Time Normalized Eye Angular Velocity Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        %         ax.YTick = [1:14];
        %         ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        plot(inpaint_nans(press)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Finger press' , 'FontSize' , 20)
        title('The Time Normalized Finger Press Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        hold on
        for i = 1:14
            [~, pt(i)] = min(abs((press - i)));
            line([pt(i) , pt(i)] , [1 14] , 'LineStyle' , ':' , 'LineWidth' , 3 , 'color' , 'r')
            line([0 pt(i)] , [i i] , 'LineStyle' , ':' , 'LineWidth' , 3 , 'color' , 'black')
        end
        
        
        
        
        
        
        subplot(2,1,2)
        plot(inpaint_nans(pressvel)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Press/Sec' , 'FontSize' , 20)
        title('The Time Normalized Press Velocity Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        %         ax.YTick = [1:14];
        %         ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        hold on
        
        
        figure('color' , [1 1 1])
        plot(press ,eye  , 'LineWidth' , 5 )
        axis square
        grid on
        xlabel('Finger presses' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title('The Time Normalized Eye Position and Finger Press Time series' , 'FontSize' , 20)
        hold on
        line([0 14], [0 14] , 'LineStyle' , ':' , 'LineWidth' , 5 , 'color' , 'r')
        ax = gca;
        ax.XTick = [1:14];
        ax.XLim = [1 14];
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        out = [];
    case 'eyePos_eyeVel'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        prompt = 'Which trail?';
        tn = input(prompt);
        window  = 10;
        ANA = getrow(ANA , ANA.TN ==  ANA.TN(tn) & ANA.BN == ANA.BN(tn));
        if ANA.AllPressIdx(14) < length(ANA.xEyePosDigit{1})
            idd   = linspace(1 , ANA.AllPressIdx(1,14)-ANA.AllPressIdx(1,1)+1 , 1000);
            press = inpaint_nans(ANA.PressTimeSeries{1}(ANA.AllPressIdx(1,1):ANA.AllPressIdx(1,14)));
            press = interp1([1 :ANA.AllPressIdx(1,14)-ANA.AllPressIdx(1,1)+1] , press , idd);
            ANA.xEyePosDigit{1} = ANA.xEyePosDigit{1}(3:end-1);
            idd = floor(linspace(ANA.AllPressIdx(1) , ANA.AllPressIdx(14) , 1001));
            for j = 1:length(idd)-1
                eye(j) = nanmean(ANA.xEyePosDigit{1}(idd(j) : idd(j+1)));
                eyevel(j) = nanmean(ANA.xEyeAngVelocity{1}(idd(j) : idd(j+1)));
                pressvel(j) = nanmean(ANA.pressVelocity{1}(idd(j) : idd(j+1)));
            end
        end
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        plot(inpaint_nans(eye)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title('The Time Normalized Eye Position Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        hold on
        for i = 1:14
            [~, pt(i)] = min(abs((press - i)));
            line([pt(i) , pt(i)] , [1 14] , 'LineStyle' , ':' , 'LineWidth' , 3 , 'color' , 'r')
        end
        
        
        subplot(2,1,2)
        plot(inpaint_nans(eyevel)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Deg/Sec' , 'FontSize' , 20)
        title('The Time Normalized Eye Angular Velocity Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        %         ax.YTick = [1:14];
        %         ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        plot(inpaint_nans(press)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Finger press' , 'FontSize' , 20)
        title('The Time Normalized Finger Press Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        hold on
        for i = 1:14
            [~, pt(i)] = min(abs((press - i)));
            line([pt(i) , pt(i)] , [1 14] , 'LineStyle' , ':' , 'LineWidth' , 3 , 'color' , 'r')
            line([0 pt(i)] , [i i] , 'LineStyle' , ':' , 'LineWidth' , 3 , 'color' , 'black')
        end
        
        
        
        
        
        
        subplot(2,1,2)
        plot(inpaint_nans(pressvel)  , 'LineWidth' , 5 )
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Press/Sec' , 'FontSize' , 20)
        title('The Time Normalized Press Velocity Time series' , 'FontSize' , 20)
        hold on
        ax = gca;
        %         ax.YTick = [1:14];
        %         ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        hold on
        
        
        figure('color' , [1 1 1])
        plot(press , eye , 'LineWidth' , 5 )
        axis square
        grid on
        xlabel('Finger presses' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title('The Time Normalized Eye Position and Finger Press Time series' , 'FontSize' , 20)
        hold on
        line([0 14], [0 14] , 'LineStyle' , ':' , 'LineWidth' , 5 , 'color' , 'r')
        ax = gca;
        ax.XTick = [1:14];
        ax.XLim = [1 14];
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        out = [];
    case 'IPI_dist'
        
        %% IPI distribution in CLA and Random seqs
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        id  = ismember(ANA.seqNumb , [1 2 3 4 5 6]) & ~ANA.isError & ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day});
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        
        A = nanmean(ANA.IPI(id,1:3) , 2);
        B = nanmean(ANA.IPI(id,4:6) , 2);
        C = nanmean(ANA.IPI(id,7:10) , 2);
        D = nanmean(ANA.IPI(id,11:13) , 2);
        figure('color' , [1 1 1])
        subplot(2,2,1)
        plot(A);
        hold on
        plot(B)
        plot(C)
        plot(D)
        
        
        legend({'First quarter' , 'Second quarter' , 'Third quarter' , 'Forth quarter'})
        grid on
        xlabel('Trial')
        ylabel('Normalized time')
        title('Average IPI over all trials in 4 quarters in CLA sequences')
        
        
        subplot(2,2,2)
        histogram(A , 'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        histogram(B , 'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        histogram(C , 'EdgeColor' , 'c' , 'FaceColor' , 'c' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        histogram(D , 'EdgeColor' , 'g' , 'FaceColor' , 'g' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        
        legend({'First quarter' , 'Second quarter' , 'Third quarter' , 'Forth quarter'})
        
        line([nanmean(A) nanmean(A)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(B) nanmean(B)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        line([nanmean(C) nanmean(C)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'c');
        line([nanmean(D) nanmean(D)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'g');
        
        title('Distribution of average IPI over all trials in 4 quarters in CLA sequences')
        
        
        % =========Random
        clear A B C D
        id  = (~ANA.seqNumb & ~ANA.isError & ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day}));
        A = nanmean(ANA.IPI(id,1:3) , 2);
        B = nanmean(ANA.IPI(id,4:6) , 2);
        C = nanmean(ANA.IPI(id,7:10) , 2);
        D = nanmean(ANA.IPI(id,11:13) , 2);
        
        
        subplot(2,2,3)
        plot(A);
        hold on
        plot(B)
        plot(C)
        plot(D)
        legend({'First quarter' , 'Second quarter' , 'Third quarter' , 'Forth quarter'})
        xlabel('Trial')
        ylabel('sec')
        title('Average IPI over all trials in 4 quarters in Random sequences')
        grid on
        
        subplot(2,2,4)
        histogram(A , 'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        histogram(B , 'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        histogram(C , 'EdgeColor' , 'c' , 'FaceColor' , 'c' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        histogram(D , 'EdgeColor' , 'g' , 'FaceColor' , 'g' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        
        legend({'First quarter' , 'Second quarter' , 'Third quarter' , 'Forth quarter'})
        
        line([nanmean(A) nanmean(A)] , [1 100] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(B) nanmean(B)] , [1 100] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        line([nanmean(C) nanmean(C)] , [1 100] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'c');
        line([nanmean(D) nanmean(D)] , [1 100] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'g');
        
        title('Distribution of average IPI over all trials in 4 Random in CLA sequences')
        
        %% Visualize IPIs based on chunk placement
        
        clear A B C D
        A = ANA.IPIChnkPlcmnt(~isnan(ANA.IPIChnkPlcmnt(:,1))  &  ~ANA.isError & ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day}),1);
        B = ANA.IPIChnkPlcmnt(~isnan(ANA.IPIChnkPlcmnt(:,2))  &  ~ANA.isError & ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day}),2);
        C = ANA.IPIChnkPlcmnt(~isnan(ANA.IPIChnkPlcmnt(:,3))  &  ~ANA.isError & ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day}),3);
        D = ANA.IPIChnkPlcmnt(~isnan(ANA.IPIChnkPlcmnt(:,4))  &  ~ANA.isError & ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day}),4);
        
        figure('color' , [1 1 1])
        histogram(A ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        histogram(B, 'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        histogram(C, 'EdgeColor' , 'c' , 'FaceColor' , 'c' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        histogram(D, 'EdgeColor' , 'g' , 'FaceColor' , 'g' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Between IPIs' , '1st chunk IPIs' , '2nd chunk IPIs' , '3rd chunk IPIs'})
        title('Distribution of IPIs per chunk placement')
        line([nanmean(A) nanmean(A)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(B) nanmean(B)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        line([nanmean(C) nanmean(C)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'c');
        line([nanmean(D) nanmean(D)] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'g');
        xlabel('IPI')
        out = [];
        
        %% IPI distances
        clear meanIPI
        for d = 2:4
            for s = 0:6
                A = getrow(ANA , ANA.seqNumb == s & ~ANA.isError & ANA.Day == d & ismember(ANA.Rep , rep));
                meanIPI{d-1 , s+1} = A.IPI;
            end
        end
        
        
        for d= 1:3
            for s = 1:7
                for s1 = 1:7
                    ipi_dist(d,s,s1) = nanmedian(nanmedian(pdist2(meanIPI{d,s} , meanIPI{d,s1} , distance)));
                end
            end
        end
        
        figure('color' , [1 1 1])
        for d = 1:3
            
            subplot(2,3,d)
            imagesc(squeeze(ipi_dist(d,:,:)) , [100 400]);
            title(['Normalized IPI distances - day ' , num2str(d)])
            hold on
            ax = gca;
            ax.FontSize = 20;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.XTickLabelRotation = 45;
            axis square
            colorbar
        end
        
        subplot(2,3,[4:6])
        imagesc(squeeze(mean(ipi_dist , 1)));
        title('all days')
        hold on
        ax = gca;
        ax.FontSize = 20;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        axis square
        colorbar
        out = [];
    case 'IPIs'
        %% within and between IPI ttest with randomization
        %% IPI per day
        nn= input('norm (1) or non-norm (2)?');
        structNumb = [1:6];
        %         plotfcn = input('nanmean or nanmedian?' , 's');
        %% IPIs vs horizon
        % this is the output of the case: 'transitions_All' that is saved to disc
        
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0 , structNumb])  & ~Dall.isError & ismember(Dall.Day , days{day}));
        
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        ANA.seqNumb(ANA.seqNumb>1) = 1;
        for tn  = 1:length(ANA.TN)
            ANA.ChunkBndry(tn , :) = diff(ANA.ChnkArrang(tn,:));
            a = find(ANA.ChunkBndry(tn , :));
            ANA.ChunkBndry(tn , a-1) = 3;
            ANA.ChunkBndry(tn , end) = 3;
            ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
            ANA.ChunkBndry(tn , 1:2) = [-1 -1];  % dont account for the first and last sseqeuce presses
            ANA.ChunkBndry(tn , end-1:end) = [-2 -2];% dont account for the first and last sseqeuce presses
            ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,13);
            ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,13);
            ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,13);
            ANA.IPI_prsnumb(tn , :) = [1 :13];
            ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,13);
        end
        
        newANA.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        newANA.IPI_norm = reshape(ANA.IPI_norm , numel(ANA.IPI) , 1);
        newANA.ChunkBndry = reshape(ANA.ChunkBndry , numel(ANA.IPI) , 1);
        newANA.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
        newANA.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
        newANA.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
        newANA.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
        newANA.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
        newANA.ChunkBndry(newANA.ChunkBndry>2) = 2;
        IPItable  = tapply(newANA , {'Day' ,'SN' , 'seqNumb' , 'ChunkBndry' , 'prsnumb'} , {'IPI' , 'nanmedian(x)'},{'IPI_norm' , 'nanmedian(x)'});
        
        figure('color' , 'white')
        Days = unique(IPItable.Day);
        % gives a scatter plot of the normalized vs non-normalized IPIs in h = 1
        for d = 1:length(Days)
            subplot(1,3,d)
            id = IPItable.seqNumb == 1 & IPItable.ChunkBndry==1 & IPItable.Day==Days(d);
            plot(IPItable.IPI(id) , IPItable.IPI_norm(id) , '*')
            hold on
            id = IPItable.seqNumb == 1 & IPItable.ChunkBndry==2 & IPItable.Day==Days(d);
            plot(IPItable.IPI(id) , IPItable.IPI_norm(id) , '*')
            xlabel('msec')
            ylabel('Norm time')
            legend({'Between Chunk' , 'Within Chunk'})
            title(['Day ' , num2str(Days(d))])
            set(gca , 'FontSize' , 20 , 'XLim' , [50 1200] , 'YLim' , [10 180]);
            grid on
        end
        
        
        if nn == 2
            ylim  = [100 650];
        elseif nn == 1
            IPItable.IPI  = IPItable.IPI_norm;
            ylim  = [70 95];
        end
        
        h1 = figure;
        for h  =1  % spill over from Se2 coding
            for d = 1:5
                [cooIPI_sb(:,h) ,pltIPI_sb(:,h) , errIPI_sb(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 1)  &  ~ismember(IPItable.ChunkBndry, [-1 , -2]) &  ismember(IPItable.ChunkBndry, [1]));
                [cooIPI_sw(:,h) ,pltIPI_sw(:,h) , errIPI_sw(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 1)  &  ~ismember(IPItable.ChunkBndry, [-1 , -2]) &  ismember(IPItable.ChunkBndry, [2]));
                [cooIPI_r(:,h) ,pltIPI_r(:,h) , errIPI_r(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 0)  &  ~ismember(IPItable.ChunkBndry, [-1 , -2]));
            end
        end
        close(h1)
        figure('color', 'white')
        h  =1;
        %             errorbar(cooIPI_sb(:,h) ,pltIPI_sb(:,h) , errIPI_sb(:,h) , 'LineWidth' , 3)
        h1 = plotshade(cooIPI_sb(:,h)' ,pltIPI_sb(:,h)' , errIPI_sb(:,h)','transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 5 , 'linestyle' , ':');
        hold on
%         errorbar(cooIPI_sw(:,h) ,pltIPI_sw(:,h) , errIPI_sw(:,h), 'LineWidth' , 3)
        h2 = plotshade(cooIPI_sw(:,h)' ,pltIPI_sw(:,h)' , errIPI_sw(:,h)','transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 5 , 'linestyle' , ':');
%         errorbar(cooIPI_r(:,h) ,pltIPI_r(:,h) , errIPI_r(:,h), 'LineWidth' , 3)
        h3 = plotshade(cooIPI_r(:,h)' ,pltIPI_r(:,h)' , errIPI_r(:,h)','transp' , .2 , 'patchcolor' , 'g' , 'linecolor' , 'g' , 'linewidth' , 5 , 'linestyle' , ':');
        title(['Inter-Press Intervals'])
        ylabel('msec' )
        xlabel('Days')
        set(gca ,'XTick' , [1:5], 'FontSize' , 40 , 'Box' , 'off' , 'XLim' , [1.5 4.5], 'GridAlpha' , 1)
        grid on
        
        legend([h1 h2 h3] , {'Between Segments' 'Within Segments' 'Random'} , 'Box' , 'off')
        
        
        for d = 2:4
            hcount = 1;
            for h  = 1
                temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.ChunkBndry,{'within/between'},'intercept',1 ,'subset' ,  ismember(IPItable.seqNumb , 1) & ismember(IPItable.Day , d))  ;
                out.chunkEffect_Chunkedseq(hcount , d-1) = temp.eff(2).p;
                temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.ChunkBndry,{'within/between'},'intercept',1 ,'subset' ,  ismember(IPItable.seqNumb , 0) & ismember(IPItable.Day , d))  ;
                out.chunkEffect_Randomseq(hcount , d-1) = temp.eff(2).p;
            end
        end
        hcount = 1;
        for h  = 1
            temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.Day,{'day'},'intercept',1 ,'subset' ,  ismember(IPItable.seqNumb , 1) & ismember(IPItable.ChunkBndry , 1))  ;
            out.dayEffect_bet_Chunkedseq(hcount) = temp.eff(2).p;
            temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.Day,{'day'},'intercept',1 ,'subset' ,  ismember(IPItable.seqNumb , 1) & ismember(IPItable.ChunkBndry , 2))  ;
            out.dayEffect_wit_Chunkedseq(hcount) = temp.eff(2).p;
        end
    case 'MTs'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0:6]) & ~Dall.isError);
        ANA.seqNumb(ANA.seqNumb>1) = 1;
        MT  = tapply(ANA , {'Day' ,'SN' , 'seqNumb','BN'} , {'MT' , 'nanmedian(x)'});
        
        MT = getrow(MT , MT.MT <= 9000 );
        
        for d=  2:4
            temp = anovaMixed(MT.MT  , MT.SN ,'within',MT.seqNumb ,{'Random/Chunked'},'intercept',1 ,'subset' , ismember(MT.Day , d))  ;
            out.ChVsRand(d-1)  = temp.eff(2).p;
        end
        h1 = figure('color' , 'white');
        [xcoords,PLOTs,ERRORs] = lineplot(MT.Day, MT.MT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 1 );
        hold on
        [xcoordr,PLOTr,ERRORr] = lineplot(MT.Day, MT.MT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 0 & ~ismember(MT.Day , 1));
        close(h1);
        
        
        
        figure('color' , 'white');
        
        h1 = plotshade(xcoords',PLOTs , ERRORs,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 5 , 'linestyle' , ':');
        hold on
        h2 = plotshade(xcoordr',PLOTr , ERRORr,'transp' , .2 , 'patchcolor' , 'g' , 'linecolor' , 'g' , 'linewidth' , 5 , 'linestyle' , ':');
        set(gca,'FontSize' , 40 , 'XLim' , [1 5] , 'XTick' , [2:4] , 'Box' , 'off' , 'GridAlpha' , 1);
        title('Execution Time')
        ylabel('msec' )
        xlabel('Days' )
        legend([h1 h2] ,{'Structured Sequences' , 'Random Sequences'})
        grid on
    case 'Errors'
        errorrate = [];
        for d = 2:4
            for sub = subjnum
                ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [1:6]) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{d}));
                errorrate = [errorrate ; [100 * sum(ANA.isError)/length(ANA.TN) , d , sub ,1]];
                ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [0]) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{d}));
                errorrate = [errorrate ; [100 * sum(ANA.isError)/length(ANA.TN) , d , sub ,0]];
            end
        end
        
        out.Err = anovaMixed(errorrate(:,1) , errorrate(:,3),'within',errorrate(:,4),{'Random/Chunked'},'intercept',1)  ;
        h1 = figure('color' , 'white');
        [xcoord,PLOTs,ERRORs] = lineplot(errorrate(:,2),errorrate(:,1) , 'subset' , errorrate(:,4) == 1);
        hold on
        [xcoord,PLOTr,ERRORr] = lineplot(errorrate(:,2),errorrate(:,1) , 'subset' , errorrate(:,4) == 0);
        close(h1);
        
        figure('color' , 'white')
%         errorbar(PLOTs , ERRORs, 'LineWidth' , 3);
%         hold on
%         errorbar(PLOTr , ERRORr, 'LineWidth' , 3);
        h1 = plotshade(xcoord',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 5 , 'linestyle' , ':');
        hold on
        h2 = plotshade(xcoord',PLOTr,ERRORr,'transp' , .2 , 'patchcolor' , 'g' , 'linecolor' , 'g' , 'linewidth' , 5 , 'linestyle' , ':');
        set(gca,'FontSize' , 40 , 'XLim' , [1 5] , 'XTick' , [2:4] , 'Box' , 'off' , 'GridAlpha' , 1);
        grid on
        ylabel('Percent')
        title(['Error rate'])
        legend([h1 , h2] , {'Structured Sequences' , 'Random Sequences'});
    case 'chunk_dist' % CityBlock
        %% chunk distances
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}));
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            IPI(tn , :) = detrend(diff(nIdx(tn ,:) , 1 , 2) , 'linear' , 2);
        end
        
        for tn = 1:length(ANA.TN)
            thresh = .3 * std(IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if IPI(tn , estChnkBndry(cb)) < nanmean(IPI(tn  ,:)) + thresh
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
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
        end
        
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
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
        end
        
        
        
        %% the boundry distance
        temp = diff(ChnkPlcmnt,1,2);
        temp(temp<0) = 0;
        chbndry = [ones(6,1) ~temp]; % all the first presses are one
        
        clear meanbnd allmeanbnd CB est_dist act_est_dist
        
        for s = 0:6
            A = getrow(ANA , ANA.seqNumb == s);
            CBD{s+1} = A.estChnkBndry;
        end
        
        for s = 1:7
            for s1 = 1:7
                est_dist(s,s1) = nanmean(nanmean(pdist2(CBD{s} , CBD{s1} , distance)));
            end
            for s2 = 1:6
                act_est_dist(s,s2) = nanmean(nanmean(pdist2(CBD{s} , chbndry(s2 , :) , distance)));
            end
        end
        
        figure('color' , [1 1 1])
        subplot(1,3,1)
        
        imagesc(est_dist );
        title('Estimated Chunking Structure Dissimilarity Matrix' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        hold on
        line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        axis square
        colorbar
        
        subplot(1,3,2)
        
        imagesc(act_est_dist);
        title('Estimated vs. Expected Chunking Structure Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        line([.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        ax.XTick = [1:6];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        xlabel('Expected')
        ylabel('Estimated')
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        
        subplot(1,3,3)
        imagesc(squareform(pdist(chbndry)));
        title('Chunking Structure Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:6];
        ax.YTick = [1:6];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        out = [];
    case 'Mychunk_dist'
        %% chunk distances
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        ANA = getrow(ANA ,  ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day}));
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            IPI(tn , :) = detrend(diff(nIdx(tn ,:) , 1 , 2) , 'linear' , 2);
        end
        
        for tn = 1:length(ANA.TN)
            thresh = .3 * std(IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if IPI(tn , estChnkBndry(cb)) < nanmean(IPI(tn  ,:)) + thresh
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
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
        end
        
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
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
        end
        
        
        for d = 2:4
            for s = 0:6
                A = getrow(ANA , ANA.seqNumb == s & ~ANA.isError & ismember(ANA.Day , d) & ismember(ANA.Rep , rep) & ANA.isgood);
                CPD{d-1 , s+1} = A.estChnkPlcmnt;
            end
        end
        for d= 1:3
            for s = 1:7
                for s1 = 1:7
                    cp_est_dist(d,s,s1) = nanmean(nanmean(pdist2(CPD{d,s} , CPD{d,s1} , distance)));
                end
                for s2 = 1:6
                    act_cp_est_dist(d,s,s2) = nanmean(nanmean(pdist2(CPD{d,s} , ChnkPlcmnt(s2 , :) , distance)));
                end
            end
        end
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        
        imagesc(squeeze(mean(cp_est_dist , 1)));
        title('Estimated Chunking Placement Dissimilarity Matrix' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        hold on
        line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        axis square
        colorbar
        
        subplot(1,2,2)
        imagesc(squareform(pdist(ChnkPlcmnt)));
        title('Chunking Placement Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:6];
        ax.YTick = [1:6];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
        ax.YTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        out = [];
        
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        
        imagesc(squeeze(mean(act_cp_est_dist , 1)));
        title('Estimated vs. Expected Chunking Placement Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        line([.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        ax.XTick = [1:6];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        xlabel('Expected')
        ylabel('Estimated')
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        
        subplot(1,2,2)
        imagesc(squareform(pdist(ChnkPlcmnt)));
        title('Chunking Placement Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:6];
        ax.YTick = [1:6];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        
        
        out = [];
    case 'Avg_pattern_sh3'
        D = load('/Users/nedakordjazi/Documents/SeqEye/SequenceHierarchical/Analysis/sh3_avrgPattern.mat');
        MT = D.MT(126:end);
        sq = D.Sequence(126:end , 1:2);
        for i = 1:length(sq)
            AvgMT(sq(i,1),sq(i,2)) = D.MT(i);
        end
        figure('color' , [1 1 1])
        imagesc(AvgMT)
        axis square
        hold on
        ax = gca;
        ax.XTick = [1:5];
        ax.YTick = [1:5];
        colorbar
        title('Double Average Patterns');
        xlabel('First press')
        ylabel('Second press')
        min(MT)/max(MT);
        
        out = [];
    case 'eye_vel_chunkplace'
        % eye velocity as a function of chunk placement
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & ~Dall.isError & Dall.isgood & ismember(Dall.Rep ,rep) & ismember(Dall.Day , days{day}));
        figure('color' , [1 1 1])
        subplot(1,2,1)
        col = {'b' , 'r' ,'c' ,'g'};
        for cp = 1:4
            A{cp} = ANA.xEyeVelCnkPlcmnt(~isnan(ANA.xEyeVelCnkPlcmnt(:,cp)),cp);
            %histogram(A{cp} ,'EdgeColor' , col{cp} , 'FaceColor' , col{cp} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            plot(cp*ones(size(A{cp})) ,A{cp} , '*' , 'color' , col{cp})
            hold on
        end
        legend({'First chunk place' , 'Second chunk place' , 'Third chunk place' , 'Forth chunk place'})
        %         for cp = 1:4
        %             line([nanmean(A{cp}) nanmean(A{cp})] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , col{cp});
        %         end
        xlabel('degree/s')
        grid on
        
        
        All = NaN*ones(max([length(A{1}) , length(A{2}) , length(A{3}) , length(A{4})]) , 4);
        for cp  = 1:4
            All(1:length(A{cp}) , cp) = A{cp};
        end
        
        subplot(1,2,2)
        hold on
        boxplot(All)
        title('Eye velocity per chunk placement')
        xlabel('Chunk place')
        
        
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        col = {'b' , 'r' ,'c' ,'g'};
        for cp = 1:4
            A{cp} = ANA.xEyeVelEstCnkPlcmnt(~isnan(ANA.xEyeVelEstCnkPlcmnt(:,cp)),cp);
            histogram(A{cp} ,'EdgeColor' , col{cp} , 'FaceColor' , col{cp} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            hold on
        end
        legend({'First chunk place' , 'Second chunk place' , 'Third chunk place' , 'Forth chunk place'})
        for cp = 1:4
            line([nanmean(A{cp}) nanmean(A{cp})] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , col{cp});
        end
        xlabel('degree/s')
        grid on
        
        
        All = NaN*ones(max([length(A{1}) , length(A{2}) , length(A{3}) , length(A{4})]) , 4);
        for cp  = 1:4
            All(1:length(A{cp}) , cp) = A{cp};
        end
        
        subplot(1,2,2)
        hold on
        boxplot(All)
        title('Eye Velocity Per Estimated Chunk Placement')
        xlabel('Chunk place')
        out =[];
    case 'eye_vel_seqplace'
        %  velocity per sequence quarter --------CLA
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        clear A
        
        id  = ismember(ANA.seqNumb , [1:6]) & ~ANA.isError;
        A{1} = nanmean(ANA.xEyePressAngVelocity(id,1:4) , 2);
        A{2} = nanmean(ANA.xEyePressAngVelocity(id,5:7) , 2);
        A{3} = nanmean(ANA.xEyePressAngVelocity(id,8:11) , 2);
        A{4} = nanmean(ANA.xEyePressAngVelocity(id,12:14) , 2);
        
        out.seqQuarter_vel = A;
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        
        col = {'b' , 'r' ,'c' ,'g'};
        for cp = 1:4
            histogram(A{cp} ,'EdgeColor' , col{cp} , 'FaceColor' , col{cp} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            hold on
        end
        
        legend({'First quarter' , 'Second quarter' , 'Third quarter' , 'Forth quarter'})
        for cp = 1:4
            line([nanmean(A{cp}) nanmean(A{cp})] , [0 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , col{cp});
        end
        grid on
        xlabel('deg/sec')
        title('Average eye velocity over all trials in 4 quarters in CLA sequences')
        
        
        All = NaN*ones(max([length(A{1}) , length(A{2}) , length(A{3}) , length(A{4})]) , 4);
        for sq  = 1:4
            All(1:length(A{sq}) , sq) = A{sq};
        end
        
        subplot(1,2,2)
        hold on
        boxplot(All)
        title('Eye velocity per sequence quarter -- CLA sequences')
        xlabel('Seq Quarter')
    case 'eye_pos_seqplace'
        %  velocity per sequence quarter --------Random
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        clear A
        ANAr   = getrow(ANA , ismember(ANA.seqNumb , [0]) & ~ANA.isError & ismember(ANA.Day , [1:4]) & ismember(ANA.Rep , rep));
        
        for seqnumb = 1:6
            subplot(3,3,seqnumb)
            ANAcla = getrow(ANA , ismember(ANA.seqNumb , [seqnumb]) & ~ANA.isError & ismember(ANA.Rep , rep));
            subplot(2,1,1)
            boxplot(ANAcla.EyePressTimePos);
            title(['Median eye position in a 50 ms vicinity of the press times - CLA ' , num2str(seqnumb)])
            ylabel('cm')
            ylim([1 14])
            grid on
        end
        
        subplot(3,3,[7:9])
        boxplot(ANAr.EyePressTimePos);
        title('Median eye position in a 50 ms vicinity of the press times - Random')
        ylabel('cm')
        grid on
        ylim([1 14])
    case 'eyepress_pos_traces'
        clear len EyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1:length(subj_name)-1])
            subjnum = 13;
        end
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.EyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            
            
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            
            out.EyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
            
        end
        
        figure('color' , [1 1 1])
        for s = 2:7
            plot(out.pressPattern(s,:) , out.EyePattern(s,:) , 'LineWidth' , 3);
            hold on
        end
        plot(out.pressPattern(1,:) , out.EyePattern(1,:)  , 'LineWidth' , 5 , 'color' , 'black');
        grid on
        xlabel('Finger Press' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average Time Normalized Eye Position vs. Finger Press Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.XTick = [1:14];
        ax.XLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        axis square
        
        
        figure('color' , [1 1 1])
        subplot(2,1,2)
        plot(out.pressPattern(2:end , :)' , 'LineWidth' , 3)
        grid on
        hold on
        plot(out.pressPattern(1 , :)' , 'LineWidth' , 5 , 'color' , 'black')
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Press Position' , 'FontSize' , 20)
        title(['The Average Time Normalized Finger Press Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        %         legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
        
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        
        subplot(2,1,1)
        plot(out.EyePattern(2:end , :)' , 'LineWidth' , 3)
        hold on
        plot(out.EyePattern(1 , :)' , 'LineWidth' , 5 , 'color' , 'black')
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average Time Normalized Eye Position Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        %         legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
        
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 max(max(out.EyePattern))];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        plot(out.EyePattern(2:end , 901:end)' , 'LineWidth' , 3)
        hold on
        plot(out.EyePattern(1 , 901:end)' , 'LineWidth' , 5 , 'color' , 'black')
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average Last 10% Eye Position Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        ax = gca;
        ax.FontSize = 20;
        ax.Box = 'off';
        ax.XLim = [1 100];
        ax.YLim = [min(min(out.EyePattern(: , 901:end))) 1.5+min(min(out.EyePattern(: , 901:end))) ]
        
        subplot(2,1,2)
        plot(out.EyePattern(2:end , 1:100)' , 'LineWidth' , 3)
        hold on
        plot(out.EyePattern(1 , 1:100)' , 'LineWidth' , 5 , 'color' , 'black')
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average First 10% Eye Position Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
        ax = gca;
        ax.XLim = [1 100];
        ax.YLim = [min(min(out.EyePattern(: , 1:100))) 1.5+min(min(out.EyePattern(: , 1:100)))];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        out = [];
    case 'eyepress_pos_distances'
        clear len EyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        %         calc = 0;
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10] , 'rows')
            subjnum = 11;
        end
        
        titleAdd = {'First Third' , 'Middle Third' , 'Last Third' , 'The Entire Time Series'};
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.EyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            % Use the random sequences in day 1 as the baseline. On day one we only have random sequences and  no CLAs
            
            
            
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
            
        end
        
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.prsDist{t}(s,s1) = nanmean(nanmean(pdist2(out.press{s}(:,tid{t}) , out.press{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.prsDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7]
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.XTickLabelRotation = 45;
            axis square
            title(['Time-Normalized Press Position Time-Series Distances -  ' , titleAdd{t} , ' - Subject ' , num2str(subjnum)])
        end
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.eyeDist{t}(s,s1) = nanmean(nanmean(pdist2(out.eye{s}(:,tid{t}) , out.eye{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.eyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7]
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.XTickLabelRotation = 45;
            axis square
            title(['Time-Normalized Eye Position Time-Series Distances -  ' , titleAdd{t} , ' - Subject ' , num2str(subjnum)])
        end
    case 'eyepress_pos_avg_distances'
        
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10] , 'rows')
            subjnum = 11;
        end
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
        end
        
        titleAdd = {'First Third' , 'Middle Third' , 'Last Third' , 'The Entire Time Series'};
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            out.prseyeDist{t} = pdist2(out.pressPattern(:,tid{t}) , out.eyePattern(:,tid{t}) , distance);
            subplot(1,length(tid),t)
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            ax.FontSize = 20;
            title(['Finger Press/Eye Position Dissimilarity Matrix - ' ,titleAdd{t}] , 'FontSize' , 20)
        end
        
        figure('color' , [1 1 1])
        for t = length(tid)
            out.prseyeDist{t} = pdist2(out.pressPattern(:,tid{t}) , out.eyePattern(:,tid{t}) , distance);
            subplot(1,length(tid),t)
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            ax.FontSize = 20;
            title('Finger Press/Eye Position Dissimilarity Matrix' , 'FontSize' , 20)
        end
    case 'eyepress_vel_traces'
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
            subjnum = 11;
        end
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            out.eyePattern  = N.norm.eyePattern;
            out.pressPattern = N.norm.pressPattern;
            out.eyevelocPattern = N.norm.eyevelocPattern;
            out.prsvelocPattern = N.norm.prsvelocPattern;
            
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
        end
        
        figure('color' , [1 1 1])
        
        hold on
        for s = 2:7
            plot(smooth(out.prsvelocPattern(s,:) , 10) ,smooth(out.eyevelocPattern(s,:) ,10), 'LineWidth' , 1);
            hold on
        end
        plot(smooth(out.prsvelocPattern(1,:) ,10),smooth(out.eyevelocPattern(1,:) ,10), 'LineWidth' , 3);
        legend({'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'})
        axis square
        xlabel('Press/sec')
        ylabel('deg/sec')
        title(['Press Velocity vs. Eye Position Time Series, Subject ' ,num2str(subjnum)])
        ax = gca;
        ax.FontSize = 20;
        grid on
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        for m  = 2:7
            plot(smooth(out.prsvelocPattern(m,:) ,10) , 'LineWidth' , 1);
            hold on
        end
        plot(smooth(out.prsvelocPattern(1,:)' ,10)   , 'LineWidth' , 3);
        xlabel('Normalized Time')
        ylabel('Press/sec')
        title(['Finger Press Velocity, Subject  ' , num2str(subjnum)])
        grid on
        ax = gca;
        ax.FontSize = 20;
        legend({'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random' })
        
        
        subplot(2,1,2)
        for m = 2:7
            plot(smooth(out.eyevelocPattern(m , :) , 10), 'LineWidth' , 1)
            hold on
        end
        plot(smooth(out.eyevelocPattern(1, :)' , 10), 'LineWidth' , 3)
        xlabel('Normalized Time')
        ylabel('deg/sec')
        title(['Eye Angular Velocity, Subject ' , num2str(subjnum)])
        grid on
        ax = gca;
        ax.FontSize = 20;
        legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'})
    case 'eyepress_vel_distances'
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
            subjnum = 11;
        end
        
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
        end
        
        eyeid = [1 2 3 4];
        prsid = [5 6 7 8];
        titleAdd = {'First Quarter' , 'Second Quarter' , 'Third Quarter' , 'Fourth Quarter' , 'The Entire Time Series'};
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.prsDist{t}(s,s1) = nanmedian(nanmedian(pdist2(out.prsveloc{s}(:,tid{t}),out.prsveloc{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.prsDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title(['Dissimilarity Between Average Time-Normalized Eye Angular Velocities - ' , titleAdd{t}])
        end
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.eyeDist{t}(s,s1) = nanmedian(nanmedian(pdist2(out.eyeveloc{s}(:,tid{t}),out.eyeveloc{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.eyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title(['Dissimilarity Between Average Time-Normalized Eye Angular Velocities - ' , titleAdd{t}])
        end
    case 'eyepress_vel_avg_distances'
        %%
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
            subjnum = 11;
        end
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye_bl{1}      = 1 + N.norm(1,subjnum).eye{1};
            out.press_bl{1}    = N.norm(1,subjnum).press{1};
            out.eyeveloc_bl{1} = N.norm(1,subjnum).eyeveloc{1};
            out.prsveloc_bl{1} = N.norm(1,subjnum).prsveloc{1};
            
            
            out.eyePattern_bl      = nanmean(out.eye_bl{1});
            out.pressPattern_bl    = nanmean(out.press_bl{1});
            out.eyevelocPattern_bl = nanmean(out.eyeveloc_bl{1});
            out.prsvelocPattern_bl = nanmean(out.prsveloc_bl{1});
            
            
            out.eye = 1 + N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        end
        
        
        eyeid = [1 2 3 4];
        prsid = [5 6 7 8];
        
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            out.prseyeDist{t} = pdist2(out.prsvelocPattern(2:7,tid{t}) , out.eyevelocPattern(2:7,tid{t}) , distance);
            subplot(1,length(tid)-1,t)
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title(['Dissimilarity between time-normalized press and eye velocity - quarter ' , num2str(t)])
        end
        
        figure('color' , [1 1 1])
        for t = length(tid)
            out.prseyeDist{t} = pdist2(out.prsvelocPattern(2:7,tid{t}) , out.eyevelocPattern(2:7,tid{t}) , distance);
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title('Dissimilarity between time-normalized press and eye velocity')
        end
    case 'run_f-tests_mt'
        % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\movement time of Structres vs. Random/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        %//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
        %\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
        ANA1 = getrow(Dall , ismember(Dall.seqNumb,[0:6]) &  ~Dall.isError & Dall.isgood);
        ANA1.tempSeqNum = ANA1.seqNumb;
        ANA1.tempSeqNum(ismember(ANA1.tempSeqNum , [1:6])) = 1;
        ANA1.mt = ANA1.AllPressTimes(:,14) - ANA1.AllPressTimes(:,1);
        
        id1 = ANA1.tempSeqNum == 1;
        id2 = ANA1.tempSeqNum == 0;
        out.CLAvsRand_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.tempSeqNum],{'Sequence type'},'subset',~isnan(ANA1.mt),'intercept',1)  ;
        out.CLAvsRand_rt = anovaMixed(ANA1.AllPressTimes(:,1) , ANA1.SN,'within',[ANA1.tempSeqNum],{'Sequence type'},'subset',~isnan(ANA1.AllPressTimes(:,1)),'intercept',1)  ;
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Structured sequences' , 'Random sequences'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 300] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 300] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of movement time, p = ' , num2str( out.CLAvsRand_mt.eff(2).p)])
        grid on
        
        subplot(1,2,2)
        histogram(ANA1.AllPressTimes(id1,1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        id2 = ANA1.tempSeqNum == 0;
        histogram(ANA1.AllPressTimes(id2,1) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Structured sequences' , 'Random sequences'})
        
        line([nanmean(ANA1.AllPressTimes(id1,1)) nanmean(ANA1.AllPressTimes(id1,1))] , [1 600] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.AllPressTimes(id2,1)) nanmean(ANA1.AllPressTimes(id2,1))] , [1 600] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of reaction time, p = ' , num2str(out.CLAvsRand_mt.eff(2).p)])
        grid on
        
        
        out.StrucvsStruct = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.seqNumb],{'Sertucture number'},'subset',ismember(ANA1.seqNumb ,[1:6]) & ~isnan(ANA1.mt),'intercept',1);
        color = {'b' , 'r' ,'g' ,'c' ,'y' ,'black'};
        figure
        for s= 1:6
            histogram(ANA1.mt(ANA1.seqNumb == s) ,'EdgeColor' , color{s} , 'FaceColor' , color{s} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            hold on
        end
        legend({'Structured 1' , 'Structured 2' ,'Structured 3' , 'Structured 4' , 'Structured 5' , 'Structured 6'})
        for s= 1:6
            line([nanmean(ANA1.mt(ANA1.seqNumb == s)) nanmean(ANA1.mt(ANA1.seqNumb == s))] , [0 150] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' ,  color{s});
            hold on
        end
        grid on
        title(['Distribution of movement time, p = ' , num2str( out.StrucvsStruct.eff(2).p)])
        
        
        
        
        % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\movement time of repetition 1 vs repetition 2/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        %//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
        %\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
        
        id1 = ANA1.Rep == 1 & ANA1.tempSeqNum == 1;
        id2 = ANA1.Rep == 2 & ANA1.tempSeqNum == 1;
        out.S_R1vsR2_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.Rep],{'Repetition'},'subset',~isnan(ANA1.mt) & ANA1.tempSeqNum == 1,'intercept',1);
        out.R_R1vsR2_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.Rep],{'Repetition'},'subset',~isnan(ANA1.mt) & ANA1.tempSeqNum == 0,'intercept',1);
        
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Repetition 1' , 'Repetition 2'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of movement time over across all structured sequences, p = ' , num2str(out.S_R1vsR2_mt.eff(2).p)])
        grid on
        
        
        id1 = ANA1.Rep == 1 & ANA1.tempSeqNum == 0;
        id2 = ANA1.Rep == 2 & ANA1.tempSeqNum == 0;
        subplot(1,2,2)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        id2 = ANA1.tempSeqNum == 0;
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Repetition 1' , 'Repetition 2'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of reaction time over across all random sequences, p = ' , num2str(out.R_R1vsR2_mt.eff(2).p)])
        grid on
        
        % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\only structure blocks vs intermixed/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        %//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
        %\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
        
        % intermixed --> BT = 2
        % CLAT --> BT = 3
        id1 = ANA1.BT == 2;
        id2 = ANA1.BT == 3;
        out.SvsI_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.BT],{'Block Type'},'subset',~isnan(ANA1.mt) & ANA1.tempSeqNum == 1,'intercept',1);
        out.SvsI_rt = anovaMixed(ANA1.AllPressTimes(:,1) , ANA1.SN,'within',[ANA1.BT],{'Block Type'},'subset',~isnan(ANA1.AllPressTimes(:,1)) & ANA1.tempSeqNum == 1,'intercept',1);
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Intermixed block' , 'Structured block'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of movement time over across all structured trails, p = ' , num2str(out.SvsI_mt.eff(2).p)])
        grid on
        
        subplot(1,2,2)
        histogram(ANA1.AllPressTimes(id1,1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        id2 = ANA1.tempSeqNum == 0;
        histogram(ANA1.AllPressTimes(id2,1) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Intermixed block' , 'Structured block'})
        
        line([nanmean(ANA1.AllPressTimes(id1,1)) nanmean(ANA1.AllPressTimes(id1,1))] , [1 500] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.AllPressTimes(id2,1)) nanmean(ANA1.AllPressTimes(id2,1))] , [1 500] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of reaction time over across all structured trails, p = ' , num2str(out.SvsI_rt.eff(2).p)])
        grid on
    case 'dtw'
        
        %% Dynamic time warping
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        ANA1 = getrow(ANA , ismember(ANA.seqNumb , [0:6]) & ~ANA.isError);
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye_bl{1}      = 1 + N.norm(1,subjnum).eye{1};
            out.press_bl{1}    = N.norm(1,subjnum).press{1};
            out.eyeveloc_bl{1} = N.norm(1,subjnum).eyeveloc{1};
            out.prsveloc_bl{1} = N.norm(1,subjnum).prsveloc{1};
            
            
            out.eyePattern_bl      = nanmean(out.eye_bl{1});
            out.pressPattern_bl    = nanmean(out.press_bl{1});
            out.eyevelocPattern_bl = nanmean(out.eyeveloc_bl{1});
            out.prsvelocPattern_bl = nanmean(out.prsveloc_bl{1});
            
            
            out.eye = 1 + N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        end
        
        for tn = 1:length(ANA1.TN)
            if ANA1.AllPressIdx(tn,14) < length(ANA1.xEyePosDigit{tn})
                ANA1.xEyePosDigit{tn}          = ANA1.xEyePosDigit{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
                ANA1.PressTimeSeries{tn}    = ANA1.PressTimeSeries{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
                %                 ANA1.ChunkTimeSeries{tn}    = ANA1.ChunkTimeSeries{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
                %                 ANA1.estChunkTimeSeries{tn} = ANA1.estChunkTimeSeries{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
            end
        end
        
        dxeye  = zeros(length(ANA1.TN) , length(ANA1.TN));
        dpress = zeros(length(ANA1.TN) , length(ANA1.TN));
        figure('color' , [1 1 1])
        for tn1 = 1:length(ANA1.TN)
            for tn2 = tn1+1:length(ANA1.TN)
                [dxeye(tn1 , tn2) , ~] = dtw(ANA1.xEyePosDigit{tn1},ANA1.xEyePosDigit{tn2});
                dxeye(tn2 , tn1)  = dxeye(tn1 , tn2);
                
                [dpress(tn1 , tn2) , ~] = dtw(ANA1.PressTimeSeries{tn1},ANA1.PressTimeSeries{tn2});
                dpress(tn2 , tn1)  = dpress(tn1 , tn2);
                [tn1 tn2]
                %                 imagesc(dxeye)
                %                 drawnow()
            end
            
        end
        out.seqNumb = ANA1.seqNumb;
        out.dxeye = dxeye;
        out.dpress = dpress;
    case 'eyeEndPos_ChunkLength'
        %%
        ANA = getrow(Dall , ismember(Dall.seqNumb , [1:6]) & Dall.Group == 1 & ismember(Dall.Rep , rep) & ~Dall.isError & Dall.isgood & ismember(Dall.Day , days{day}));
        ANA.EndChunkLength = NaN*ones(length(ANA.TN) , 1);
        ANA.EyeLastPos = NaN*ones(length(ANA.TN) , 1);
        ANA.numChunks = NaN*ones(length(ANA.TN) , 1);
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
        end
        IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
        IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        id = find(ismember(ANA.seqNumb , [1:6]));
        
        for tn = 1:length(id)
            ANA.EndChunkLength(id(tn)) = ChnkPlcmnt(ANA.seqNumb(id(tn)) , end);
            ANA.EyeLastPos(id(tn) ,1) = ANA.EyePressTimePos(id(tn) , end);
            ANA.numChunks(id(tn) ,1) = length(find(ANA.ChnkPlcmnt(tn,:) == 1));
        end
        
        out.endEffect1 = anovaMixed(ANA.EyeLastPos , ANA.SN,'within',[ANA.EndChunkLength],{'Last Chunk Length'},'subset',~isnan(ANA.EyeLastPos),'intercept',1)  ;
        out.endEffect2 = anovaMixed(ANA.EyeLastPos , ANA.SN,'within',[ANA.numChunks],{'Number of Chunks'},'subset',~isnan(ANA.EyeLastPos),'intercept',1)  ;
    case 'eyeFrstPos_ChunkLength'
        ANA = getrow(Dall , Dall.Group == 1 & ismember(Dall.Rep , rep) & ~Dall.isError & Dall.isgood & ismember(Dall.Day , days{day}));
        
        ANA.FrstChunkLength = NaN*ones(length(ANA.TN) , 1);
        ANA.EyeFrstPos = NaN*ones(length(ANA.TN) , 1);
        ANA.numChunks = NaN*ones(length(ANA.TN) , 1);
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
        end
        IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
        IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        id = find(ismember(ANA.seqNumb , [1:6]));
        for tn = 1:length(id)
            [a,b] = unique(ChnkArrng(ANA.seqNumb(id(tn)),:));
            ANA.FrstChunkLength(id(tn)) = b(2)-1;
            ANA.EyeFrstPos(id(tn) ,1) = ANA.EyePressTimePos(id(tn) , 1);
            ANA.numChunks(id(tn) ,1) = length(find(ANA.ChnkPlcmnt(tn,:) == 1));
        end
        A1 = ANA;
        
        
        out.strtEffect1 = anovaMixed(ANA.EyeFrstPos , ANA.SN,'within',[ANA.FrstChunkLength],{'Last Chunk Length'},'subset',~isnan(ANA.EyeFrstPos),'intercept',1)  ;
        out.strtEffect2 = anovaMixed(ANA.EyeFrstPos , ANA.SN,'within',[ANA.numChunks],{'Last Chunk Length'},'subset',~isnan(ANA.EyeFrstPos),'intercept',1)  ;
    case 'eye10p_crossvaldist_pos'
        %         %         if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
        %             subjnum = 11;
        %         end
        calcDist = 1;
        if calcDist
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [1:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).first10p_eye= N.norm(d,subjnum(SubN)).eye;
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).first10p_eye{seqs(s)}(:,1:100)];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; seqs(s)*ones(size(out(SubN, d-1).first10p_eye{seqs(s)}(:,1:100) , 1) ,1)];
                    end
                    lab(lab>1) = 2;
                    seqs = unique(lab);
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).first10_D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eyfirst10_D_eyee(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                first10_D_eyeSig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,f10_ePLOTw(SubN , :),f10_eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,f10_ePLOTb(SubN , :),f10_eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,f10_ePLOTw(SubN+1 , :),f10_eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,f10_ePLOTb(SubN+1 , :),f10_eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            seqs = [1:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).last10p_eye= N.norm(d,subjnum(SubN)).eye;
                    
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).last10p_eye{seqs(s)}(:,901:1000)];
                        lab = [lab ; seqs(s)*ones(size(out(SubN, d-1).last10p_eye{seqs(s)}(:,1:100) , 1) ,1)];
                    end
                    lab(lab>1) = 2;
                    seqs = unique(lab);
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).last10_D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eylast10_D_eyee(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                last10_D_eyeSig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,l10_ePLOTw(SubN , :),l10_eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,l10_ePLOTb(SubN , :),l10_eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,l10_ePLOTw(SubN+1 , :),l10_eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,l10_ePLOTb(SubN+1 , :),l10_eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        else
            seqs = [2:7];
            load([baseDir , '/eyepos_dbetclass.mat'])
            load([baseDir , '/eyepos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_eyepos_dbetclass.mat'])
                load([baseDir , '/EUC_eyepos_dwithclass.mat'])
            end
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prspos_dbetclass.mat'])
            load([baseDir , '/prspos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_prspos_dbetclass.mat'])
                load([baseDir , '/EUC_prspos_dwithclass.mat'])
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
        end
        
        
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(f10_ePLOTw(SubN , :),f10_eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(f10_ePLOTb(SubN , :),f10_eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['F10 Eye Position, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
            %             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(f10_ePLOTw(SubN+1 , :),f10_eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(f10_ePLOTb(SubN+1 , :),f10_eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['f10 Subject Average Dissimilarity in Eye position,  p =' , num2str(first10_D_eyeSig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(l10_ePLOTw(SubN , :),l10_eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(l10_ePLOTb(SubN , :),l10_eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['l10 eye Position, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
            %             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(l10_ePLOTw(SubN+1 , :),l10_eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(l10_ePLOTb(SubN+1 , :),l10_eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['l10 eye Velocity,  p =' , num2str(last10_D_eyeSig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        out = [];
    case 'crossvaldist_pos'
        %         %         if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
        %             subjnum = 11;
        %         end
        calcDist = 1;
        if calcDist
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [2:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).eye= N.norm(d,subjnum(SubN)).eye;
                    
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).eye{seqs(s)}];
                        lab = [lab ; s*ones(size(out(SubN, d-1).eye{seqs(s)} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                eyesig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).press= N.norm(d,subjnum(SubN)).press;
                    
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(out(SubN, d-1).press)
                        X = [X ; out(SubN, d-1).press{s}];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; s*ones(size(out(SubN, d-1).press{s} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                prssig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        else
            seqs = [2:7];
            load([baseDir , '/eyepos_dbetclass.mat'])
            load([baseDir , '/eyepos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_eyepos_dbetclass.mat'])
                load([baseDir , '/EUC_eyepos_dwithclass.mat'])
            end
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prspos_dbetclass.mat'])
            load([baseDir , '/prspos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_prspos_dbetclass.mat'])
                load([baseDir , '/EUC_prspos_dwithclass.mat'])
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            
        end
        
        
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(ePLOTw(SubN , :),eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(ePLOTb(SubN , :),eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Eye Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
            %             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Eye Velocity,  p =' , num2str(eyesig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(pPLOTw(SubN , :),pERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(pPLOTb(SubN , :),pERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Press Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
            %             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Press Velocity,  p =' , num2str(prssig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        out = [];
    case 'crossval_IPI_dist'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0:6]) & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}));
        % time normalizing the IPIs (sp the press indecies) to 1 : 1000 normalized samples
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        X = ANA.IPI;clear d
        
        lab = ANA.seqNumb +1;
        
        for i = 1:size(X,1)
            id = ones(size(X , 1) , 1);
            id(i) = 0;
            Y = inpaint_nans(X(~id , :));
            X1 = X(id==1 , :);
            lab1 = lab(id==1);
            for s = 1:7
                m = nanmean(X1(lab1==s , :));
                d(i , s) = pdist([Y;m], distance);%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
            end
            
        end
        seqs =  unique(lab);
        for l =  1:length(seqs)
            for l1 =  1:length(seqs)
                id = lab == l;s
                out.D_IPI(l,l1) = nanmean(d(id , l1));
            end
        end
        
        figure('color' , [1 1 1])
        
        imagesc(out.D_IPI);
        title('Crossvalidated IPI Dissimilarity Matrix' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        hold on
        line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        axis square
        colorbar
        
        count = 1;
        sn = unique(ANA.SN);
        for sub = 1:length(sn)
            ANAs = getrow(ANA , ANA.SN == sn(sub));
            X = ANAs.IPI;clear d
            
            lab = ANAs.seqNumb +1;
            
            for i = 1:size(X,1)
                id = ones(size(X , 1) , 1);
                id(i) = 0;
                Y = inpaint_nans(X(~id , :));
                X1 = X(id==1 , :);
                lab1 = lab(id==1);
                for s = 1:7
                    m = nanmean(X1(lab1==s , :));
                    d(sub , i , s) = pdist([Y;m], distance);%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                end
                
            end
            seqs =  unique(lab);
            for l =  1:length(seqs)
                for l1 =  1:length(seqs)
                    id = lab == l;
                    out.D_IPI_sub.sub(count,1) = sub;
                    out.D_IPI_sub.d(count,1) = nanmean(d(sub , id , l1));
                    out.D_IPI_sub.s1(count,1) = l1;
                    out.D_IPI_sub.s2(count,1) = l;
                    out.D_IPI_sub.isWithin(count , 1) = (l == l1);
                    count = count + 1;
                end
            end
        end
        figure('color' , 'white')
        [xw , pw , ew] = lineplot(out.D_IPI_sub.isWithin , out.D_IPI_sub.d,'style_thickline');
        set(gca , 'FontSize' , 20 , 'XTick' , [0 1] , 'XTickLabels' , {'Between Structure' , 'Within Structure'},...
            'Box' , 'off');
        ylabel('Correlation Distance')
        title('Within Vs. Between Structure Crossvalidated Dissimilarity')
        grid on
        out = [];
    case 'crossvaldist_vel'
        %         if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
        %             subjnum = 11;
        %         end
        calcDist = 0;
        if calcDist
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [2:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).eyeveloc= N.norm(d,subjnum(SubN)).eyeveloc;
                    
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).eyeveloc{seqs(s)}];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; seqs(s)*ones(size(out(SubN, d-1).eyeveloc{seqs(s)} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                eyesig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).prsveloc= N.norm(d,subjnum(SubN)).prsveloc;
                    
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).prsveloc{seqs(s)}];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; seqs(s)*ones(size(out(SubN, d-1).prsveloc{seqs(s)} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                prssig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        else
            seqs = [2:7];
            load([baseDir , '/eye_dbetclass.mat'])
            load([baseDir , '/eye_dwithclass.mat'])
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            
            
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prs_dbetclass.mat'])
            load([baseDir , '/prs_dwithclass.mat'])
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            
        end
        
        
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(ePLOTw(SubN , :),eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(ePLOTb(SubN , :),eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Eye Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
            %             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Eye Velocity,  p =' , num2str(eyesig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(pPLOTw(SubN , :),pERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(pPLOTb(SubN , :),pERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Press Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
            %             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Press Velocity,  p =' , num2str(prssig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        out = [];
    case 'crossvaldist_chunk'
        %% chunk distances
        ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:6]) & ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.Rep , rep) & ~Dall.isError);
        % time normalizing the IPIs (sp the press indecies) to 1 : 1000 normalized samples
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , ANA.seqlength(tn))  - ANA.AllPressIdx(tn , 1)) / 500;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            IPI(tn , :) = detrend(diff(nIdx(tn ,:) , 1 , 2) , 'linear' , 2);
        end
        
        for tn = 1:length(ANA.TN)
            thresh = .3 * std(IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if IPI(tn , estChnkBndry(cb)) < nanmean(IPI(tn  ,:)) + thresh
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
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
        end
        
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
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
        end
        
        
        
        %% the boundry distance
        temp = diff(ChnkPlcmnt,1,2);
        temp(temp<0) = 0;
        chbndry = [ones(6,1) ~temp]; % all the first presses are one
        
        
        
        for s = 0:6
            A = getrow(ANA , ANA.seqNumb == s);
            CBD{s+1} = A.estChnkBndry;
        end
        diag_CnkBndry_es_es = [];
        offdiag_CnkBndry_es_es = [];
        
        diag_CnkBndry_es_ac = [];
        offdiag_CnkBndry_es_ac = [];
        figure('color' , 'white');
        
        for d = 2:5
            
            clear meanbnd allmeanbnd CB est_dist lab lab1
            ANA1 = getrow(ANA , ismember(ANA.Day , days{d}));
            
            X{d-1} = ANA1.estChnkBndry;
            lab = ANA1.seqNumb +1;
            
            for i = 1:size(X{d-1},1)
                id = ones(size(X{d-1} , 1) , 1);
                id(i) = 0;
                Y = X{d-1}(~id , :);
                X1 = X{d-1}(id==1 , :);
                lab1 = lab(id==1);
                for s = 1:7
                    m{d-1}(s,:) = nanmean(X1(lab1==s , :));
                    dis{d-1}(i , s) = pdist([Y;m{d-1}(s,:)], 'cityblock');%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                end
            end
            
            seqs =  unique(lab);
            for l =  1:length(seqs)
                for l1 =  1:length(seqs)
                    id = lab == l;
                    D_CnkBndry{d-1}(l,l1) = nanmean(dis{d-1}(id , l1));
                    Dist_est_est(d-1).A{l,l1} = dis{d-1}(id , l1);
                    if l1~=l
                        offdiag_CnkBndry_es_es = [offdiag_CnkBndry_es_es  ; [ Dist_est_est(d-1).A{l,l1}, d*ones(length(Dist_est_est(d-1).A{l,l1}) , 1)]];
                    else
                        diag_CnkBndry_es_es = [diag_CnkBndry_es_es  ; [ Dist_est_est(d-1).A{l,l1}, d*ones(length(Dist_est_est(d-1).A{l,l1}) , 1)]];
                    end
                    
                end
            end
            % between estimated and actual
            for s = 1:7
                for s2 = 1:6
                    act_est_dist{d-1}(s,s2) = nanmean(nanmean(pdist2(CBD{s} , chbndry(s2 , :) , 'cityblock')));
                    Dist_est_act(d-1).A{s,s2} = pdist2(CBD{s} , chbndry(s2 , :) , 'cityblock');
                    if s~=s2+1
                        offdiag_CnkBndry_es_ac = [offdiag_CnkBndry_es_ac  ; [Dist_est_act(d-1).A{s,s2}, d*ones(length(Dist_est_act(d-1).A{s,s2}) , 1)]];
                    else
                        diag_CnkBndry_es_ac = [diag_CnkBndry_es_ac  ; [Dist_est_act(d-1).A{s,s2}, d*ones(length(Dist_est_act(d-1).A{s,s2}) , 1)]];
                    end
                end
                
            end
            
            %             offdiag_CnkBndry = [offdiag_CnkBndry  ; [reshape(D_CnkBndry{d-1}(2:end ,2:end)-diag(NaN*ones(length(D_CnkBndry{d-1}(2:end ,2:end)),1)) , (length(seqs)-1)^2 , 1) d*ones((length(seqs)-1)^2 , 1)]];
            %             diag_CnkBndry = [diag_CnkBndry  ; [diag(D_CnkBndry{d-1}(2:end ,2:end))  d*ones((length(seqs)-1) , 1)]];
            %
            
            
            figure('color' , [1 1 1])
            
            imagesc(D_CnkBndry{d-1});
            title(['Crossvalidated RDM for Est. Chunk boundries - day ' , num2str(days{d})] , 'FontSize' , 20)
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
            ax.XTickLabelRotation = 45;
            ax.FontSize = 20;
            hold on
            line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            axis square
            colorbar
            
            figure('color' , [1 1 1])
            
            imagesc(act_est_dist{d-1});
            title(['Est. vs. Expected Chunking RDM - day ' , num2str(days{d})]  , 'FontSize' , 20)
            hold on
            ax = gca;
            line([.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0]);
            ax.XTick = [1:6];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
            xlabel('Expected');
            ylabel('Estimated');
            ax.XTickLabelRotation = 45;
            ax.FontSize = 20;
            axis square
            colorbar
        end
        
        
        figure('color' , 'white')
        imcounter = {[7:9] 1 2 3 4 5 6 };
        titl = {'Random'};
        for s = 1:size(ChnkPlcmnt , 1)
            titl = [titl num2str(ChnkPlcmnt(s , :))];
        end
        cnt = 1;
        for s = 1:7
            subplot(3,3,imcounter{cnt})
            for d = 2:5
                plot(m{d-1}(s,:), 'LineWidth' , 3)
                hold on
            end
            grid on
            set(gca , 'XLim' , [1 14] )
            if s>=2
                ids  = find(ChnkPlcmnt(s-1 , :) == 1);
                for j = 1:length(ids)
                    line([ids(j) ids(j)] , [0 1] , 'LineWidth' , 3 , 'color' , 'c' , 'LineStyle' , ':')
                end
            end
            legend({'Day2' , 'Day3' , 'Day4' , 'AllDays'})
            title(['Average estimated Chunk boundries - structure ' , titl{s}])
            cnt = cnt+1;
        end
        
        
        h1 = figure;
        [CO_off , PL_off , ER_off] = lineplot(offdiag_CnkBndry_es_es(:,2) , offdiag_CnkBndry_es_es(:,1) , 'plotfcn' , 'nanmean');
        hold on
        [CO_on , PL_on , ER_on] = lineplot(diag_CnkBndry_es_es(:,2) , diag_CnkBndry_es_es(:,1) , 'plotfcn' , 'nanmean');
        close(h1)
        
        h1 = figure;
        [CO1_off , PL1_off , ER1_off] = lineplot(offdiag_CnkBndry_es_ac(:,2) , offdiag_CnkBndry_es_ac(:,1) , 'plotfcn' , 'nanmean');
        hold on
        [CO1_on , PL1_on , ER1_on] = lineplot(diag_CnkBndry_es_ac(:,2) , diag_CnkBndry_es_ac(:,1) , 'plotfcn' , 'nanmean');
        close(h1)
        
        figure('color' , 'white')
        
        subplot(2,1,1)
        h1 = plotshade(CO_off' , PL_off , ER_off,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':');
        hold on
        h2 = plotshade(CO_on' , PL_on , ER_on,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':');
        legend([h1 h2] ,{'Between structures' , 'Within structures'})
        grid on
        title('Average Dissimilarity Between and Within estimated Chunking Structures in City Block Chunk boundary distance')
        ax = gca;
        ax.XLim = [1.5 5.5];
        ax.XTick = [2:5];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 20;
        ylabel('Average distance');
        
        subplot(2,1,2)
        h1 = plotshade(CO1_off' , PL1_off , ER1_off,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':');
        hold on
        h2 = plotshade(CO1_on' , PL1_on , ER1_on,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':');
        legend([h1 h2] ,{'Between structures' , 'Within structures'})
        grid on
        title('Average Dissimilarity Between and Within estimated/actual Chunking Structures in City Block Chunk boundary distance')
        ax = gca;
        ax.XLim = [1.5 5.5];
        ax.XTick = [2:5];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 20;
        ylabel('Average distance');
        
        
        
        
        
        out = [];
    case 'Sacc_All'
        horz  = 13;
        if calc
            isSymmetric = 1;
            %         subjnum = subjnum(~ismember(subjnum , [2 4 6])); % Chao's eye data is fucked up!!!!
            eyeinfo.PB         = [];
            eyeinfo.CB         = [];
            eyeinfo.Hor        = [];
            eyeinfo.sn         = [];
            eyeinfo.day        = [];
            eyeinfo.BN         = [];
            eyeinfo.sacPerSec  = [];
            eyeinfo.sacDur     = [];
            eyeinfo.sacPeakVel = [];
            eyeinfo.sacAmp     = [];
            eyeinfo.seqNumb    = [];
            eyeinfo.DigFixDur  = [];
            eyeinfo.prsnumb    = [];
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:6]) & ismember(Dall.SN , subjnum) & ismember(Dall.Day , [2 3 4]) & Dall.isgood & ~Dall.isError & cellfun(@length , Dall.xEyePosDigit)>1 & Dall.Group == 1);
            ignorDig = 2;
            for tn  = 1:length(ANA.TN)
                ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                a = find(ANA.ChunkBndry(tn , :));
                ANA.ChunkBndry(tn , a(ignorDig:end)-1) = 3;
                ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                ANA.ChunkBndry(tn , 1:ignorDig+1) = [-1 -1 -1];  % dont account for the first and last sseqeuce presses
                ANA.ChunkBndry(tn , end-ignorDig:end) = [-1 -1 -1];% dont account for the first and last sseqeuce presses
                ANA.DigFixWeighted(tn , :) = zeros(1 ,14);
                
                
                if isSymmetric
                    window = 50; %samples = 100ms
                    for p = 1:14
                        id = ANA.xEyePosDigit{tn , 1}<=p+.5 & ANA.xEyePosDigit{tn , 1}>p-.5;
                        if sum(id)
                            ANA.DigFixWeighted(tn , p) = mean(abs(ANA.xEyePosDigit{tn , 1}(id) - p))*(sum(id)/500)*1000;
                        else
                            ANA.DigFixWeighted(tn , p) = 0;
                        end
                        id = [ANA.AllPressIdx(tn , p) - window :ANA.AllPressIdx(tn , p) + window];
                        if id(1) > length(ANA.xEyePosDigit{tn}) | sum(id<0)>0
                             ANA.EyePressTimePos(tn , p) = NaN;
                        elseif id(end)>length(ANA.xEyePosDigit{tn})
                            ANA.EyePressTimePos(tn , p) = nanmedian(ANA.xEyePosDigit{tn}(id(1):end));
                        else
                            ANA.EyePressTimePos(tn , p) = nanmedian(ANA.xEyePosDigit{tn}(id));
                        end
                    end
                end
                perv_Ben           = [1:14] - ANA.EyePressTimePos(tn , :);
                goodid             = (perv_Ben>=-3.5 & perv_Ben<-.9);
                prsnumb            = find(goodid);
                count              = sum(goodid);
                eyeinfo.PB         = [eyeinfo.PB ;perv_Ben(goodid)'];
                eyeinfo.prsnumb    = [eyeinfo.prsnumb ;find(goodid')];
                eyeinfo.CB         = [eyeinfo.CB ;ANA.ChunkBndry(tn ,goodid)'];
                eyeinfo.Hor        = [eyeinfo.Hor ; ANA.Horizon(tn)*ones(count , 1)];
                eyeinfo.sn         = [eyeinfo.sn ; ANA.SN(tn)*ones(count , 1)];
                eyeinfo.day        = [eyeinfo.day ; ANA.Day(tn)*ones(count , 1)];
                eyeinfo.BN         = [eyeinfo.BN ; ANA.BN(tn)*ones(count , 1)];
                eyeinfo.sacPerSec  = [eyeinfo.sacPerSec ; ANA.SaccPerSec(tn)*ones(count , 1)];
                eyeinfo.sacDur     = [eyeinfo.sacDur ; mean(ANA.SaccDuration{tn})*ones(count , 1)];
                eyeinfo.sacPeakVel = [eyeinfo.sacPeakVel ; mean(ANA.SaccPeakVel{tn})*ones(count , 1)];
                eyeinfo.sacAmp     = [eyeinfo.sacAmp ; mean(ANA.SaccAmplitude{tn})*ones(count , 1)];
                eyeinfo.seqNumb    = [eyeinfo.seqNumb ; ANA.seqNumb(tn)*ones(count , 1)];
                eyeinfo.DigFixDur  = [eyeinfo.DigFixDur ;ANA.DigFixWeighted(tn ,goodid)'];
            end
            K_withinChunk = tapply(eyeinfo , {'day' , 'sn','CB'} , {'sacDur' , 'nanmedian'}, {'DigFixDur' , 'nanmedian'},...
                {'sacAmp' , 'nanmedian'} , {'sacPerSec' , 'nanmedian'}, {'PB' , 'nanmedian'} , 'subset' ,...
                eyeinfo.seqNumb ~= 0 & eyeinfo.CB~=-1);
            
            K_sqnum = tapply(eyeinfo , {'day' , 'sn','seqNumb' , 'prsnumb'} , {'sacDur' , 'nanmedian'}, {'sacPeakVel' , 'nanmedian'},...
                {'sacAmp' , 'nanmedian'} , {'sacPerSec' , 'nanmedian'}, {'PB' , 'nanmedian'},{'DigFixDur' , 'nanmedian'});
            
            K_rand_chnk = tapply(eyeinfo , {'day' , 'sn','seqNumb'} , {'sacDur' , 'nanmedian'}, {'sacPeakVel' , 'nanmedian'},...
                {'sacAmp' , 'nanmedian'} , {'sacPerSec' , 'nanmedian'}, {'PB' , 'nanmedian'},{'DigFixDur' , 'nanmedian'},...
                'subset' , eyeinfo.CB~=-1);
            
            save([baseDir , '/se1_eyeInfo.mat'] , 'eyeinfo' ,'K_withinChunk' , 'K_sqnum' , 'K_rand_chnk' , '-v7.3')
        else
            load([baseDir , '/se1_eyeInfo.mat'])
        end
         ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:6]) & ismember(Dall.SN , subjnum) & ismember(Dall.Day , [2 3 4]) & Dall.isgood & ~Dall.isError & cellfun(@length , Dall.xEyePosDigit)>1 & Dall.Group == 1);
        h = figure;
        dayz = {[2] , [3] , [4]};
        for sqn = 0:6
            for dd = 1:length(dayz)
            [xDur{sqn+1,dd},pDur{sqn+1,dd},eDur{sqn+1,dd}] = lineplot([K_sqnum.prsnumb] ,  K_sqnum.DigFixDur , 'plotfcn' , 'nanmean','subset' ,...
               ismember(K_sqnum.seqNumb , sqn) & ismember(K_sqnum.day , dayz{dd}), 'style_shade');
            end
        end
        close(h)
        
        h = figure('color' , 'white');
        
        for sqn = 0:6
            chnkpl = unique(ANA.ChnkPlcmnt(ANA.seqNumb == sqn , :) , 'rows');
            subplot(7 , 1 , sqn+1);
            for d = 1:length(dayz)
                errorbar(xDur{sqn+1,d}' , pDur{sqn+1,d},eDur{sqn+1,d} , 'LineWidth' , 3);
                hold on
                
            end
            firsts = find(chnkpl == 1);
            for j = 1:length(firsts)
                line([firsts(j) firsts(j)] , [0 300] , 'LineWidth'  ,2 , 'color' , 'black')
            end
            xlabel('Digits')
            
            grid on
            
            set(gca , 'FontSize' , 16 , 'YLim' , [0 300] , 'XTick' , [1:13] , 'XLim' , [0 14])
            legend({'day 2' , 'days 3' , 'days 4'})
            title(['Fixation Duration on Digits on stucure ' , num2str(sqn)])
        end
        
       

        
        K1 = K_rand_chnk;
        h = figure('color' , 'white');
        subplot(411)
        id = K1.day == 2;
        K0.day = K1.day(id) ;
        K0.DigFixDur = K1.DigFixDur(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        templab = K0.seqNumb >0; % just test Random vs all structured
        temp = anovaMixed(K0.DigFixDur , K0.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
        out.rand_chnk(1) = temp.eff(2).p;
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'DigFixDur' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}));
            K00.day = (d+1)*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
            templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.DigFixDur , K00.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  K0.DigFixDur , 'plotfcn' , 'nanmean','style_thickline');
        
        grid on
        set(gca , 'FontSize' , 20)
        title(['Weighted Average Fixation Duration  - Days 2 , 3 , 4'])
        xlabel('Structure number')
        ylabel('weighted ms')
        legend({['Day 1 P = ' , num2str(out.rand_chnk(1))] , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        %========================================= Preview
        subplot(412)
        id = K1.day == 2;
        K0.day = K1.day(id) ;
        K0.PB = K1.PB(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        templab = K0.seqNumb >0; % just test Random vs all structured
        temp = anovaMixed(K0.PB , K0.sn,'within',templab,{'seqNumb'},'intercept',1)  ;
        out.rand_chnk(1) = temp.eff(2).p;
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'PB' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}) );
            K00.day = (d+1)*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
            templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.PB , K00.sn,'within',templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  -K0.PB , 'plotfcn' , 'nanmean','style_thickline');
        ylabel('Digits')
        grid on
        set(gca , 'FontSize' , 20)
        title(['Amount Looking Ahead - Days 2 , 3 , 4'])
        xlabel('Structure number')
        legend({['Day 1 P = ' , num2str(out.rand_chnk(1))] , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        %========================================= Saccade per second
        subplot(413)
        id = K1.day ==2;
        K0.day = K1.day(id) ;
        K0.sacPerSec = K1.sacPerSec(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        templab = K0.seqNumb >0; % just test Random vs all structured
        temp = anovaMixed(K0.sacPerSec , K0.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
        out.rand_chnk(1) = temp.eff(2).p;
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'sacPerSec' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}));
            K00.day = (d+1)*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
            templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.sacPerSec , K00.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  K0.sacPerSec , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        ylabel('Saccade per sec')
        title(['Saccade rate  - Days 2 , 3 , 4'])
        xlabel('Structure number')
        legend({['Day 1 P = ' , num2str(out.rand_chnk(1))] , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        %========================================= Saccade amplitude
        subplot(414)
        id = K1.day == 2;
        K0.day = K1.day(id) ;
        K0.sacAmp = K1.sacAmp(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        templab = K0.seqNumb >0; % just test Random vs all structured
        temp = anovaMixed(K0.sacAmp , K0.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
        out.rand_chnk(1) = temp.eff(2).p;
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'sacAmp' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}));
            K00.day = (d+1)*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
            templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.sacAmp , K00.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  K0.sacAmp , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Sccade Amplitude - Days 2 , 3 , 4'])
        xlabel('Structure number')
        legend({['Day 1 P = ' , num2str(out.rand_chnk(1))] , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        ylabel('Digits')
       
        
        
        K1 = K_withinChunk;
        % within chunked
        h = figure('color' , 'white');
        subplot(211)
        id = K1.day == 2;
        K0.day = K1.day(id) ;
        K0.DigFixDur = K1.DigFixDur(id);
        K0.CB  = K1.CB(id);
        K0.sn = K1.sn(id);
        templab = K0.CB ~= 2; % just test middle vs rest
        temp = anovaMixed(K0.DigFixDur , K0.sn,'within', templab ,{'seqNumb'},'intercept',1)  ;
        out.withinChunk(1) = temp.eff(2).p;
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'CB' , 'sn'} , {'DigFixDur' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}));
            K00.day = (d+1)*ones(size(K00.CB));
            K0 = addstruct(K0 , K00);
            templab = K00.CB ~= 2; % just test middle vs rest
            temp = anovaMixed(K00.DigFixDur , K00.sn,'within', templab ,{'seqNumb'},'intercept',1)  ;
            out.withinChunk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.CB] ,  K0.DigFixDur , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Average Fixation Duration in Horizon(s)  - Days 2 , 3 , 4'])
        xlabel('Segment Placement - 1: first      2: within    3: Last')
        legend({['Day 1 P = ' , num2str(out.withinChunk(1))] , ['Day 2 P = ' , num2str(out.withinChunk(2))], ['Day 3 P = ' , num2str(out.withinChunk(3))]})
        ylabel('Weighted ms')
        %========================================= Preview
         subplot(212)
        id = K1.day == 2;
        K0.day = K1.day(id) ;
        K0.PB = K1.PB(id);
        K0.CB  = K1.CB(id);
        K0.sn = K1.sn(id);
        templab = K0.CB ~= 3; % just test middle vs rest
        temp = anovaMixed(K0.PB , K0.sn,'within', templab ,{'seqNumb'},'intercept',1)  ;
        out.withinChunk(1) = temp.eff(2).p;
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'CB' , 'sn'} , {'PB' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}));
            K00.day = (d+1)*ones(size(K00.CB));
            K0 = addstruct(K0 , K00);
            ttemplab = K00.CB ~= 3; % just test middle vs rest
            temp = anovaMixed(K00.PB , K00.sn,'within', templab ,{'seqNumb'},'intercept',1)  ;
            out.withinChunk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.CB] ,  -K0.PB , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Looking ahead - Days 2 , 3 , 4'])
        xlabel('Segment Placement - 1: first      2: within    3: Last')
        ylabel('Digits')
        legend({['Day 1 P = ' , num2str(out.withinChunk(1))] , ['Day 2 P = ' , num2str(out.withinChunk(2))], ['Day 3 P = ' , num2str(out.withinChunk(3))]})
    case 'saccs_Singlesubj'
        SNu = [1 3 4 5 6 7 8 9 10];
        for subjnum = 1:length(SNu)
            for d = 2:4
                ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:6]) & ismember(Dall.SN , SNu(subjnum)) & Dall.isgood & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , d));
                ANA.seqNumb(ANA.seqNumb > 0) = 1;
                sp(subjnum , d-1,1)   = nanmean(ANA.SaccPerSec(ismember(ANA.seqNumb , 0)));
                sp(subjnum , d-1,2)   = nanmean(ANA.SaccPerSec(ismember(ANA.seqNumb , [1:6])));
                sp_s(subjnum , d-1,1) = nanstd(ANA.SaccPerSec(ismember(ANA.seqNumb , 0)));
                sp_s(subjnum , d-1,2) = nanstd(ANA.SaccPerSec(ismember(ANA.seqNumb , [1:6])));
                
                [~ , out.SP(subjnum , d-1)] = ttest2(ANA.SaccPerSec(ANA.seqNumb == 1) , ANA.SaccPerSec(ANA.seqNumb == 0));
                [subjnum , d]
                %             out.SP(d-1) = anovaMixed(ANA.SaccPerSec , ANA.SN ,'within',[ANA.seqNumb],{'Seq type'},'subset' , ~isnan(ANA.SaccPerSec) , 'intercept',1);
                
                
                ns(subjnum , d-1,1)   = nanmean(ANA.NumSaccs(ismember(ANA.seqNumb , 0)));
                ns(subjnum , d-1,2)   = nanmean(ANA.NumSaccs(ismember(ANA.seqNumb , [1:6])));
                ns_s(subjnum , d-1,1) = nanstd(ANA.NumSaccs(ismember(ANA.seqNumb , 0)));
                ns_s(subjnum , d-1,2) = nanstd(ANA.NumSaccs(ismember(ANA.seqNumb , [1:6])));
                
                [~ , out.NS(subjnum , d-1)] = ttest2(ANA.NumSaccs(ANA.seqNumb == 1) , ANA.NumSaccs(ANA.seqNumb == 0));
                %             out.NS(d-1) = anovaMixed(ANA.NumSaccs , ANA.SN ,'within',[ANA.seqNumb],{'Seq type'},'subset' , ~isnan(ANA.NumSaccs) , 'intercept',1);
                
                FD = [];
                seq = [];
                subs = [];
                for tn = 1:length(ANA.TN)
                    FD  = [FD ; ANA.EyeFixDuration{tn}];
                    seq = [seq ; ANA.seqNumb(tn) *ones(length(ANA.EyeFixDuration{tn}) ,1)];
                    subs = [subs ; ANA.SN(tn)*ones(length(ANA.EyeFixDuration{tn}) ,1)];
                end
                [~ , out.FD(subjnum , d-1)] = ttest2(FD(seq == 1) , FD(seq == 0));
                %             out.FD(d-1) = anovaMixed(FD , subs ,'within',[seq],{'Seq type'},'subset' , ~isnan(FD) , 'intercept',1);
                
                
                fd(subjnum , d-1,1)   = nanmean(FD(ismember(seq , 0)));
                fd(subjnum , d-1,2)   = nanmean(FD(ismember(seq , [1:6])));
                fd_s(subjnum , d-1,1) = nanstd(FD(ismember(seq , 0)));
                fd_s(subjnum , d-1,2) = nanstd(FD(ismember(seq , [1:6])));
                SD = [];
                seq = [];
                subs = [];
                for tn = 1:length(ANA.TN)
                    SD  = [SD ; ANA.SaccDuration{tn}];
                    seq = [seq ; ANA.seqNumb(tn) *ones(length(ANA.SaccDuration{tn}) ,1)];
                    subs = [subs ; ANA.SN(tn)*ones(length(ANA.SaccDuration{tn}) ,1)];
                end
                %             out.SD(d-1) = anovaMixed(SD , subs ,'within',[seq],{'Seq type'},'subset' , ~isnan(SD) , 'intercept',1);
                [~ , out.SD(subjnum , d-1)] = ttest2(SD(seq == 1) , SD(seq == 0));
                
                sd(subjnum , d-1,1)   = nanmean(SD(ismember(seq , 0)));
                sd(subjnum , d-1,2)   = nanmean(SD(ismember(seq , [1:6])));
                sd_s(subjnum , d-1,1) = nanstd(SD(ismember(seq , 0)));
                sd_s(subjnum , d-1,2) = nanstd(SD(ismember(seq , [1:6])));
                
                for tn = 1:length(ANA.TN)
                    SPV(tn,1) = nanmean(ANA.SaccPeakVel{tn});
                    subs1(tn ,1) = ANA.SN(tn);
                    seq1(tn ,1) = ANA.seqNumb(tn);
                end
                [~ , out.SPV(subjnum , d-1)] = ttest2(SPV(seq1 == 1) , SPV(seq1 == 0));
                %             out.SPV(d-1) = anovaMixed(SPV , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(SPV) , 'intercept',1);
                spv(subjnum , d-1,1)   = nanmean(SPV(ismember(ANA.seqNumb , 0)));
                spv(subjnum , d-1,2)   = nanmean(SPV(ismember(ANA.seqNumb , [1:6])));
                spv_s(subjnum , d-1,1) = nanstd(SPV(ismember(ANA.seqNumb , 0)));
                spv_s(subjnum , d-1,2) = nanstd(SPV(ismember(ANA.seqNumb , [1:6])));
                
                
                for tn = 1:length(ANA.TN)
                    SA(tn,1) = nanmean(ANA.SaccAmplitude{tn});
                end
                [~ , out.SA(subjnum , d-1)] = ttest2(SA(seq1 == 1) , SA(seq1 == 0));
                %             out.SA(d-1) = anovaMixed(SA , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(SA) , 'intercept',1);
                sa(subjnum , d-1,1)   = nanmean(SA(ismember(ANA.seqNumb , 0)));
                sa(subjnum , d-1,2)   = nanmean(SA(ismember(ANA.seqNumb , [1:6])));
                sa_s(subjnum , d-1,1) = nanstd(SA(ismember(ANA.seqNumb , 0)));
                sa_s(subjnum , d-1,2) = nanstd(SA(ismember(ANA.seqNumb , [1:6])));
                
                for tn = 1:length(ANA.TN)
                    perv_Ben(tn , :)        = nanmean([1:14] - ANA.EyePressTimePos(tn , :));
                end
                %             out.PB(d-1) = anovaMixed(perv_Ben , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(perv_Ben) , 'intercept',1);
                [~ , out.PB(subjnum , d-1)] = ttest2(perv_Ben(seq1 == 1) , perv_Ben(seq1 == 0));
                pb(subjnum , d-1,1)   = nanmean(perv_Ben(ismember(ANA.seqNumb , 0)));
                pb(subjnum , d-1,2)   = nanmean(perv_Ben(ismember(ANA.seqNumb , [1:6])));
                pb_s(subjnum , d-1,1) = nanstd(perv_Ben(ismember(ANA.seqNumb , 0)));
                pb_s(subjnum , d-1,2) = nanstd(perv_Ben(ismember(ANA.seqNumb , [1:6])));
            end
        end
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(sp(subjnum,:,2) , sp_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(sp(subjnum,:,1) , sp_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            title(['Sub ' , num2str(SNu(subjnum)), ' Saccade rate - LDp = ',num2str(out.SP(subjnum , 3))])
            ylabel('Saccades per Second')
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SaccPerSec p-val = ' , num2str(out.SP(subjnum , 1))])
            disp(['Day 2 SaccPerSec p-val = ' , num2str(out.SP(subjnum , 2))])
            disp(['Day 3 SaccPerSec p-val = ' , num2str(out.SP(subjnum , 3))])
        end
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(ns(subjnum,:,2) , ns_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(ns(subjnum,:,1) , ns_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            title(['Sub ' , num2str(subjnum), ' Number of Saccades per trial - LDp = ' , num2str(out.NS(subjnum , 3))])
            ylabel('Saccades')
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 NumSaccades p-val = ' , num2str(out.NS(subjnum , 1))])
            disp(['Day 2 NumSaccades p-val = ' , num2str(out.NS(subjnum , 2))])
            disp(['Day 3 NumSaccades p-val = ' , num2str(out.NS(subjnum , 3))])
        end
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(fd(subjnum,:,2) , fd_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(fd(subjnum,:,1) , fd_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('msec')
            title(['Average Fixation Duration - LDp = ' ,num2str(out.FD(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 FixDuration p-val = ' , num2str(out.FD(subjnum , 1))])
            disp(['Day 2 FixDuration p-val = ' , num2str(out.FD(subjnum , 2))])
            disp(['Day 3 FixDuration p-val = ' , num2str(out.FD(subjnum , 3))])
        end
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(sd(subjnum,:,2) , sd_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(sd(subjnum,:,1) , sd_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('msec')
            title(['Average Saccade Duration - LDp = ',num2str(out.SD(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SacDuration p-val = ' , num2str(out.SD(subjnum , 1))])
            disp(['Day 2 SacDuration p-val = ' , num2str(out.SD(subjnum , 2))])
            disp(['Day 3 SacDuration p-val = ' , num2str(out.SD(subjnum , 3))])
        end
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(sa(subjnum,:,2) , sa_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(sa(subjnum,:,1) , sa_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('deg')
            title(['Average Saccade Amplitude - LDp = ' , num2str(out.SA(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SacAmplitude p-val = ' , num2str(out.SA(subjnum , 1))])
            disp(['Day 2 SacAmplitude p-val = ' , num2str(out.SA(subjnum , 2))])
            disp(['Day 3 SacAmplitude p-val = ' , num2str(out.SA(subjnum , 3))])
        end
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(spv(subjnum,:,2) , spv_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(spv(subjnum,:,1) , spv_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('deg/sec')
            title(['Average Saccade Peak Velocity - LDp = ' , num2str(out.SPV(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SacPeakVelocity p-val = ' , num2str(out.SPV(subjnum , 1))])
            disp(['Day 2 SacPeakVelocity p-val = ' , num2str(out.SPV(subjnum , 2))])
            disp(['Day 3 SacPeakVelocity p-val = ' , num2str(out.SPV(subjnum , 3))])
        end
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(-pb(subjnum,:,2) , pb_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(-pb(subjnum,:,1) , pb_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('Digits')
            title(['Average Preview Benefit in - LDp = ' , num2str(out.PB(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 PervBenefit p-val = ' , num2str(out.PB(subjnum , 1))])
            disp(['Day 2 PervBenefit p-val = ' , num2str(out.PB(subjnum , 2))])
            disp(['Day 3 PervBenefit p-val = ' , num2str(out.PB(subjnum , 3))])
        end
    case 'sacc_Chunk_all'
        % The alalysis routine allocated the space to digits asymmetrically,
        % so that all the space to the left of a digit is allocated to that digit
        % p-1<x<p ---> allocated to digit p
        % if you want symmetrical digit allocation
        %         prompt = 'Which Horizon?';
        %         h = input(prompt);
        isSymmetric = 1;
        %         subjnum = subjnum(~ismember(subjnum , [2 4 6])); % Chao's eye data is fucked up!!!!
        days = {2 3 4};
        for d = 1:length(days)
            DigFix{d}  = [];
            PervBen{d} = [];
            isSacc{d}  = [];
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [1:6]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & ismember(Dall.Day , days{d}) & cellfun(@length , Dall.xEyePosDigit)>1 & ismember(Dall.Rep , rep) );
            for tn  = 1:length(ANA.TN)
                ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                a = find(ANA.ChunkBndry(tn , :));
                ANA.ChunkBndry(tn , a(2:end)-1) = 3;
                ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                ANA.ChunkBndry(tn , 1:3) = [-1 -1 -1];  % dont account for the first and last sseqeuce presses
                ANA.ChunkBndry(tn , end-2:end) = [-1 -1 -1];% dont account for the first and last sseqeuce presses
                ANA.DigFixWeight(tn , :) = zeros(1 ,14);
                
                
                if isSymmetric
                    window = 12;
                    ANA.EyeFixDigit{tn , 1} = ANA.xEyePosDigit{tn}((ANA.AllPressIdx(tn , 1)+2)-window :ANA.AllPressIdx(tn , ANA.seqlength(tn)) + 2 + window) .* ANA.SaccFlag{tn};
                    for p = 1:14
                        id = ANA.EyeFixDigit{tn , 1}<=p+.5 & ANA.EyeFixDigit{tn , 1}>p-.5;
                        %                         ANA.EyeFixDigit{tn , 1}(id) = p;
                        if sum(id)
                            ANA.DigFixWeight(tn , p) = (sum(abs(1- (ANA.EyeFixDigit{tn , 1}(id) - p)))/sum(id))*(sum(id)/500);
                        else
                            ANA.DigFixWeight(tn , p) = 0;
                        end
                    end
                end
                perv_Ben        = [1:14] - ANA.EyePressTimePos(tn , :);
                perv_Ben(perv_Ben<=-3.5 | perv_Ben>-.9) = NaN;
                PervBen{d} = [PervBen{d} ;[perv_Ben(ANA.ChunkBndry(tn , :) == 1)' ...
                    (d)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                    ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))]];
                
                PervBen{d} = [PervBen{d} ;[perv_Ben(ANA.ChunkBndry(tn , :) == 2)' ...
                    (d)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                    2*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))]];
                
                PervBen{d} = [PervBen{d} ;[perv_Ben(ANA.ChunkBndry(tn , :) == 3)' ...
                    (d)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))...
                    3*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))...
                    ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))...
                    ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))]];
                
                
                
                
                ANA.DigFixWeight(tn , perv_Ben<=-3.5 | perv_Ben>-.9) = NaN;
                DigFix{d} = [DigFix{d} ;[ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)' ...
                    (d)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                    ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.SN(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.seqNumb(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))]];
                
                DigFix{d} = [DigFix{d} ;[ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)' ...
                    (d)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                    2* ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.SN(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.seqNumb(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))]];
                
                DigFix{d} = [DigFix{d} ;[ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)' ...
                    (d)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))...
                    3*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))...
                    ANA.SN(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))...
                    ANA.seqNumb(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))]];
                
                
                
                
                
                
                isSacc{d} = [isSacc{d} ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 1))/sum(ANA.ChunkBndry(tn , :) == 1) ...
                    d 1 ANA.SN(tn)  ANA.seqNumb(tn)]];
                
                isSacc{d} = [isSacc{d} ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 2))/sum(ANA.ChunkBndry(tn , :) == 2) ...
                    d 2 ANA.SN(tn)  ANA.seqNumb(tn)]];
                
                isSacc{d} = [isSacc{d} ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 3))/sum(ANA.ChunkBndry(tn , :) == 3) ...
                    d 3 ANA.SN(tn)  ANA.seqNumb(tn)]];
                
            end
            
            out.FD(d) = anovaMixed(DigFix{d}(:,1) , DigFix{d}(:,4) ,'within',[DigFix{d}(:,3)],{'Chunk place'},'subset' , DigFix{d}(:,2) == d & ~isnan(DigFix{d}(:,1)), 'intercept',1);
            out.PB(d) = anovaMixed(PervBen{d}(:,1) , PervBen{d}(:,4) ,'within',[PervBen{d}(:,3)],{'Chunk place'},'subset' , PervBen{d}(:,2) == d & ~isnan(PervBen{d}(:,1)), 'intercept',1);
            out.iS(d) = anovaMixed(isSacc{d}(:,1) , isSacc{d}(:,4) ,'within',[isSacc{d}(:,3)],{'Chunk place'},'subset' , isSacc{d}(:,2) == d &~isnan(isSacc{d}(:,1)), 'intercept',1);
            if d == 1
                all_ANA = ANA;
            else
                all_ANA = addstruct(all_ANA , ANA);
            end
        end
        clear TEMP lengtemp leng
        FD = [];
        cp1 = [];
        for chnkprs = 1:3
            idoi = DigFix{end}(:,3) == chnkprs;
            leng(chnkprs) = length(DigFix{end}(idoi , 1));
            TEMP{chnkprs} = NaN * zeros(leng(chnkprs) ,d);
            TEMP{chnkprs}(:,d) =  DigFix{end}(idoi , 1);
            for dd =  1:length(days)
                idtemp = DigFix{dd}(:,3) == chnkprs;
                lengtemp(chnkprs) = length(DigFix{dd}(idtemp , 1));
                TEMP{chnkprs}(1:lengtemp(chnkprs),dd) = DigFix{dd}(idtemp,1);
            end
            FD = [FD;TEMP{chnkprs}];
            cp1 = [cp1;chnkprs*ones(length(TEMP{chnkprs}) , 1)];
        end
        
        PB = [];
        cp2 = [];
        for chnkprs = 1:3
            idoi = PervBen{end}(:,3) == chnkprs;
            leng(chnkprs) = length(PervBen{end}(idoi , 1));
            TEMP{chnkprs} = NaN * zeros(leng(chnkprs) ,d);
            TEMP{chnkprs}(:,d) =  PervBen{end}(idoi , 1);
            for dd =  1:length(days)
                idtemp = PervBen{dd}(:,3) == chnkprs;
                lengtemp(chnkprs) = length(PervBen{dd}(idtemp , 1));
                TEMP{chnkprs}(1:lengtemp(chnkprs),dd) = PervBen{dd}(idtemp,1);
            end
            PB = [PB;TEMP{chnkprs}];
            cp2 = [cp2;chnkprs*ones(length(TEMP{chnkprs}) , 1)];
        end
        
        iS = [];
        cp3 = [];
        for chnkprs = 1:3
            idoi = isSacc{end}(:,3) == chnkprs;
            leng(chnkprs) = length(isSacc{end}(idoi , 1));
            TEMP{chnkprs} = NaN * zeros(leng(chnkprs) ,d);
            TEMP{chnkprs}(:,d) =  isSacc{end}(idoi , 1);
            for dd =  1:length(days)
                idtemp = isSacc{dd}(:,3) == chnkprs;
                lengtemp(chnkprs) = length(isSacc{dd}(idtemp , 1));
                TEMP{chnkprs}(1:lengtemp(chnkprs),dd) = isSacc{dd}(idtemp,1);
            end
            iS = [iS;TEMP{chnkprs}];
            cp3 = [cp3;chnkprs*ones(length(TEMP{chnkprs}) , 1)];
        end
        
        %% Random
        for d = 1:length(days)
            DigFix_rnd{d}  = [];
            PervBen_rnd{d} = [];
            isSacc_rnd{d}  = [];
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [0]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & ismember(Dall.Day , days{d}) & cellfun(@length , Dall.xEyePosDigit)>1 & ismember(Dall.Rep , rep));
            for tn  = 1:length(ANA.TN)
                ANA.ChunkBndry(tn , :) = zeros(1 ,14);
                ANA.ChunkBndry(tn , 1:3) = [-1 -1 -1];  % dont account for the first and last sseqeuce presses
                ANA.ChunkBndry(tn , end-2:end) = [-1 -1 -1];% dont account for the first and last sseqeuce presses
                ANA.DigFix_rndWeight(tn , :) = zeros(1 ,14);
                
                
                if isSymmetric
                    window = 12;
                    [d tn]
                    ANA.EyeFixDigit{tn , 1} = ANA.xEyePosDigit{tn}((ANA.AllPressIdx(tn , 1)+2)-window :ANA.AllPressIdx(tn , ANA.seqlength(tn)) + 2 + window) .* ANA.SaccFlag{tn};
                    for p = 1:14
                        id = ANA.EyeFixDigit{tn , 1}<=p+.5 & ANA.EyeFixDigit{tn , 1}>p-.5;
                        %                         ANA.EyeFixDigit{tn , 1}(id) = p;
                        if sum(id)
                            ANA.DigFix_rndWeight(tn , p) = (sum(abs(1- (ANA.EyeFixDigit{tn , 1}(id) - p)))/sum(id))*(sum(id)/500);
                        else
                            ANA.DigFix_rndWeight(tn , p) = 0;
                        end
                    end
                end
                perv_Ben_rnd        = [1:14] - ANA.EyePressTimePos(tn , :);
                perv_Ben_rnd(perv_Ben_rnd<=-3.5 | perv_Ben_rnd>-.9) = NaN;
                PervBen_rnd{d} = [PervBen_rnd{d} ;[perv_Ben_rnd(ANA.ChunkBndry(tn , :) == 0)' ...
                    (d)*ones(size(perv_Ben_rnd(ANA.ChunkBndry(tn , :) == 0)'))...
                    zeros(size(perv_Ben_rnd(ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.SN(tn)*ones(size(perv_Ben_rnd(ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.seqNumb(tn)*ones(size(perv_Ben_rnd(ANA.ChunkBndry(tn , :) == 0)'))]];
                
                
                
                ANA.DigFix_rndWeight(tn , perv_Ben_rnd<=-3.5 | perv_Ben_rnd>-.9) = NaN;
                DigFix_rnd{d} = [DigFix_rnd{d} ;[ANA.DigFix_rndWeight(tn , ANA.ChunkBndry(tn , :) == 0)' ...
                    (d)*ones(size(ANA.DigFix_rndWeight(tn , ANA.ChunkBndry(tn , :) == 0)'))...
                    zeros(size(ANA.DigFix_rndWeight(tn , ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.SN(tn)*ones(size(ANA.DigFix_rndWeight(tn , ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.seqNumb(tn)*ones(size(ANA.DigFix_rndWeight(tn , ANA.ChunkBndry(tn , :) == 0)'))]];
                
                
                isSacc_rnd{d} = [isSacc_rnd{d} ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 0))/sum(ANA.ChunkBndry(tn , :) == 0) d 0 ANA.SN(tn)  ANA.seqNumb(tn)]];
                
                
            end
            
            if d == 1
                all_ANA = ANA;
            else
                all_ANA = addstruct(all_ANA , ANA);
            end
        end
        clear TEMP lengtemp leng
        FD_rnd = [];
        cp1_rnd = [];
        idoi = DigFix_rnd{end}(:,3) == 0;
        leng = length(DigFix_rnd{end}(idoi , 1));
        TEMP = NaN * zeros(leng ,d);
        TEMP(:,d) =  DigFix_rnd{end}(idoi , 1);
        for dd =  1:length(days)
            idtemp = DigFix_rnd{dd}(:,3) == 0;
            lengtemp = length(DigFix_rnd{dd}(idtemp , 1));
            TEMP(1:lengtemp,dd) = DigFix_rnd{dd}(idtemp,1);
        end
        FD_rnd = [FD_rnd;TEMP];
        cp1_rnd = [cp1_rnd;zeros(length(TEMP) , 1)];
        
        
        PB_rnd = [];
        cp2_rnd = [];
        idoi = PervBen_rnd{end}(:,3) == 0;
        leng = length(PervBen_rnd{end}(idoi , 1));
        TEMP = NaN * zeros(leng ,d);
        TEMP(:,d) =  PervBen_rnd{end}(idoi , 1);
        for dd =  1:length(days)
            idtemp = PervBen_rnd{dd}(:,3) == 0;
            lengtemp = length(PervBen_rnd{dd}(idtemp , 1));
            TEMP(1:lengtemp,dd) = PervBen_rnd{dd}(idtemp,1);
        end
        PB_rnd = [PB_rnd;TEMP];
        cp2_rnd = [cp2_rnd;zeros(length(TEMP) , 1)];
        
        
        iS_rnd = [];
        cp3_rnd = [];
        idoi = isSacc_rnd{end}(:,3) == 0;
        leng = length(isSacc_rnd{end}(idoi , 1));
        TEMP = NaN * zeros(leng ,d);
        TEMP(:,d) =  isSacc_rnd{end}(idoi , 1);
        for dd =  1:length(days)
            idtemp = isSacc_rnd{dd}(:,3) == 0;
            lengtemp = length(isSacc_rnd{dd}(idtemp , 1));
            TEMP(1:lengtemp,dd) = isSacc_rnd{dd}(idtemp,1);
        end
        iS_rnd = [iS_rnd;TEMP];
        cp3_rnd = [cp3_rnd;zeros(length(TEMP) , 1)];
        
        justFix = 0;
        
        cols = {'blue' , 'red' , 'c'};
        if justFix
            figure('color' , 'white');
            htemp = figure;
            [xcoords,PLOT_df,ERROR_df] = lineplot([cp1 ;cp1_rnd] , [FD ;FD_rnd],'plotfcn','nanmean','leg' ,{'Day1' , 'Day2,3' , 'Day4,5'},'markersize',3,'linewidth',3,'markertype' , {'o' , '*' , '+'},'markercolor' , {'blue' , 'red' , 'c'},...
                'errorcolor' ,{'blue' , 'red' , 'c'} , 'linecolor' ,{'blue' , 'red' , 'c'} );
            close(htemp);
            
            day = 1;
            h1 = plotshade(xcoords',PLOT_df(day,:) , ERROR_df(day,:),'transp' , .2 , 'patchcolor' , cols{day} , 'linecolor' , cols{day} , 'linewidth' , 3 , 'linestyle' , ':');
            hold on
            day = 2;
            h2 = plotshade(xcoords',PLOT_df(day,:) , ERROR_df(day,:),'transp' , .2 , 'patchcolor' , cols{day} , 'linecolor' , cols{day} , 'linewidth' , 3 , 'linestyle' , ':');
            hold on
            day = 3;
            h3 = plotshade(xcoords',PLOT_df(day,:) , ERROR_df(day,:),'transp' , .2 , 'patchcolor' , cols{day} , 'linecolor' , cols{day} , 'linewidth' , 3 , 'linestyle' , ':');
            
            legend([h1 h2 h3] ,{'Day 2' , 'Day 3' , 'Day 4'})
            
            hold on
            set(gca ,'XTick', [0:3] ,'XTickLabel' ,{'Random' 'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'} ,'XTickLabelRotation' , 45 ,'FontSize' , 16)
            grid on
            ylabel('Sec')
            title(['Digit Fix Duration'])
            
        else
            figure('color' , 'white');
            htemp = figure;
            [xcoords,PLOT_df,ERROR_df] = lineplot([cp1 ;cp1_rnd] , [FD ;FD_rnd],'plotfcn','nanmean','leg' ,{'Day1' , 'Day2,3' , 'Day4,5'},'markersize',3,'linewidth',3,'markertype' , {'o' , '*' , '+'},'markercolor' , {'blue' , 'red' , 'c'},...
                'errorcolor' ,{'blue' , 'red' , 'c'} , 'linecolor' ,{'blue' , 'red' , 'c'} );
            close(htemp);
            
            day = 1;
            h1 = plotshade(xcoords',PLOT_df(day,:) , ERROR_df(day,:),'transp' , .2 , 'patchcolor' , cols{day} , 'linecolor' , cols{day} , 'linewidth' , 3 , 'linestyle' , ':');
            hold on
            day = 2;
            h2 = plotshade(xcoords',PLOT_df(day,:) , ERROR_df(day,:),'transp' , .2 , 'patchcolor' , cols{day} , 'linecolor' , cols{day} , 'linewidth' , 3 , 'linestyle' , ':');
            hold on
            day = 3;
            h3 = plotshade(xcoords',PLOT_df(day,:) , ERROR_df(day,:),'transp' , .2 , 'patchcolor' , cols{day} , 'linecolor' , cols{day} , 'linewidth' , 3 , 'linestyle' , ':');
            
            legend([h1 h2 h3] ,{'Day 2' , 'Day 3' , 'Day 4'})
            
            hold on
            set(gca ,'XTick', [] ,'FontSize' , 16)
            set(gca ,'XTick', [0:3] ,'XTickLabel' ,{'Random','First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'} ,'XTickLabelRotation' , 45)
            grid on
            ylabel('Sec')
            title(['Digit Fix Duration'])
            figure('color' , 'white');
            htemp = figure;
            [xcoords,PLOT_df,ERROR_df] = lineplot([cp2; cp2_rnd] , [-PB ;-PB_rnd],'plotfcn','nanmean','leg' ,{'Day1' , 'Day2,3' , 'Day4,5' },'markersize',3,'linewidth',3,'markertype' , {'o' , '*' , '+'},'markercolor' , {'blue' , 'red' , 'c'},...
                'errorcolor' ,{'blue' , 'red' , 'c'} , 'linecolor' ,{'blue' , 'red' , 'c'});
            close(htemp);
            
            chnkprs = 1;
            h1 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
            hold on
            chnkprs = 2;
            h2 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
            hold on
            chnkprs = 3;
            h3 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
            
            legend([h1 h2 h3] ,{'Day 2' , 'Day 3' , 'Day 4'})
            
            hold on
            set(gca ,'XTick', [] ,'FontSize' , 16)
            set(gca ,'XTick', [0:3] ,'XTickLabel' ,{'Random','First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'} ,'XTickLabelRotation' , 45)
            grid on
            ylabel('Sec')
            title('Fixation/Press distance')
            
        end
        
        disp(['Day 1 FixDuration p-val = ' , num2str(out.FD(1).eff(2).p)])
        disp(['Day 2 FixDuration p-val = ' , num2str(out.FD(2).eff(2).p)])
        disp(['Day 3 FixDuration p-val = ' , num2str(out.FD(3).eff(2).p)])
        
        
        disp(['Day 1 isSacc p-val = ' , num2str(out.iS(1).eff(2).p)])
        disp(['Day 2 isSacc p-val = ' , num2str(out.iS(2).eff(2).p)])
        disp(['Day 3 isSacc p-val = ' , num2str(out.iS(3).eff(2).p)])
        
        
        disp(['Day 1 PervBenefit p-val = ' , num2str(out.PB(1).eff(2).p)])
        disp(['Day 2 PervBenefit p-val = ' , num2str(out.PB(2).eff(2).p)])
        disp(['Day 3 PervBenefit p-val = ' , num2str(out.PB(3).eff(2).p)])
    case 'Random_error'
        for sub = subjnum
            clear C chunks
            ANAChunk = getrow(Dall ,Dall.seqNumb > 7 & ismember(Dall.SN , sub));
            ANAChunk.AllPress(isnan(ANAChunk.AllPress)) = 0;
            C = unique(ANAChunk.AllPress , 'rows');
            DuoCounter  = 1;
            TrpCounter  = 1;
            QudCounter  = 1;
            for c = 1:size(C , 1)
                chunks{c} = C(c, C(c,:)~=0);
                switch length(chunks{c})
                    case 2
                        doubles(DuoCounter , :) = chunks{c};
                        doubleInfo(DuoCounter , 1) = 2; % comes form double chunk
                        doubleInfo(DuoCounter , 2) = 1;
                        DuoCounter = DuoCounter + 1;
                    case 3
                        for pcntr = 1:2
                            doubles(DuoCounter , :) = chunks{c}(pcntr : pcntr+1);
                            doubleInfo(DuoCounter , 1) = 3; % comes form triple chunk
                            doubleInfo(DuoCounter , 2) = pcntr; % to know if it's the first or second duo in the reiplet
                            DuoCounter = DuoCounter + 1;
                        end
                        triples(TrpCounter , :) = chunks{c};
                        tripletInfo(TrpCounter , 1) = 3;
                        tripletInfo(TrpCounter , 2) = 1;
                        TrpCounter = TrpCounter + 1;
                    case 4
                        for pcntr = 1:3
                            doubles(DuoCounter , :) = chunks{c}(pcntr : pcntr+1);
                            doubleInfo(DuoCounter , 1) = 4; % comes form quadruple chunk
                            doubleInfo(DuoCounter , 2) = pcntr; % to know if it's the first, second or third duo in the reiplet
                            DuoCounter = DuoCounter + 1;
                        end
                        for pcntr = 1:2
                            triples(TrpCounter , :) = chunks{c}(pcntr : pcntr+2);
                            tripletInfo(TrpCounter , 1) = 4;
                            tripletInfo(TrpCounter , 2)  = pcntr;
                            TrpCounter = TrpCounter + 1;
                        end
                        quadrpl(QudCounter , :) = chunks{c};
                        quadrplInfo(QudCounter , 1) = 4;
                        quadrplInfo(QudCounter , 2) = 1;
                        QudCounter = QudCounter + 1;
                end
            end
            
            out.DoubInducedError{sub} = [];
            out.TripInducedError{sub} = [];
            out.QuadInducedError{sub} = [];
            out.NotInducedByChunk(sub) = 0;
            out.InducedByChunk(sub) = 0;
            ANARand = getrow(Dall , Dall.seqNumb == 0 & ismember(Dall.SN , sub) & Dall.isgood & Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , [2:4]));
            out.SumofError(sub) = sum(sum(ANARand.badPress));
            out.WholeChunkErr(sub , :) = [0 0 0];
            for tn = 1:length(ANARand.TN)
                for CL = 2:4
                    BP = find(ANARand.badPress(tn , CL:end))+(CL-1);
                    if ~isempty(BP)
                        for bp = BP
                            switch CL
                                case 2
                                    if ismember(ANARand.AllPress(tn , bp-(CL-1):bp) , doubles , 'rows')
                                        ind = ismember(doubles , ANARand.AllPress(tn , bp-(CL-1):bp) , 'rows');
                                        out.DoubInducedError{sub} = [out.DoubInducedError{sub} ; doubleInfo(ind, :)];
                                        out.InducedByChunk(sub) = out.InducedByChunk(sub)  +1;
                                    else
                                        out.NotInducedByChunk(sub) = out.NotInducedByChunk(sub)  +1;
                                    end
                                    
                                case 3
                                    if ismember(ANARand.AllPress(tn , bp-(CL-1):bp) , triples , 'rows')
                                        ind = ismember(triples , ANARand.AllPress(tn , bp-(CL-1):bp) , 'rows');
                                        out.TripInducedError{sub} = [out.TripInducedError{sub} ; tripletInfo(ind, :)];
                                    end
                                case 4
                                    if ismember(ANARand.AllPress(tn , bp-(CL-1):bp) , quadrpl , 'rows')
                                        ind = ismember(quadrpl , ANARand.AllPress(tn , bp-(CL-1):bp) , 'rows');
                                        out.QuadInducedError{sub} = [out.QuadInducedError{sub} ; quadrplInfo(ind, :)];
                                    end
                            end
                        end
                    end
                end
            end
            
            if length(out.DoubInducedError{sub})>0
                out.WholeChunkErr(sub,1) = out.WholeChunkErr(sub,1) + sum(out.DoubInducedError{sub}(:,1) == 2);
            end
            if length(out.TripInducedError{sub})>0
                out.WholeChunkErr(sub,2) = out.WholeChunkErr(sub,2) + sum(out.TripInducedError{sub}(:,1) == 3);
            end
            if length(out.QuadInducedError{sub})>0
                out.WholeChunkErr(sub,3) =  out.WholeChunkErr(sub,3) + sum(out.QuadInducedError{sub}(:,1) == 4);
            end
        end
        out.InducedByChunk = 100*(out.InducedByChunk./out.SumofError);
        out.NotInducedByChunk = 100*(out.NotInducedByChunk./out.SumofError);
        figure('color' , 'white')
        errorbar([nanmean(out.InducedByChunk) nanmean(out.NotInducedByChunk)] , [nanstd(out.InducedByChunk) nanstd(out.NotInducedByChunk)] , 'LineWidth' , 3)
        grid
        ylabel('percent')
        title('Errors in Random Sequences Induces by Chunks vs not Induced by Chunks')
        hold on
        ax = gca;
        ax.XTick = [1:2];
        ax.XTickLabel = {'Induced By Chunk' , 'Not Induced By Chunk'};
        ax.XTickLabelRotation = 45;
    case 'Pupil'
        % only look at the intermixed blocks
        ANA1 = getrow(Dall , Dall.BT == 2 & ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [1:6])  & ~Dall.isError &  ismember(Dall.Rep , rep));
        ANA0 = getrow(Dall , Dall.BT == 2 & ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0])  & ~Dall.isError &  ismember(Dall.Rep , rep));
        for tn = 1:length(ANA1.TN)
            ANA1.Pupil{tn}(ANA1.Pupil{tn}<= .95 *nanmedian(ANA1.Pupil{tn})) = NaN;
            ANA1.meanPup(tn,1) = nanmean(ANA1.Pupil{tn});
            idx = (ANA1.state{tn} == 7);
            l = ceil(length(find(idx))/10);
            idx(1:end - l) = 0;
            ANA1.EndFixPup(tn , 1) = nanmean(ANA1.Pupil{tn}(idx));
        end
        for tn = 1:length(ANA0.TN)
            ANA0.Pupil{tn}(ANA0.Pupil{tn}<= .95 *nanmedian(ANA0.Pupil{tn})) = NaN;
            ANA0.meanPup(tn,1) = nanmean(ANA0.Pupil{tn});
            idx = (ANA0.state{tn} == 7);
            l = ceil(length(find(idx))/10);
            idx(1:end - l) = 0;
            ANA0.EndFixPup(tn , 1) = nanmean(ANA0.Pupil{tn}(idx));
        end
        
        out.BlPupEffect = anovaMixed([ANA1.meanPup ; ANA0.meanPup] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.meanPup ; 1+0*ANA0.meanPup],{'Random/Chunked'},'intercept',1,'subset' , ~isnan([0*ANA1.meanPup ; 1+0*ANA0.meanPup]))  ;
        out.TrPupEffect = anovaMixed([ANA1.EndFixPup ; ANA0.EndFixPup] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.EndFixPup ; 1+0*ANA0.EndFixPup],{'Random/Chunked'},'intercept',1,'subset' , ~isnan([0*ANA1.EndFixPup ; 1+0*ANA0.EndFixPup]))  ;
        BN = unique([ANA0.BN ; ANA1.BN]);
        allPup = [];
        allBN = [];
        D = [];
        for bn = 1:length(BN)
            l1 = sum(ANA1.BN == BN(bn) & ~isnan(ANA1.meanPup));
            l0 = sum(ANA0.BN == BN(bn) & ~isnan(ANA0.meanPup));
            temp = NaN*ones(max(l1 , l0) , 2);
            temp(1:l1,1) = ANA1.meanPup(ANA1.BN == BN(bn) & ~isnan(ANA1.meanPup));
            temp(1:l0,2) = ANA0.meanPup(ANA0.BN == BN(bn) & ~isnan(ANA0.meanPup));
            allPup = [allPup ; temp];
            allBN = [allBN ; BN(bn)*ones(length(temp) , 1)];
            D = [D ; unique(Dall.Day(Dall.BN == BN(bn)))];
        end
        %         [xcoord,ePLOT0,ERROR0] = lineplot(allBN, allPup ,'plotfcn','nanmean','leg' ,{'Random' , 'Chunked'},'markersize',3,'linewidth',3);
        h1 = figure('color' , 'white');
        [xcoord1,ePLOT1,ERROR1] = lineplot(allBN, allPup(:,1) ,'markersize',3,'linewidth',3 , 'subset' , ~isnan(allPup(:,1)));
        [xcoord0,ePLOT0,ERROR0] = lineplot(allBN, allPup(:,2) ,'markersize',3,'linewidth',3 , 'subset' , ~isnan(allPup(:,2)));
        close(h1);
        
        D = find(diff(D));
        figure('color' , 'white')
        errorbar(xcoord0,ePLOT0 , ERROR0 , 'LineWidth' , 3);
        hold on
        errorbar(xcoord1, ePLOT1 , ERROR1 , 'LineWidth' , 3);
        ax = gca;
        ax.FontSize = 20;
        for i = 1:length(D)
            line([BN(D(i)+1) , BN(D(i)+1)] , [min(min(ePLOT1)) max(max(ePLOT0))] ,  'LineWidth' , 3 , 'LineStyle' , ':' , 'color' , 'g');
        end
        legend({'Random' , 'Chunked'} , 'FontSize' , 20)
        title(['Pupil size during Random vs. chunked sequences , p = ' , num2str(out.BlPupEffect.eff(2).p)],'FontSize' , 20)
        xlabel('Training Blocks', 'FontSize' , 20)
        grid on
        
        
        
        BN = unique([ANA0.BN ; ANA1.BN]);
        allPup = [];
        allBN = [];
        D = [];
        for bn = 1:length(BN)
            l1 = sum(ANA1.BN == BN(bn) & ~isnan(ANA1.EndFixPup));
            l0 = sum(ANA0.BN == BN(bn) & ~isnan(ANA0.EndFixPup));
            temp = NaN*ones(max(l1 , l0) , 2);
            temp(1:l1,1) = ANA1.EndFixPup(ANA1.BN == BN(bn) & ~isnan(ANA1.EndFixPup));
            temp(1:l0,2) = ANA0.EndFixPup(ANA0.BN == BN(bn) & ~isnan(ANA0.EndFixPup));
            allPup = [allPup ; temp];
            allBN = [allBN ; BN(bn)*ones(length(temp) , 1)];
            D = [D ; unique(Dall.Day(Dall.BN == BN(bn)))];
        end
        %         [xcoord,ePLOT0,ERROR0] = lineplot(allBN, allPup ,'plotfcn','nanmean','leg' ,{'Random' , 'Chunked'},'markersize',3,'linewidth',3);
        h1 = figure('color' , 'white');
        [xcoord1,ePLOT1,ERROR1] = lineplot(allBN, allPup(:,1) ,'markersize',3,'linewidth',3 , 'subset' , ~isnan(allPup(:,1)));
        [xcoord0,ePLOT0,ERROR0] = lineplot(allBN, allPup(:,2) ,'markersize',3,'linewidth',3 , 'subset' , ~isnan(allPup(:,2)));
        close(h1);
        
        D = find(diff(D));
        figure('color' , 'white')
        errorbar(xcoord0,ePLOT0 , ERROR0 , 'LineWidth' , 3);
        hold on
        errorbar(xcoord1, ePLOT1 , ERROR1 , 'LineWidth' , 3);
        ax = gca;
        ax.FontSize = 20;
        for i = 1:length(D)
            line([BN(D(i)+1) , BN(D(i)+1)] , [min(min(ePLOT1)) max(max(ePLOT0))] ,  'LineWidth' , 3 , 'LineStyle' , ':' , 'color' , 'g');
        end
        legend({'Random' , 'Chunked'} , 'FontSize' , 20)
        title(['Pupil size at the end of Random vs. chunked sequences , p = ' , num2str(out.TrPupEffect.eff(2).p)],'FontSize' , 20)
        xlabel('Training Blocks', 'FontSize' , 20)
        grid on
    case 'Points'
        for tn = 1:length(Dall.TN)
            Dall.MT(tn , 1) = Dall.AllPressTimes(tn , Dall.seqlength(tn)) - Dall.AllPressTimes(tn , 1);
        end
        points_rate = [];
        MT = [];
        MTall = [];
        h = 13; % Horizon
        for sub = subjnum
            ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [1:6]) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
            points_rate = [points_rate ; [100 * sum(ANA.points == 3)/length(ANA.TN) ,...
                100 * sum(ANA.points == 2)/length(ANA.TN) ,...
                100 * sum(ANA.points == 1)/length(ANA.TN) ,...
                100 * sum(ANA.points == 0)/length(ANA.TN), h , sub ,1]];
            MTall = [MTall ; [ANA.MT ANA.Horizon]];
            MT = [MT ; [nanmean(ANA.MT(ANA.points == 3)) ,...
                nanmean(ANA.MT(ANA.points == 2)) ,...
                nanmean(ANA.MT(ANA.points == 1)) ,...
                nanmean(ANA.MT(ANA.points == 0)) , h , sub ,1]];
            ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [0]) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
            points_rate = [points_rate ; [100 * sum(ANA.points == 3)/length(ANA.TN) ,...
                100 * sum(ANA.points == 2)/length(ANA.TN) ,...
                100 * sum(ANA.points == 1)/length(ANA.TN) ,...
                100 * sum(ANA.points == 0)/length(ANA.TN) , h , sub ,0]];
            MT = [MT ; [nanmean(ANA.MT(ANA.points == 3)) ,...
                nanmean(ANA.MT(ANA.points == 2)) ,...
                nanmean(ANA.MT(ANA.points == 1)) ,...
                nanmean(ANA.MT(ANA.points == 0)) , h , sub ,0 ]];
            MTall = [MTall ; [ANA.MT ANA.Horizon]];
        end
        MTall = MTall(MTall(:,1) <= 9000 , :);
        
        P3_temp = pivottable(points_rate(:,5) ,points_rate(:,7), points_rate(:,1) , 'nanmedian');
        P2_temp = pivottable(points_rate(:,5) ,points_rate(:,7), points_rate(:,2) , 'nanmedian');
        P1_temp = pivottable(points_rate(:,5) ,points_rate(:,7), points_rate(:,3) , 'nanmedian');
        E_temp  = pivottable(points_rate(:,5) ,points_rate(:,7), points_rate(:,4) , 'nanmedian');
        
        MT3_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,1) , 'nanmedian');
        MT2_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,2) , 'nanmedian');
        MT1_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,3) , 'nanmedian');
        MTe_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,4) , 'nanmedian');
        
        out.Err = anovan(points_rate(:,4) , points_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        out.Po3 = anovan(points_rate(:,1) , points_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        out.Po2 = anovan(points_rate(:,2) , points_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        out.Po1 = anovan(points_rate(:,3) , points_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        
        figure('color' , 'white')
        for sub = 1:12
            subplot(4,3,sub)
            b = bar([points_rate(points_rate(:,6) == sub,1:4)]);
            hold on
            ax = gca;
            title(['Subject ' ,num2str(sub)] )
            ax.XTickLabel = {'Chunked' , 'Random'};
        end
        
        figure('color' , 'white')
        barplot(points_rate(:,7) , points_rate(:,1:4));
        hold on
        ax = gca;
        title('Points rate')
        ax.XTickLabel = {'Random 3' ,'Random 2' ,'Random 1' ,'Random 0' , 'Chunked 3', 'Chunked 2', 'Chunked 1', 'Chunked 0'};
        ax.XTickLabelRotation = 45;
    case 'transitions_All'
        %         prompt = 'Which Horizon?';
        %         h = input(prompt);
        lastIPI = 0;
        load([baseDir , '/CMB.mat'])
        if calc
            for subjnum = 1:length(subj_name)-1
                Dall.isWrong = Dall.AllPress ~=Dall.AllResponse;
                ANA_allh = getrow(Dall , ismember(Dall.SN , subjnum));
                ANA1_allh = getrow(Dall ,ismember(Dall.seqNumb , [1:6]) & ismember(Dall.SN , subjnum));
                ANA0_allh = getrow(Dall ,ismember(Dall.seqNumb , 0) & ismember(Dall.SN , subjnum));
                
                
                
                allT2 = [length(ANA1_allh.AllPress)*(size(ANA1_allh.AllPress , 2)-1) length(ANA0_allh.AllPress)*(size(ANA0_allh.AllPress , 2)-1) ...
                    sum(ismember(ANA_allh.seqNumb , [0:6]))*(size(ANA_allh.AllPress , 2)-1) + sum(ismember(ANA_allh.seqNumb , [102 202])) + sum(ismember(ANA_allh.seqNumb , [103 203 ]))*2 + sum(ismember(ANA_allh.seqNumb , [104 204 ]))*3];
                t2_Nums_allh(subjnum).Chunked = zeros(length(CMB.comb2) , 1);
                t2_Nums_allh(subjnum).Rand    = zeros(length(CMB.comb2) , 1);
                t2_Nums_allh(subjnum).All     = zeros(length(CMB.comb2) , 1);
                
                ANA1_allh.t2_Nums = zeros(length(ANA1_allh.AllPress) , size(ANA1_allh.AllPress , 2) -1);
                ANA0_allh.t2_Nums = zeros(length(ANA0_allh.AllPress) , size(ANA0_allh.AllPress , 2) -1);
                ANA_allh.t2_Nums = zeros(length(ANA_allh.AllPress) , size(ANA_allh.AllPress , 2) -1);
                for t2 = 1:length(CMB.comb2)
                    
                    t2_Nums_allh(subjnum).TranNumb(t2 , 1) = t2;
                    t2_Nums_allh(subjnum).Transition(t2 , 1:2) = CMB.comb2(t2,:);
                    for p = 1:size(ANA1_allh.AllPress , 2) -1
                        t2_Nums_allh(subjnum).Chunked(t2,1) =  t2_Nums_allh(subjnum).Chunked(t2,1) + sum(ismember(ANA1_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                        t2_Nums_allh(subjnum).Rand(t2,1) =  t2_Nums_allh(subjnum).Rand(t2,1) + sum(ismember(ANA0_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                        t2_Nums_allh(subjnum).All(t2,1) =  t2_Nums_allh(subjnum).All(t2,1) + sum(ismember(ANA_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                        
                        ANA1_allh.t2_Nums(ismember(ANA1_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                        ANA0_allh.t2_Nums(ismember(ANA0_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                        ANA_allh.t2_Nums(ismember(ANA_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                    end
                end
                t2_Nums_allh(subjnum).Chunked = t2_Nums_allh(subjnum).Chunked/allT2(1);
                t2_Nums_allh(subjnum).Rand    = t2_Nums_allh(subjnum).Rand/allT2(2);
                t2_Nums_allh(subjnum).All    = t2_Nums_allh(subjnum).All/allT2(3);
                
                for p1 = 1:5
                    for p2 = 1:5
                        i = find(ismember(CMB.comb2 , [p1 p2] , 'rows'));
                        PoC2(p1 , p2) = t2_Nums_allh(subjnum).All(i);
                    end
                    PoC2_n(p1, :) = PoC2(p1, :)/sum(PoC2(p1, :)); % sets the sum of every row to 1
                    for p2 = 1:5
                        i = find(ismember(CMB.comb2 , [p1 p2] , 'rows'));
                        t2_Nums_allh(subjnum).All_normalized(i) = PoC2_n(p1, p2);
                    end
                end
                
                allT3 = [length(ANA1_allh.AllPress)*(size(ANA1_allh.AllPress , 2)-2) length(ANA0_allh.AllPress)*(size(ANA0_allh.AllPress , 2)-2) ...
                    sum(ismember(ANA_allh.seqNumb , [0:6]))*(size(ANA_allh.AllPress , 2)-2) + sum(ismember(ANA_allh.seqNumb , [103 203 ])) + sum(ismember(ANA_allh.seqNumb , [104 204 ]))*2];
                t3_Nums_allh(subjnum).Chunked = zeros(length(CMB.comb3) , 1);
                t3_Nums_allh(subjnum).Rand    = zeros(length(CMB.comb3) , 1);
                t3_Nums_allh(subjnum).All    = zeros(length(CMB.comb3) , 1);
                
                ANA1_allh.t3_Nums = zeros(length(ANA1_allh.AllPress) , size(ANA1_allh.AllPress , 2) -2);
                ANA0_allh.t3_Nums = zeros(length(ANA0_allh.AllPress) , size(ANA0_allh.AllPress , 2) -2);
                ANA_allh.t3_Nums = zeros(length(ANA_allh.AllPress) , size(ANA_allh.AllPress , 2) -2);
                for t3 = 1:length(CMB.comb3)
                    
                    
                    t3_Nums_allh(subjnum).TranNumb(t3 , 1) = t3;
                    t3_Nums_allh(subjnum).Transition(t3 , :) = CMB.comb3(t3,:);
                    for p = 1:size(ANA1_allh.AllPress , 2) -2
                        t3_Nums_allh(subjnum).Chunked(t3,1) =  t3_Nums_allh(subjnum).Chunked(t3,1) + sum(ismember(ANA1_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                        t3_Nums_allh(subjnum).Rand(t3,1) =  t3_Nums_allh(subjnum).Rand(t3,1) + sum(ismember(ANA0_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                        t3_Nums_allh(subjnum).All(t3,1) =  t3_Nums_allh(subjnum).All(t3,1) + sum(ismember(ANA_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                        
                        ANA1_allh.t3_Nums(ismember(ANA1_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                        ANA0_allh.t3_Nums(ismember(ANA0_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                        ANA_allh.t3_Nums(ismember(ANA_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                    end
                end
                t3_Nums_allh(subjnum).Chunked = t3_Nums_allh(subjnum).Chunked/allT3(1);
                t3_Nums_allh(subjnum).Rand = t3_Nums_allh(subjnum).Rand/allT3(2);
                t3_Nums_allh(subjnum).All = t3_Nums_allh(subjnum).All/allT3(3);
                for p12 = 1:25
                    for p3 = 1:5
                        i = find(ismember(CMB.comb3 , [CMB.comb2(p12 , :) , p3] , 'rows'));
                        PoC3(p12 , p3) = t3_Nums_allh(subjnum).All(i);
                    end
                    PoC3_n(p12, :) = PoC3(p12, :)/sum(PoC3(p12, :)); % sets the sum of every row to 1
                    PoC3_n(isnan(PoC3_n)) = 0;
                    for p3 = 1:5
                        i = find(ismember(CMB.comb3 , [CMB.comb2(p12 , :) p3] , 'rows'));
                        t3_Nums_allh(subjnum).All_normalized(i) = PoC3_n(p12, p3);
                    end
                end
                
                allT4 = [length(ANA1_allh.AllPress)*(size(ANA1_allh.AllPress , 2)-3) length(ANA0_allh.AllPress)*(size(ANA0_allh.AllPress , 2)-3) ...
                    sum(ismember(ANA_allh.seqNumb , [0:6]))*(size(ANA_allh.AllPress , 2)-3)+sum(ismember(ANA_allh.seqNumb , [104 204]))];
                t4_Nums_allh(subjnum).Chunked = zeros(length(CMB.comb4) , 1);
                t4_Nums_allh(subjnum).Rand    = zeros(length(CMB.comb4) , 1);
                t4_Nums_allh(subjnum).All    = zeros(length(CMB.comb4) , 1);
                
                ANA1_allh.t4_Nums = zeros(length(ANA1_allh.AllPress) , size(ANA1_allh.AllPress , 2) -3);
                ANA0_allh.t4_Nums = zeros(length(ANA0_allh.AllPress) , size(ANA0_allh.AllPress , 2) -3);
                ANA_allh.t4_Nums = zeros(length(ANA_allh.AllPress) , size(ANA_allh.AllPress , 2) -3);
                for t4 = 1:length(CMB.comb4)
                    
                    t4_Nums_allh(subjnum).TranNumb(t4 , 1) = t4;
                    t4_Nums_allh(subjnum).Transition(t4 , :) = CMB.comb4(t4,:);
                    for p = 1:size(ANA1_allh.AllPress , 2) -3
                        t4_Nums_allh(subjnum).Chunked(t4,1) =  t4_Nums_allh(subjnum).Chunked(t4,1) + sum(ismember(ANA1_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                        t4_Nums_allh(subjnum).Rand(t4,1) =  t4_Nums_allh(subjnum).Rand(t4,1) + sum(ismember(ANA0_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                        t4_Nums_allh(subjnum).All(t4,1) =  t4_Nums_allh(subjnum).All(t4,1) + sum(ismember(ANA_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                        
                        ANA1_allh.t4_Nums(ismember(ANA1_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                        ANA0_allh.t4_Nums(ismember(ANA0_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                        ANA_allh.t4_Nums(ismember(ANA_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                    end
                end
                t4_Nums_allh(subjnum).Chunked = t4_Nums_allh(subjnum).Chunked/allT4(1);
                t4_Nums_allh(subjnum).Rand = t4_Nums_allh(subjnum).Rand/allT4(2);
                t4_Nums_allh(subjnum).All = t4_Nums_allh(subjnum).All/allT4(3);
                for p123 = 1:125
                    for p4 = 1:5
                        i = find(ismember(CMB.comb4 , [CMB.comb3(p123 , :) p4] , 'rows'));
                        PoC4(p123 , p4) = t4_Nums_allh(subjnum).All(i);
                    end
                    PoC4_n(p123, :) = PoC4(p123, :)/sum(PoC4(p123, :)); % sets the sum of every row to 1
                    PoC4_n(isnan(PoC4_n)) = 0;
                    for p4 = 1:5
                        i = find(ismember(CMB.comb4 , [CMB.comb3(p123 , :) p4] , 'rows'));
                        t4_Nums_allh(subjnum).All_normalized(i) = PoC4_n(p123, p4);
                    end
                end
                
            end
            
            
            
            %%
            for subjnum  = 1:length(subj_name)-1
                id = ANA_allh.seqlength <=4;
                
                A = [ANA_allh.AllPress(id , :) ANA_allh.seqlength(id)];
                A(isnan(A)) = 0;
                A = unique(A , 'rows');
                for sl = 2:4
                    id = A(:,end) == sl;
                    CMB.Chunks{sl} = A(id,1:4);
                end
                
                ANA1 = getrow(Dall ,ismember(Dall.seqNumb , [1:6]) & ismember(Dall.SN , subjnum));
                ANA0 = getrow(Dall ,ismember(Dall.seqNumb , 0) & ismember(Dall.SN , subjnum));
                ANA = getrow(Dall ,ismember(Dall.SN , subjnum));
                
                allT2 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-1) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-1) length(ANA.AllPress)*(size(ANA.AllPress , 2)-1)];
                t2_Nums(subjnum).Chunked = zeros(length(CMB.comb2) , 1);
                t2_Nums(subjnum).Rand    = zeros(length(CMB.comb2) , 1);
                t2_Nums(subjnum).All    = zeros(length(CMB.comb2) , 1);
                
                ANA1.t2_Nums = zeros(length(ANA1.AllPress) , size(ANA1.AllPress , 2) -1);
                ANA0.t2_Nums = zeros(length(ANA0.AllPress) , size(ANA0.AllPress , 2) -1);
                ANA.t2_Nums = zeros(length(ANA.AllPress) , size(ANA.AllPress , 2) -1);
                for t2 = 1:length(CMB.comb2)
                    t2_Nums(subjnum).Chunked_IPI{t2,1} = [];
                    t2_Nums(subjnum).Rand_IPI{t2,1} = [];
                    t2_Nums(subjnum).All_IPI{t2,1} = [];
                    
                    t2_Nums(subjnum).TranNumb(t2 , 1) = t2;
                    t2_Nums(subjnum).Transition(t2 , 1:2) = CMB.comb2(t2,:);
                    for p = 1:size(ANA1.AllPress , 2) -1
                        t2_Nums(subjnum).Chunked(t2,1) =  t2_Nums(subjnum).Chunked(t2,1) + sum(ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                        t2_Nums(subjnum).Rand(t2,1)    =  t2_Nums(subjnum).Rand(t2,1) + sum(ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                        t2_Nums(subjnum).All(t2,1)     =  t2_Nums(subjnum).Rand(t2,1) + sum(ismember(ANA.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                        
                        ANA1.t2_Nums(ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                        ANA0.t2_Nums(ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                        ANA.t2_Nums(ismember(ANA.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                        
                        CorID = ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') & ~sum(ANA1.isWrong(:,p:p+1) , 2);
                        t2_Nums(subjnum).Chunked_IPI{t2} = [t2_Nums(subjnum).Chunked_IPI{t2} ; [ANA1.IPI(CorID , p) t2*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                        
                        CorID = ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') & ~sum(ANA0.isWrong(:,p:p+1) , 2);
                        t2_Nums(subjnum).Rand_IPI{t2} = [t2_Nums(subjnum).Rand_IPI{t2} ; [ANA0.IPI(CorID , p) t2*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                        
                        CorID = ismember(ANA.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') & ~sum(ANA.isWrong(:,p:p+1) , 2);
                        t2_Nums(subjnum).All_IPI{t2} = [t2_Nums(subjnum).All_IPI{t2} ; [ANA.IPI(CorID , p) t2*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                    end
                    t2_Nums(subjnum).MeanChunked_IPI(t2,1) = nanmean(t2_Nums(subjnum).Chunked_IPI{t2}(:,1));
                    t2_Nums(subjnum).MeanRand_IPI(t2,1) = nanmean(t2_Nums(subjnum).Rand_IPI{t2}(:,1));
                    t2_Nums(subjnum).MeanAll_IPI(t2,1) = nanmean(t2_Nums(subjnum).All_IPI{t2}(:,1));
                end
                
                
                t2_Nums(subjnum).Chunked = t2_Nums_allh(subjnum).Chunked;%t2_Nums(h).Chunked/allT2(1);
                t2_Nums(subjnum).Rand = t2_Nums_allh(subjnum).Rand;%t2_Nums(h).Rand/allT2(2);
                t2_Nums(subjnum).All = t2_Nums_allh(subjnum).All;%t2_Nums(h).Rand/allT2(2);
                [~ , t2_Nums(subjnum).sort_ID] = sort(t2_Nums(subjnum).All , 'descend');
                t2_Nums(subjnum).All_normalized = t2_Nums_allh(subjnum).All_normalized;%t2_Nums(h).Rand/allT2(2);
                [~ , t2_Nums(subjnum).sort_norm_ID] = sort(t2_Nums(subjnum).All_normalized , 'descend');
                
                
                
                t2_Nums(subjnum).SortedMeanIPI_Chunked = t2_Nums(subjnum).MeanChunked_IPI(t2_Nums(subjnum).sort_ID);
                t2_Nums(subjnum).SortedIPI_Chunked = t2_Nums(subjnum).Chunked_IPI(t2_Nums(subjnum).sort_ID);
                t2_Nums(subjnum).ReadyToPlot_Chunked = cell2mat(t2_Nums(subjnum).SortedIPI_Chunked);
                clear xtick_r xtick_c xticklab_r xticklab_c
                counter = 1;
                counter_n = 1;
                ChunkNum2 = [];
                ChunkNum2_34 = [];
                for i = 1:length(CMB.comb2)
                    idd = t2_Nums(subjnum).ReadyToPlot_Chunked(:,2) == t2_Nums(subjnum).sort_ID(i);
                    if sum(ismember(CMB.Chunks{2}(:,1:2) , CMB.comb2(i,:) , 'rows'))
                        ChunkNum2 = [ChunkNum2 ; [i , find(t2_Nums(subjnum).sort_ID == i)]];
                    end
                    if sum(ismember(CMB.Chunks{3}(:,1:2) , CMB.comb2(i,:) , 'rows')) | sum(ismember(CMB.Chunks{3}(:,2:3) , CMB.comb2(i,:) , 'rows'))
                        ChunkNum2_34 = [ChunkNum2_34 ; [i , find(t2_Nums(subjnum).sort_ID == i)]];
                    end
                    if sum(ismember(CMB.Chunks{4}(:,1:2) , CMB.comb2(i,:) , 'rows')) | sum(ismember(CMB.Chunks{4}(:,2:3) , CMB.comb2(i,:) , 'rows')) | sum(ismember(CMB.Chunks{4}(:,3:4) , CMB.comb2(i,:) , 'rows'))
                        ChunkNum2_34 = [ChunkNum2_34 ; [i , find(t2_Nums(subjnum).sort_ID == i)]];
                    end
                    if sum(idd)
                        T2(subjnum).xticklab_c{counter} = num2str(unique(t2_Nums(subjnum).ReadyToPlot_Chunked(idd,2)));
                        t2_Nums(subjnum).ReadyToPlot_Chunked(idd,7) = t2_Nums(subjnum).All(t2_Nums(subjnum).sort_ID(i));
                        T2(subjnum).xtick_c(counter) = i;
                        counter = counter + 1;
                    end
                    t2_Nums(subjnum).ReadyToPlot_Chunked(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    
                    idd_norm = t2_Nums(subjnum).ReadyToPlot_Chunked(:,2) == t2_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd_norm)
                        T2(subjnum).xticklab_a{counter_n} = num2str(unique(t2_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,2)));
                        t2_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,9) = t2_Nums(subjnum).All_normalized(t2_Nums(subjnum).sort_norm_ID(i));
                        T2(subjnum).xtick_a(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t2_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                end
                
                t2_Nums(subjnum).SortedMeanIPI_Rand = t2_Nums(subjnum).MeanRand_IPI(t2_Nums(subjnum).sort_ID);
                t2_Nums(subjnum).SortedIPI_Rand = t2_Nums(subjnum).Rand_IPI(t2_Nums(subjnum).sort_ID);
                t2_Nums(subjnum).ReadyToPlot_Rand = cell2mat(t2_Nums(subjnum).SortedIPI_Rand);
                counter = 1;
                counter_n = 1;
                for i = 1:length(CMB.comb2)
                    idd = t2_Nums(subjnum).ReadyToPlot_Rand(:,2) == t2_Nums(subjnum).sort_ID(i);
                    if sum(idd)
                        T2(subjnum).xticklab_r{counter} = num2str(unique(t2_Nums(subjnum).ReadyToPlot_Rand(idd,2)));
                        t2_Nums(subjnum).ReadyToPlot_Rand(idd,7) = t2_Nums(subjnum).All(t2_Nums(subjnum).sort_ID(i));
                        T2(subjnum).xtick_r(counter) = i;
                        counter = counter + 1;
                    end
                    t2_Nums(subjnum).ReadyToPlot_Rand(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    
                    idd_norm = t2_Nums(subjnum).ReadyToPlot_Rand(:,2) == t2_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd_norm)
                        T2(subjnum).xticklab_a{counter_n} = num2str(unique(t2_Nums(subjnum).ReadyToPlot_Rand(idd_norm,2)));
                        t2_Nums(subjnum).ReadyToPlot_Rand(idd_norm,9) = t2_Nums(subjnum).All_normalized(t2_Nums(subjnum).sort_norm_ID(i));
                        T2(subjnum).xtick_a(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t2_Nums(subjnum).ReadyToPlot_Rand(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                    
                end
                
                t2_Nums(subjnum).SortedMeanIPI_All= t2_Nums(subjnum).MeanAll_IPI(t2_Nums(subjnum).sort_ID);
                t2_Nums(subjnum).SortedIPI_All = t2_Nums(subjnum).All_IPI(t2_Nums(subjnum).sort_ID);
                t2_Nums(subjnum).ReadyToPlot_All = cell2mat(t2_Nums(subjnum).SortedIPI_All);
                counter = 1;
                counter_n = 1;
                for i = 1:length(CMB.comb2)
                    idd = t2_Nums(subjnum).ReadyToPlot_All(:,2) == t2_Nums(subjnum).sort_ID(i);
                    if sum(idd)
                        T2(subjnum).xticklab_a{counter} = num2str(unique(t2_Nums(subjnum).ReadyToPlot_All(idd,2)));
                        t2_Nums(subjnum).ReadyToPlot_All(idd,7) = t2_Nums(subjnum).All(t2_Nums(subjnum).sort_ID(i));
                        T2(subjnum).xtick_a(counter) = i;
                        counter = counter + 1;
                    end
                    t2_Nums(subjnum).ReadyToPlot_All(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    
                    idd_norm = t2_Nums(subjnum).ReadyToPlot_All(:,2) == t2_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd_norm)
                        T2(subjnum).xticklab_a{counter_n} = num2str(unique(t2_Nums(subjnum).ReadyToPlot_All(idd_norm,2)));
                        t2_Nums(subjnum).ReadyToPlot_All(idd_norm,9) = t2_Nums(subjnum).All_normalized(t2_Nums(subjnum).sort_norm_ID(i));
                        T2(subjnum).xtick_a(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t2_Nums(subjnum).ReadyToPlot_All(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                end
                
                
                h1 = figure;
                hold on
                subjnum
                [xcoordC_2{subjnum},PLOTC_2{subjnum},ERRORC_2{subjnum}]  = lineplot(t2_Nums(subjnum).ReadyToPlot_Chunked(:,6) , t2_Nums(subjnum).ReadyToPlot_Chunked(:,1), 'subset' , ismember(t2_Nums(subjnum).ReadyToPlot_Chunked(:,4) , [2 4]));
                [xcoordR_2{subjnum},PLOTR_2{subjnum},ERRORR_2{subjnum}]  = lineplot(t2_Nums(subjnum).ReadyToPlot_Rand(:,6) , t2_Nums(subjnum).ReadyToPlot_Rand(:,1), 'subset' , ismember(t2_Nums(subjnum).ReadyToPlot_Rand(:,4) , [2 4]));
                [xcoordA_2{subjnum},PLOTA_2{subjnum},ERRORA_2{subjnum}]  = lineplot(t2_Nums(subjnum).ReadyToPlot_All(:,6) , t2_Nums(subjnum).ReadyToPlot_All(:,1), 'subset' , ismember(t2_Nums(subjnum).ReadyToPlot_All(:,4) , [2 4]));
                close(h1)
                
                temp = corrcoef(t2_Nums(subjnum).ReadyToPlot_Chunked(:,1) , t2_Nums(subjnum).ReadyToPlot_Chunked(:,7));
                C2_chunked(subjnum) = temp(2);
                temp  = corrcoef(t2_Nums(subjnum).ReadyToPlot_Rand(:,1) , t2_Nums(subjnum).ReadyToPlot_Rand(:,7));
                C2_random(subjnum) = temp(2);
                temp  = corrcoef(t2_Nums(subjnum).ReadyToPlot_All(:,1) , t2_Nums(subjnum).ReadyToPlot_All(:,7));
                C2_all(subjnum) = temp(2);
                
                
                
                
                clear xtick_r xtick_c xticklab_r xticklab_c
                allT3 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-2) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-2) length(ANA.AllPress)*(size(ANA.AllPress , 2)-2)];
                t3_Nums(subjnum).Chunked = zeros(length(CMB.comb3) , 1);
                t3_Nums(subjnum).Rand    = zeros(length(CMB.comb3) , 1);
                t3_Nums(subjnum).All     = zeros(length(CMB.comb3) , 1);
                
                ANA1.t3_Nums = zeros(length(ANA1.AllPress) , size(ANA1.AllPress , 2) -2);
                ANA0.t3_Nums = zeros(length(ANA0.AllPress) , size(ANA0.AllPress , 2) -2);
                ANA.t3_Nums = zeros(length(ANA.AllPress) , size(ANA.AllPress , 2) -2);
                for t3 = 1:length(CMB.comb3)
                    t3_Nums(subjnum).Chunked_IPI{t3,1} = [];
                    t3_Nums(subjnum).Rand_IPI{t3,1} = [];
                    t3_Nums(subjnum).All_IPI{t3,1} = [];
                    
                    t3_Nums(subjnum).TranNumb(t3 , 1) = t3;
                    t3_Nums(subjnum).Transition(t3 , :) = CMB.comb3(t3,:);
                    for p = 1:size(ANA1.AllPress , 2) -2
                        t3_Nums(subjnum).Chunked(t3,1) =  t3_Nums(subjnum).Chunked(t3,1) + sum(ismember(ANA1.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                        t3_Nums(subjnum).Rand(t3,1) =  t3_Nums(subjnum).Rand(t3,1) + sum(ismember(ANA0.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                        t3_Nums(subjnum).All(t3,1) =  t3_Nums(subjnum).All(t3,1) + sum(ismember(ANA.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                        
                        ANA1.t3_Nums(ismember(ANA1.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                        ANA0.t3_Nums(ismember(ANA0.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                        ANA.t3_Nums(ismember(ANA.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                        
                        CorID = ismember(ANA1.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') & ~sum(ANA1.isWrong(:,p:p+2) , 2);
                        if ~lastIPI
                            t3_Nums(subjnum).Chunked_IPI{t3} = [t3_Nums(subjnum).Chunked_IPI{t3} ; [sum(ANA1.IPI(CorID , p:p+1),2) t3*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                        else
                            t3_Nums(subjnum).Chunked_IPI{t3} = [t3_Nums(subjnum).Chunked_IPI{t3} ; [sum(ANA1.IPI(CorID , p+1),2) t3*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                        end
                        
                        CorID = ismember(ANA0.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') & ~sum(ANA0.isWrong(:,p:p+2) , 2);
                        if ~lastIPI
                            t3_Nums(subjnum).Rand_IPI{t3} = [t3_Nums(subjnum).Rand_IPI{t3} ; [sum(ANA0.IPI(CorID , p:p+1),2) t3*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                        else
                            t3_Nums(subjnum).Rand_IPI{t3} = [t3_Nums(subjnum).Rand_IPI{t3} ; [sum(ANA0.IPI(CorID , p+1),2) t3*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                        end
                        
                        CorID = ismember(ANA.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') & ~sum(ANA.isWrong(:,p:p+2) , 2);
                        if ~lastIPI
                            t3_Nums(subjnum).All_IPI{t3} = [t3_Nums(subjnum).All_IPI{t3} ; [sum(ANA.IPI(CorID , p:p+1),2) t3*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                        else
                            t3_Nums(subjnum).All_IPI{t3} = [t3_Nums(subjnum).All_IPI{t3} ; [sum(ANA.IPI(CorID , p+1),2) t3*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                        end
                    end
                    t3_Nums(subjnum).MeanChunked_IPI(t3,1) = nanmean(t3_Nums(subjnum).Chunked_IPI{t3}(:,1));
                    t3_Nums(subjnum).MeanRand_IPI(t3,1) = nanmean(t3_Nums(subjnum).Rand_IPI{t3}(:,1));
                    t3_Nums(subjnum).MeanAll_IPI(t3,1) = nanmean(t3_Nums(subjnum).All_IPI{t3}(:,1));
                end
                t3_Nums(subjnum).Chunked = t3_Nums_allh(subjnum).Chunked;%t3_Nums(h).Chunked/allT3(1);
                t3_Nums(subjnum).Rand = t3_Nums_allh(subjnum).Rand;%t3_Nums(h).Rand/allT3(2);
                t3_Nums(subjnum).All = t3_Nums_allh(subjnum).All;%t3_Nums(h).Rand/allT3(2);
                [~ , t3_Nums(subjnum).sort_ID] = sort(t3_Nums(subjnum).All , 'descend');
                t3_Nums(subjnum).All_normalized = t3_Nums_allh(subjnum).All_normalized;%t3_Nums(h).Rand/allT3(2);
                [~ , t3_Nums(subjnum).sort_norm_ID] = sort(t3_Nums(subjnum).All_normalized , 'descend');
                
                t3_Nums(subjnum).SortedMeanIPI_Chunked = t3_Nums(subjnum).MeanChunked_IPI(t3_Nums(subjnum).sort_ID);
                t3_Nums(subjnum).SortedIPI_Chunked = t3_Nums(subjnum).Chunked_IPI(t3_Nums(subjnum).sort_ID);
                t3_Nums(subjnum).ReadyToPlot_Chunked = cell2mat(t3_Nums(subjnum).SortedIPI_Chunked);
                counter = 1;
                counter_n = 1;
                ChunkNum3 = [];
                ChunkNum3_4 = [];
                for i = 1:length(CMB.comb3)
                    
                    idd = t3_Nums(subjnum).ReadyToPlot_Chunked(:,2) == t3_Nums(subjnum).sort_ID(i);
                    if sum(ismember(CMB.Chunks{3}(:,1:3) , CMB.comb3(i,:) , 'rows'))
                        ChunkNum3 = [ChunkNum3 ; [i , find(t3_Nums(subjnum).sort_ID == i)]];
                    end
                    if sum(ismember(CMB.Chunks{4}(:,1:3) , CMB.comb3(i,:) , 'rows')) | sum(ismember(CMB.Chunks{4}(:,2:4) , CMB.comb3(i,:) , 'rows'))
                        ChunkNum3_4 = [ChunkNum3_4 ; [i , find(t3_Nums(subjnum).sort_ID == i)]];
                    end
                    if sum(idd)
                        T3(subjnum).xticklab_c{counter} = num2str(unique(t3_Nums(subjnum).ReadyToPlot_Chunked(idd,2)));
                        t3_Nums(subjnum).ReadyToPlot_Chunked(idd,7) = t3_Nums(subjnum).All(t3_Nums(subjnum).sort_ID(i));
                        T3(subjnum).xtick_c(counter) = i;
                        counter = counter + 1;
                    end
                    t3_Nums(subjnum).ReadyToPlot_Chunked(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    
                    idd_norm = t3_Nums(subjnum).ReadyToPlot_Chunked(:,2) == t3_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd_norm)
                        T3(subjnum).xticklab_a{counter_n} = num2str(unique(t3_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,2)));
                        t3_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,9) = t3_Nums(subjnum).All_normalized(t3_Nums(subjnum).sort_norm_ID(i));
                        T3(subjnum).xtick_a(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t3_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                    
                end
                
                
                t3_Nums(subjnum).SortedMeanIPI_Rand = t3_Nums(subjnum).MeanRand_IPI(t3_Nums(subjnum).sort_ID);
                t3_Nums(subjnum).SortedIPI_Rand = t3_Nums(subjnum).Rand_IPI(t3_Nums(subjnum).sort_ID);
                t3_Nums(subjnum).ReadyToPlot_Rand = cell2mat(t3_Nums(subjnum).SortedIPI_Rand);
                counter = 1;
                counter_n = 1;
                for i = 1:length(CMB.comb3)
                    idd = t3_Nums(subjnum).ReadyToPlot_Rand(:,2) == t3_Nums(subjnum).sort_ID(i);
                    if sum(idd)
                        T3(subjnum).xticklab_r{counter} = num2str(unique(t3_Nums(subjnum).ReadyToPlot_Rand(idd,2)));
                        t3_Nums(subjnum).ReadyToPlot_Rand(idd,7) = t3_Nums(subjnum).All(t3_Nums(subjnum).sort_ID(i));
                        T3(subjnum).xtick_r(counter) = i;
                        counter = counter + 1;
                    end
                    t3_Nums(subjnum).ReadyToPlot_Rand(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    
                    idd_norm = t3_Nums(subjnum).ReadyToPlot_Rand(:,2) == t3_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd_norm)
                        T3(subjnum).xticklab_a{counter_n} = num2str(unique(t3_Nums(subjnum).ReadyToPlot_Rand(idd_norm,2)));
                        t3_Nums(subjnum).ReadyToPlot_Rand(idd_norm,9) = t3_Nums(subjnum).All_normalized(t3_Nums(subjnum).sort_norm_ID(i));
                        T3(subjnum).xtick_a(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t3_Nums(subjnum).ReadyToPlot_Rand(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                    
                end
                
                
                t3_Nums(subjnum).SortedMeanIPI_All = t3_Nums(subjnum).MeanAll_IPI(t3_Nums(subjnum).sort_ID);
                t3_Nums(subjnum).SortedIPI_All = t3_Nums(subjnum).All_IPI(t3_Nums(subjnum).sort_ID);
                t3_Nums(subjnum).ReadyToPlot_All = cell2mat(t3_Nums(subjnum).SortedIPI_All);
                counter = 1;
                counter_n = 1;
                for i = 1:length(CMB.comb3)
                    idd = t3_Nums(subjnum).ReadyToPlot_All(:,2) == t3_Nums(subjnum).sort_ID(i);
                    if sum(idd)
                        T3(subjnum).xticklab_a{counter} = num2str(unique(t3_Nums(subjnum).ReadyToPlot_All(idd,2)));
                        t3_Nums(subjnum).ReadyToPlot_All(idd,7) = t3_Nums(subjnum).All(t3_Nums(subjnum).sort_ID(i));
                        T3(subjnum).xtick_a(counter) = i;
                        counter = counter + 1;
                    end
                    t3_Nums(subjnum).ReadyToPlot_All(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    
                    idd_norm = t3_Nums(subjnum).ReadyToPlot_All(:,2) == t3_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd_norm)
                        T3(subjnum).xticklab_a{counter_n} = num2str(unique(t3_Nums(subjnum).ReadyToPlot_All(idd_norm,2)));
                        t3_Nums(subjnum).ReadyToPlot_All(idd_norm,9) = t3_Nums(subjnum).All_normalized(t3_Nums(subjnum).sort_norm_ID(i));
                        T3(subjnum).xtick_a(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t3_Nums(subjnum).ReadyToPlot_All(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                end
                
                
                h1 = figure;
                hold on
                [xcoordC_3{subjnum},PLOTC_3{subjnum},ERRORC_3{subjnum}]  = lineplot(t3_Nums(subjnum).ReadyToPlot_Chunked(:,6) , t3_Nums(subjnum).ReadyToPlot_Chunked(:,1) , 'subset' , ismember(t3_Nums(subjnum).ReadyToPlot_Chunked(:,4) , [2 4]));
                [xcoordR_3{subjnum},PLOTR_3{subjnum},ERRORR_3{subjnum}]  = lineplot(t3_Nums(subjnum).ReadyToPlot_Rand(:,6) , t3_Nums(subjnum).ReadyToPlot_Rand(:,1) , 'subset' , ismember(t3_Nums(subjnum).ReadyToPlot_Rand(:,4) , [2 4]));
                [xcoordA_3{subjnum},PLOTA_3{subjnum},ERRORA_3{subjnum}]  = lineplot(t3_Nums(subjnum).ReadyToPlot_All(:,6) , t3_Nums(subjnum).ReadyToPlot_All(:,1) , 'subset' , ismember(t3_Nums(subjnum).ReadyToPlot_All(:,4) , [2 4]));
                close(h1)
                
                %
                %             anovan(t3_Nums(h).ReadyToPlot_Rand(:,1) , t3_Nums(h).ReadyToPlot_Rand(:,2))
                temp = corrcoef(t3_Nums(subjnum).ReadyToPlot_Chunked(:,1) , t3_Nums(subjnum).ReadyToPlot_Chunked(:,7));
                C3_chunked(subjnum) = temp(2);
                temp  = corrcoef(t3_Nums(subjnum).ReadyToPlot_Rand(:,1) , t3_Nums(subjnum).ReadyToPlot_Rand(:,7));
                C3_random(subjnum) = temp(2);
                temp  = corrcoef(t3_Nums(subjnum).ReadyToPlot_All(:,1) , t3_Nums(subjnum).ReadyToPlot_All(:,7));
                C3_all(subjnum) = temp(2);
                
                
                
                clear xtick_r xtick_c xticklab_r xticklab_c
                allt4 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-3) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-3) length(ANA.AllPress)*(size(ANA.AllPress , 2)-3)];
                t4_Nums(subjnum).Chunked = zeros(length(CMB.comb4) , 1);
                t4_Nums(subjnum).Rand    = zeros(length(CMB.comb4) , 1);
                t4_Nums(subjnum).All     = zeros(length(CMB.comb4) , 1);
                
                ANA1.t4_Nums = zeros(length(ANA1.AllPress) , size(ANA1.AllPress , 2) -3);
                ANA0.t4_Nums = zeros(length(ANA0.AllPress) , size(ANA0.AllPress , 2) -3);
                ANA.t4_Nums  = zeros(length(ANA.AllPress) , size(ANA.AllPress , 2) -3);
                for t4 = 1:length(CMB.comb4)
                    t4_Nums(subjnum).Chunked_IPI{t4,1} = [];
                    t4_Nums(subjnum).Rand_IPI{t4,1} = [];
                    t4_Nums(subjnum).All_IPI{t4,1} = [];
                    
                    t4_Nums(subjnum).TranNumb(t4 , 1) = t4;
                    t4_Nums(subjnum).Transition(t4 , :) = CMB.comb4(t4,:);
                    for p = 1:size(ANA1.AllPress , 2) -3
                        t4_Nums(subjnum).Chunked(t4,1) =  t4_Nums(subjnum).Chunked(t4,1) + sum(ismember(ANA1.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                        t4_Nums(subjnum).Rand(t4,1) =  t4_Nums(subjnum).Rand(t4,1) + sum(ismember(ANA0.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                        t4_Nums(subjnum).All(t4,1) =  t4_Nums(subjnum).All(t4,1) + sum(ismember(ANA.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                        
                        ANA1.t4_Nums(ismember(ANA1.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                        ANA0.t4_Nums(ismember(ANA0.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                        ANA.t4_Nums(ismember(ANA.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                        
                        CorID = ismember(ANA1.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') & ~sum(ANA1.isWrong(:,p:p+3) , 2);
                        if ~lastIPI
                            t4_Nums(subjnum).Chunked_IPI{t4} = [t4_Nums(subjnum).Chunked_IPI{t4} ; [sum(ANA1.IPI(CorID , p:p+2),2) t4*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                        else
                            t4_Nums(subjnum).Chunked_IPI{t4} = [t4_Nums(subjnum).Chunked_IPI{t4} ; [sum(ANA1.IPI(CorID , p+2),2) t4*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                        end
                        
                        CorID = ismember(ANA0.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') & ~sum(ANA0.isWrong(:,p:p+3) , 2);
                        if ~lastIPI
                            t4_Nums(subjnum).Rand_IPI{t4} = [t4_Nums(subjnum).Rand_IPI{t4} ; [sum(ANA0.IPI(CorID , p:p+2),2) t4*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                        else
                            t4_Nums(subjnum).Rand_IPI{t4} = [t4_Nums(subjnum).Rand_IPI{t4} ; [sum(ANA0.IPI(CorID , p+2),2) t4*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                        end
                        
                        CorID = ismember(ANA.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') & ~sum(ANA.isWrong(:,p:p+3) , 2);
                        if ~lastIPI
                            t4_Nums(subjnum).All_IPI{t4} = [t4_Nums(subjnum).All_IPI{t4} ; [sum(ANA.IPI(CorID , p:p+2),2) t4*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                        else
                            t4_Nums(subjnum).All_IPI{t4} = [t4_Nums(subjnum).All_IPI{t4} ; [sum(ANA.IPI(CorID , p+2),2) t4*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                        end
                    end
                    t4_Nums(subjnum).MeanChunked_IPI(t4,1) = nanmean(t4_Nums(subjnum).Chunked_IPI{t4}(:,1));
                    t4_Nums(subjnum).MeanRand_IPI(t4,1) = nanmean(t4_Nums(subjnum).Rand_IPI{t4}(:,1));
                    t4_Nums(subjnum).MeanAll_IPI(t4,1) = nanmean(t4_Nums(subjnum).All_IPI{t4}(:,1));
                end
                t4_Nums(subjnum).Chunked = t4_Nums_allh(subjnum).Chunked;%t4_Nums(h).Chunked/allt4(1);
                t4_Nums(subjnum).Rand = t4_Nums_allh(subjnum).Rand;%t4_Nums(h).Rand/allt4(2);
                t4_Nums(subjnum).All = t4_Nums_allh(subjnum).All;%t4_Nums(h).Rand/allt4(2);
                t4_Nums(subjnum).All_normalized = t4_Nums_allh(subjnum).All_normalized;%t4_Nums(h).Rand/allt4(2);
                [~ , t4_Nums(subjnum).sort_ID] = sort(t4_Nums(subjnum).All , 'descend');
                [~ , t4_Nums(subjnum).sort_norm_ID] = sort(t4_Nums(subjnum).All_normalized , 'descend');
                
                t4_Nums(subjnum).SortedMeanIPI_Chunked = t4_Nums(subjnum).MeanChunked_IPI(t4_Nums(subjnum).sort_ID);
                t4_Nums(subjnum).SortedIPI_Chunked = t4_Nums(subjnum).Chunked_IPI(t4_Nums(subjnum).sort_ID);
                t4_Nums(subjnum).ReadyToPlot_Chunked = cell2mat(t4_Nums(subjnum).SortedIPI_Chunked);
                counter = 1;
                counter_n = 1;
                ChunkNum4 = [];
                for i = 1:length(CMB.comb4)
                    idd = t4_Nums(subjnum).ReadyToPlot_Chunked(:,2) == t4_Nums(subjnum).sort_ID(i);
                    
                    if sum(ismember(CMB.Chunks{4} , CMB.comb4(i,:)  , 'rows'))
                        ChunkNum4 = [ChunkNum4 ; [i , find(t4_Nums(subjnum).sort_ID == i)]];
                    end
                    
                    if sum(idd)
                        T4(subjnum).xticklab_c{counter} = num2str(unique(t4_Nums(subjnum).ReadyToPlot_Chunked(idd,2)));
                        t4_Nums(subjnum).ReadyToPlot_Chunked(idd,7) = t4_Nums(subjnum).All(t4_Nums(subjnum).sort_ID(i));
                        T4(subjnum).xtick_c(counter) = i;
                        counter = counter + 1;
                    end
                    idd_norm = t4_Nums(subjnum).ReadyToPlot_Chunked(:,2) == t4_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd_norm)
                        T4(subjnum).xticklab_c{counter_n} = num2str(unique(t4_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,2)));
                        t4_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,9) = t4_Nums(subjnum).All_normalized(t4_Nums(subjnum).sort_norm_ID(i));
                        T4(subjnum).xtick_c(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t4_Nums(subjnum).ReadyToPlot_Chunked(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    t4_Nums(subjnum).ReadyToPlot_Chunked(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                    
                end
                
                t4_Nums(subjnum).SortedMeanIPI_Rand = t4_Nums(subjnum).MeanRand_IPI(t4_Nums(subjnum).sort_ID);
                t4_Nums(subjnum).SortedIPI_Rand = t4_Nums(subjnum).Rand_IPI(t4_Nums(subjnum).sort_ID);
                t4_Nums(subjnum).ReadyToPlot_Rand = cell2mat(t4_Nums(subjnum).SortedIPI_Rand);
                counter = 1;
                counter_n = 1;
                for i = 1:length(CMB.comb4)
                    idd = t4_Nums(subjnum).ReadyToPlot_Rand(:,2) == t4_Nums(subjnum).sort_ID(i);
                    idd_norm = t4_Nums(subjnum).ReadyToPlot_Rand(:,2) == t4_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd)
                        T4(subjnum).xticklab_r{counter} = num2str(unique(t4_Nums(subjnum).ReadyToPlot_Rand(idd,2)));
                        t4_Nums(subjnum).ReadyToPlot_Rand(idd,7) = t4_Nums(subjnum).All(t4_Nums(subjnum).sort_ID(i));
                        T4(subjnum).xtick_r(counter) = i;
                        counter = counter + 1;
                    end
                    if sum(idd_norm)
                        T4(subjnum).xticklab_r{counter_n} = num2str(unique(t4_Nums(subjnum).ReadyToPlot_Rand(idd_norm,2)));
                        t4_Nums(subjnum).ReadyToPlot_Rand(idd_norm,9) = t4_Nums(subjnum).All_normalized(t4_Nums(subjnum).sort_norm_ID(i));
                        T4(subjnum).xtick_r(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t4_Nums(subjnum).ReadyToPlot_Rand(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    t4_Nums(subjnum).ReadyToPlot_Rand(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                end
                
                t4_Nums(subjnum).SortedMeanIPI_All = t4_Nums(subjnum).MeanAll_IPI(t4_Nums(subjnum).sort_ID);
                t4_Nums(subjnum).SortedIPI_All = t4_Nums(subjnum).All_IPI(t4_Nums(subjnum).sort_ID);
                t4_Nums(subjnum).ReadyToPlot_All = cell2mat(t4_Nums(subjnum).SortedIPI_All);
                counter = 1;
                counter_n = 1;
                for i = 1:length(CMB.comb4)
                    idd = t4_Nums(subjnum).ReadyToPlot_All(:,2) == t4_Nums(subjnum).sort_ID(i);
                    idd_norm = t4_Nums(subjnum).ReadyToPlot_All(:,2) == t4_Nums(subjnum).sort_norm_ID(i);
                    if sum(idd)
                        T4(subjnum).xticklab_a{counter} = num2str(unique(t4_Nums(subjnum).ReadyToPlot_All(idd,2)));
                        t4_Nums(subjnum).ReadyToPlot_All(idd,7) = t4_Nums(subjnum).All(t4_Nums(subjnum).sort_ID(i));
                        T4(subjnum).xtick_a(counter) = i;
                        counter = counter + 1;
                    end
                    if sum(idd_norm)
                        T4(subjnum).xticklab_a{counter_n} = num2str(unique(t4_Nums(subjnum).ReadyToPlot_All(idd,2)));
                        t4_Nums(subjnum).ReadyToPlot_All(idd_norm,9) = t4_Nums(subjnum).All_normalized(t4_Nums(subjnum).sort_norm_ID(i));
                        T4(subjnum).xtick_a(counter_n) = i;
                        counter_n = counter_n + 1;
                    end
                    t4_Nums(subjnum).ReadyToPlot_All(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
                    t4_Nums(subjnum).ReadyToPlot_All(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
                    
                end
                
                h1 = figure;
                hold on
                [xcoordC_4{subjnum},PLOTC_4{subjnum},ERRORC_4{subjnum}]  = lineplot(t4_Nums(subjnum).ReadyToPlot_Chunked(:,6) , t4_Nums(subjnum).ReadyToPlot_Chunked(:,1) , 'subset' , ismember(t4_Nums(subjnum).ReadyToPlot_Chunked(:,4) , [2 4]));
                [xcoordR_4{subjnum},PLOTR_4{subjnum},ERRORR_4{subjnum}]  = lineplot(t4_Nums(subjnum).ReadyToPlot_Rand(:,6) , t4_Nums(subjnum).ReadyToPlot_Rand(:,1) , 'subset' , ismember(t4_Nums(subjnum).ReadyToPlot_Rand(:,4) , [2 4]));
                [xcoordA_4{subjnum},PLOTA_4{subjnum},ERRORA_4{subjnum}]  = lineplot(t4_Nums(subjnum).ReadyToPlot_All(:,6) , t4_Nums(subjnum).ReadyToPlot_All(:,1) , 'subset' , ismember(t4_Nums(subjnum).ReadyToPlot_All(:,4) , [2 4]));
                close(h1)
                
                %
                %             anovan(t4_Nums(h).ReadyToPlot_Rand(:,1) , t4_Nums(h).ReadyToPlot_Rand(:,2))
                temp = corrcoef(t4_Nums(subjnum).ReadyToPlot_Chunked(:,1) , t4_Nums(subjnum).ReadyToPlot_Chunked(:,7));
                C4_chunked(subjnum) = temp(2);
                temp  = corrcoef(t4_Nums(subjnum).ReadyToPlot_Rand(:,1) , t4_Nums(subjnum).ReadyToPlot_Rand(:,7));
                C4_random(subjnum) = temp(2);
                temp  = corrcoef(t4_Nums(subjnum).ReadyToPlot_All(:,1) , t4_Nums(subjnum).ReadyToPlot_All(:,7));
                C4_all(subjnum) = temp(2);
                
            end
            
            %%
            
            
            C.IPI = [];
            C.IPI_norm = [];
            C.t2 = [];
            C.SN =[];
            C.BN =[];
            C.Day = [];
            C.Horizon  = [];
            C.IPIarrangement = [];
            C.estIPIarrangement = [];
            C.t2Rank = [];
            C.t2Prob = [];
            C.t2Rank_n = [];
            C.t2Prob_n = [];
            C.t3 = [];
            C.t4 = [];
            C.T3Rank = [];
            C.t3Prob = [];
            C.T4Rank = [];
            C.t4Prob = [];
            C.T3Rank_n = [];
            C.t3Prob_n = [];
            C.T4Rank_n = [];
            C.t4Prob_n = [];
            R = C;
            All = C;
            Dall.estChnkBndry(Dall.estChnkBndry == 0) = 2;
            for subjnum = 1:length(subj_name)-1
                ANA1 = getrow(Dall ,ismember(Dall.seqNumb , [1:6]) & ismember(Dall.SN , subjnum) & ~Dall.isError);
                for tn = 1:length(ANA1.TN)
                    n = (ANA1.AllPressIdx(tn , sum(~isnan(ANA1.AllPressIdx(tn , :))))  - ANA1.AllPressIdx(tn , 1)) / 1000;
                    nIdx(tn , :) = (ANA1.AllPressIdx(tn , :) - ANA1.AllPressIdx(tn , 1))/n;
                    ANA1.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
                end
                
                
                [a,d]=find(isnan(ANA1.IPI) | (ANA1.IPI> 2000));
                tns = ones(length(ANA1.TN),1);
                tns(a) = 0;
                ANA1 = getrow(ANA1 , logical(tns));
                ANA0 = getrow(Dall ,ismember(Dall.seqNumb , 0) & ismember(Dall.SN , subjnum)  & ~Dall.isError);
                for tn = 1:length(ANA0.TN)
                    n = (ANA0.AllPressIdx(tn , sum(~isnan(ANA0.AllPressIdx(tn , :))))  - ANA0.AllPressIdx(tn , 1)) / 1000;
                    nIdx(tn , :) = (ANA0.AllPressIdx(tn , :) - ANA0.AllPressIdx(tn , 1))/n;
                    ANA0.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
                end
                
                [a,d]=find(isnan(ANA0.IPI) | (ANA0.IPI> 2000));
                tns = ones(length(ANA0.TN),1);
                tns(a) = 0;
                ANA0 = getrow(ANA0 , logical(tns));
                ANA  = addstruct(ANA1 , ANA0);
                
                allT2 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-1) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-1) length(ANA.AllPress)*(size(ANA.AllPress , 2)-1)];
                
                [t2sort_ID_n(:,1) , t2sort_ID_n(:,2)] = sort(t2_Nums_allh(subjnum).All_normalized , 'descend');
                [t3sort_ID_n(:,1) , t3sort_ID_n(:,2)] = sort(t3_Nums_allh(subjnum).All_normalized , 'descend');
                [t4sort_ID_n(:,1) , t4sort_ID_n(:,2)] = sort(t4_Nums_allh(subjnum).All_normalized , 'descend');
                [t2sort_ID(:,1) , t2sort_ID(:,2)] = sort(t2_Nums_allh(subjnum).All , 'descend');
                [t3sort_ID(:,1) , t3sort_ID(:,2)] = sort(t3_Nums_allh(subjnum).All , 'descend');
                [t4sort_ID(:,1) , t4sort_ID(:,2)] = sort(t4_Nums_allh(subjnum).All , 'descend');
                
                for p = 3:11
                    for t2 = 1:length(CMB.comb2)
                        
                        t2id = find(t2sort_ID(:,2) == t2);
                        t2id_n = find(t2sort_ID_n(:,2) == t2);
                        
                        clear t3  t4 T3Rank T4Rank t4Prob t3Prob T3Rank_n T4Rank_n t4Prob_n t3Prob_n
                        CorID = ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows');
                        tempID = find(CorID);
                        if ~isempty(tempID)
                            C.IPI = [C.IPI ; ANA1.IPI(CorID , p)];
                            C.IPI_norm = [C.IPI_norm ; ANA1.IPI_norm(CorID , p)];
                            C.t2 = [C.t2 ; t2*ones(length(ANA1.IPI(CorID , p)) , 1)];
                            C.SN = [C.SN ; ANA1.SN(CorID)];
                            C.BN = [C.BN ; ANA1.BN(CorID)];
                            C.Day = [C.Day ;ANA1.Day(CorID)];
                            C.Horizon  = [C.Horizon ;ANA1.Horizon(CorID)];
                            C.IPIarrangement = [C.IPIarrangement ;ANA1.IPIarrangement(CorID , p)];
                            C.estIPIarrangement = [C.estIPIarrangement ;ANA1.estChnkBndry(CorID , p)];
                            C.t2Rank = [C.t2Rank ; t2id*ones(length(ANA1.IPI(CorID , p)) , 1)];
                            C.t2Prob = [C.t2Prob ; t2sort_ID(t2id,1)*ones(length(ANA1.IPI(CorID , p)) , 1)];
                            
                            C.t2Rank_n = [C.t2Rank_n ; t2id_n*ones(length(ANA1.IPI(CorID , p)) , 1)];
                            C.t2Prob_n = [C.t2Prob_n ; t2sort_ID_n(t2id_n,1)*ones(length(ANA1.IPI(CorID , p)) , 1)];
                            
                            for jj = 1:length(tempID)
                                t3(jj,:) = [find(ismember(CMB.comb3 , ANA1.AllPress(tempID(jj),p-1:p+1), 'rows'))   ,   find(ismember(CMB.comb3 , ANA1.AllPress(tempID(jj),p:p+2) ,'rows'))];
                                T3Rank(jj , :) = [find(t3sort_ID(:,2) == t3(jj,1)) , find(t3sort_ID(:,2) == t3(jj,2))];
                                t3Prob(jj,:) = [t3sort_ID(T3Rank(jj , 1),1) , t3sort_ID(T3Rank(jj , 2),1)];
                                
                                T3Rank_n(jj , :) = [find(t3sort_ID_n(:,2) == t3(jj,1)) , find(t3sort_ID_n(:,2) == t3(jj,2))];
                                t3Prob_n(jj,:) = [t3sort_ID_n(T3Rank_n(jj , 1),1) , t3sort_ID_n(T3Rank(jj , 2),1)];
                                
                                
                                t4(jj,:) = [find(ismember(CMB.comb4 , ANA1.AllPress(tempID(jj),p-2:p+1), 'rows'))   ,   find(ismember(CMB.comb4 , ANA1.AllPress(tempID(jj),p-1:p+2) ,'rows'))   ,   find(ismember(CMB.comb4 , ANA1.AllPress(tempID(jj),p:p+3) ,'rows'))];
                                T4Rank(jj , :) = [find(t4sort_ID(:,2) == t4(jj,1))  , find(t4sort_ID(:,2) == t4(jj,2))  , find(t4sort_ID(:,2) == t4(jj,3))];
                                t4Prob(jj,:) = [t4sort_ID(T4Rank(jj , 1),1) , t4sort_ID(T4Rank(jj , 2),1)  ,  t4sort_ID(T4Rank(jj , 3),1)];
                                
                                T4Rank_n(jj , :) = [find(t4sort_ID_n(:,2) == t4(jj,1))  , find(t4sort_ID_n(:,2) == t4(jj,2))  , find(t4sort_ID_n(:,2) == t4(jj,3))];
                                t4Prob_n(jj,:) = [t4sort_ID_n(T4Rank_n(jj , 1),1) , t4sort_ID_n(T4Rank_n(jj , 2),1)  ,  t4sort_ID_n(T4Rank_n(jj , 3),1)];
                            end
                            C.t3 = [C.t3 ; t3];
                            C.t4 = [C.t4 ; t4];
                            C.T3Rank = [C.T3Rank ; T3Rank];
                            C.t3Prob = [C.t3Prob ; t3Prob];
                            C.T4Rank = [C.T4Rank ; T4Rank];
                            C.t4Prob = [C.t4Prob ; t4Prob];
                            
                            C.T3Rank_n = [C.T3Rank_n ; T3Rank_n];
                            C.t3Prob_n = [C.t3Prob_n ; t3Prob_n];
                            C.T4Rank_n = [C.T4Rank_n ; T4Rank_n];
                            C.t4Prob_n = [C.t4Prob_n ; t4Prob_n];
                        end
                        
                        
                        
                        
                        clear t3  t4 T3Rank T4Rank t4Prob t3Prob T3Rank_n T4Rank_n t4Prob_n t3Prob_n
                        CorID = ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows');
                        tempID = find(CorID);
                        if ~isempty(tempID)
                            R.IPI = [R.IPI ; ANA0.IPI(CorID , p)];
                            R.IPI_norm = [R.IPI_norm ; ANA0.IPI_norm(CorID , p)];
                            R.t2 = [R.t2 ; t2*ones(length(ANA0.IPI(CorID , p)) , 1)];
                            R.SN = [R.SN ; ANA0.SN(CorID)];
                            R.BN  = [R.BN ; ANA0.BN(CorID)];
                            R.Day = [R.Day ;ANA0.Day(CorID)];
                            R.Horizon  = [R.Horizon; ANA0.Horizon(CorID)];
                            R.IPIarrangement = [R.IPIarrangement ;ANA0.IPIarrangement(CorID , p)];
                            R.estIPIarrangement = [R.estIPIarrangement ;ANA0.estChnkBndry(CorID , p)];
                            R.t2Rank = [R.t2Rank ; t2id*ones(length(ANA0.IPI(CorID , p)) , 1)];
                            R.t2Prob = [R.t2Prob ; t2sort_ID(t2id,1)*ones(length(ANA0.IPI(CorID , p)) , 1)];
                            
                            R.t2Rank_n = [R.t2Rank_n ; t2id_n*ones(length(ANA0.IPI(CorID , p)) , 1)];
                            R.t2Prob_n = [R.t2Prob_n ; t2sort_ID_n(t2id_n,1)*ones(length(ANA0.IPI(CorID , p)) , 1)];
                            for jj = 1:length(tempID)
                                t3(jj,:) = [find(ismember(CMB.comb3 , ANA0.AllPress(tempID(jj),p-1:p+1), 'rows'))   ,   find(ismember(CMB.comb3 , ANA0.AllPress(tempID(jj),p:p+2) ,'rows'))];
                                T3Rank(jj , :) = [find(t3sort_ID(:,2) == t3(jj,1)) , find(t3sort_ID(:,2) == t3(jj,2))];
                                t3Prob(jj,:) = [t3sort_ID(T3Rank(jj , 1),1) , t3sort_ID(T3Rank(jj , 2),1)];
                                
                                T3Rank_n(jj , :) = [find(t3sort_ID_n(:,2) == t3(jj,1)) , find(t3sort_ID_n(:,2) == t3(jj,2))];
                                t3Prob_n(jj,:) = [t3sort_ID_n(T3Rank_n(jj , 1),1) , t3sort_ID_n(T3Rank(jj , 2),1)];
                                
                                t4(jj,:) = [find(ismember(CMB.comb4 , ANA0.AllPress(tempID(jj),p-2:p+1), 'rows'))   ,   find(ismember(CMB.comb4 , ANA0.AllPress(tempID(jj),p-1:p+2) ,'rows'))   ,   find(ismember(CMB.comb4 , ANA0.AllPress(tempID(jj),p:p+3) ,'rows'))];
                                T4Rank(jj , :) = [find(t4sort_ID(:,2) == t4(jj,1))  , find(t4sort_ID(:,2) == t4(jj,2))  , find(t4sort_ID(:,2) == t4(jj,3))];
                                t4Prob(jj,:) = [t4sort_ID(T4Rank(jj , 1),1) , t4sort_ID(T4Rank(jj , 2),1)  ,  t4sort_ID(T4Rank(jj , 3),1)];
                                
                                T4Rank_n(jj , :) = [find(t4sort_ID_n(:,2) == t4(jj,1))  , find(t4sort_ID_n(:,2) == t4(jj,2))  , find(t4sort_ID_n(:,2) == t4(jj,3))];
                                t4Prob_n(jj,:) = [t4sort_ID_n(T4Rank_n(jj , 1),1) , t4sort_ID_n(T4Rank_n(jj , 2),1)  ,  t4sort_ID_n(T4Rank_n(jj , 3),1)];
                            end
                            R.t3 = [R.t3 ; t3];
                            R.t4 = [R.t4 ; t4];
                            R.T3Rank = [R.T3Rank ; T3Rank];
                            R.t3Prob = [R.t3Prob ; t3Prob];
                            R.T4Rank = [R.T4Rank ; T4Rank];
                            R.t4Prob = [R.t4Prob ; t4Prob];
                            
                            R.T3Rank_n = [R.T3Rank_n ; T3Rank_n];
                            R.t3Prob_n = [R.t3Prob_n ; t3Prob_n];
                            R.T4Rank_n = [R.T4Rank_n ; T4Rank_n];
                            R.t4Prob_n = [R.t4Prob_n ; t4Prob_n];
                        end
                        
                        
                        clear t3  t4 T3Rank T4Rank t4Prob t3Prob T3Rank_n T4Rank_n t4Prob_n t3Prob_n
                        CorID = ismember(ANA.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows');
                        tempID = find(CorID);
                        if ~isempty(tempID)
                            All.IPI = [All.IPI ; ANA.IPI(CorID , p)];
                            All.IPI_norm = [All.IPI_norm ; ANA.IPI_norm(CorID , p)];
                            All.t2 = [All.t2 ; t2*ones(length(ANA.IPI(CorID , p)) , 1)];
                            All.SN = [All.SN ; ANA.SN(CorID)];
                            All.BN = [All.BN ; ANA.BN(CorID)];
                            All.Day = [All.Day ;ANA.Day(CorID)];
                            All.Horizon  = [All.Horizon ;ANA.Horizon(CorID)];
                            All.IPIarrangement = [All.IPIarrangement ;ANA.IPIarrangement(CorID , p)];
                            All.estIPIarrangement = [All.estIPIarrangement ;ANA.estChnkBndry(CorID , p)];
                            All.t2Rank = [All.t2Rank ; t2id*ones(length(ANA.IPI(CorID , p)) , 1)];
                            All.t2Prob = [All.t2Prob ; t2sort_ID(t2id,1)*ones(length(ANA.IPI(CorID , p)) , 1)];
                            
                            All.t2Rank_n = [All.t2Rank_n ; t2id_n*ones(length(ANA.IPI(CorID , p)) , 1)];
                            All.t2Prob_n = [All.t2Prob_n ; t2sort_ID_n(t2id_n,1)*ones(length(ANA.IPI(CorID , p)) , 1)];
                            for jj = 1:length(tempID)
                                t3(jj,:) = [find(ismember(CMB.comb3 , ANA.AllPress(tempID(jj),p-1:p+1), 'rows'))   ,   find(ismember(CMB.comb3 , ANA.AllPress(tempID(jj),p:p+2) ,'rows'))];
                                T3Rank(jj , :) = [find(t3sort_ID(:,2) == t3(jj,1)) , find(t3sort_ID(:,2) == t3(jj,2))];
                                t3Prob(jj,:) = [t3sort_ID(T3Rank(jj , 1),1) , t3sort_ID(T3Rank(jj , 2),1)];
                                
                                T3Rank_n(jj , :) = [find(t3sort_ID_n(:,2) == t3(jj,1)) , find(t3sort_ID_n(:,2) == t3(jj,2))];
                                t3Prob_n(jj,:) = [t3sort_ID_n(T3Rank_n(jj , 1),1) , t3sort_ID_n(T3Rank(jj , 2),1)];
                                
                                t4(jj,:) = [find(ismember(CMB.comb4 , ANA.AllPress(tempID(jj),p-2:p+1), 'rows'))   ,   find(ismember(CMB.comb4 , ANA.AllPress(tempID(jj),p-1:p+2) ,'rows'))   ,   find(ismember(CMB.comb4 , ANA.AllPress(tempID(jj),p:p+3) ,'rows'))];
                                T4Rank(jj , :) = [find(t4sort_ID(:,2) == t4(jj,1))  , find(t4sort_ID(:,2) == t4(jj,2))  , find(t4sort_ID(:,2) == t4(jj,3))];
                                t4Prob(jj,:) = [t4sort_ID(T4Rank(jj , 1),1) , t4sort_ID(T4Rank(jj , 2),1)  ,  t4sort_ID(T4Rank(jj , 3),1)];
                                
                                T4Rank_n(jj , :) = [find(t4sort_ID_n(:,2) == t4(jj,1))  , find(t4sort_ID_n(:,2) == t4(jj,2))  , find(t4sort_ID_n(:,2) == t4(jj,3))];
                                t4Prob_n(jj,:) = [t4sort_ID_n(T4Rank_n(jj , 1),1) , t4sort_ID_n(T4Rank_n(jj , 2),1)  ,  t4sort_ID_n(T4Rank_n(jj , 3),1)];
                            end
                            All.t3 = [All.t3 ; t3];
                            All.t4 = [All.t4 ; t4];
                            All.T3Rank = [All.T3Rank ; T3Rank];
                            All.t3Prob = [All.t3Prob ; t3Prob];
                            All.T4Rank = [All.T4Rank ; T4Rank];
                            All.t4Prob = [All.t4Prob ; t4Prob];
                            
                            All.T3Rank_n = [All.T3Rank_n ; T3Rank_n];
                            All.t3Prob_n = [All.t3Prob_n ; t3Prob_n];
                            All.T4Rank_n = [All.T4Rank_n ; T4Rank_n];
                            All.t4Prob_n = [All.t4Prob_n ; t4Prob_n];
                        end
                    end
                    
                    
                end
            end
        else
            load([baseDir , '/se1_TranProb.mat']);
        end
        
        %%
        ReadytobePlotted_Chunked_AllSubj_2 = [];
        ReadytobePlotted_Chunked_AllSubj_3 = [];
        ReadytobePlotted_Chunked_AllSubj_4 = [];
        
        ReadytobePlotted_All_AllSubj_2 = [];
        ReadytobePlotted_All_AllSubj_3 = [];
        ReadytobePlotted_All_AllSubj_4 = [];
        
        ReadytobePlotted_Rand_AllSubj_2 = [];
        ReadytobePlotted_Rand_AllSubj_3 = [];
        ReadytobePlotted_Rand_AllSubj_4 = [];
        
        for subjnum = 1:length(subj_name) -1
            ReadytobePlotted_Chunked_AllSubj_2 = [ReadytobePlotted_Chunked_AllSubj_2 ; t2_Nums(subjnum).ReadyToPlot_Chunked];
            ReadytobePlotted_Chunked_AllSubj_3 = [ReadytobePlotted_Chunked_AllSubj_3 ; t3_Nums(subjnum).ReadyToPlot_Chunked];
            ReadytobePlotted_Chunked_AllSubj_4 = [ReadytobePlotted_Chunked_AllSubj_4 ; t4_Nums(subjnum).ReadyToPlot_Chunked];
            
            ReadytobePlotted_All_AllSubj_2 = [ReadytobePlotted_All_AllSubj_2 ; t2_Nums(subjnum).ReadyToPlot_All];
            ReadytobePlotted_All_AllSubj_3 = [ReadytobePlotted_All_AllSubj_3 ; t3_Nums(subjnum).ReadyToPlot_All];
            ReadytobePlotted_All_AllSubj_4 = [ReadytobePlotted_All_AllSubj_4 ; t4_Nums(subjnum).ReadyToPlot_All];
            
            ReadytobePlotted_Rand_AllSubj_2 = [ReadytobePlotted_Rand_AllSubj_2 ; t2_Nums(subjnum).ReadyToPlot_Rand];
            ReadytobePlotted_Rand_AllSubj_3 = [ReadytobePlotted_Rand_AllSubj_3 ; t3_Nums(subjnum).ReadyToPlot_Rand];
            ReadytobePlotted_Rand_AllSubj_4 = [ReadytobePlotted_Rand_AllSubj_4 ; t4_Nums(subjnum).ReadyToPlot_Rand];
        end
        %%
        % column order for ready to be plotted... :IPI/Transition number/Subject Number/Day/Horizon/non-normSort Rank/non-normProbability/normSort Rank/normProbability
        h1 = figure;
        hold on
        [XC2,PC2,EC2] = lineplot(ReadytobePlotted_Chunked_AllSubj_2(:,6) , ReadytobePlotted_Chunked_AllSubj_2(:,1));
        [XC3,PC3,EC3] = lineplot(ReadytobePlotted_Chunked_AllSubj_3(:,6) , ReadytobePlotted_Chunked_AllSubj_3(:,1));
        [XC4,PC4,EC4] = lineplot(ReadytobePlotted_Chunked_AllSubj_4(:,6) , ReadytobePlotted_Chunked_AllSubj_4(:,1));
        
        [XR2,PR2,ER2]= lineplot(ReadytobePlotted_Rand_AllSubj_2(:,6) , ReadytobePlotted_Rand_AllSubj_2(:,1));
        [XR3,PR3,ER3] = lineplot(ReadytobePlotted_Rand_AllSubj_3(:,6) , ReadytobePlotted_Rand_AllSubj_3(:,1));
        [XR4,PR4,ER4] = lineplot(ReadytobePlotted_Rand_AllSubj_4(:,6) , ReadytobePlotted_Rand_AllSubj_4(:,1));
        
        [XA2,PA2,EA2] = lineplot(ReadytobePlotted_All_AllSubj_2(:,6) , ReadytobePlotted_All_AllSubj_2(:,1));
        [XA3,PA3,EA3] = lineplot(ReadytobePlotted_All_AllSubj_3(:,6) , ReadytobePlotted_All_AllSubj_3(:,1));
        [XA4,PA4,EA4] = lineplot(ReadytobePlotted_All_AllSubj_4(:,6) , ReadytobePlotted_All_AllSubj_4(:,1));
        close(h1)
        
        
        h1 = figure;
        hold on
        [XC2,PC2,EC2] = lineplot(ReadytobePlotted_Chunked_AllSubj_2(:,6) , ReadytobePlotted_Chunked_AllSubj_2(:,1));
        [XC3,PC3,EC3] = lineplot(ReadytobePlotted_Chunked_AllSubj_3(:,6) , ReadytobePlotted_Chunked_AllSubj_3(:,1));
        [XC4,PC4,EC4] = lineplot(ReadytobePlotted_Chunked_AllSubj_4(:,6) , ReadytobePlotted_Chunked_AllSubj_4(:,1));
        
        [XR2,PR2,ER2]= lineplot(ReadytobePlotted_Rand_AllSubj_2(:,6) , ReadytobePlotted_Rand_AllSubj_2(:,1));
        [XR3,PR3,ER3] = lineplot(ReadytobePlotted_Rand_AllSubj_3(:,6) , ReadytobePlotted_Rand_AllSubj_3(:,1));
        [XR4,PR4,ER4] = lineplot(ReadytobePlotted_Rand_AllSubj_4(:,6) , ReadytobePlotted_Rand_AllSubj_4(:,1));
        
        [XA2,PA2,EA2] = lineplot(ReadytobePlotted_All_AllSubj_2(:,6) , ReadytobePlotted_All_AllSubj_2(:,1));
        [XA3,PA3,EA3] = lineplot(ReadytobePlotted_All_AllSubj_3(:,6) , ReadytobePlotted_All_AllSubj_3(:,1));
        [XA4,PA4,EA4] = lineplot(ReadytobePlotted_All_AllSubj_4(:,6) , ReadytobePlotted_All_AllSubj_4(:,1));
        close(h1)
        
        figCount = 1;
        figure('color' , 'white')
        subplot(3,3,figCount)
        errorbar(XC2,PC2,EC2, 'LineWidth' , 3);
        hold on
        grid on
        title('Double trans chunked seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        subplot(3,3,figCount)
        errorbar(XR2,PR2,ER2, 'LineWidth' , 3);
        hold on
        
        grid on
        title('Double trans rand seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        subplot(3,3,figCount)
        errorbar(XA2,PA2,EA2, 'LineWidth' , 3);
        
        grid on
        title('Double trans all seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        
        subplot(3,3,figCount)
        errorbar(XC3,PC3,EC3 , 'LineWidth' , 3);
        hold on
        grid on
        title('Triplet trans chunked seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        subplot(3,3,figCount)
        errorbar(XR3,PR3,ER3, 'LineWidth' , 3);
        hold on
        grid on
        title('Triplet trans rand seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        subplot(3,3,figCount)
        errorbar(XA3,PA3,EA3, 'LineWidth' , 3);
        hold on
        grid on
        title('Triplet trans all seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        
        
        subplot(3,3,figCount)
        errorbar(XC4,PC4,EC4 , 'LineWidth' , 3);
        hold on
        grid on
        title('Quad trans chunked seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        subplot(3,3,figCount)
        errorbar(XR4,PR4,ER4, 'LineWidth' , 3);
        hold on
        grid on
        title('Quad trans rand seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        subplot(3,3,figCount)
        errorbar(XA4,PA4,EA4, 'LineWidth' , 3);
        hold on
        grid on
        title('Quad trans all seqs sorted by non-normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        %%
        h1 = figure;
        hold on
        [XC2,PC2,EC2] = lineplot(ReadytobePlotted_Chunked_AllSubj_2(:,8) , ReadytobePlotted_Chunked_AllSubj_2(:,1));
        [XC3,PC3,EC3] = lineplot(ReadytobePlotted_Chunked_AllSubj_3(:,8) , ReadytobePlotted_Chunked_AllSubj_3(:,1));
        [XC4,PC4,EC4] = lineplot(ReadytobePlotted_Chunked_AllSubj_4(:,8) , ReadytobePlotted_Chunked_AllSubj_4(:,1));
        
        [XR2,PR2,ER2]= lineplot(ReadytobePlotted_Rand_AllSubj_2(:,8) , ReadytobePlotted_Rand_AllSubj_2(:,1));
        [XR3,PR3,ER3] = lineplot(ReadytobePlotted_Rand_AllSubj_3(:,8) , ReadytobePlotted_Rand_AllSubj_3(:,1));
        [XR4,PR4,ER4] = lineplot(ReadytobePlotted_Rand_AllSubj_4(:,8) , ReadytobePlotted_Rand_AllSubj_4(:,1));
        
        [XA2,PA2,EA2] = lineplot(ReadytobePlotted_All_AllSubj_2(:,8) , ReadytobePlotted_All_AllSubj_2(:,1));
        [XA3,PA3,EA3] = lineplot(ReadytobePlotted_All_AllSubj_3(:,8) , ReadytobePlotted_All_AllSubj_3(:,1));
        [XA4,PA4,EA4] = lineplot(ReadytobePlotted_All_AllSubj_4(:,8) , ReadytobePlotted_All_AllSubj_4(:,1));
        close(h1)
        
        figCount = 1;
        figure('color' , 'white')
        subplot(3,3,figCount)
        errorbar(XC2,PC2,EC2, 'LineWidth' , 3);
        hold on
        grid on
        title('Double trans chunked seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        subplot(3,3,figCount)
        errorbar(XR2,PR2,ER2, 'LineWidth' , 3);
        hold on
        
        grid on
        title('Double trans rand seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        subplot(3,3,figCount)
        errorbar(XA2,PA2,EA2, 'LineWidth' , 3);
        
        grid on
        title('Double trans all seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        
        subplot(3,3,figCount)
        errorbar(XC3,PC3,EC3 , 'LineWidth' , 3);
        hold on
        grid on
        title('Triplet trans chunked seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        subplot(3,3,figCount)
        errorbar(XR3,PR3,ER3, 'LineWidth' , 3);
        hold on
        grid on
        title('Triplet trans rand seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        subplot(3,3,figCount)
        errorbar(XA3,PA3,EA3, 'LineWidth' , 3);
        hold on
        grid on
        title('Triplet trans all seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        
        
        subplot(3,3,figCount)
        errorbar(XC4,PC4,EC4 , 'LineWidth' , 3);
        hold on
        grid on
        title('Quad trans chunked seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        subplot(3,3,figCount)
        errorbar(XR4,PR4,ER4, 'LineWidth' , 3);
        hold on
        grid on
        title('Quad trans rand seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        
        subplot(3,3,figCount)
        errorbar(XA4,PA4,EA4, 'LineWidth' , 3);
        hold on
        grid on
        title('Quad trans all seqs sorted by normP in all seqs','FontSize' , 14)
        figCount = figCount +1;
        
        %% %%%%%%%%%%%%%************************************* 1st order
        subjnum= 1;
        for p1 = 1:5
            for p2 = 1:5
                i = find(ismember(CMB.comb2 , [p1 p2] , 'rows'));
                PoC2(p1 , p2) = t2_Nums_allh(subjnum).All(i);
                for subjnum = 1
                    MT2_C{subjnum}(p1, p2) = t2_Nums(subjnum).MeanChunked_IPI(i);
                    MT2_R{subjnum}(p1, p2) = t2_Nums(subjnum).MeanRand_IPI(i);
                    MT2_A{subjnum}(p1, p2) = t2_Nums(subjnum).MeanAll_IPI(i);
                end
            end
            PoC2(p1, :) = PoC2(p1, :)/sum(PoC2(p1, :)); % sets the sum of every row to 1
        end
        figure('color' , 'white')
        subplot(3,4,1)
        imagesc(PoC2);
        ylabel('Press 1')
        xlabel('Press 2')
        axis square
        colorbar
        set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
            'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
        title('Probability of Occurence for fist-order transitions in All sequences')
        for j = 1:length(ChunkNum2)
            y = CMB.comb2(ChunkNum2(j,1) , 1);
            x = CMB.comb2(ChunkNum2(j,1) , 2);
            rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
        end
        for j = 1:length(ChunkNum2_34)
            y = CMB.comb2(ChunkNum2_34(j,1) , 1);
            x = CMB.comb2(ChunkNum2_34(j,1) , 2);
            rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
        end
        for subjnum = 1
            subplot(3,4,2)
            imagesc(MT2_C{subjnum} , [200 600]);
            ylabel('Press 1')
            xlabel('Press 2')
            axis square
            set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
                'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
            title(['Average Movement time (Chunked) H = ' , num2str(subjnum)])
            for j = 1:length(ChunkNum2)
                y = CMB.comb2(ChunkNum2(j,1) , 1);
                x = CMB.comb2(ChunkNum2(j,1) , 2);
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            colorbar
            
        end
        
        
        for subjnum = 1
            subplot(3,4,3)
            imagesc(MT2_R{subjnum} , [200 600]);
            ylabel('Press 1')
            xlabel('Press 2')
            axis square
            set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
                'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
            title(['Average Movement time (Random) H = ' , num2str(subjnum)])
            for j = 1:length(ChunkNum2)
                y = CMB.comb2(ChunkNum2(j,1) , 1);
                x = CMB.comb2(ChunkNum2(j,1) , 2);
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            colorbar
        end
        
        
        
        for subjnum = 1
            subplot(3,4,4)
            imagesc(MT2_A{subjnum} , [200 600]);
            ylabel('Press 1')
            xlabel('Press 2')
            axis square
            set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
                'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
            title(['Average Movement time (All) H = ' , num2str(subjnum)])
            for j = 1:length(ChunkNum2)
                y = CMB.comb2(ChunkNum2(j,1) , 1);
                x = CMB.comb2(ChunkNum2(j,1) , 2);
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            colorbar
        end
        
        %%%%%%%%%%%%%%%************************************* 2nd order
        
        for p1 = 1:5
            for p23 = 1:25
                i = find(ismember(CMB.comb3 , [p1 CMB.comb2(p23 , :)] , 'rows'));
                PoC3(p1 , p23) = t3_Nums_allh(subjnum).All(i);
                for subjnum = 1
                    MT3_C{subjnum}(p1, p23) = t3_Nums(subjnum).MeanChunked_IPI(i);
                    MT3_R{subjnum}(p1, p23) = t3_Nums(subjnum).MeanRand_IPI(i);
                    MT3_A{subjnum}(p1, p23) = t3_Nums(subjnum).MeanAll_IPI(i);
                end
            end
            PoC3(p1, :) = PoC3(p1, :)/sum(PoC3(p1, :)); % sets the sum of every row to 1
        end
        
        subplot(345)
        imagesc(PoC3);
        ylabel('Press 1')
        xlabel('Press 2 3')
        %         axis square
        colorbar
        set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
            'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
            'XTickLabelRotation' , 45)
        title('Probability of Occurence for 2nd-order transitions in All sequences')
        for j = 1:length(ChunkNum3)
            y = CMB.comb3(ChunkNum3(j,1) , 1);
            x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
            rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
        end
        for j = 1:length(ChunkNum3_4)
            y = CMB.comb3(ChunkNum3_4(j,1) , 1);
            x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
            rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
        end
        for subjnum = 1
            subplot(346)
            if ~lastIPI
                imagesc(MT3_C{subjnum} , [500 , 1300]); % for sum of IPIS
            else
                imagesc(MT3_C{subjnum},[200 600]);  % for the last IPIS
            end
            ylabel('Press 1')
            xlabel('Press 2 3')
            %             axis square
            set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
                'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
                'XTickLabelRotation' , 45)
            title(['Average Movement time (Chunked) H = ' , num2str(subjnum)])
            
            for j = 1:length(ChunkNum3)
                y = CMB.comb3(ChunkNum3(j,1) , 1);
                x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            for j = 1:length(ChunkNum3_4)
                y = CMB.comb3(ChunkNum3_4(j,1) , 1);
                x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
            end
            
            colorbar
        end
        
        
        
        for subjnum = 1
            subplot(347)
            if ~lastIPI
                imagesc(MT3_R{subjnum} , [500 , 1300]);% for sum of IPIS
            else
                imagesc(MT3_R{subjnum},[200 600]);  % for the last IPIS
            end
            ylabel('Press 1')
            xlabel('Press 2 3')
            %             axis square
            set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
                'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
                'XTickLabelRotation' , 45)
            title(['Average Movement time (Random) H = ' , num2str(subjnum)])
            
            for j = 1:length(ChunkNum3)
                y = CMB.comb3(ChunkNum3(j,1) , 1);
                x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            for j = 1:length(ChunkNum3_4)
                y = CMB.comb3(ChunkNum3_4(j,1) , 1);
                x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
            end
            
            colorbar
            
        end
        
        
        
        for subjnum = 1
            subplot(348)
            if ~lastIPI
                imagesc(MT3_A{subjnum} , [500 , 1300]);% for sum of IPIS
            else
                imagesc(MT3_A{subjnum},[200 600]);  % for the last IPIS
            end
            ylabel('Press 1')
            xlabel('Press 2 3')
            %             axis square
            set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
                'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
                'XTickLabelRotation' , 45)
            title(['Average Movement time (All) H = ' , num2str(subjnum)])
            
            for j = 1:length(ChunkNum3)
                y = CMB.comb3(ChunkNum3(j,1) , 1);
                x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            for j = 1:length(ChunkNum3_4)
                y = CMB.comb3(ChunkNum3_4(j,1) , 1);
                x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
            end
            
            colorbar
        end
        
        %%%%%%%%%%%%%%%************************************* 3rd order
        
        for p12 = 1:25
            for p34 = 1:25
                i = find(ismember(CMB.comb4 , [CMB.comb2(p12 , :) CMB.comb2(p34 , :)] , 'rows'));
                PoC4(p12 , p34) = t4_Nums_allh(subjnum).All(i);
                for subjnum = 1
                    MT4_C{subjnum}(p12, p34) = t4_Nums(subjnum).MeanChunked_IPI(i);
                    MT4_R{subjnum}(p12, p34) = t4_Nums(subjnum).MeanRand_IPI(i);
                    MT4_A{subjnum}(p12, p34) = t4_Nums(subjnum).MeanAll_IPI(i);
                end
            end
            PoC4(p12, :) = PoC4(p12, :)/sum(PoC4(p12, :)); % sets the sum of every row to 1
        end
        
        figure('color' , 'white')
        subplot(349)
        imagesc(PoC4);
        ylabel('Press 1 2')
        xlabel('Press 3 4')
        axis square
        colorbar
        set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
            'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
            'XTickLabelRotation' , 45)
        title('Probability of Occurence for 3rd-order transitions in All sequences')
        for j = 1:length(ChunkNum4)
            y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
            x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
            rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
        end
        for subjnum = 1
            subplot(3,4,10)
            if ~lastIPI
                imagesc(MT4_C{subjnum},[800 2700]);  % for sum of IPIS
            else
                imagesc(MT4_C{subjnum},[200 600]);  % for the last IPIS
            end
            ylabel('Press1')
            xlabel('Press2')
            axis square
            set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' ,cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
                'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
                'XTickLabelRotation' , 45)
            title(['Average Movement time (Chunked) H = ' , num2str(subjnum)])
            
            for j = 1:length(ChunkNum4)
                y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
                x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            
            colorbar
        end
        
        
        for subjnum = 1
            subplot(3,4,11)
            if ~lastIPI
                imagesc(MT4_R{subjnum},[800 2700]);% for sum of IPIS
            else
                imagesc(MT4_R{subjnum},[200 600]);  % for the last IPIS
            end
            
            ylabel('Press1')
            xlabel('Press2')
            axis square
            set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' ,cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
                'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
                'XTickLabelRotation' , 45)
            title(['Average Movement time (Random) H = ' , num2str(subjnum)])
            
            for j = 1:length(ChunkNum4)
                y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
                x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            
            colorbar
            
            imcount = imcount+1;
        end
        
        
        for subjnum = 1
            subplot(3,4,12)
            if ~lastIPI
                imagesc(MT4_A{subjnum} , [800 2700]);% for sum of IPIS
            else
                imagesc(MT4_A{subjnum},[200 600]);  % for the last IPIS
            end
            ylabel('Press1')
            xlabel('Press2')
            axis square
            set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' ,cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
                'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
                'XTickLabelRotation' , 45)
            title(['Average Movement time (All) H = ' , num2str(subjnum)])
            
            for j = 1:length(ChunkNum4)
                y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
                x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
                rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
            end
            
            colorbar
            
        end
        
        
        
        out = [];
    case 'glm_IPIs'
        
        %% Significance F-tests
        %-------- 1st order transitions
        
        temp = anovaMixed(C.IPI_norm ,C.SN ,'within',C.t2Rank_n_binned,{'DoubleTrans'},'intercept',1,'subset' , ismember(C.Day , [2:4]));% & C.t2Rank_n(:,1)<=25);
        out.chunked.probEffect_N(1,1) = temp.eff(2).p;
        
        %-------- 2nd order transitions
        
        temp= anovaMixed(C.IPI_norm ,C.SN ,'within',C.T3Rank_n_binned(:,1),{'tripTrans'},'intercept',1,'subset' , ismember(C.Day , [2:4]));% & C.T3Rank_n(:,1)<=10);
        out.chunked.probEffect_N(1,2) = temp.eff(2).p;
        
        %-------- 3rd order transitions
        
        
        temp = anovaMixed(C.IPI_norm ,C.SN ,'within',C.T4Rank_n_binned,{'QuadTrans'},'intercept',1,'subset' , ismember(C.Day , [2:4]));% & C.T4Rank_n(:,1)<=20);
        out.chunked.probEffect_N(1,3)=temp.eff(2).p;
        
        %-------- chunk placement
        
        temp = anovaMixed(C.IPI_norm ,C.SN ,'within',C.IPIarrangement,{'DoubleTrans'},'intercept',1,'subset' , ismember(C.Day , [2:4]));% & C.t2Rank_n(:,1)<=25);
        out.chunked.chnkEffect_N = temp.eff(2).p;
        %-------- 1st order transitions
        
        temp = anovaMixed(R.IPI_norm ,R.SN ,'within',R.t2Rank_n_binned,{'DoubleTrans'},'intercept',1,'subset' , ismember(R.Day , [2:4]));
        out.Random.probEffect_N(1,1)=temp.eff(2).p;
        
        %-------- 2nd order transitions
        
        
        temp = anovaMixed(R.IPI_norm ,R.SN ,'within',R.T3Rank_n_binned(:,1),{'tripTrans'},'intercept',1,'subset' , ismember(R.Day , [2:4]));
        out.Random.probEffect_N(1,2) = temp.eff(2).p;
        
        %-------- 3rd order transitions
        
        temp = anovaMixed(R.IPI_norm ,R.SN ,'within',R.T4Rank_n_binned(:,1),{'QuadTrans'},'intercept',1,'subset' , ismember(R.Day , [2:4]) );
        out.Random.probEffect_N(1,3)=temp.eff(2).p;
        
        %-------- chunk placement
        
        
        temp = anovaMixed(R.IPI_norm ,R.SN ,'within',R.IPIarrangement,{'DoubleTrans'},'intercept',1,'subset' , ismember(R.Day , [2:4]));
        out.Random.chnkEffect_N=temp.eff(2).p;
        
        %% Visualize by rank and not by transitions cz most probable ones are different for every subject
        for h =1
            
            h1 = figure;
            hold on
            [XC2_n{h},PC2_n{h},EC2_n{h}] = lineplot(C.t2Rank_n_binned , C.IPI_norm ,'plotfcn' , 'nanmedian', 'subset' , ismember(C.Day , [2:5]));
            [XC3_n{h},PC3_n{h},EC3_n{h}] = lineplot(C.T3Rank_n_binned(:,1) , C.IPI_norm ,'plotfcn' , 'nanmedian', 'subset' , ismember(C.Day , [2:5]));
            [XC4_n{h},PC4_n{h},EC4_n{h}] = lineplot(C.T4Rank_n_binned(:,1) , C.IPI_norm ,'plotfcn' , 'nanmedian', 'subset' , ismember(C.Day , [2:5]) );
            [XCA{h},PCA{h},ECA{h}] = lineplot(C.IPIarrangement , C.IPI_norm ,'plotfcn' , 'nanmedian', 'subset' , ismember(C.Day , [2:5]));
            
            
            [XR2_n{h},PR2_n{h},ER2_n{h}] = lineplot(R.t2Rank_n_binned , R.IPI_norm , 'plotfcn' , 'nanmedian','subset' , ismember(R.Day , [2:5]));
            [XR3_n{h},PR3_n{h},ER3_n{h}] = lineplot(R.T3Rank_n_binned(:,1) , R.IPI_norm , 'plotfcn' , 'nanmedian','subset' , ismember(R.Day , [2:5]));
            [XR4_n{h},PR4_n{h},ER4_n{h}] = lineplot(R.T4Rank_n_binned(:,1) , R.IPI_norm ,'plotfcn' , 'nanmedian','subset' , ismember(R.Day , [2:5]) );
            [XRA{h},PRA{h},ERA{h}] = lineplot(R.IPIarrangement , R.IPI_norm ,'plotfcn' , 'nanmedian', 'subset' , ismember(R.Day , [2:5]));
            
        end
        
        %%
        cCount = 1;
        figure('color' , 'white')
        for h = 1
            subplot(121)
            errorbar(XCA{h},PCA{h},ECA{h}, 'LineWidth' , 3);
            hold on
            grid on
            title('chunk placement in chunked seqs','FontSize' , 14)
            cCount = cCount+1;
        end
        set(gca,'XTick' , [1 2] , 'FontSize' , 16 , 'XTickLabel' , {'Between' , 'Within'})
        xlabel('Chunk placement')
        ylabel('Median Normalized IPI')
        cCount = 1;
        for h = 1
            subplot(122)
            errorbar(XRA{h},PRA{h},ERA{h}, 'LineWidth' , 3);
            hold on
            grid on
            title('chunk placement in random seqs','FontSize' , 14)
            cCount = cCount+1;
        end
        set(gca,'XTick' , [0 2] , 'FontSize' , 16 , 'XTickLabel' , {'Random' , 'Within'})
        xlabel('Probability Class')
        ylabel('Median Normalized IPI')
        
        figure('color' , 'white')
        cCount = 1;
        for h = 1
            subplot(3,2,1)
            errorbar(XC2_n{h},PC2_n{h},EC2_n{h}, 'LineWidth' , 3 );
            hold on
            grid on
            title('Double trans chunked seqs sorted by norm-prob','FontSize' , 14)
            cCount = cCount+1;
        end
        set(gca,'XTick' , [1:5] , 'FontSize' , 16)
        xlabel('Probability Class')
        ylabel('Median Normalized IPI')
        cCount = 1;
        for h = 1
            subplot(3,2,2)
            errorbar(XR2_n{h},PR2_n{h},ER2_n{h}, 'LineWidth' , 3 );
            hold on
            grid on
            title('Double trans rand seqs sorted by norm-prob','FontSize' , 14)
            cCount = cCount+1;
        end
        set(gca,'XTick' , [1:5] , 'FontSize' , 16)
        xlabel('Probability Class')
        ylabel('Median Normalized IPI')
        
        
        cCount = 1;
        %         figure('color' , 'white')
        for h = 1
            subplot(323)
            errorbar(XC3_n{h},PC3_n{h},EC3_n{h}, 'LineWidth' , 3);
            hold on
            grid on
            title('Trip trans chunked seqs sorted by norm-prob','FontSize' , 14)
            cCount = cCount+1;
        end
        set(gca,'XTick' , [1:5] , 'FontSize' , 16)
        xlabel('Probability Class')
        ylabel('Median Normalized IPI')
        cCount = 1;
        for h = 1
            subplot(324)
            errorbar(XR3_n{h},PR3_n{h},ER3_n{h}, 'LineWidth' , 3 );
            hold on
            grid on
            title('Double trans rand seqs sorted by norm-prob','FontSize' , 14)
            cCount = cCount+1;
        end
        set(gca,'XTick' , [1:5] , 'FontSize' , 16)
        xlabel('Probability Class')
        ylabel('Median Normalized IPI')
        
        
        
        cCount = 1;
        %         figure('color' , 'white')
        for h = 1
            subplot(325)
            errorbar(XC4_n{h},PC4_n{h},EC4_n{h}, 'LineWidth' , 3);
            hold on
            grid on
            title('Quad trans chunked seqs sorted by norm-prob','FontSize' , 14)
            cCount = cCount +1;
        end
        set(gca,'XTick' , [1:5] , 'FontSize' , 16)
        xlabel('Probability Class')
        ylabel('Median Normalized IPI')
        cCount = 1;
        for h = 1
            subplot(326)
            errorbar(XR4_n{h},PR4_n{h},ER4_n{h}, 'LineWidth' , 3);
            hold on
            grid on
            title('Quad trans rand seqs sorted by norm-prob','FontSize' , 14)
            cCount = cCount +1;
        end
        set(gca,'XTick' , [1:5] , 'FontSize' , 16)
        xlabel('Probability Class')
        ylabel('Median Normalized IPI')
        
        
        
        %% The GLM
        switch N3
            case 'c'
                M = C;
                titleSuffix = 'Chunked Sequences';
            case 'r'
                M = R;
                titleSuffix = 'Random Sequences';
            case 'a'
                M = All;
                titleSuffix = 'All Sequences';
        end
        for dd = 2:4
            for sn = 1:length(subj_name) - 1
                clear T TSS RSS0 RSS1 FSS1 FSS0 L X X1 X2 X3 X4 X5 X6
                T = getrow(M , ismember(M.SN , sn) & ismember(M.Day , dd));
                
                
                L = length(T.IPI);
                X1 = ones(L , 1); % intercept
                X2 = T.IPIarrangement;
                %             X2 = T.estIPIarrangement;
                
                X2(X2==1) = -1; % between
                X2(X2==0) = -1; % Random = between
                X2(X2==2) = 1;  % within
                
                
                
                switch norma
                    case 1
                        X3 = T.t2Rank_n_binned;
                        switch LastIPI
                            case 1
                                %                                 X4 = T.t3Prob_n(:,1);
                                %                                 X5 = T.t4Prob_n(:,1);
                                X4 = T.T3Rank_n_binned;
                                X5 = T.T4Rank_n_binned;
                            otherwise
                                X4 = mean(T.t3Prob_n , 2);
                                X5 = mean(T.t4Prob_n , 2);
                        end
                    case 0
                        X3 = T.t2Prob;
                        switch LastIPI
                            case 1
                                X4 = T.t3Prob(:,1);
                                X5 = T.t4Prob(:,1);
                            otherwise
                                X4 = mean(T.t3Prob , 2);
                                X5 = mean(T.t4Prob , 2);
                        end
                end
                
                X6 = max(T.BN) - [T.BN] +1;
                switch N4
                    case {'n' 'N'}
                        X = [X1 X2 X3 X4 X5 X6];
                        xx    = {[1 6] , [1 2 6] ,  [1 3 6]  , [1,3:4,6]           ,[1 3:6] ,      [1:3,6]  ,        [1:4,6]          [1:6]};
                        Y = T.IPI;
                        label = {'I+L' '  I+C+L', 'I+1st+L' ,'I+1st+2nd+L'    'I+1st+2nd+3rd+L' , 'I+C+1st+L'  ,  'I+C+1st+2nd+L'  ,  'Full'};
                    otherwise
                        X = [X1 X2 X3 X4 X5 X6];
                        xx    = {[1 6] , [1 2 6] ,  [1 3 6]  , [1,3:4,6]           ,[1 3:6] ,      [1:3,6]  ,        [1:4,6]          [1:6]};
                        Y = T.IPI_norm;
                        label = {'I+L' '  I+C+L', 'I+1st+L' ,'I+1st+2nd+L'    'I+1st+2nd+3rd+L' , 'I+C+1st+L'  ,  'I+C+1st+2nd+L'  ,  'Full'};
                end
                Xnew = [X1];
                B = pinv(Xnew'*Xnew)*Xnew'*Y;
                Ypred = Xnew*B;
                Res  = Y-Ypred;
                TSS = sum(Y.^2); % Total Variance
                FSS1 = sum(Ypred.^2); % Fitted Variance of the Null Model (just the intercept)
                RSS1 = sum((Y-Ypred).^2); % Residual Variance of the Null Model
                for k = 1:length(xx)
                    Xnew = X(:,xx{k});
                    Bnew = pinv(Xnew'*Xnew)*Xnew'*Y;
                    Ypred_new = Xnew*Bnew;
                    FSS0(k) = sum(Ypred_new.^2);  % Fitted Variance of the partial model
                    RSS0(k) = sum((Y-Ypred_new).^2); % Residual Variance of the partial Model
                end
                % ___________________________      R_squared = 1 - (Residual variance of the partial model/Residual variance of the null model)
                R2{dd-1}(sn , :) = 1 - (RSS0./RSS1);
                xlab{dd-1}(sn , :) = 1:length(xx);
                % ___________________________
            end
        end
        figure('color' , 'white')
        for dd = 1:3
            subplot(1,3,dd)
            cCount = 1;
            f = figure;
            [xcoord,ePLOT,ERROR] = lineplot(reshape(xlab{dd} , numel(xlab{dd}) , 1) , reshape(R2{dd} , numel(R2{dd}) , 1) , 'plotfcn' , 'nanmean', 'linewidth' , 3);
            close(f)
            errorbar(xcoord,ePLOT,ERROR  , 'LineWidth' , 3)
            hold on
            plotshade(xcoord',ePLOT , ERROR,'transp' , .2 ,'linewidth' , 3 , 'linestyle' , ':');
            ylabel('R^2')
            title([titleSuffix, ' , Days ' , num2str(dd)])
            set(gca , 'XLim' , [0 length(xx)+1] , 'XTick' , [1: length(xx)] , 'XTickLabels' , label , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'YLim' , [0 .4],'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' Day ' , num2str(dd+1)])
            grid on
        end
        
        
        out = [];
    case 'Crossval_GLM_ols'
        
        
        %  N1 = input('Use Conditional Transition Probabilities? (y/n)' , 's');
        %         switch N1
        %             case 'y'
        %                 norm = 1;
        %             otherwise
        %                 norm = 0;
        %         end
        %         N2 = input('Use the last IPI for 2nd/3rd order transitions? (y/n)' , 's');
        %         switch N2
        %             case 'y'
        %                 LastIPI = 1;
        %             otherwise
        %                 LastIPI = 0;
        %         end
        cond = 1;
        LastIPI = 1;
        N3 = input('look into Chunked/Random/All sequences? (c/r/a)'  , 's');
        N4 = input('Use normalized IPIs? (Y/N)'  , 's');
        %         N5 = input('Analyse Medians, or raw IPIs? (M/R)' , 's');
        %         N3 = 'c';
        N5 = 'r';
        N6 = input('Bar or Shade? (b/s)' , 's');
        load([baseDir , '/CMB.mat'])
        
        
        load([baseDir , '/se1_TranProb.mat'])
        for count = 1:length(C.BN)
            C.Presses(count , 1:2) = CMB.comb2(C.t2(count),:);
        end
        for count = 1:length(R.BN)
            R.Presses(count , 1:2) = CMB.comb2(R.t2(count),:);
        end
        All = addstruct(C,R);
        
        switch N5
            case {'M' , 'm'}
                C = C_Summarized;
                R = R_Summarized;
        end
        
        
        
        %% develope the GLM
        
        switch N3
            case 'c'
                M = C;
                titleSuffix = 'Chunked';
                
            case 'r'
                M = R;
                titleSuffix = 'Random';
                
            case 'a'
                M = All;
                titleSuffix = 'All Sequences';
        end
        % =================== Form Design Matrix
        L = length(M.IPI);
        M.IPIarrangement(M.IPIarrangement == 0) = 1;
        M.X1 = ones(L , 1); % intercept
        M.X2 = M.IPIarrangement;
        %             T.X2 = T.estIPIarrangement;
        
        M.X2(M.X2==1) = 1; % between
        M.X2(M.X2==0) = 1; % Random = between
        M.X2(M.X2==2) = -1;  % within
        
        switch cond
            case 1
                temp = M.t2Prob_n;
                M.X3 = ceil(1.5*(1 + mapminmax(temp')))';
                switch LastIPI
                    case 1
                        temp = M.t3Prob_n(:,1);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = M.t4Prob_n(:,1);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                    otherwise
                        temp = mean(M.t3Prob_n , 2);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = mean(M.t4Prob_n , 2);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                end
            case 0
                temp = M.t2Prob;
                M.X3 = ceil(1.5*(1 + mapminmax(temp')))';
                switch LastIPI
                    case 1
                        temp = M.t3Prob(:,1);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = M.t4Prob(:,1);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                    otherwise
                        temp = mean(M.t3Prob , 2);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = mean(M.t4Prob , 2);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                end
        end
        nB = fliplr(round(linspace(1 , 20  ,max(M.BN))));
        M.X6 = M.BN;
        for q = 1:length(nB)
            M.X6(M.X6 == q) = nB(q);
        end
        
        M.X7 = ismember(M.t2 , [21:25])  +1;
        M.X = [M.X1 M.X2 M.X3 M.X4 M.X5 M.X6 M.X7];
        switch N4
            case {'n' 'N'}
                xx =     {[1 6]   [1 6:7]    [1 2 6:7] ,   [1 3 6:7] ,   [1 4 6:7]        [1 3:5,6:7]          ,[1:3  , 6:7] ,      [1:4,6:7]  ,    [1:7]};
                label = {'L'    'L+R'  'C+L+R',     '1st+L+R' ,'1st+2nd+L+R'    '1st+2nd+3rd+L+R' , 'C+1st+L+R'  ,  'C+1st+2nd+L+R'  ,  'Full'};
                 plotIND = [2:5 , 8];
            otherwise
                M.IPI = M.IPI_norm;
                titleSuffix = [titleSuffix , '-norm'];
                xx =     {[7]    [2 7] ,   [3 7] ,   [4 7]        [3:5,7]          ,[2:3  , 7] ,      [2:4,7]  ,    [2:5 , 7]};
                label = {'R'  'C+R',     '1st+R' ,'1st+2nd+R'    '1st+2nd+3rd+R' , 'C+1st+R'  ,  'C+1st+2nd+R'  ,  'Full'};
                 plotIND = [2:5 , 8];
        end
        cleanLabel = {'within/between Chunk', '1st order probability' ,'1st + 2nd order probability' ,'1st + 2nd + 3rd order probability' ,'Full Model'};
         % =================== % =================== % =================== % =================== Make Modelzzzzz BITCH!
        
        CVfol = 3;
        if calc
            count = 1;
            h = 1;
            for dd = 2:4
                for sn = 1:length(subj_name) - 1
                    T = getrow(M , ismember(M.SN , sn)  & ismember(M.Day , dd));
                    L = length(T.IPI);
                    CVI = crossvalind('Kfold', L, CVfol);
                    for cvl = 1:CVfol
                        Test = getrow(T , CVI==cvl);
                        Train = getrow(T , CVI~=cvl);
                        params  = xx{1};   % Null Model
                        X = Train.X(:,params);
                        X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                        Y = Train.IPI - mean(Train.IPI);  % mean subtract the output
                        Mdl = fitglm(X , Y ,'Intercept',false);
                        
                        %Test
                        X = Test.X(:,params);
                        X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                        Y = Test.IPI - mean(Test.IPI);  % mean subtract the output
                        [Ypred0,Posterior] = predict(Mdl,X);
                        for ml = 1:length(xx)
                            params  = xx{ml};
                            X = Train.X(:,params);
                            X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                            Y = Train.IPI - mean(Train.IPI);  % mean subtract the output
                            Mdl = fitglm(X , Y,'Intercept',false);
                            
                            %test
                            X = Test.X(:,params);
                            X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                            Y = Test.IPI - mean(Test.IPI);  % mean subtract the output
                            [Ypred,Posterior] = predict(Mdl,X);
                            
                            [cv_Dev , lh_comp] = se2_crossval_DEVandLH(Y , Ypred , Ypred0);
                            out.R2(count,:)   = se2_R2ModelComp(Y , Ypred0 , Ypred);
                            out.R2_adjusted(count , 1) = se2_R2Adjusted(Y, Ypred , length(params) + 1);
                            
                            out.B{count,:} = Mdl.Coefficients.Estimate;
                            % Deviation  = -2*ln[(likelihood of fitted model)/(likelihood of saturated model)]
                            out.Dev(count , 1) = Mdl.Deviance;
                            out.cv_Dev(count , 1) = cv_Dev;
                            out.lh_comp(count , 1) = lh_comp;
                            out.hor(count , 1) = h;
                            out.subj(count , 1) = sn;
                            out.day(count , 1) = dd;
                            out.xx(count , 1) = ml;
                            out.cv (count , 1) = cvl;
                            temp = corrcoef(Y , Ypred);
                            out.corYY (count , 1) = temp(2);
                            count = count+1;
                            disp(['Model ' , num2str(ml) , ' - Day ' , num2str(dd) , ' - Subject ' , num2str(sn) , ' - Horizon ' , num2str(h)])
                        end
                    end
                end
            end
            Mdl = out;
            save([baseDir , '/se1_CrossvalIPI_',titleSuffix,'_OLS.mat'] , 'Mdl' , '-v7.3')
        else
            load([baseDir , '/se1_CrossvalIPI_',titleSuffix,'_OLS.mat']);
            out =  Mdl;
        end
        
        
        clear Mdl
        % =================== % =================== % =================== % =================== Compare Modelzzzzz BITCH!
        dayz = {[2] [3] [4]};
        K = tapply(out , {'hor' , 'subj' , 'day' , 'xx'} , {'R2' , 'nanmean'} ,{'R2_adjusted' , 'nanmean'} ,  {'corYY' , 'nanmean'} ,  {'lh_comp' , 'nanmean'}); % average over cv loops
        
        clear xp_dev pp_dev ep_dev xp_r2 pp_r2 ep_r2 xp_r2a pp_r2a ep_r2a dev_image R2_image xp_dn  pp_dn ep_dn dvn_image
        h1 = figure;
        h = 1
        for dd = 1:length(dayz)
            [xp_cor{dd}(h, :) , pp_cor{dd}(h,:) , ep_cor{dd}(h,:)]  = lineplot(K.xx, K.corYY , 'plotfcn','nanmean' ,  'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            hold on
            [xp_r2{dd}(h, :) , pp_r2{dd}(h,:) , ep_r2{dd}(h,:)]  = lineplot(K.xx, K.R2 , 'plotfcn','nanmean' , 'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            [xp_r2a{dd}(h, :) , pp_r2a{dd}(h,:) , ep_r2a{dd}(h,:)]  = lineplot(K.xx, K.R2_adjusted , 'plotfcn','nanmean' , 'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            [xp_ml{dd}(h, :) , pp_ml{dd}(h,:) , ep_ml{dd}(h,:)]  = lineplot(K.xx, K.lh_comp ,  'plotfcn','nanmean' ,'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            hold on
        end
        
        close(h1)
        % =================== % =================== % =================== % =================== Visualize model R2 comparisons!
        
         switch N6
            case{'s'}
                
                
                figure('color' , 'white')
                for dd = 1:length(dayz)
                    subplot(1,length(dayz),dd)
                    for h = 1
                        %   errorbar(xp{g}(h, :) , pp{g}(h,:) , ep{g}(h,:) , 'color' , colors(cCount,:) , 'LineWidth' , 3)
                        hold on
                        eval(['h' , num2str(h) , ' = plotshade(xp_r2{dd}(h, 1:end) , pp_r2{dd}(h,1:end) , ep_r2{dd}(h,1:end),''transp'' , .2 , ''patchcolor'' , ''b'', ''linecolor'' , ''b'' , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                        hold on
                    end
                    ylabel('R^2 nomarlized to the Null model')
                    set(gca , 'XLim' , [0 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , cleanLabel , 'FontSize' , 20 ,...
                        'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
                    title([titleSuffix , ' , Days ' , num2str(dayz{dd})])
                    grid on
                end
                out = [];
                
                
                
                
                figure('color' , 'white')
                for dd = 1:length(dayz)
                    subplot(1,length(dayz),dd)
                    for h = 1
                        %   errorbar(xp{g}(h, :) , pp{g}(h,:) , ep{g}(h,:) , 'color' , colors(cCount,:) , 'LineWidth' , 3)
                        hold on
                        eval(['h' , num2str(h) , ' = plotshade(xp_cor{dd}(h, 1:end) , pp_cor{dd}(h,1:end) , ep_cor{dd}(h,1:end),''transp'' , .2 , ''patchcolor'' , ''b'', ''linecolor'' , ''b'' , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                        hold on
                    end
                    ylabel('Prediction-output correlation to the Null model')
                    set(gca , 'XLim' , [1 length(xx)],'XTick' , [1: length(xx)] , 'XTickLabels' , label(1:end) , 'FontSize' , 20 ,...
                        'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
                    title([titleSuffix , ' , Days ' , num2str(dayz{dd})])
                    grid on
                end
                out = [];
            otherwise
                xp_ = reshape(cell2mat(xp_r2) , size(xp_r2{1} , 1) , size(xp_r2{1} , 2) , length(xp_r2));
                pp_ = reshape(cell2mat(pp_r2) , size(pp_r2{1} , 1) , size(pp_r2{1} , 2) , length(pp_r2));
                ep_ = reshape(cell2mat(ep_r2) , size(ep_r2{1} , 1) , size(ep_r2{1} , 2) , length(ep_r2));
                
                
                lo = min(min(min(pp_)));
                hi = max(max(max(pp_)));
                ylim = [1.1*lo   1.1*hi];
                figure('color' , 'white')
                bar(squeeze(pp_(:,plotIND , :))');
                grid on
                set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] , 'XTickLabel' , {'Day 1' , 'Day 2' , 'Day 3'},'YLim' , ylim)
                title([titleSuffix ,' - Model Crossvalidated R^2'])
                legend(cleanLabel)
                
                xp_ = reshape(cell2mat(xp_cor) , size(xp_cor{1} , 1) , size(xp_cor{1} , 2) , length(xp_cor));
                pp_ = reshape(cell2mat(pp_cor) , size(pp_cor{1} , 1) , size(pp_cor{1} , 2) , length(pp_cor));
                ep_ = reshape(cell2mat(ep_cor) , size(ep_cor{1} , 1) , size(ep_cor{1} , 2) , length(ep_cor));
                
                lo = min(min(min(pp_)));
                hi = max(max(max(pp_)));
                ylim = [1.1*lo   1.1*hi];
                figure('color' , 'white')
                
                bar(squeeze(pp_(:,plotIND , :))');
                grid on
                set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] , 'XTickLabel' , {'Day 1' , 'Day 2' , 'Day 3'},...
                    'YLim' , ylim)
                title([titleSuffix , ' - Model Prediction - Output Correlation'])
                legend(cleanLabel)
        end
        
    case 'Crossval_GLM_ridge' 
        % for ridge the trick is to mean subtract each column and dont
        % include intecept - the mean of the intercept would be too high
        % and it will push the week regreesors towards zero oo much
        % but do mean subtract the Y
        
        
        
        %             N1 = input('Use Conditional Transition Probabilities? (y/n)' , 's');
        %         switch N1
        %             case 'y'
        %                 norm = 1;
        %             otherwise
        %                 norm = 0;
        %         end
        %         N2 = input('Use the last IPI for 2nd/3rd order transitions? (y/n)' , 's');
        %         switch N2
        %             case 'y'
        %                 LastIPI = 1;
        %             otherwise
        %                 LastIPI = 0;
        %         end
        cond = 1;
        LastIPI = 1;
        N3 = input('look into Chunked/Random/All sequences? (c/r/a)'  , 's');
        N4 = input('Use normalized IPIs? (Y/N)'  , 's');
        %         N5 = input('Analyse Medians, or raw IPIs? (M/R)' , 's');
        N6 = input('Bar or Shade? (b/s)' , 's');
        %
        N5 = 'r';
        
        
            load([baseDir , '/CMB.mat']);
     
        
        load([baseDir , '/se1_TranProb.mat'])
        for count = 1:length(C.BN)
            C.Presses(count , 1:2) = CMB.comb2(C.t2(count),:);
        end
        for count = 1:length(R.BN)
            R.Presses(count , 1:2) = CMB.comb2(R.t2(count),:);
        end
        All = addstruct(C,R);
        
        switch N5
            case {'M' , 'm'}
                C = C_Summarized;
                R = R_Summarized;
        end
        
        
        
        %% develope the GLM
        
        switch N3
            case 'c'
                M = C;
                titleSuffix = 'Chunked';                
            case 'r'
                M = R;
                titleSuffix = 'Random';
                
            case 'a'
                M = All;
                titleSuffix = 'All Sequences';
        end
        % =================== Form Design Matrix
        L = length(M.IPI);
        M.IPIarrangement(M.IPIarrangement == 0) = 1;
        M.X1 = ones(L , 1); % intercept
        M.X2 = M.IPIarrangement;
        %             T.X2 = T.estIPIarrangement;
        
        M.X2(M.X2==1) = 1; % between
        M.X2(M.X2==0) = 1; % Random = between
        M.X2(M.X2==2) = -1;  % within
        
        switch cond
            case 1
                temp = M.t2Prob_n;
                M.X3 = ceil(1.5*(1 + mapminmax(temp')))';
                switch LastIPI
                    case 1
                        temp = M.t3Prob_n(:,1);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = M.t4Prob_n(:,1);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                    otherwise
                        temp = mean(M.t3Prob_n , 2);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = mean(M.t4Prob_n , 2);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                end
            case 0
                temp = M.t2Prob;
                M.X3 = ceil(1.5*(1 + mapminmax(temp')))';
                switch LastIPI
                    case 1
                        temp = M.t3Prob(:,1);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = M.t4Prob(:,1);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                    otherwise
                        temp = mean(M.t3Prob , 2);
                        M.X4 = ceil(1.5*(1 + mapminmax(temp')))';
                        temp = mean(M.t4Prob , 2);
                        M.X5 = ceil(1.5*(1 + mapminmax(temp')))';
                end
        end
        nB = fliplr(round(linspace(1 , 20  ,max(M.BN))));
        M.X6 = M.BN;
        for q = 1:length(nB)
            M.X6(M.X6 == q) = nB(q);
        end
        
        M.X7 = ismember(M.t2 , [21:25])  +1;
        M.X = [M.X1 M.X2 M.X3 M.X4 M.X5 M.X6 M.X7];
        
        switch N4
            case {'n' 'N'}
                xx =     {[6]   [6:7]    [2 6:7] ,   [3 6:7] ,   [4 6:7]        [3:5,6:7]          ,[2:3  , 6:7] ,      [2:4,6:7]  ,    [2:7]};
                label = {'L'    'L+R'  'C+L+R',     '1st+L+R' ,'1st+2nd+L+R'    '1st+2nd+3rd+L+R' , 'C+1st+L+R'  ,  'C+1st+2nd+L+R'  ,  'Full'};
                ylim = [-0.001 .01];
                plotIND = [3:6 , 9];
            otherwise
                M.IPI = M.IPI_norm;
                titleSuffix = [titleSuffix , '-norm'];
                xx =     {[7]    [2 7] ,   [3 7] ,   [4 7]        [3:5,7]          ,[2:3  , 7] ,      [2:4,7]  ,    [2:5 , 7]};
                label = {'R'  'C+R',     '1st+R' ,'1st+2nd+R'    '1st+2nd+3rd+R' , 'C+1st+R'  ,  'C+1st+2nd+R'  ,  'Full'};
                 ylim  =[-0.0055 .06];
                plotIND = [2:5 , 8];
        end
        cleanLabel = {'within/between Chunk', '1st order probability' ,'1st + 2nd order probability' ,'1st + 2nd + 3rd order probability' ,'Full Model'};
%         cat_cleanLabel = categorical(repmat(cleanLabel , length(hor)-1 , 1));
        % =================== % =================== % =================== % =================== Make Modelzzzzz BITCH!
        
        CVfol = 3;
        if calc
            count = 1;
            h = 1;
            for dd = 2:4
                for sn = 1:length(subj_name) - 1
                    T = getrow(M , ismember(M.SN , sn) & ismember(M.Day , dd));
                    L = length(T.IPI);
                    CVI = crossvalind('Kfold', L, CVfol);
                    for cvl = 1:CVfol
                        Test = getrow(T , CVI==cvl);
                        Train = getrow(T , CVI~=cvl);
                        params  = xx{1};   % Null Model
                        X = Train.X(:,params);
                        X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                        Y = Train.IPI - mean(Train.IPI);  % mean subtract the output
                        Mdl = fitrlinear(X , Y,'Lambda',.1,'Regularization','ridge' , 'FitBias' , false);
                        
                        X = Test.X(:,params);
                        X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                        Y = Test.IPI - mean(Test.IPI);  % mean subtract the output
                        Ypred0 = predict(Mdl,X);
                        for ml = 1:length(xx)
                            params  = xx{ml};
                            X = Train.X(:,params);
                            X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                            Y = Train.IPI - mean(Train.IPI);  % mean subtract the output
                            Mdl = fitrlinear(X , Y,'Lambda',.1,'Regularization','ridge' , 'FitBias' , false);
                            
                            X = Test.X(:,params);
                            X = X - repmat(mean(X) , length(X) , 1); % mean subtract the design matrix
                            Y = Test.IPI - mean(Test.IPI);  % mean subtract the output
                            
                            
                            Ypred = predict(Mdl,X);
                            [cv_Dev , lh_comp] = se2_crossval_DEVandLH(Y , Ypred , Ypred0);
                            out.R2(count,:)   = se2_R2ModelComp(Y , Ypred0 , Ypred);
                            out.R2_adjusted(count , 1) = se2_R2Adjusted(Y , Ypred , length(params) + 1);
                            out.B{count,:} =  Mdl.Beta;
                            % Deviation  = -2*ln[(likelihood of fitted model)/(likelihood of saturated model)]
                            %                                 out.Dev(count , 1) = Mdl.Deviance;
                            out.cv_Dev(count , 1) = cv_Dev;
                            out.lh_comp(count , 1) = lh_comp;
                            out.hor(count , 1) = h;
                            out.subj(count , 1) = sn;
                            out.day(count , 1) = dd;
                            temp = corrcoef(Y , Ypred);
                            out.corYY (count , 1) = temp(2);
                            out.xx(count , 1) = ml;
                            out.cv (count , 1) = cvl;
                            count = count+1;
                            disp(['Model ' , num2str(ml) , ' - Day ' , num2str(dd) , ' - Subject ' , num2str(sn) , ' - Horizon ' , num2str(h)])
                        end
                    end
                end
            end
            Mdl = out;
            save([baseDir , '/se1_CrossvalIPI_',titleSuffix,'_RR.mat'] , 'Mdl' , '-v7.3')
        else
            load([baseDir , '/se1_CrossvalIPI_',titleSuffix,'_RR.mat']);
            out =  Mdl;
        end
        
        clear Mdl
        % =================== % =================== % =================== % =================== Compare Modelzzzzz BITCH!
        dayz = {[2] [3] [4]};
        K = tapply(out , {'hor' , 'subj' , 'day' , 'xx'} , {'R2' , 'nanmean'} ,{'R2_adjusted' , 'nanmean'} ,  {'corYY' , 'nanmean'} ,  {'lh_comp' , 'nanmean'}); % average over cv loops
        
        clear xp_dev pp_dev ep_dev xp_r2 pp_r2 ep_r2 xp_r2a pp_r2a ep_r2a dev_image R2_image xp_dn  pp_dn ep_dn dvn_image
        h1 = figure;
        h = 1
        for dd = 1:length(dayz)
            [xp_cor{dd}(h, :) , pp_cor{dd}(h,:) , ep_cor{dd}(h,:)]  = lineplot(K.xx, K.corYY , 'plotfcn','nanmean' ,  'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            hold on
            [xp_r2{dd}(h, :) , pp_r2{dd}(h,:) , ep_r2{dd}(h,:)]  = lineplot(K.xx, K.R2 , 'plotfcn','nanmean' , 'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            [xp_r2a{dd}(h, :) , pp_r2a{dd}(h,:) , ep_r2a{dd}(h,:)]  = lineplot(K.xx, K.R2_adjusted , 'plotfcn','nanmean' , 'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            [xp_ml{dd}(h, :) , pp_ml{dd}(h,:) , ep_ml{dd}(h,:)]  = lineplot(K.xx, K.lh_comp ,  'plotfcn','nanmean' ,'subset' , K.hor == h & ismember(K.day , dayz{dd}));
            hold on
        end
        
        close(h1)
        % =================== % =================== % =================== % =================== Visualize model R2 comparisons!
        
        switch N6
            case{'s'}
                
                
                figure('color' , 'white')
                for dd = 1:length(dayz)
                    subplot(1,length(dayz),dd)
                    for h = 1
                        %   errorbar(xp{g}(h, :) , pp{g}(h,:) , ep{g}(h,:) , 'color' , colors(cCount,:) , 'LineWidth' , 3)
                        hold on
                        eval(['h' , num2str(h) , ' = plotshade(xp_r2{dd}(h, 1:end) , pp_r2{dd}(h,1:end) , ep_r2{dd}(h,1:end),''transp'' , .2 , ''patchcolor'' , ''b'', ''linecolor'' , ''b'' , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                        hold on
                    end
                    ylabel('R^2 nomarlized to the Null model')
                    set(gca , 'XLim' , [0 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , cleanLabel , 'FontSize' , 20 ,...
                        'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
                    title([titleSuffix , ' , Days ' , num2str(dayz{dd})])
                    grid on
                end
                out = [];
                
                
                
                
                figure('color' , 'white')
                for dd = 1:length(dayz)
                    subplot(1,length(dayz),dd)
                    for h = 1
                        %   errorbar(xp{g}(h, :) , pp{g}(h,:) , ep{g}(h,:) , 'color' , colors(cCount,:) , 'LineWidth' , 3)
                        hold on
                        eval(['h' , num2str(h) , ' = plotshade(xp_cor{dd}(h, 1:end) , pp_cor{dd}(h,1:end) , ep_cor{dd}(h,1:end),''transp'' , .2 , ''patchcolor'' , ''b'', ''linecolor'' , ''b'' , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                        hold on
                    end
                    ylabel('Prediction-output correlation to the Null model')
                    set(gca , 'XLim' , [1 length(xx)],'XTick' , [1: length(xx)] , 'XTickLabels' , label(1:end) , 'FontSize' , 20 ,...
                        'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
                    title([titleSuffix , ' , Days ' , num2str(dayz{dd})])
                    grid on
                end
                out = [];
            otherwise
                xp_ = reshape(cell2mat(xp_r2) , size(xp_r2{1} , 1) , size(xp_r2{1} , 2) , length(xp_r2));
                pp_ = reshape(cell2mat(pp_r2) , size(pp_r2{1} , 1) , size(pp_r2{1} , 2) , length(pp_r2));
                ep_ = reshape(cell2mat(ep_r2) , size(ep_r2{1} , 1) , size(ep_r2{1} , 2) , length(ep_r2));
                
                
                
                figure('color' , 'white')
                    bar(squeeze(pp_(:,plotIND , :))');
                    grid on
                    set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] , 'XTickLabel' , {'Day 1' , 'Day 2' , 'Day 3'},'YLim' , ylim)
                    title([titleSuffix ,' - Model Crossvalidated R^2'])
                legend(cleanLabel)
                
                xp_ = reshape(cell2mat(xp_cor) , size(xp_cor{1} , 1) , size(xp_cor{1} , 2) , length(xp_cor));
                pp_ = reshape(cell2mat(pp_cor) , size(pp_cor{1} , 1) , size(pp_cor{1} , 2) , length(pp_cor));
                ep_ = reshape(cell2mat(ep_cor) , size(ep_cor{1} , 1) , size(ep_cor{1} , 2) , length(ep_cor));
                
                
                figure('color' , 'white')
                
                    bar(squeeze(pp_(:,plotIND , :))');
                    grid on
                    set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] , 'XTickLabel' , {'Day 1' , 'Day 2' , 'Day 3'},...
                        'YLim' , [-.1 .4])
                    title([titleSuffix , ' - Model Prediction - Output Correlation'])
                legend(cleanLabel)
        end
        
        
        
end








