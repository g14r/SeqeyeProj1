
function out  = se1_runGLM(subjnum , what , distance)
%%  distances:
% 'euclidean'	      Euclidean distance (default).
% 'squaredeuclidean'  Squared Euclidean distance. (This option is provided for efficiency only. It does not satisfy the triangle inequality.)
% 'seuclidean'	      Standardized Euclidean distance. Each coordinate difference between rows in out.X and Y is scaled by dividing by the corresponding element of the standard deviation computed from out.X, S=nanstd(out.X). To specify another value for S, use D = PDIST2(out.X,Y,'seuclidean',S).
% 'cityblock'	      City block metric.
% 'minkowski'         Minkowski distance. The default exponent is 2. To compute the distance with a different exponent, use D = pdist2(out.X,Y,'minkowski',P), where the exponent P is a scalar positive value.
% 'chebychev'	      Chebychev distance (maximum coordinate difference).
% 'mahalanobis'	      Mahalanobis distance, using the sample covariance of out.X as computed by nancov. To compute the distance with a different covariance, use D = pdist2(out.X,Y,'mahalanobis',C) where the matrix C is symmetric and positive definite.
% 'cosine'	          One minus the cosine of the included angle between points (treated as vectors).
% 'correlation'	      One minus the sample correlation between points (treated as sequences of values).
% 'spearman'	      One minus the sample Spearman's rank correlation between observations, treated as sequences of values.
% 'hamming'           Hamming distance, the percentage of coordinates that differ.
% 'jaccard'           One minus the Jaccard coefficient, the percentage of nonzero coordinates that differ.%     case 'chunk_est_instance'
%%  Cases
%     case 'eye_vel_instance'
%     case 'IPI_dist'
%     case 'IPI_ttest_rand'
%     case 'chunk_dist'
%     case 'Avg_pattern_sh3'
%     case 'eye_vel_seqplace'
%     case 'eye_pos_seqplace'
%     case 'eyepress_windowed_distances'
%     case 'eyepress_windowed_avg_distances'
%     case 'eyepress_avg_traces'
%     case 'diff_eyepress_avg_traces'
%     case 'diff_eyepress_distances'          Distances between velocities---trial wise
%     case 'diff_eyepress_avg_distances'      Distances between velocities---Average pattern
%     case 'dtw'
%%
prefix = 'se1_';
baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye1/analyze';
%baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/se1_data/analyze';
subj_name = {'SZ1' , 'JN2' , 'SP1'};
load([baseDir , '/CMB.mat'])
N1 = load([baseDir , '/se1_tnorm_r1.mat']);
N(1) = N1.N;

N2 = load([baseDir , '/se1_tnorm_r2.mat']);
N(2) = N2.N;
%load([baseDir , '/se1_all.mat'])
switch what
    case 'press'
        out.Y.press    = [];
        out.Y.prsveloc = [];
        out.X = [];
        for rep = 1:2
            for d = 1:4
                for seq = 1:7
                    out.Y.press = [out.Y.press ; N(rep).norm(d , subjnum).press{seq}];
                    out.Y.prsveloc = [out.Y.prsveloc ; N(rep).norm(d , subjnum).prsveloc{seq}];
                    num = size(N(rep).norm(d , subjnum).prsveloc{seq} , 1);
                    tempx = [ones(num , 1) zeros(num , 13)];
                    tempx(:,1+seq)    = ones(num , 1);
                    tempx(:,8+d)    = ones(num , 1);
                    tempx(:,12+rep) = ones(num , 1);
                    out.X = [out.X ; tempx];
                end
            end
        end
        out.beta_press = pinv(out.X)*out.Y.press;
        out.beta_prsveloc = pinv(out.X)*out.Y.prsveloc;
       
        
        figure
        imagesc(out.X)
        title('Design matrix')
        hold on
        ax = gca;
        ax.XTick = [1:14];
        ax.XTickLabel = {'Intercept' 'Random' 'Structure1' 'Structur2' , 'Structure3' ,'Structure4' , 'Structure5' ,'Structure6' ,'Day1','Day2','Day3','Day4','Repetition1','Repetition2'};
        ax.XTickLabelRotation = 45;
        ylabel('Trials')
        
        figure
        subplot(3,1,1)
        
        plot(out.beta_press(3:8,:)' , 'LineWidth' , 1);
        hold on
        plot(out.beta_press(2,:)' , 'LineWidth' , 4);
        title('Betas obtained from the press timeseries')
        grid on
        xlabel('Normalized time')
        legend({ 'Structure1' 'Structur2' , 'Structure3' ,'Structure4' , 'Structure5' ,'Structure6' 'Random'})
        
        subplot(3,1,2)
        plot(out.beta_press(9:12,:)' , 'LineWidth' , 2);
        title('Betas obtained from the days of training')
        grid on
        xlabel('Normalized time')
        legend({'Day1','Day2','Day3','Day4'})
        
        subplot(3,1,3)
        plot(out.beta_press(13:14,:)' , 'LineWidth' , 2);
        title('Betas obtained from the repetition')
        grid on
        xlabel('Normalized time')
        legend({'Repetition1','Repetition2'})
        
        
        figure
        subplot(3,1,1)
        
        plot(out.beta_prsveloc(3:8,:)' , 'LineWidth' , 1);
        hold on
        plot(out.beta_prsveloc(2,:)' , 'LineWidth' , 4);
        title('Betas obtained from the press velocity timeseries')
        grid on
        xlabel('Normalized time')
        legend({'Structure1' 'Structur2' , 'Structure3' ,'Structure4' , 'Structure5' ,'Structure6' 'Random' })
        
        subplot(3,1,2)
        plot(out.beta_prsveloc(9:12,:)' , 'LineWidth' , 2);
        title('Betas obtained from the days of training')
        grid on
        xlabel('Normalized time')
        legend({'Day1','Day2','Day3','Day4'})
        
        subplot(3,1,3)
        plot(out.beta_prsveloc(13:14,:)' , 'LineWidth' , 2);
        title('Betas obtained from the repetition')
        grid on
        xlabel('Normalized time')
        legend({'Repetition1','Repetition2'})
        
        figure
        subplot(1,2,1)
        imagesc(squareform(pdist(out.beta_press(2:8,:) , distance)));
        title('Dissimilarity matrix of average patterns of press positions')
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.YTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.XTickLabelRotation = 45;
        axis square
        colorbar
        
        subplot(1,2,2)
        imagesc(squareform(pdist(out.beta_prsveloc(2:8,:) , distance)));
        title('Dissimilarity matrix of average patterns of press velocities')
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.YTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.XTickLabelRotation = 45;
        axis square
        colorbar
        
    case 'eye'
        out.Y.eye      = [];
        out.Y.eyebntn  = [];
        out.Y.eyeveloc = [];
        out.X = [];
        for rep = 1:2
            for d = 1:4
                for seq = 1:7
                    out.Y.eye = [out.Y.eye ; N(rep).norm(d , subjnum).eye{seq}];
                    out.Y.eyeveloc = [out.Y.eyeveloc ; N(rep).norm(d , subjnum).eyeveloc{seq}];
                    num = size(N(rep).norm(d , subjnum).eyeveloc{seq} , 1);
                    tempx = [ones(num , 1) zeros(num , 13)];
                    tempx(:,1+seq)    = ones(num , 1);
                    tempx(:,8+d)    = ones(num , 1);
                    tempx(:,12+rep) = ones(num , 1);
                    out.X = [out.X ; tempx];
                end
            end
        end
        out.beta_eyeveloc = pinv(out.X)*out.Y.eyeveloc;
        out.beta_eye = pinv(out.X)*out.Y.eye;
        
        figure  
 
        imagesc(out.X)
        title('Design matrix')
        hold on
        ax = gca;
        ax.XTick = [1:14];
        ax.XTickLabel = {'Intercept' 'Random' 'Structure1' 'Structur2' , 'Structure3' ,'Structure4' , 'Structure5' ,'Structure6' ,'Day1','Day2','Day3','Day4','Repetition1','Repetition2'};
        ax.XTickLabelRotation = 45;
        ylabel('Trials')

        
        figure
        subplot(3,1,1)
        
        plot(out.beta_eye(3:8,:)' , 'LineWidth' , 1);
        hold on
        plot(out.beta_eye(2,:)' , 'LineWidth' , 4);
        title('Betas obtained from the eye position timeseries')
        grid on
        xlabel('Normalized time')
        legend({'Structure1' 'Structur2' , 'Structure3' ,'Structure4' , 'Structure5' ,'Structure6' 'Random' })
        
        subplot(3,1,2)
        plot(out.beta_eye(9:12,:)' , 'LineWidth' , 2);
        title('Betas obtained from the days of training')
        grid on
        xlabel('Normalized time')
        legend({'Day1','Day2','Day3','Day4'})
        
        subplot(3,1,3)
        plot(out.beta_eye(13:14,:)' , 'LineWidth' , 2);
        title('Betas obtained from the repetition')
        grid on
        xlabel('Normalized time')
        legend({'Repetition1','Repetition2'})
        
        
        figure
        subplot(3,1,1)
        
        plot(out.beta_eyeveloc(3:8,:)' , 'LineWidth' , 1);
        hold on 
        plot(out.beta_eyeveloc(2,:)' , 'LineWidth' , 4);
        title('Betas obtained from the eye velocity timeseries')
        grid on
        xlabel('Normalized time')
        legend({'Structure1' 'Structur2' , 'Structure3' ,'Structure4' , 'Structure5' ,'Structure6' 'Random' })
        
        subplot(3,1,2)
        plot(out.beta_eyeveloc(9:12,:)' , 'LineWidth' , 2);
        title('Betas obtained from the days of training')
        grid on
        xlabel('Normalized time')
        legend({'Day1','Day2','Day3','Day4'})
        
        subplot(3,1,3)
        plot(out.beta_eyeveloc(13:14,:)' , 'LineWidth' , 2);
        title('Betas obtained from the repetition')
        grid on
        xlabel('Normalized time')
        legend({'Repetition1','Repetition2'})
        
        figure
        subplot(1,2,1)
        imagesc(squareform(pdist(out.beta_eye(2:8,:) , distance)));
        title('Dissimilarity matrix of average patterns of eye positions')
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.YTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.XTickLabelRotation = 45;
        axis square
        colorbar
        
        subplot(1,2,2)
        imagesc(squareform(pdist(out.beta_eyeveloc(2:8,:) , distance)));
        title('Dissimilarity matrix of average patterns of eye velocities')
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.YTickLabel = {'Random' 'CLA1' 'CLA2' , 'CLA3' ,'CLA4' , 'CLA5' ,'CLA6'}
        ax.XTickLabelRotation = 45;
        axis square
        colorbar
end


