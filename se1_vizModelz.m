function sig = se1_vizModelz(what)

% what  is 's' for shadeplot and 'b' for barplot



% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';     %macbook
% baseDir = '/Volumes/MotorControl/data/SeqEye1/analyze';  % server
baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye1/analyze';          %iMac


iN1 = input('Visualize Chunked or Random? (c/r)'  , 's');
iN3 = input('Fit or Crossvalidated? (f/c)' , 's');
switch iN3
    case 'c'
        iN2 = input('Ridge or OLS? (r/o)' , 's');
end
iN4 = input('Within, between, all? (w/b/a)'  , 's');
colors = [0 0 1;...
    0 1 0;...
    1 0 0;...
    0 1 1;...
    1 0 1;...
    1 0.69 0.39;...
    0.6 0.2 0;...
    0 0.75 0.75;...
    0.22 0.44 0.34;...
    0.32 0.19 0.19];
%% set up legends and labels
clear xp_dev pp_dev ep_dev xp_r2 pp_r2 ep_r2 xp_r2a pp_r2a ep_r2a dev_image R2_image xp_dn  pp_dn ep_dn dvn_image xp_cor  pp_cor  ep_cor sig
% x labels for shadeplot
label = {'R'  'C+R',     '1st+R' ,'1st+2nd+R'    '1st+2nd+3rd+R' , 'C+1st+R'  ,  'C+1st+2nd+R'  ,  'Full'};
dayz = {[2] [3] [4]};
plotIND = [2 5 , 8]; % the model indices to include in bar plot
%% load proper model

switch iN1
    case 'c'
        titleSuffix = 'Chunked';
    case 'r'
        titleSuffix = 'Random';
end
switch iN4
    case 'w'
        titleSuffix = [titleSuffix , '_within'];
        plotIND = [5];
    case 'b'
        titleSuffix = [titleSuffix , '_Between'];
        plotIND = [5];
end
switch iN3
    case 'c'
        switch iN2
            case 'o'
                load([baseDir , '/se1_CrossvalIPI_',titleSuffix,'-norm_OLS.mat'])
            case 'r'
                load([baseDir , '/se1_CrossvalIPI_',titleSuffix,'-norm_RR.mat'])
        end
    case 'f'
        load([baseDir , '/se1_fitIPI_',titleSuffix,'-norm_OLS.mat'])
end
switch iN1
    case 'c'
        titleSuffix = 'Chunked';
    case 'r'
        titleSuffix = 'Random';
end
switch iN4
    case 'w'
        titleSuffix = [titleSuffix , ' Within Segment'];
    case 'b'
        titleSuffix = [titleSuffix , ' Between Segment'];
end


% do significance tests:

% 'right' tail testing that A > B
% 'left' tail testing that A < B
subj = unique(Mdl.SN);
count = 1;
for dd = 1:length(dayz)
    ido = ismember(Mdl.Day , dayz{dd});
    
    A = getrow(Mdl , Mdl.modelNum == 2 & ido);
    B = getrow(Mdl , Mdl.modelNum == 5 & ido);
    C = getrow(Mdl , Mdl.modelNum == 8 & ido);
    N = getrow(Mdl , Mdl.modelNum == 1 & ido);
    
    [~ , sig.R2.chunkVprob(count , 1)] = ttest2(A.R2 , B.R2);
    [~ , sig.R2.chunkVfull(count , 1)] = ttest2(A.R2 , C.R2);
    [~ , sig.R2.probVfull(count , 1)] = ttest2(B.R2 , C.R2);
    [~ , sig.R2.NVprob(count , 1)] = ttest2(N.R2 , B.R2);
    [~ , sig.R2.NVchunk(count , 1)] = ttest2(A.R2 , N.R2);
    [~ , sig.R2.NVfull(count , 1)] = ttest2(N.R2 , C.R2);
    
    [~ , sig.rel_R.chunkVprob(count , 1)] = ttest2(A.rel_R , B.rel_R);
    [~ , sig.rel_R.chunkVfull(count , 1)] = ttest2(A.rel_R , C.rel_R);
    [~ , sig.rel_R.probVfull(count , 1)] = ttest2(B.rel_R , C.rel_R);
    [~ , sig.rel_R.NVprob(count , 1)] = ttest2(N.rel_R , B.rel_R);
    [~ , sig.rel_R.NVchunk(count , 1)] = ttest2(A.rel_R , N.rel_R);
    [~ , sig.rel_R.NVfull(count , 1)] = ttest2(N.rel_R , C.rel_R);
    
    [~ , sig.R.chunkVprob(count , 1)] = ttest2(A.R , B.R);
    [~ , sig.R.chunkVfull(count , 1)] = ttest2(A.R , C.R);
    [~ , sig.R.probVfull(count , 1)] = ttest2(B.R , C.R);
    [~ , sig.R.NVprob(count , 1)] = ttest2(N.R , B.R);
    [~ , sig.R.NVchunk(count , 1)] = ttest2(A.R , N.R);
    [~ , sig.R.NVfull(count , 1)] = ttest2(N.R , C.R);
    
    [~ , sig.aic.chunkVprob(count , 1)] = ttest2(A.rel_AIC , B.rel_AIC);
    [~ , sig.aic.chunkVfull(count , 1)] = ttest2(A.rel_AIC , C.rel_AIC);
    [~ , sig.aic.probVfull(count , 1)] = ttest2(B.rel_AIC , C.rel_AIC);
    [~ , sig.aic.NVprob(count , 1)] = ttest2(N.rel_AIC , B.rel_AIC);
    [~ , sig.aic.NVchunk(count , 1)] = ttest2(A.rel_AIC , N.rel_AIC);
    [~ , sig.aic.NVfull(count , 1)] = ttest2(N.rel_AIC , C.rel_AIC);
    
    sig.day(count , 1) = dd;
    sig.h(count , 1) = 13;
    count  = count +1;
end



% Model descriptions
% cleanLabel = {'within/between Chunk', '1st order probability' ,'1st + 2nd order probability' ,'1st + 2nd + 3rd order probability' ,'Full Model'};
cleanLabel = {'Within/Between Segments' ,'1st + 2nd + 3rd Order Transition Probabilities' ,'Full Model'};
% =================== % =================== % =================== % =================== summarize

K = tapply(Mdl , {'SN' , 'Day' , 'modelNum'} , {'R2' , 'nanmean'} , {'rel_R' , 'nanmean'}, {'rel_AIC' , 'nanmean'},{'R' , 'nanmean'}); % average over cv loops


h1 = figure;
% I find color setting and legending  in lineplot a pain, so I generate the data using lineplot, close that figure and pop them into shadeplot and barplot
h = 1;
for dd = 1:length(dayz)
    ido = ismember(K.Day , dayz{dd});
    [xp_cor_rel{dd}(h, :) , pp_cor_rel{dd}(h,:) , ep_cor_rel{dd}(h,:)]  = lineplot(K.modelNum, K.rel_R , 'plotfcn','nanmean' ,  'subset' , ido);
    hold on
    [xp_r2{dd}(h, :) , pp_r2{dd}(h,:) , ep_r2{dd}(h,:)]  = lineplot(K.modelNum, K.R2 , 'plotfcn','nanmean' , 'subset', ido);
    [xp_aic{dd}(h, :) , pp_aic{dd}(h,:) , ep_aic{dd}(h,:)]  = lineplot(K.modelNum, K.rel_AIC , 'plotfcn','nanmean' , 'subset', ido);
    [xp_cor{dd}(h, :) , pp_cor{dd}(h,:) , ep_cor{dd}(h,:)]  = lineplot(K.modelNum, K.R , 'plotfcn','nanmean' ,  'subset' , ido);
    hold on
end

close(h1)


h1 = figure;
h = 1;
ido = ismember(K.modelNum , plotIND);
subplot(311)
lineplot([K.Day    K.modelNum], -K.rel_R , 'plotfcn','nanmean' ,  'subset' , ido ,'style_thickline');
subplot(312)
lineplot([K.Day    K.modelNum], K.rel_AIC , 'plotfcn','nanmean' ,  'subset' , ido ,'style_thickline');
subplot(313)
lineplot([K.Day    K.modelNum], K.R2 , 'plotfcn','nanmean' ,  'subset' , ido ,'style_thickline');
close(h1)


%% do the plotting
switch what
    
    case 's'
        
        % =================== % =================== % =================== % =================== Visualize model R2 comparisons
        figure('color' , 'white')
        
        for dd = 1:length(dayz)
            subplot(1,3,dd)
            h = 1;
            hold on
            eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_r2{dd}(h,plotIND) , ep_r2{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
            hold on
            
            ylabel('R^2')
            set(gca , 'XLim' , [0 length(plotIND)+1] ,'XTick' , [1: length(plotIND)] , 'XTickLabels' , cleanLabel , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' R^2, Days ' , num2str(dayz{dd})])
             
            grid on
        end
        
       
        % =================== % =================== % =================== % =================== Visualize model correlation comparisons
        
        figure('color' , 'white')
        
        for dd = 1:length(dayz)
            subplot(1,3,dd)
            h = 1
            hold on
            eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_cor_rel{dd}(h,plotIND) , ep_cor_rel{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
            hold on
            
            ylabel('correlation')
            set(gca , 'XLim' , [1 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , cleanLabel , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' Model Correlations, Days ' , num2str(dayz{dd})])
            
            grid on
        end
       
        % =================== % =================== % =================== % =================== Visualize model AIC comparisons
       
        figure('color' , 'white')
        
        for dd = 1:length(dayz)
            subplot(1,3,dd)
            h = 1;
            hold on
            eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_aic{dd}(h,plotIND) , ep_aic{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
            hold on
            set(gca , 'XLim' , [1 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , cleanLabel , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1 , 'YLim' , [-2 16])
            title([titleSuffix , ' Model Relative AIC, Days ' , num2str(dayz{dd})])
            
            grid on
        end
        
        
    case 'b'
        % *******************barplot
        % =================== % =================== % =================== % =================== Visualize model R2 comparisons
        xp_ = squeeze(reshape(cell2mat(xp_r2) , size(xp_r2{1} , 1) , size(xp_r2{1} , 2) , length(xp_r2)));
        pp_ = squeeze(reshape(cell2mat(pp_r2) , size(pp_r2{1} , 1) , size(pp_r2{1} , 2) , length(pp_r2)));
        ep_ = squeeze(reshape(cell2mat(ep_r2) , size(ep_r2{1} , 1) , size(ep_r2{1} , 2) , length(ep_r2)));
        
        
        lo = min(min(min(pp_ - ep_)));
        hi = max(max(max(pp_+ep_)));
        ylim = [1.1*lo   1.1*hi];
        figure('color' , 'white')
        barwitherr(ep_(plotIND ,:) , pp_(plotIND ,:))
        set(gca , 'FontSize' , 40 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] ,...
            'YLim' , [-0.05 .2],'XTickLabel' , {'Day 2' , 'Day 3' , 'Day 4'})
        grid on
        title([titleSuffix , ' - Model Relative R^2']);
        
        legend(cleanLabel)
        % =================== % =================== % =================== % =================== Visualize model correlation comparisons
        xp_ = squeeze(reshape(cell2mat(xp_cor_rel) , size(xp_cor_rel{1} , 1) , size(xp_cor_rel{1} , 2) , length(xp_cor_rel)));
        pp_ = squeeze(reshape(cell2mat(pp_cor_rel) , size(pp_cor_rel{1} , 1) , size(pp_cor_rel{1} , 2) , length(pp_cor_rel)));
        ep_ = squeeze(reshape(cell2mat(ep_cor_rel) , size(ep_cor_rel{1} , 1) , size(ep_cor_rel{1} , 2) , length(ep_cor_rel)));
        

        lo = min(min(min(pp_ - ep_)));
        hi = max(max(max(pp_ + ep_)));
        ylim = [1.1*lo   1.1*hi];
        figure('color' , 'white')
%         bar(ep_);
        barwitherr(ep_(plotIND ,:) , pp_(plotIND ,:))
        set(gca , 'FontSize' , 40 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] ,...
            'YLim' , [-.05 .26],'XTickLabel' , {'Day 2' , 'Day 3' , 'Day 4'})
        grid on
        title([titleSuffix , ' - Model Prediction - Output Correlation']);
        
        legend(cleanLabel)
        
       
        
        % =================== % =================== % =================== %
        % =================== Visualize model AIC
        xp_ = squeeze(reshape(cell2mat(xp_aic) , size(xp_aic{1} , 1) , size(xp_aic{1} , 2) , length(xp_aic)));
        pp_ = squeeze(reshape(cell2mat(pp_aic) , size(pp_aic{1} , 1) , size(pp_aic{1} , 2) , length(pp_aic)));
        ep_ = squeeze(reshape(cell2mat(ep_aic) , size(ep_aic{1} , 1) , size(ep_aic{1} , 2) , length(ep_aic)));
        
        lo = min(min(min(-pp_ - ep_)));
        hi = max(max(max(-pp_ + ep_)));
        ylim = [1.1*lo   1.1*hi];
        figure('color' , 'white')
        barwitherr(ep_(plotIND ,:) , -pp_(plotIND ,:))
        set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] ,...
            'YLim' ,ylim,'XTickLabel' , {'Day 2' , 'Day 3' , 'Day 4'})
        grid on
        title([titleSuffix , ' - Model Relative AIC']);
                
        
      
        
end


