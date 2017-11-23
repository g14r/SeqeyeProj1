function [Out,varout] = DigDeeper(what , plt , subj , BT , sequence)



if isempty(plt)
    plt = 1;
end
if isempty(subj)
    subj = 1;
end
if isempty(BT)
    BT = 7;
end
if isempty(sequence)
    sequence = 1;
end


D3 = load('sh3_alldata.mat');
D3 = D3.D;
ST = load('sequenceType.mat');
D3.IPI = diff(D3.peakForceTime , 1 , 2);
%D3.IPI(ismember(D3.ChunkNumb , [126:150]) , 2) = 0;
%----- for every double check if the ipi is different if it's followed by a different finger
id2  = find(ST.Sequence(:,3) == 0);
S2 = ST.Sequence(id2,:);
id3  = find(ST.Sequence(:,3) ~= 0);
fcomb = nchoosek(1:5 , 2);
subs = unique(D3.Subj);

switch what
    case 'double'
        
        for f = 1:5
            for i = 1:5
                id = (ST.Sequence(:,1) == f & ST.Sequence(:,2) == i & ST.Sequence(:,3) == 0 );
                Out.mt_f(:,f,i) = D3.MT(ismember(D3.seqType , find(id)) );
            end
        end
        
        
        for f = 1:5
            for i = 1:5
                id = (ST.Sequence(:,2) == f & ST.Sequence(:,1) == i & ST.Sequence(:,3) == 0 );
                Out.mt_p(:,f,i) = D3.MT(ismember(D3.seqType , find(id)));
            end
            for j = 1:length(fcomb)
                a = squeeze(Out.mt_f(:, f, fcomb(j,1)));
                b = squeeze(Out.mt_f(:, f, fcomb(j,2)));
                [h , Out.p_mt_f(f, j)] = ttest2(a,b);
                
                
                a = squeeze(Out.mt_p(:, f, fcomb(j,1)));
                b = squeeze(Out.mt_p(:, f, fcomb(j,2)));
                [h , Out.p_mt_p(f, j)] = ttest2(a,b);
            end
        end
        
        
        Out.mt_f_average  = squeeze(mean(Out.mt_f , 1));
        Out.mt_p_average  = squeeze(mean(Out.mt_p , 1));
        
        if plt
            figure(1)
            
            subplot(1,3,1);
            imagesc(Out.mt_f_average)
            hold on
            ax = gca;
            title('Average double MT')
            ax.XTick = [1:5];
            ax.YTick = [1:5];
            xlabel('First press')
            ylabel('Second press')
            axis square
            colorbar
            
            subplot(1,3,2);
            temp = squareform(squeeze(mean(Out.p_mt_f , 1)))  + eye(5);
            imagesc(temp);
            hold on
            ax = gca;
            axis square
            title(['Average p-value for IPIs when followed by different fingers, mean = ' num2str(mean(mean(Out.p_mt_f )))])
            ax.XTick = [1:5];
            ax.YTick = [1:5];
            axis square
            colorbar
            
            
            subplot(1,3,3);
            temp = squareform(squeeze(mean(Out.p_mt_p , 1)))  + eye(5);
            imagesc(temp);
            hold on
            ax = gca;
            axis square
            title(['Average p-value for IPIs when preceded by different fingers, mean = ' num2str(mean(mean(Out.p_mt_p )))])
            ax.XTick = [1:5];
            ax.YTick = [1:5];
            axis square
            colorbar
            
        end
        
        
        
    case 'triple'
        
        
        for f = 1:length(id2)
            Out.mt2(:,f) = D3.MT(ismember(D3.seqType ,id2(f)) );
            for i = 1:5
                id = (ismember(ST.Sequence(id3,1:2) , ST.Sequence(id2(f) , 1:2) , 'rows') & ST.Sequence(id3,3) == i);
                Out.mt_po(:,f,i) = D3.MT((D3.seqType == find(id) ));
                Out.ipi1(:,f,i)  = D3.IPI((D3.seqType == find(id)) , 1);
                
                id = (ismember(ST.Sequence(id3,2:3) , ST.Sequence(id2(f) , 1:2) , 'rows') & ST.Sequence(id3,1) == i);
                Out.mt_pre(:,f,i) = D3.MT((D3.seqType == find(id) ));
                Out.ipi2(:,f,i)  = D3.IPI((D3.seqType == find(id) ) , 2);
                
                a = squeeze(Out.ipi1(:,f, i));
                b = squeeze(Out.mt2(:,f));
                [h , Out.p_ipi1(f , i)] = ttest2(a,b);
                
                a = squeeze(Out.ipi2(:,f, i));
                b = squeeze(Out.mt2(:,f));
                [h , Out.p_ipi2(f , i)] = ttest2(a,b);
            end
            
            
            for j = 1:length(fcomb)
                a = squeeze(Out.mt_po(:,f,fcomb(j,1)));
                b = squeeze(Out.mt_po(:,f ,fcomb(j,2)));
                [h , Out.p_mt_po(f , j)] = ttest2(a,b);
                
                a = squeeze(Out.mt_pre(:,f,fcomb(j,1)));
                b = squeeze(Out.mt_pre(:,f ,fcomb(j,2)));
                [h , Out.p_mt_pre(f , j)] = ttest2(a,b);
            end
        end
        
        
        
        if plt
            
            figure
            subplot(2,2,[1 3])
            plot(mean(Out.p_ipi1,2) , 'LineWidth' , 3);
            hold on
            plot(mean(Out.p_ipi2,2) , 'LineWidth' , 3);
            ax = gca;
            title('Average effect on the pair IPIs when followed/preceded by any finger')
            ylabel('Average p-value')
            for p  = 1:length(S2)
                pairLabel{p} = num2str(S2(p,1:2));
            end
            ax.XTick = [1:25];
            ax.XTickLabel = pairLabel;
            ax.XTickLabelRotation=45;
            legend({'Followed by' , 'Preceded by'})
            
            
            subplot(2,2,2)
            imagesc(Out.p_ipi1)
            hold on
            ax = gca;
            ax.YTick = [1:25];
            ax.YTickLabel = pairLabel;
            ax.XTick = [1:5];
            title('Effect on the pair IPIs when followed by a finger')
            colorbar
            
            subplot(2,2,4)
            imagesc(Out.p_ipi2)
            hold on
            ax = gca;
            ax.YTick = [1:25];
            ax.YTickLabel = pairLabel;
            ax.XTick = [1:5];
            title('Effect on the pair IPIs when preceded by a finger')
            colorbar
            
            figure('color' , [1 1 1])
            subplot(3,1,1)
            histogram(D3.IPI(ismember(D3.ChunkNumb , [126:150]) , 1) , 'BinWidth' , 10 , 'FaceColor' , 'b' , 'FaceAlpha' , .5,'EdgeColor' , 'b' , 'EdgeAlpha' , .6)
            hold on
            ax = gca;
            ax.XLim = [0 3000];
            a=  mean(D3.IPI(ismember(D3.ChunkNumb , [126:150]) , 1));
            line([a a] , [0 400] , 'LineWidth' , 3 , 'color' , 'b' , 'LineStyle' , ':')
            title(['Distribution of the IPI in all Double chunks, Average = ' , num2str(a)])
            
            subplot(3,1,2)
            hold on
            histogram(D3.IPI(ismember(D3.ChunkNumb , [1:125]) , 1), 'BinWidth' , 10 , 'FaceColor' , 'r' , 'FaceAlpha' , .5,'EdgeColor' , 'r' , 'EdgeAlpha' , .6)
            title('')
            hold on
            ax = gca;
            ax.XLim = [0 3000];
            a=  mean(D3.IPI(ismember(D3.ChunkNumb , [1:125]) , 1));
            line([a a] , [0 1400] , 'LineWidth' , 3 , 'color' , 'r' , 'LineStyle' , ':')
            title(['Distribution of the first IPI in all Triple chunks, Average = ' , num2str(a)])
            
            subplot(3,1,3)
            hold on
            histogram(D3.IPI(ismember(D3.ChunkNumb , [1:125]) , 2), 'BinWidth' , 10 , 'FaceColor' , 'g' , 'FaceAlpha' , .5,'EdgeColor' , 'g' , 'EdgeAlpha' , .6)
            hold on
            ax = gca;
            ax.XLim = [0 3000];
            a=  mean(D3.IPI(ismember(D3.ChunkNumb , [1:125]) , 2));
            line([a a] , [0 1000] , 'LineWidth' , 3 , 'color' , 'g' , 'LineStyle' , ':')
            title(['Distribution of the second IPI in all Triple chunks, Average = ' , num2str(a)])
            
            % following effect
            figure('color' , 'white')
            temp_f = ST.Sequence(1:150,3);
            fign = 1;
            for f = 0:5
                idf = (temp_f == f);
                id  = ismember(D3.ChunkNumb ,  find(idf)) & D3.good == 1;
                Out.IPI_f{f+1} = D3.IPI(id, 1);
                subplot(6,2,fign)
                histogram(Out.IPI_f{f+1});
                title(['IPI1 when followed by finger ' , num2str(f) , ', Average = ' , num2str(nanmean(Out.IPI_f{f+1}))])
                hold on
                ax = gca;
                ax.XLim = [0 800];
                a=  nanmean(Out.IPI_f{f+1});
                line([a a] , [0 1000] , 'LineWidth' , 3 , 'LineStyle' , ':')
                ax.YGrid = 'on'
                fign = fign + 2;
            end
            
            
            % precedence effect
            temp_p = [ST.Sequence(1:125,1); zeros(25,1)];
            fign = 2;
            for f = 0:5
                idp = (temp_p == f);
                id  = ismember(D3.ChunkNumb ,  find(idp)) & D3.good == 1;
                switch f
                    case 0
                        Out.IPI_p{f+1} = D3.IPI(id, 1);
                    otherwise
                        Out.IPI_p{f+1} = D3.IPI(id, 2);
                end
                subplot(6,2,fign)
                histogram(Out.IPI_p{f+1});
                title(['IPI1 when preceded by finger ' , num2str(f) , ', Average = ' , num2str(nanmean(Out.IPI_p{f+1}))])
                hold on
                ax = gca;
                ax.XLim = [0 800];
                a=  nanmean(Out.IPI_p{f+1});
                line([a a] , [0 1000] , 'LineWidth' , 3 , 'LineStyle' , ':')
                ax.YGrid = 'on'
                fign = fign + 2;
            end
            
        end
        
        
        
        
        %%
    case 'tCompare';          % Generate the preprocessed residual MT and error data for a single subject / sequence
        
        % Generate the right sub-structure
        load('sh4_alldata.mat','D');
        
        
        D.SupposedPress = [D.press0 D.press1 D.press2 D.press3 D.press4 D.press5 D.press6 D.press7 D.press8 D.press9 D.press10];
        D.MadePress     =  [D.response0 D.response1 D.response2 D.response3 D.response4 D.response5 D.response6 D.response7 D.response8 D.response9 D.response10];
        D.Error         = D.SupposedPress~=D.MadePress;
        D.IPI = D.PeakTimeNew; % has ipis for errors anf afterwrds to zero
        
        % BT is Block Type
        % 7   --> Random sequences
        % 1,2 --> Chunk training
        % 3,4 --> Sequence training...
        % loads only one sequence form only one subject
        % model the doubles only in Random sequences - Crossvalidated
        E=load(fullfile('sh3_avrgPattern.mat'));
        E = getrow(E,126:150);
        T = getrow(D,ismember(D.Subj , subj) & ismember(D.BT , BT) & ~D.isError);
        T.DouMadeUpIPI = zeros(size(T.IPI));
        for i  = 1 : length(T.TN)
            for j = 1:size(T.IPI , 2)
                T.DouMadeUpIPI(i , j)  = E.MT(ismember(E.Sequence(:,1:2) , [T.MadePress(i,j) T.MadePress(i,j+1)] , 'rows')) ;
            end
        end
        % If made-up<actual then the t-statistic will be negative. If it is larger, the t-statistic will be positive.
        [h, p_duo      ,ci,  stat_duo]     = ttest2(T.DouMadeUpIPI ,T.IPI ,'tail' , 'both');
        Out. stat_duo = mean(stat_duo.tstat);
        Out.p_duo = mean(p_duo);
        
        % model the triples only in Random sequences - Crossvalidated
        
        ChunkAverageIPI = tapply(D3,{'ChunkNumb'} ,{'IPI', 'mean(x)'},'subset', D3.good == 1);
        T.TriMadeUpIPI = zeros(size(T.IPI));
        for i  = 1 : length(T.TN)
            tempipi = zeros(size(T.IPI , 2)-1 , size(T.IPI , 2));
            for j = 1:size(T.IPI , 2)-1
                Out.tempipi(j,j:j+1) = ChunkAverageIPI.IPI(ismember(ST.Sequence , T.MadePress(i,j:j+2) , 'rows') , :);
            end
            T.TriMadeUpIPI(i , :) = sum(tempipi)./sum((tempipi~=0));
        end
        
        % If made-up<actual then the t-statistic will be negative. If it is larger, the t-statistic will be positive.
        [h, p_tri      ,ci,  stat_tri]     = ttest2(T.TriMadeUpIPI ,T.IPI ,'tail' , 'both');
        Out. tstat_tri = mean(stat_tri.tstat);
        Out.p_tri = mean(p_tri);
        
        %%
    case 'CrossVal reg sh3 tri'
        Out = getrow(D3 , logical(D3.good));
        includeDoub = 0;
        if includeDoub
            varout.X{1} = [ones(275 , 1) eye(275)];
            
            % Model accounting for the mean of firsts and seconds and the single doubles
            varout.X{2} = [ones(275 , 1) [ones(125,1) ; zeros(150,1)] [zeros(125,1) ; ones(25,1); zeros(125,1)] [zeros(150,1) ; ones(125,1)]];
            
            % the Model accounting for the double effects and means
            varout.X{3} = varout.X{2};
            temp = [ST.Sequence(1:150,1:2);ST.Sequence(1:125,2:3)];
            for i = 126:150
                id = ismember(temp , [ST.Sequence(i,1:2)] , 'rows');
                varout.X{3} = [varout.X{3} id];
            end
            
            % saturated Model acounting for the first and second
            varout.X{4} = [ones(275 , 1) [eye(125); zeros(25,125) ; -eye(125)] [zeros(125,1) ; ones(25, 1) ; zeros(125,1)]];
        else
            subs = unique(D3.Subj);
            % the full saturated Model = Noise celing
            varout.X{1} = [ones(250 , 1) eye(250)];
            
            % Model accounting for the mean of firsts and seconds
            varout.X{2} = [ones(250 , 1) [ones(125,1) ; zeros(125,1)] [zeros(125,1) ; ones(125,1)]];
            
            % the Model accounting for the double effects and means
            varout.X{3} = varout.X{2};
            temp = [ST.Sequence(1:125,1:2);ST.Sequence(1:125,2:3)];
            for i = 126:150
                id = ismember(temp , [ST.Sequence(i,1:2)] , 'rows');
                varout.X{3} = [varout.X{3} id];
            end
            %saturated Model acounting for the first ans second
            varout.X{4} = [ones(250 , 1) [eye(125); -eye(125)]];
        end
        
        clear temp
        if includeDoub
            for i = 1:125
                temp(i,1:2) = 'a ';
                temp(i+150,1:2) = 'b ';
            end
            for i = 126:150
                temp(i,1:2) = 'o ';
            end
            ylab = [temp [num2str(ST.Sequence(1:150,1:2));num2str(ST.Sequence(1:125,2:3))]];
        else
            
            for i = 1:125
                temp(i,1:2) = 'a ';
                temp(i+125,1:2) = 'b ';
            end
            ylab = [temp [num2str(ST.Sequence(1:125,1:2));num2str(ST.Sequence(1:125,2:3))]];
        end
        
        
        
        clear cr
        for x = 1:3 ...length(varout.X)
            for s = 1:length(subs)
                Yall1  = pivottable(D3.ChunkNumb, [D3.Day], D3.IPI(:,1) ,'nanmean(x)','subset', D3.good == 1 & D3.Subj == s);
                Yall2  = pivottable(D3.ChunkNumb, [D3.Day], D3.IPI(:,2) ,'nanmean(x)','subset', D3.good == 1 & D3.Subj == s);
                for d = 1:3
                    if includeDoub
                        Y = [Yall1(1:150,:) ; Yall2(1:125,:)];
                    else
                        Y = [Yall1(1:125,:) ; Yall2(1:125,:)];
                    end
                    Y(:,d) = [];
                    y = mean(Y,2);
                    varout.b{x}(:,s,d) = pinv(varout.X{x})*y;
                    ypred = varout.X{x} * varout.b{x}(:,s,d);
                    if includeDoub
                        Y = [Yall1(1:150,:) ; Yall2(1:125,:)];
                    else
                        Y = [Yall1(1:125,:) ; Yall2(1:125,:)];
                    end
                    A = Y(:,d);
                    varout.cr(d,s,x) = corr(A,ypred);
                end
                
            end
            varout.b{x} = squeeze(mean(mean(varout.b{x} , 2) , 3));
        end
        
        varout.Cr = squeeze((nanmean(varout.cr , 1)));
        
        
        for x = 1:length(varout.X)
            figure
            imagesc(varout.X{x});
            hold on
            ax= gca;
            if includeDoub
                ax.YTick = [1:275];
            else
                ax.YTick = [1:250];
            end
            ax.YTickLabel = ylab;
        end
        
        figure('color' , 'white')
        %boxplot(varout.Cr)
        xlab = ones(size(varout.Cr));
        for i = 1:length(varout.b)
            xlab(:,i) = i;
        end
        myboxplot(xlab , varout.Cr,'style_tukey','plotall' , 2);
        ylabel('Correlation')
        hold on
        ax = gca;
        ax.YLim = [0 .8];
        ax.XTickLabel = {'Full Model(250 reg)' 'Mean Model(2 reg)' 'Mean + Double Model(2+25 reg)' 'First or Second IPI(125 reg)'...
            'following' 'preceding' 'following + preceding' 'following +preceding +doubles'};
        ax.XTickLabelRotation=45;
        ax.YGrid = 'on';
        %%
    case 'CrossVal reg sh3 fol-pre'
        
        %following
        varout.X{1} = ones(150 , 1);
        temp_f = ST.Sequence(1:150,3);
        for f = 1:5
            idf = (temp_f == f);
            %idp = find(temp_p == p);
            varout.X{1} = [varout.X{1} idf];
        end
        
        %preceding
        varout.X{2} = ones(150 , 1);
        %temp_f = ST.Sequence(1:125,3);
        temp_p = [ST.Sequence(1:125,1) ; zeros(25,1)];
        for f = 1:5
            %idf = (temp_f == f);
            idp  = (temp_p == f);
            varout.X{2} = [varout.X{2}  idp];
        end
        
        clear temp
        for i = 1:125
            temp(i,1:2) = 'a ';
            temp(i+125,1:2) = 'b ';
        end
        ylab = [temp [num2str(ST.Sequence(1:125,1:2));num2str(ST.Sequence(1:125,2:3))]];
        
        clear cr
        for x = 1
            for s = 1:length(subs)
                Yall1  = pivottable(D3.ChunkNumb, [D3.Day], D3.IPI(:,1) ,'nanmean(x)','subset', D3.good == 1 & D3.Subj == s);
                for d = 1:3
                    Y = Yall1(1:150,:);
                    Y(:,d) = [];
                    y = mean(Y,2);
                    varout.b{x}(:,s,d) = pinv(varout.X{x})*y;
                    ypred = varout.X{x} * varout.b{x}(:,s,d);
                    Y = Yall1(1:150,:);
                    A = Y(:,d);
                    varout.cr(d,s,x) = corr(A,ypred);
                end
                
            end
            varout.b{x} = squeeze(mean(mean(varout.b{x} , 2) , 3));
        end
        
        for x = 2
            for s = 1:length(subs)
                Yall1  = pivottable(D3.ChunkNumb, [D3.Day], D3.IPI(:,1) ,'nanmean(x)','subset', D3.good == 1 & D3.Subj == s);
                Yall2  = pivottable(D3.ChunkNumb, [D3.Day], D3.IPI(:,2) ,'nanmean(x)','subset', D3.good == 1 & D3.Subj == s);
                for d = 1:3
                    Y = [Yall1(126:150,:) ; Yall2(1:125,:)];
                    Y(:,d) = [];
                    y = mean(Y,2);
                    varout.b{x}(:,s,d) = pinv(varout.X{x})*y;
                    ypred = varout.X{x} * varout.b{x}(:,s,d);
                    Y = [Yall1(126:150,:) ; Yall2(1:125,:)];
                    A = Y(:,d);
                    varout.cr(d,s,x) = corr(A,ypred);
                end
                
            end
            varout.b{x} = squeeze(mean(mean(varout.b{x} , 2) , 3));
        end
        
        varout.Cr = squeeze(nanmean(varout.cr , 1));
        
        
        for x = 1:length(varout.X)
            figure
            imagesc(varout.X{x});
            hold on
            ax= gca;
            ax.YTick = [1:250];
            ax.YTickLabel = ylab;
        end
        
        figure('color' , 'white')
        bar(varout.Cr)
        ylabel('Correlation')
        hold on
        ax = gca;
        ax.YLim = [0 .8];
        ax.XTickLabel = {'following' 'preceding' 'following + preceding' 'following +preceding +doubles'};
        ax.XTickLabelRotation=45;
        ax.YGrid = 'on';
        %%            
    case 'MixedAnova'
        Out = getrow(D3 , logical(D3.good));
        Out.firstTriDoub = zeros(size(Out.TN));
        Out.lastTriDoub  = zeros(size(Out.TN));
        %tripels
        for i  = 126:150
            id = ismember(Out.ChunkNumb , [1:125]) & ismember(Out.Chunk(: , 1:2) , ST.Sequence(i,1:2) , 'rows');
            Out.firstTriDoub(id , :) = i;
            id = ismember(Out.ChunkNumb , [1:125]) & ismember(Out.Chunk(: , 2:3) , ST.Sequence(i,1:2) , 'rows');
            Out.lastTriDoub(id , :) = i;           
        end
        % important notes about Anova mixed. In the results in will give
        % you a df1 and a df2. The test statistic is calculated by the
        % division of two khi squared distributoins. so the df1 is the
        % degrees of freedome for the numerator distibution, so always one
        % less than the number of conditions in each factor (in this case 1 for the intercept 4
        % for fingers and 24 for duos and 4*24 for the interactions).
        % However df2 is the DoF for the denominator distribution, and is always one less than the 
        % number of subjects for the intercept (6) and 6 * df1 for each
        % factor : 6*4 and 6*24. also remember you can compare the f-values
        % for different degrees of freedome. you can only compare f-values
        % if the degrees of freedome are equal.

        varout.Follow_Anova = anovaMixed(Out.IPI(:,1),Out.Subj,'within',[Out.response3,Out.firstTriDoub],{'FollowingFinger','DuoNumber'},'subset',ismember(Out.ChunkNumb ,[1:125]) & ~isnan(Out.IPI(:,1)),'intercept',1)  ;      
        varout.Preced_ANova = anovaMixed(Out.IPI(:,2),Out.Subj,'within',[Out.response1,Out.lastTriDoub],{'PrecedingFinger','DuoNumber'},'subset',ismember(Out.ChunkNumb ,[1:125]) & ~isnan(Out.IPI(:,2)),'intercept',1);        
end



