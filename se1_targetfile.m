function se1_targetfile(SubCode, GroupCode, WantTheStars,genchunks, CMB)


% genchunks = 0;
% WantTheStars = 0;
% SubCode = 'SZ';
% GroupCode= 1;
% provode CMB from the subject's target file folder only if you already have the chunks and need to  produce more sequences
% baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/TargetFiles';
baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye1/TargetFiles';
cd(baseDir)
Fname  = [SubCode , num2str(GroupCode)  , '_tgtFiles'];
mkdir(Fname)


if genchunks
    CMB = se1_getChunks(4,2,2,2, 0);   
end
if GroupCode == 1
    CMB = rmfield(CMB, 'Group2');
    CMB.Group = CMB.Group1;   % just to have consistent filed name across groups
    CMB = rmfield(CMB, 'Group1');
    for i  = 1:length(CMB.Group)
        temp = CMB.DesiredCnhkargmnt(CMB.Group(i),:);
        temp1 = [];
        for j = 1:length(find(temp))
            temp1 = [temp1 1:temp(j)];
        end
        CMB.Seq2Chunk(i , :) = temp1;
    end
elseif GroupCode == 2
    CMB = rmfield(CMB, 'Group1');
    CMB.Group = CMB.Group2;   % just to have consistent filed name across groups
    CMB = rmfield(CMB, 'Group2');
    for i  = 1:length(CMB.Group)
        temp = CMB.DesiredCnhkargmnt(CMB.Group(i),:);
        temp1 = [];
        for j = 1:length(find(temp))
            temp1 = [temp1 1:temp(j)];
        end
        CMB.Seq2Chunk(i , :) = temp1;
    end
end
save([Fname , '/' , SubCode , num2str(GroupCode) , '_CMB'] , 'CMB')


% make target file
% switch what
%     case 'targetfile'

rng('shuffle');


NumofChunks = 6;
SequenceLength = size(CMB.Seq2Chunk , 2);
Chunks = cell2mat(CMB.Chunks);
ChunkNumb = [102 202 103 203 104 204];


Cstar = {'' , '**' , '***' , '****'};
repeating  = 1; % determins that every sequence in the CLAT and Intermixed blocks happens twice

FT = 1:6;

OrderFields = {'seqNumb','FT','press1','press2','press3','press4','press5','press6','press7','press8','press9','press10','press11','press12','press13','press14','hand','cueS','cueC','cueP','iti','sounds' , 'Horizon' , 'StimTimeLim'};



%% Chunk Training Block
%WantTheStars = 0;
Trials = 1:66;

if WantTheStars
    NumStarTrials = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 9 9 9];
    
    clear ChunkTrain
    for e = 1:length(NumStarTrials)
        ChunkTrain.cueS = cellstr(repmat('£',length(Trials),1));
        ChunkTrain.FT = 2 * ones(length(Trials),1);
        ChunkTrain.iti(1:length(Trials),:) = 500;
        ChunkTrain.hand(1:length(Trials),:) = 2;
        ChunkTrain.sounds(1:length(Trials),:) = 1;
     
        % making sure tht all the chunks are ocurring equal times in each chunk trianing block
        X = [];
        for rep = 1: ceil((length(Trials)/ (NumStarTrials(e) + 1))/length(NumofChunks))
            X = [X ; sample_wor(NumofChunks,length(NumofChunks))];
        end
        X = kron(X, ones(NumStarTrials(e) + 1 , 1));
        X = X(1:length(Trials));
        % Indeices in where the stras should apear as 0
        starInd = repmat([1 ; zeros(NumStarTrials(e) ,1)] , length(Trials)/(NumStarTrials(e)+1) ,1);
        for i = 1:length(starInd)
            if ~starInd(i)
                NumStars = sum((Chunks(X(i),:) ~= 0 ));
                ChunkTrain.cueP{i,1} = Cstar{NumStars};
                ChunkTrain.Horizon(i,:) = NumStars - 1;
            else
                NumStars = sum((Chunks(X(i),:) ~= 0 ));
                ChunkTrain.cueP{i,1} = char(regexprep(cellstr(num2str(Chunks(X(i),1:NumStars))),'\s',''));
            end
            for press= 1:NumStars
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = Chunks(X(i),press);'];
                eval(comnd);
            end
            for press= NumStars + 1 : SequenceLength
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = 0;'];
                eval(comnd);
            end
            % Chunk length arrangement number
            ChunkTrain.seqNumb(i,1) = ChunkNumb(X(i));
            ChunkTrain.cueC(i,:) = {'£'} ;
        end
        name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CTs' , num2str(NumStarTrials(e)) , '_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkTrain,OrderFields));
        
        clear x name
        clear ChunkTrain
    end
else
    %             NumStarTrials = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 9 9 9];
    %             OrderFields = {'seqNumb','FT','press1','press2','press3','press4','press5','press6','press7','press8','press9','press10','press11','press12','press13','press14','hand','cueS','cueC','cueP','iti','sounds'};
    
    clear ChunkTrain
    StimTimeLim = zeros(length(Trials) , 1);
    for e = 1:20
        ChunkTrain.cueS = cellstr(repmat('£',length(Trials),1));
        ChunkTrain.FT = 2 * ones(length(Trials),1);
        ChunkTrain.iti(1:length(Trials),:) = 500;
        ChunkTrain.hand(1:length(Trials),:) = 2;
        ChunkTrain.sounds(1:length(Trials),:) = 1;
        ChunkTrain.StimTimeLim = StimTimeLim;
        
        % making sure that all the chunks are ocurring equal times in each chunk trianing block
        X = [];
        for rep = 1: length(Trials)/NumofChunks
            X = [X  randperm(NumofChunks)];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible chunks in one run is not detected
        %                 X = kron(X, ones(NumStarTrials(e) + 1 , 1));
        %                 X = X(1:length(Trials));
        % Indeices in where the stras should apear as 0
        %                 starInd = repmat([1 ; zeros(NumStarTrials(e) ,1)] , length(Trials)/(NumStarTrials(e)+1) ,1);
        for i = 1:length(Trials)
            Numdigs = sum((Chunks(X(i),:) ~= 0 ));
            ChunkTrain.cueP{i,1} = char(regexprep(cellstr(num2str(Chunks(X(i),1:Numdigs))),'\s',''));
            
            for press= 1:Numdigs
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = Chunks(X(i),press);'];
                eval(comnd);
            end
            ChunkTrain.Horizon(i,:) = Numdigs - 1;
            for press= Numdigs + 1 : SequenceLength
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = 0;'];
                eval(comnd);
            end
            % Chunk length arrangement number
            ChunkTrain.seqNumb(i,1) = ChunkNumb(X(i));
            ChunkTrain.cueC(i,:) = {'£'} ;
        end
        name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CT_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkTrain,OrderFields));
        
        clear x name
        clear ChunkTrain
    end
    
end

%% Chunk Arrangemnet Training Block

if ~repeating

    CLA = CMB.Group;

    h = [1:12]';
    h =  kron(h, ones(4 , 1));
    % Generate 30 Full Horizons, and then 4 of each for 3 : 12
    % Horizon = [(SequenceLength - 1) * ones(30 , 1) ; h];
    
    % just generate full horizons for now
    Horizon = [(SequenceLength - 1) * ones(30 , 1)];
    Trials = 1:30;
    StimTimeLim = zeros(length(Trials) , 1);
    
    clear ChunkArrangeLearn
    
    for e= 1:length(Horizon)
        ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        for rep = 1: length(Trials) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        ChunkArrangeLearn.seqNumb = X;
        ChunkArrangeLearn.FT(1:length(Trials),1) = 2;
        for i = 1:length(X)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                chnkind = randi(2);
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            ChunkArrangeLearn.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(i,:) ={'£'};
            ChunkArrangeLearn.cueS(i,:) ={'£'};
            ChunkArrangeLearn.iti(1:i,:) = 500;
            ChunkArrangeLearn.hand(i,:) = 2;
            ChunkArrangeLearn.sounds(i,:) = 1;
            ChunkArrangeLearn.Horizon(i,:) = Horizon(e);
            
            
            for press= 1:14
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
        
        name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CLAT_h' ,num2str(Horizon(e)), '_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkArrangeLearn,OrderFields));
        
        clear x
        clear ChunkArrangeLearn
    end
end
%% Test blocks - Random Sequences
clear RandomSeq
ElimChunkStart  = 0; % 1 means that random seqs that start with a known chunk will be elimminated
Trials = 1:30;
StimTimeLim = zeros(length(Trials) , 1); 
h = [1:12]';
h =  kron(h, ones(4 , 1));
% Generate 30 Full Horizons, and then 4 of each for 3 : 12
% Horizon = [(SequenceLength - 1) * ones(30 , 1) ; h];

% just generate full horizons for now
Horizon = [(SequenceLength - 1) * ones(10 , 1)];



for e= 1:length(Horizon)
    RandomSeq.StimTimeLim = StimTimeLim; 
    RandomSeq.seqNumb = zeros(length(Trials) , 1);
    RandomSeq.FT(1:length(Trials),1) = 2;
    RandomSeq.Horizon(1:length(Trials),:) = Horizon(e);
    for i = 1:length(Trials)
        %                 ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
        if ElimChunkStart
            seq = Chunks(1,1:2);
            while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks) 
                seq = sample_wor([1:5],1,14);
            end
        else
            seq = sample_wor([1:5],1,14);
            if length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) < length(Chunks) 
                RandomSeq.seqNumb(i) = 7; % represents random sequences that start with a known chunk
            end
        end
        RandomSeq.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
        RandomSeq.cueC(i,:) ={'£'};
        RandomSeq.cueS(i,:) ={'£'};
        RandomSeq.iti(1:i,:) = 500;
        RandomSeq.hand(i,:) = 2;
        RandomSeq.sounds(i,:) = 1;
        
        for press= 1:14
            comnd  = [' RandomSeq.press' , num2str(press) , '(i,1) = seq(press);'];
            eval(comnd);
        end
    end
    
    name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_RAND_h' ,num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
    dsave(name,orderfields(RandomSeq,OrderFields));
    
    clear x
    clear RandomSeq
end

%% The inrtermixed CLAT and Random blocks 
if ~repeating

    CLA = CMB.Group;
    Trials = 1:36;
    StimTimeLim = zeros(length(Trials) , 1);
    RandProportion  = 1/2; % Can be 1/3 , 1/2 or 2/3 - portion of the sequences that are going to be random
    ElimChunkStart  = 1; % 1 means that random seqs that start with a known chunk will be elimminated
    h = [1:12]';
    h =  kron(h, ones(4 , 1));
    % Generate 30 Full Horizons, and then 4 of each for 3 : 12
    % Horizon = [(SequenceLength - 1) * ones(20 , 1) ; h];
    
    % just generate full horizons for now
    Horizon = [(SequenceLength - 1) * ones(20 , 1)];
    clear Intermixed
    for e= 1:length(Horizon)
        Intermixed.StimTimeLim = StimTimeLim;
        Intermixed.Horizon(1:length(Trials),1) = Horizon(e);
        RandIndex = randperm(length(Trials));
        RandIndex = RandIndex(1: RandProportion * length(Trials));
        
        CLATIndex = Trials(~ismember(Trials , RandIndex));
        % first create the CLAT sequecnces
        X = [];
        for rep = 1: (length(Trials) - length(RandIndex)) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        Intermixed.seqNumb  = zeros(length(Trials) , 1);
        Intermixed.seqNumb(CLATIndex) = X;
        Intermixed.FT(1:length(Trials),1) = 2;
        Intermixed.iti(1:length(Trials),:) = 500;
        Intermixed.hand(1:length(Trials),:) = 2;
        Intermixed.sounds(1:length(Trials),:) = 1;
        
        for i = 1:length(CLATIndex)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                chnkind = randi(2);
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            Intermixed.cueP{CLATIndex(i),:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            Intermixed.cueC(CLATIndex(i),:) ={'£'};
            Intermixed.cueS(CLATIndex(i),:) ={'£'};
            
            for press= 1:14
                comnd  = [' Intermixed.press' , num2str(press) , '(CLATIndex(i),1) = seq(press);'];
                eval(comnd);
            end
        end
        
        for i = 1:length(RandIndex)
            %                 ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            if ElimChunkStart
                seq = Chunks(1,1:2);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <length(Chunks) + 1
                    seq = sample_wor([1:5],1,14);
                end
            else
                seq = sample_wor([1:5],1,14);
                if length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) < length(Chunks)
                    Intermixed.seqNumb(i) = 7; % represents random sequences that start with a known chunk
                end
            end
            Intermixed.cueP{RandIndex(i),:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            Intermixed.cueC(RandIndex(i),:) ={'£'};
            Intermixed.cueS(RandIndex(i),:) ={'£'};
            for press= 1:14
                comnd  = [' Intermixed.press' , num2str(press) , '(RandIndex(i),1) = seq(press);'];
                eval(comnd);
            end
        end
        switch RandProportion
            case 1/3
                name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_IM_TR_h', num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
            case 1/2
                name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_IM_HR_h', num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
            case 2/3
                name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_IM_2TR_h', num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
        end
        
        dsave(name,orderfields(Intermixed,OrderFields));
        
        clear x Intermixed
    end
end


%% Chunk Arrangemnet Training Block, with two repetitions of every sequence
if repeating
    CLA = CMB.Group;
    h = [1:12]';
    h =  kron(h, ones(4 , 1));
    % Generate 30 Full Horizons, and then 4 of each for 3 : 12
    % Horizon = [(SequenceLength - 1) * ones(30 , 1) ; h];
    
    % just generate full horizons for now
    Horizon = [(SequenceLength - 1) * ones(30 , 1)]; 
    
    Trials = 1:36;
    StimTimeLim = zeros(length(Trials) , 1);
    
    clear ChunkArrangeLearn
    
    for e= 1:length(Horizon)
        ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        for rep = 1: .5*length(Trials) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        
        ChunkArrangeLearn.FT(1:length(Trials),1) = 2;
        cn = 1;
        for i = 1:length(X)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                chnkind = randi(2);
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            ChunkArrangeLearn.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(cn,:) ={'£'};
            ChunkArrangeLearn.cueS(cn,:) ={'£'};
            ChunkArrangeLearn.iti(1:cn,:) = 500;
            ChunkArrangeLearn.hand(cn,:) = 2;
            ChunkArrangeLearn.sounds(cn,:) = 1;
            ChunkArrangeLearn.Horizon(cn,:) = Horizon(e);
            ChunkArrangeLearn.seqNumb(cn , :) = X(i);
            for press= 1:14
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn = cn + 1;
            
            ChunkArrangeLearn.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(cn,:) ={'£'};
            ChunkArrangeLearn.cueS(cn,:) ={'£'};
            ChunkArrangeLearn.iti(1:cn,:) = 500;
            ChunkArrangeLearn.hand(cn,:) = 2;
            ChunkArrangeLearn.sounds(cn,:) = 1;
            ChunkArrangeLearn.Horizon(cn,:) = Horizon(e);
            ChunkArrangeLearn.seqNumb(cn , :) = X(i);
            
            for press= 1:14
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn  = cn + 1;
        end
        
        name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CLAT2_h' ,num2str(Horizon(e)), '_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkArrangeLearn,OrderFields));
        
        clear x
        clear ChunkArrangeLearn
    end
    
end

%% The inrtermixed CLAT and Random blocks, with two repetitions of every CLAT sequence 

if repeating
    RandProportion  = 1/2; % Can be 1/3 , 1/2 or 2/3 - portion of the sequences that are going to be random
    ElimChunkStart  = 1; 
    CLA = CMB.Group;
    h = [1:12]';
    h =  kron(h, ones(4 , 1));
    % Generate 30 Full Horizons, and then 4 of each for 3 : 12
    % Horizon = [(SequenceLength - 1) * ones(30 , 1) ; h];
    
    % just generate full horizons for now
    Horizon = [(SequenceLength - 1) * ones(30 , 1)];
    
    
    Trials = 1:12;
    StimTimeLim = zeros(length(Trials) , 1);
    
    clear ChunkArrangeLearn
    clear RandomSeq
    
    
    % Test blocks - Random Sequences

    ElimChunkStart  = 1; % 1 means that random seqs that start with a known chunk will be elimminated
    StimTimeLim = zeros(length(Trials) , 1);
    h = [1:12]';
    h =  kron(h, ones(4 , 1));
    % Generate 30 Full Horizons, and then 4 of each for 3 : 12

    
    
    for e= 1:length(Horizon)
        
        
        ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        
        for rep = 1: .5*length(Trials) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        
        ChunkArrangeLearn.FT(1:length(Trials),1) = 2;
        cn = 1;
        uniqidx = 1;
        for i = 1:length(X)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                chnkind = randi(2);
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            ChunkArrangeLearn.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(cn,:) ={'£'};
            ChunkArrangeLearn.cueS(cn,:) ={'£'};
            ChunkArrangeLearn.iti(1:cn,:) = 500;
            ChunkArrangeLearn.hand(cn,:) = 2;
            ChunkArrangeLearn.sounds(cn,:) = 1;
            ChunkArrangeLearn.Horizon(cn,:) = Horizon(e);
            ChunkArrangeLearn.seqNumb(cn , :) = X(i);
            ChunkArrangeLearn.unique(cn , :) = uniqidx;
            for press= 1:14
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn = cn + 1;
            
            ChunkArrangeLearn.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(cn,:) ={'£'};
            ChunkArrangeLearn.cueS(cn,:) ={'£'};
            ChunkArrangeLearn.iti(1:cn,:) = 500;
            ChunkArrangeLearn.hand(cn,:) = 2;
            ChunkArrangeLearn.sounds(cn,:) = 1;
            ChunkArrangeLearn.Horizon(cn,:) = Horizon(e);
            ChunkArrangeLearn.seqNumb(cn , :) = X(i);
            ChunkArrangeLearn.unique(cn , :) = uniqidx;
            
            for press= 1:14
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn  = cn + 1;
            uniqidx = uniqidx +1;
        end
        
        %name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CLAT2_h' ,num2str(Horizon(e)), '_B' , num2str(e) , '.tgt'];
        %dsave(name,orderfields(ChunkArrangeLearn,OrderFields));
        
        clear x
        
        
        
        
        
        allTrials = length(Trials)/(1 - RandProportion);
        rTrials = 1:allTrials*RandProportion;
        RandomSeq.StimTimeLim = StimTimeLim;
        RandomSeq.seqNumb = zeros(length(rTrials) , 1);
        RandomSeq.FT(1:length(rTrials),1) = 2;
        RandomSeq.Horizon(1:length(rTrials),:) = Horizon(e);
        cn  = 1;
        for i = 1:length(rTrials)/2
            %                 ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            if ElimChunkStart
                seq = Chunks(1,1:2);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks)
                    seq = sample_wor([1:5],1,14);
                end
            else
                seq = sample_wor([1:5],1,14);
            end
            
            if ~ElimChunkStart
                if length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) < length(Chunks)
                    RandomSeq.seqNumb(cn) = 7; % represents random sequences that start with a known chunk
                end
            end
            RandomSeq.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            RandomSeq.cueC(cn,:) ={'£'};
            RandomSeq.cueS(cn,:) ={'£'};
            RandomSeq.iti(1:cn,:) = 500;
            RandomSeq.hand(cn,:) = 2;
            RandomSeq.sounds(cn,:) = 1;
            RandomSeq.unique(cn , :) = uniqidx;
            for press= 1:14
                comnd  = [' RandomSeq.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn = cn+1;
            
            
            if ~ElimChunkStart
                if length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) < length(Chunks)
                    RandomSeq.seqNumb(cn) = 7; % represents random sequences that start with a known chunk
                end
            end
            RandomSeq.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            RandomSeq.cueC(cn,:) ={'£'};
            RandomSeq.cueS(cn,:) ={'£'};
            RandomSeq.iti(1:cn,:) = 500;
            RandomSeq.hand(cn,:) = 2;
            RandomSeq.sounds(cn,:) = 1;
            RandomSeq.unique(cn , :) = uniqidx;
            for press= 1:14
                comnd  = [' RandomSeq.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn = cn+1;
            uniqidx = uniqidx +1;
        end
        
        %     name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_RAND_h' ,num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
        %     dsave(name,orderfields(RandomSeq,OrderFields));
        uniqidx = uniqidx -1;
        temp = randperm(uniqidx);
        temp_Intermixed = addstruct(RandomSeq , ChunkArrangeLearn);
        Intermixed = getrow(temp_Intermixed , temp_Intermixed.FT == 0);
        for m = 1:length(temp)
            T = getrow(temp_Intermixed , temp_Intermixed.unique == temp(m));
            Intermixed = addstruct(Intermixed , T);
        end
        
        
        clear x
        clear RandomSeq
        clear ChunkArrangeLearn
        
        
        switch RandProportion
            case 1/3
                name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_IM2_TR_h', num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
            case 1/2
                name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_IM2_HR_h', num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
            case 2/3
                name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_IM2_2TR_h', num2str(Horizon(e)) , '_B' , num2str(e) , '.tgt'];
        end
        
        Intermixed = rmfield(Intermixed , 'unique');
        
        dsave(name,orderfields(Intermixed,OrderFields));
        
        
        clear Intermixed
        
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Merge






