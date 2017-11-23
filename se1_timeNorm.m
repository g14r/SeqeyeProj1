function N = se1_timeNorm(Dall , what, subjnum, day , rep)

prefix = 'se1_';
baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye1/analyze';
%baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/se1_data/analyze';
subj_name = {'SZ1','JN2' ,'SP1','AT1','DW1','KL1','JG1','GP1','SK1' ,'NM1','VR1','PB1'};
load([baseDir , '/CMB.mat'])
%load([baseDir , '/se1_all.mat'])
Days  = {1 ,2 ,3 ,4 ,[2:4] ,[3:4] [1:4]};
window  = 10;
switch what
    case 'all_Sub'
        
        for sub = 1:length(subj_name)
            for d  = 1:length(Days)
                switch d
                    case {1 2 3 4}
                        ANA = getrow(Dall , Dall.SN == sub & ismember(Dall.Day , Days{d}) & ismember(Dall.Rep , rep) & ~Dall.isError & Dall.isgood & ismember(Dall.seqNumb , [0:6]));
                        %%              first calculate then time normalize the velocities
                        for s = 0:6
                            id = find(ANA.seqNumb == s);
                            N.norm(d,sub).eye{s+1}   = [];
                            N.norm(d,sub).eyebntn{s+1} = [];
                            N.norm(d,sub).press{s+1}  = [];
                            N.norm(d,sub).eyeveloc{s+1} = [];
                            N.norm(d,sub).prsveloc{s+1} = [];
                            for i = 1:length (id)
                                if sum(isnan(ANA.xEyePosDigit{id(i)})) < .4*length(ANA.xEyePosDigit{id(i)})
                                    idd   = linspace(1 , ANA.AllPressIdx(id(i),14)-ANA.AllPressIdx(id(i),1)+1 , 1000);
                                    idde = floor(linspace(ANA.AllPressIdx(id(i),1) , ANA.AllPressIdx(id(i),14) , 1001));
                                    for j = 1:length(idde)-1
                                        temp1(j) = nanmedian(ANA.xEyePosDigit{id(i)}(idde(j) : idde(j+1)));
                                        temp2(j) = nanmedian(ANA.PressTimeSeries{id(i)}(idde(j) : idde(j+1)));
                                        temp3(j) = nanmedian(ANA.xEyeAngVelocity{id(i)}(idde(j) : idde(j+1)));
                                        temp4(j) = nanmedian(ANA.pressVelocity{id(i)}(idde(j) : idde(j+1)));
                                    end
                                    N.norm(d,sub).eye{s+1}      = [N.norm(d,sub).eye{s+1}      ; inpaint_nans(temp1,3)];
                                    N.norm(d,sub).eyebntn{s+1}  = [N.norm(d,sub).eyebntn{s+1}  ; [ANA.BN(id(i)) ANA.TN(id(i))]];
                                    N.norm(d,sub).press{s+1}    = [N.norm(d,sub).press{s+1}    ; inpaint_nans(temp2)];
                                    N.norm(d,sub).eyeveloc{s+1} = [N.norm(d,sub).eyeveloc{s+1} ; inpaint_nans(temp3)];
                                    N.norm(d,sub).prsveloc{s+1} = [N.norm(d,sub).prsveloc{s+1} ; inpaint_nans(temp4)];
                                    
                                end
                                
                            end
                            N.norm(d,sub).eyePattern(s+1, :)  = nanmedian(N.norm(d,sub).eye{s+1});
                            N.norm(d,sub).pressPattern(s+1,:) = nanmedian(N.norm(d,sub).press{s+1});
                            N.norm(d,sub).eyevelocPattern(s+1,:) = nanmedian(N.norm(d,sub).eyeveloc{s+1});
                            N.norm(d,sub).prsvelocPattern(s+1,:) = nanmedian(N.norm(d,sub).prsveloc{s+1});
                            [sub d s]
                        end
                        
                        
                        
                        
                        %%                first time normalize then calculate velocity
                                        for s = 0:6
                                            id = find(ANA.seqNumb == s);
                                            N.norm(d,sub).eye{s+1}   = [];
                                            N.norm(d,sub).eyebntn{s+1} = [];
                                            N.norm(d,sub).press{s+1}  = [];
                                            N.norm(d,sub).eyeveloc{s+1} = [];
                                            N.norm(d,sub).prsveloc{s+1} = [];
                                            for i = 1:length (id)
                                                if sum(isnan(ANA.xEyePosDigit{id(i)})) < .2*length(ANA.xEyePosDigit{id(i)})
                                                    idd   = linspace(1 , ANA.AllPressIdx(id(i),14)-ANA.AllPressIdx(id(i),1)+1 , 1000);
                                                    idde = floor(linspace(ANA.AllPressIdx(id(i),1) , ANA.AllPressIdx(id(i),14) , 1001));
                                                    for j = 1:length(idde)-1
                                                        temp1(j) = nanmedian(ANA.xEyePosDigit{id(i)}(idde(j) : idde(j+1)));
                                                    end
                                                    N.norm(d,sub).eye{s+1}   = [N.norm(d,sub).eye{s+1} ;inpaint_nans(temp1,3)];
                        
                                                    temp2 = inpaint_nans(ANA.PressTimeSeries{id(i)}(ANA.AllPressIdx(id(i),1):ANA.AllPressIdx(id(i),14)));
                                                    temp2 = interp1([1 :ANA.AllPressIdx(id(i),14)-ANA.AllPressIdx(id(i),1)+1] , temp2 , idd);
                                                    N.norm(d,sub).press{s+1}   = [N.norm(d,sub).press{s+1} ; temp2];
                        
                                                    N.norm(d,sub).eyeveloc{s+1} = [N.norm(d,sub).eyeveloc{s+1} ; temp1(window+1:end) -  temp1(1:end-window)];
                                                    N.norm(d,sub).prsveloc{s+1} = [N.norm(d,sub).prsveloc{s+1} ; temp2(window+1:end) -  temp2(1:end-window)];
                                                    N.norm(d,sub).eyebntn{s+1}  = [N.norm(d,sub).eyebntn{s+1} ; [ANA.BN(id(i)) ANA.TN(id(i))]];
                                                end
                                            end
                                            N.norm(d,sub).eyePattern(s+1, :)  = nanmedian(N.norm(d,sub).eye{s+1})/2.3;
                                            N.norm(d,sub).pressPattern(s+1,:) = nanmedian(N.norm(d,sub).press{s+1});
                                            N.norm(d,sub).eyevelocPattern(s+1,:) = nanmedian(N.norm(d,sub).eyeveloc{s+1});
                                            N.norm(d,sub).prsvelocPattern(s+1,:) = nanmedian(N.norm(d,sub).prsveloc{s+1});
                                        end
                    otherwise
                        for s  = 0:6
                            N.norm(d,sub).eye{s+1}   = [];
                            N.norm(d,sub).eyebntn{s+1} = [];
                            N.norm(d,sub).press{s+1}  = [];
                            N.norm(d,sub).eyeveloc{s+1} = [];
                            N.norm(d,sub).prsveloc{s+1} = [];
                            for d1 = 1:length(Days{d})
                                N.norm(d,sub).eye{s+1}      = [N.norm(d,sub).eye{s+1}      ; N.norm(Days{d}(d1),sub).eye{s+1}];
                                N.norm(d,sub).eyebntn{s+1}  = [N.norm(d,sub).eyebntn{s+1}  ; N.norm(Days{d}(d1),sub).eyebntn{s+1}];
                                N.norm(d,sub).press{s+1}    = [N.norm(d,sub).press{s+1}    ; N.norm(Days{d}(d1),sub).press{s+1}];
                                N.norm(d,sub).eyeveloc{s+1} = [N.norm(d,sub).eyeveloc{s+1} ; N.norm(Days{d}(d1),sub).eyeveloc{s+1}];
                                N.norm(d,sub).prsveloc{s+1} = [N.norm(d,sub).prsveloc{s+1} ; N.norm(Days{d}(d1),sub).prsveloc{s+1}];
                            end
                            N.norm(d,sub).eyePattern(s+1, :)  = nanmedian(N.norm(d,sub).eye{s+1});
                            N.norm(d,sub).pressPattern(s+1,:) = nanmedian(N.norm(d,sub).press{s+1});
                            N.norm(d,sub).eyevelocPattern(s+1,:) = nanmedian(N.norm(d,sub).eyeveloc{s+1});
                            N.norm(d,sub).prsvelocPattern(s+1,:) = nanmedian(N.norm(d,sub).prsveloc{s+1});
                        end
                end
                
            end
        end
        suball = sub  +1;
        for d  = 1:length(Days)
            N.norm(d,suball).eye      = cell(1,7);
            N.norm(d,suball).eyebntn  = cell(1,7);
            N.norm(d,suball).press    = cell(1,7);
            N.norm(d,suball).eyeveloc = cell(1,7);
            N.norm(d,suball).prsveloc = cell(1,7);
        end
        for s  = 0:6
            for d  = 1:length(Days)
                N.norm(d,suball).subj{s+1} = [];
                for subj = [1 3 4 5 6 7 8 9 10 11 12]
                    N.norm(d,suball).eye{s+1}      = [N.norm(d,suball).eye{s+1}      ; N.norm(d,subj).eye{s+1}];
                    N.norm(d,suball).eyebntn{s+1}  = [N.norm(d,suball).eyebntn{s+1}  ; N.norm(d,subj).eyebntn{s+1}];
                    N.norm(d,suball).press{s+1}    = [N.norm(d,suball).press{s+1}    ; N.norm(d,subj).press{s+1}];
                    N.norm(d,suball).eyeveloc{s+1} = [N.norm(d,suball).eyeveloc{s+1} ; N.norm(d,subj).eyeveloc{s+1}];
                    N.norm(d,suball).prsveloc{s+1} = [N.norm(d,suball).prsveloc{s+1} ; N.norm(d,subj).prsveloc{s+1}];
                    N.norm(d,suball).subj{s+1}     = [N.norm(d,suball).subj{s+1}     ; subj*ones(size(N.norm(d,suball).subj{s+1} ,1) ,1)];
                end
                N.norm(d,suball).eyePattern(s+1, :)     = nanmedian(N.norm(d,suball).eye{s+1});
                N.norm(d,suball).pressPattern(s+1,:)    = nanmedian(N.norm(d,suball).press{s+1});
                N.norm(d,suball).eyevelocPattern(s+1,:) = nanmedian(N.norm(d,suball).eyeveloc{s+1});
                N.norm(d,suball).prsvelocPattern(s+1,:) = nanmedian(N.norm(d,suball).prsveloc{s+1});
            end
        end
        
        
        %% =====================================================================================================================
    case 'single_sub'
        %%
        %% =====================================================================================================================
        sub = subjnum;
        ANA = getrow(Dall , Dall.SN == sub & ismember(Dall.Day , day));
        for s = 0:6
            id = find(ANA.seqNumb == s & ~ANA.isError & ANA.isgood);
            N.norm.eye{s+1}   = [];
            N.norm.eyebntn{s+1} = [];
            N.norm.press{s+1} = [];
            N.norm.eyeveloc{s+1} = [];
            N.norm.prsveloc{s+1} = [];
            
            for i = 1:length (id)
                if ANA.AllPressIdx(id(i),14) < length(ANA.xEyePosDigit{id(i)})
                    idd = floor(linspace(ANA.AllPressIdx(id(i),1) , ANA.AllPressIdx(id(i),14) , 1001));
                    iddv = floor(linspace(ANA.AllPressIdx(id(i),1)-window , ANA.AllPressIdx(id(i),14)-window , 1001));
                    clear temp1 temp2 temp3 temp4
                    for j = 1:length(idd)-1
                        temp1(j) = nanmean(ANA.xEyePosDigit{id(i)}(idd(j) : idd(j+1)));
                        temp2(j) = nanmean(ANA.PressTimeSeries{id(i)}(idd(j) : idd(j+1)));
                        temp3(j) = nanmean(ANA.xEyeAngVelocity{id(i)}(iddv(j)-window : iddv(j+1)-window));
                        temp4(j) = nanmean(ANA.pressVelocity{id(i)}(iddv(j)-window : iddv(j+1)-window));
                    end
                    N.norm.eye{s+1}   = [N.norm.eye{s+1} ;inpaint_nans(temp1)];
                    N.norm.eyebntn{s+1} = [N.norm.eyebntn{s+1} ; [ANA.BN(id(i)) ANA.TN(id(i))]];
                    N.norm.press{s+1} = [N.norm.press{s+1} ; inpaint_nans(temp2)];
                    N.norm.eyeveloc{s+1} = [N.norm.eyeveloc{s+1} ; inpaint_nans(temp3)* 1000];
                    N.norm.prsveloc{s+1} = [N.norm(d,sub).prsveloc{s+1} inpaint_nans(temp4)* 1000];
                    
                end
            end
            N.norm.eyePattern(s+1, :)  = nanmean(N.norm.eye{s+1})/2.3;
            N.norm.pressPattern(s+1,:) = nanmean(N.norm.press{s+1});
            N.norm.eyevelocPattern(s+1,:) = nanmean(N.norm.eyeveloc{s+1});
            N.norm.prsvelocPattern(s+1,:) = nanmean(N.norm.prsveloc{s+1});
        end
end
