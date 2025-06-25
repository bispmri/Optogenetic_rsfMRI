clc; clear all; close all;
% newpath='D:\Linshan';
% userpath(newpath)
addpath(pathdef)

%% Stage 1: Train HMM model
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%        Train HMM-Gaussian for all scans          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; clear all; close all;
    cd('.....\Dataset B\HMM raw data_cut\');
    load ("data_all.mat"); % concatenated all fMRI data
    load ("T_all.mat"); % time length for each trial
    data_modality = 'fMRI' ; % one of: 'fMRI', 'M/EEG', 'M/EEG power' or 'LFP' 
    no_states = 18; % the number of states depends a lot on the question at hand
    PCA=70; 
    Hz = 1; % the frequency of the data
    stochastic_inference = 0; % set to 1 if a normal run is too computationally expensive (memory or time)
    N = length(T); % number of subjects (total trial number?)
    M = size(Whloesubjtimecourse,1);
    options = struct();

    % PCA can be used to reduce dimentionality
    e = explainedvar_PCA(Whloesubjtimecourse,T,options); % how much variance PCA explains on the data
    pca=PCA/100;
    pc1 = find(e>pca,1); % get no. of PCA components to explain pca% of variance

    % Setting the options
    options.K = no_states;
    options.detrend=1;
    options.standardise = 1; %each trial has the same SD and mean
    options.verbose = 1;
    options.Fs = Hz;
    options.plotGamma = 2; % 1 for states in selected time scale; 2 for averaged across trials
    options.pca = pc1;
    
    % Gaussian observation model fMRI
    options.order = 0;
    options.zeromean = 0; % if 1, the mean of the time series will not be used to drive the states; 0 for PCA-HMM model
    options.covtype = 'full'; 
    options.useParallel=0 ;

    [hmm, Gamma] = hmmmar(Whloesubjtimecourse,T,options);

    % save hmm and Gamma
    cd('......'); 
    writefilename=strcat('hmm_', num2str(no_states), '_PCA',num2str(PCA),'%.mat');
    save(writefilename, 'hmm');
    writefilename=strcat('Gamma_', num2str(no_states), '_PCA',num2str(PCA),'%.mat');
    save(writefilename, 'Gamma');
    %% estimation from different dataset 
    clc; clear all; close all;
    cd('.....\Dataset B\HMM raw data_cut\');
    load ("data_RS.mat");
    load ("T_RS.mat");
    load ("data_all_OG.mat");
    load ("T_all_OG.mat");
    statesnum =22;
    PCA=70;
    cd('......'); 
    filename=strcat('Gamma_', num2str(statesnum), '_PCA',num2str(PCA),'%.mat');
    load (filename);
    filename=strcat('hmm_', num2str(statesnum), '_PCA',num2str(PCA),'%.mat');
    load (filename);

    [hmm_RS,Gamma_RS,vpath_RS] = hmmdual(WhloesubjtimecourseRS,T_RS,hmm);
    [hmm_OG,Gamma_OG,vpath_OG] = hmmdual(Whloesubjtimecourse,T_all_OG,hmm);
    writefilename=strcat('hmm_', num2str(statesnum), '_PCA',num2str(PCA),'%_RS.mat');
    save(writefilename, 'hmm_RS');
    writefilename=strcat('Gamma_', num2str(statesnum), '_PCA',num2str(PCA),'%_RS.mat');
    save(writefilename, 'Gamma_RS');
    writefilename=strcat('hmm_', num2str(statesnum), '_PCA',num2str(PCA),'%_allOG.mat');
    save(writefilename, 'hmm_all_OG');
    writefilename=strcat('Gamma_', num2str(statesnum), '_PCA',num2str(PCA),'%_allOG.mat');
    save(writefilename, 'Gamma_all_OG');


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%                load previous model                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; clear all; close all;
    modelplace = '......';
    cd(modelplace);

    % basic info
    statesnum=22;
    PCA=70;


    filename=strcat('Gamma_', num2str(statesnum), '_PCA',num2str(PCA),'%.mat');
    load (filename);
    filename=strcat('hmm_', num2str(statesnum), '_PCA',num2str(PCA),'%.mat');
    load (filename);
    filename=strcat('Gamma_', num2str(statesnum), '_PCA',num2str(PCA),'%_allOG.mat');
    load (filename);
    filename=strcat('hmm_', num2str(statesnum), '_PCA',num2str(PCA),'%_allOG.mat');
    load (filename);
    filename=strcat('Gamma_', num2str(statesnum), '_PCA',num2str(PCA),'%_RS.mat');
    load (filename);
    filename=strcat('hmm_', num2str(statesnum), '_PCA',num2str(PCA),'%_RS.mat');
    load (filename);


    cd('M:\1pulse_30cycle RS Data_v2\HMM raw data\Dataset B\HMM raw data_cut\');
    load ("data_RS.mat");
    load ("T_RS.mat");
    load ("data_all_OG.mat");
    load ("T_all_OG.mat");

    %% generate state transition space
    % transition matrix
    taskstate=[13 7 18 5 17 12 6 4 10 8 1 9 3 15]; %1pulseRS
    TP2=getTransProbs(hmm_RS);
    TP=TP2;
    imagesc(TP(taskstate,taskstate));colorbar;
    colormap(jet);
    clim([0 0.5]);
    TP0 = TP(taskstate,taskstate);

    % plot transition space
    n=0;
    s=[];
    t_idx=[];
    weights=[];
    names={};
    taskstate=sort([1 3:10 12:13 15 17 18]); %similar state in 1 pulse
    TP2=TP(taskstate,taskstate);
    for i=1:length(taskstate)
        for ii=1:length(taskstate)
%             if TP2(i,ii)>=0.15
%             if TP2(i,ii)>=min(maxk(reshape(TP2,[length(taskstate)*length(taskstate),1]),fix((length(taskstate)-1)*(length(taskstate)-1)*0.1)))
            if TP2(i,ii)>=min(maxk(reshape(TP2,[length(taskstate)*length(taskstate),1]),fix(length(taskstate)*length(taskstate)*0.15))) % top% connectivity
                n=n+1;
                s(n)=i;
                t_idx(n)=ii;
                weights(n)=TP2(i,ii);
            end
        end
    end

    taskstate_rank=[11,13,8,4,7,2,10,12,9,6,1,14,5,3]; %RS+1pulse model
    for i=1:length(taskstate)
%         names(i)={' '};
        names(i)={strcat('state',num2str(taskstate_rank(i)))};
    end

    G = digraph(s,t_idx,weights,names);
    LWidths = 6*G.Edges.Weight/max(G.Edges.Weight);

    g=plot(G,'LineWidth',LWidths, 'ArrowSize',15,'NodeFontSize',15);  
    g.Marker = 'o';
    g.NodeColor = [0.6350 0.0780 0.1840];
    x=[-1 -1.5 1 0.3 0.3 0.2 -0.5 1.8 1 -0.5 0.9 -1.5 -0.5 0.3]; % og+rs model
    y=[1.7 1 -0.5 0.2 -0.7 1.3 -1.5 0.5 -1.5 -0.5 1.7 0 0.5 1];
    g.XData(:) = x;
    g.YData(:) = y;
    highlight(g,'Edges',[14,15,8,19,6,13,17,18],'EdgeColor','#FFA500'); % IUS %RS!!
    highlight(g,'Edges',[12,7,5,11,9,4],'EdgeColor','#32CD32'); %RL
    highlight(g,'Edges',[1,2,10,3,16],'EdgeColor','#7B68EE'); %BS
    highlight(g,'Edges',[1,8],'EdgeColor','#FF0000'); %red
%     highlight(g,'Edges',[14,15,19,7,13,17,18],'EdgeColor','#FFA500'); % IUS %OG!!
%     highlight(g,'Edges',[12,8,5,11,9,4],'EdgeColor','#32CD32'); %RL
%     highlight(g,'Edges',[1,3,16,10],'EdgeColor','#7B68EE'); %BS
%     highlight(g,'Edges',[2,6],'EdgeColor','#FF0000'); %red
    g.EdgeAlpha=0.8;
    axis off
     %% state decomposition and generate the voxel-wise activation maps of single slice
    scantype=3; % 1=30sog; 2=15sog; 3=rsdata;
    taskstate =1:14;
    DIM=[128 128 1 length(Gamma_RS)];
    VOX = [0.25 0.25 1.0];
    path1='......\All RSdata\All\';
    slice=1:16;
    timep=1:length(Gamma_RS);
    trialname=9901;
    thrmin=0.1;
    thrmax=1;
    % find the lower threshold
    P=0;
    for state=1:statesnum
        LocNosig=find(Gamma(:,state)<2/3);
        P(state)=prctile(Gamma(LocNosig,state),95);
        LocSig=find(Gamma(:,state)>P);
    end

    rtemp1=0;
    rtemp2=0;
    dtemp1=0;
    dtemp2=0;
    ctemp1=0;
    ctemp2=0;
    ctemp=0;
    rtemp=0;
    dtemp=0;
    n=0;
    temp={};

    for state=taskstate
    % for state=13
%         % find the rise, during and decrease time
        rtemp1=find(diff(Gamma_decomposed(:,state))>0)+1;
        rtemp2=find(P(state)<Gamma_decomposed(:,state) & Gamma_decomposed(:,state)<0.9);
        rtemp=intersect(rtemp1,rtemp2)-1; % rising time
    
        dtemp1=find(diff(Gamma_decomposed(:,state))<0);
        dtemp2=find(P(state)<Gamma_decomposed(:,state) & Gamma_decomposed(:,state)<0.9);
        dtemp=intersect(dtemp1,dtemp2); % decreasing time
    
        ctemp=find(0.9<=Gamma_decomposed(:,state)); % during time

        alltemp=find(2/3<=Gamma_decomposed(:,state)); % all time

        StateMeanActivationmap_Singleslice_Group(path1,slice, VOX,DIM, timep,ctemp,trialname,thrmin,thrmax,state,scantype); %during
        StateMeanActivationmap_Singleslice_Group_decreasing(path1,slice, VOX,DIM, timep,dtemp,trialname,thrmin,thrmax,state,scantype);
        StateMeanActivationmap_Singleslice_Group_rising(path1,slice, VOX,DIM, timep,rtemp,trialname,thrmin,thrmax,state,scantype);
        StateMeanActivationmap_Singleslice_Group_all(path1,slice, VOX,DIM, timep,alltemp,trialname,thrmin,thrmax,state,scantype);
        n=n+1;
        temp_c(n)={ctemp};
        temp_r(n)={rtemp};
        temp_d(n)={dtemp};
        temp_all(n)={alltemp};
    end


