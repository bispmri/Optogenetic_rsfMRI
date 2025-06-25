function StateMeanActivationmap_Singleslice_Group_rising(path1,iifiles, VOX,DIM, timep,statetc,trialname,thrmin,thrmax,state,scantype);
% input x is the time series from each single 4D voxel from the realigned data.
% the realigned fMRI data is first converted from a 4D data(time domain) to a 3D data(frequency domain)

cohere = zeros(DIM(1),DIM(2),DIM(3));

for abc=1:length(iifiles)    
    % Go to realigned functional data directory
    directory = strcat(path1, num2str(iifiles(abc)),'\NoFilter');
    %% Read Map Mask

    cd(path1);

    fid=fopen('lrrsra2dseq.img','r');
    temp0=fread(fid,'int8');
    scale_back=max(temp0);
    SBAmask=reshape(temp0,[DIM(1) DIM(2) 16])/scale_back;
    SBAmask=SBAmask(:,:,iifiles(abc));
    fclose(fid);
% end 

    %%
    %Read realigned data
    cd(directory);
    
  % read .mat fold
    filename=strcat('rrsra2dseq.mat');
    load(filename);
    RSdata=double(reshape(GvgEPIn,DIM(1), DIM(2), DIM(3),DIM(4)));
    RSdata0=double(RSdata);
    clear GvgEPIn
    
    %mask out   
    for kk=1:DIM(4)
             RSdata(:,:,:,kk)=SBAmask.*RSdata(:,:,:,kk);
    end
    T2 = ['M:\Atlas_Template\T2_template\mmean_T2_128.img'];

    %%
%     % convert 30s og and 15s og
    if ismember(scantype,[1,2])==1
        RSdata0(:,:,:,1:12825)=RSdata(:,:,:,35056:47880);
        RSdata0(:,:,:,12826:47880)=RSdata(:,:,:,1:35055);
        RSdata=RSdata0;
        clear RSdata0
    end
    % extract state time course
    RSstate=zeros(DIM(1),DIM(2),DIM(3),length(statetc));

    for n=1:length(statetc)
        RSstate(:,:,:,n)=RSdata(:,:,:,statetc(n));
    end
    cohere=mean(RSstate./100,4);
    cohere(isnan(cohere)==1) = 0;

    % Save state activation across time
    if ~exist(strcat(directory,'\State_MeanActivationMap\State_Activation_rising_v2\'), 'dir')
        mkdir(strcat(directory,'\State_MeanActivationMap\State_Activation_rising_v2\'));
    end
    cd(strcat(directory,'\State_MeanActivationMap\State_Activation_rising_v2\'));
    ofilename=strcat('State_Activation',num2str(state));
    nii = make_nii(RSstate(:,:,:,:),VOX,[0 0 0], 4); % maximum is 32727!
    save_nii(nii, ofilename);

    fid=fopen(strcat('Meanactivation_State',num2str(state), '_all.img'),'w');
    fwrite(fid,cohere(:,:,:), 'float32');
    spm_hwrite(strcat('Meanactivation_State',num2str(state),'_all'),[DIM(1) DIM(2) DIM(3)],VOX,1,16,0,[0 0 0],'spm compatible');
    fclose(fid); 



   % Save Overlay without phase on Realigned EPI
        % Save Overlay
        fid=fopen(T2,'r');
        img1=fread(fid,'int16');
        fclose(fid);
        img1=reshape(img1,[128 128 16]);
        structv=zeros([78 78 DIM(3)]);
        structv(:,:,:)=img1(26:103,26:103,iifiles(abc));
        ovCCmap=cohere(26:103,26:103,:);
        tmp2=zeros(78,78*10,3);
        for kkk=1:DIM(3)
            tmp=zeros(78,78,3);
            for jjj=1:78
                for iii=1:78                   
                    if(ovCCmap(iii,jjj,kkk)<-thrmin)
                    ovCCmap(iii,jjj,kkk)=double((ovCCmap(iii,jjj,kkk)+thrmin)./(thrmax-thrmin));
                    tmp(iii,jjj,1)=0;
                    tmp(iii,jjj,2)=abs(round(ovCCmap(iii,jjj,kkk)*255));
                    tmp(iii,jjj,3)=255-abs(round(ovCCmap(iii,jjj,kkk)*255));
                    %[0 1] -> [red yellow]
                elseif(ovCCmap(iii,jjj,kkk)>thrmin)
                    ovCCmap(iii,jjj,kkk)=double((ovCCmap(iii,jjj,kkk)-thrmin)./(thrmax-thrmin));
                    tmp(iii,jjj,1)=255.0;
                    tmp(iii,jjj,2)=ovCCmap(iii,jjj,kkk)*255;
                    tmp(iii,jjj,3)=0;
                else
                    tmp(iii,jjj,1)=round(structv(iii,jjj,kkk)/128*3);
%                     tmp(iii,jjj,1)=round(structv(iii,jjj,kkk)*32640/max(max(max(structv(iii,jjj,kkk))))/128);
                    tmp(iii,jjj,2)=round(structv(iii,jjj,kkk)/128*3);
                    tmp(iii,jjj,3)=round(structv(iii,jjj,kkk)/128*3);
                    end
                end
            end
            tmp(:,:,1)=flipud(tmp(:,:,1));
            tmp(:,:,2)=flipud(tmp(:,:,2));
            tmp(:,:,3)=flipud(tmp(:,:,3));
            tmp(:,:,1)=rot90(tmp(:,:,1),-1);
            tmp(:,:,2)=rot90(tmp(:,:,2),-1);
            tmp(:,:,3)=rot90(tmp(:,:,3),-1);
            tmp2(:,((kkk-1)*78+1):kkk*78,:)=tmp;
            
        end
        writefilename=strcat('MeanActivation_State',num2str(state),'_rising_',num2str(thrmin),'-',num2str(thrmax),'.tif');
        imwrite(uint8(tmp2), writefilename, 'tif');
end   
end