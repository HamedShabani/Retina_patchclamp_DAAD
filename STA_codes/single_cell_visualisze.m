function [PSTH_all SPK_all]=single_cell_visualisze(fpath,fname,jnk,Stim_dur,bwid,trls,meafig,cell_numb)
% this function plots the soike raster of data. in addition It makes a
% structure as outpt in which you can find the PSTH and the OGB-PSTH of
% data.
% iputs:
%       fpath= location of spike timing matlab file
%       fname= name of file contains spike timing
%       jnk= nimber of junk files in the data set(visual and electrical Triggers)
%           that shoud not be included in plotting
%       trgs= here you enter the type of stimulus you want to see the
%       results. the trigger file should be saved as a structure in data
%       file. e.g(trgs=RawDat.trgss.Flash).
%       stim_dur= stimulus duration sec
%       bwid= binwidth in sec used for PSTH calculation
%       trls= Moving bar stimulus trial number
%       cell_numb= number of cell (a single number between one to the last channel)
%Output:
%       
% example:
% close all
% set (0, 'DefaultTextInterpreter' , 'none' )
% FolderName='D:\Hame2\Data\2018_12_12\merged\2018_12_12_right_nasal\';cd(FolderName)
% bwid=.02;%binwidth sec
% jnk=3;
% trls=2
% meafig='HT_2018_12_12_right_nasal';
% HT.manual_corr=[-.1,0,.2,.1,0,-.2,.3,.1];
% stim_dur=[4,32,12];
% cell_numb=1;
% single_cell_visualisze(FolderName,'2018_12_12_right_nasal.mat',jnk,stim_dur,bwid,trls,meafig,cell_numb)

%%
set (0, 'DefaultTextInterpreter' , 'none' )

    %%
        % Moving bar
    %cell_numb=i2;
    bar_PSTH=single_cell_movingbar(fpath,fname,jnk,'trgs=RawDat.trgss.MovingBar',bwid,trls,meafig,cell_numb);


baseline_time=1;% duration of baseline activity selected for subtracted from evoked activity in sec
Fs=50000; bs_t=[0:1/Fs:baseline_time];

bs_l=length(bs_t);%length base line
trgname(1).trigger='RawDat.trgss.Flash';
trgname(2).trigger='RawDat.trgss.Chirp';
trgname(3).trigger='RawDat.trgss.BG';

for trig=1:3
    stim_dur=Stim_dur(trig);
   t=[0:1/Fs:stim_dur]; 
l=length(t);% length evoked
trgs=trgname(trig).trigger;
k = strfind(trgs,'Flash');


if isempty(k)
    k = strfind(trgs,'Noise');

    
end
if isempty(k)
    k = strfind(trgs,'BG');

    
end
if isempty(k)
    k = strfind(trgs,'Chirp');

end
stm_name=trgs(k:end);

RawDat=load([fpath fname]); %load the workspace
vnames=sort(fieldnames(RawDat)); %get unit names
trgrs=eval(trgs); %assign the triggers
%trgs=[trgs+xlim(1),trgs+xlim(2)]; %add beginning and end times for response binning and plotting


trgs=[];
if strcmp(stm_name,'BG') % because both blue and green colors have a separate ttls I need to remove one of them
    true_trgrs=[1:2:length(trgrs)];% select only first color's ttl(Green)
    trgrs=trgrs(true_trgrs);
end

between_trials =find(diff(trgrs)>2*stim_dur);% detect start of different recording session in one experiment using the timing between trials
trgs=[trgrs(1:end-1),trgrs(2:end)];% modify triggers timing
trgs(end+1,:)=[trgrs(end),trgrs(end)+stim_dur];% add the last trg stim_dur=4 sec for flash
if strcmp(stm_name,'BG') 
trgs(:,2)=trgs(:,1)+12;
end

trgs(between_trials,2)= trgrs(between_trials,1)+stim_dur;% correct the wrong value of last trigger of each session

ndata=[];
trg_data=[];trg_label=[];trg_label2=[];c=0;ytk=[];
count=[];
for i1=1:length(vnames)-jnk
    Data(i1).name=vnames(i1+jnk);
    eval(['Data(i1).spks=RawDat.',Data(i1).name{:},';']);
end
%% raster plot
for i2=cell_numb%1:length(Data)-1
    

%%
    i2;
    ndata=Data(i2).spks;
    c2=0;
    psth=[];psth_ogb=[];nonzero_trls=[];
    psth_bs_ogb=[];psth_bs=[];
    for trg=1:length(trgrs)
        
        trg_data=[trg_data ;ndata(ndata>trgs(trg,1)-baseline_time & ndata<trgs(trg,2))-trgs(trg,1)];
        trg_label2=[trg_label2 ;c+ones(length(ndata(ndata>trgs(trg,1)-baseline_time & ndata<trgs(trg,2))),1)];
        temp=zeros(1,l);
        temp_bs=zeros(1,bs_l);
        
        spk_tmp=[ndata(ndata>trgs(trg,1) & ndata<trgs(trg,2))-trgs(trg,1)];
        spk_tmp(spk_tmp>stim_dur)=[];
        
        bsline_spk_tmp=[ndata(ndata>trgs(trg,1)-baseline_time & ndata<trgs(trg,1))-trgs(trg,1)+baseline_time];
        
        if ~isempty(spk_tmp) % remove trials with no data
        
            psth=[psth; (1/bwid)*histcounts(ndata(ndata>trgs(trg,1)...
                -baseline_time & ndata<trgs(trg,2)),[[trgs(trg,1)-baseline_time]: bwid:trgs(trg,1)+stim_dur])];
            

        end
        c=c+1;% nuber of trg of all channels
        c2=c2+1;% number of trg of one channel
    end

    
    
    
  %%%%%%%%%%%% normalization based on baden's paper

  count =[count;c];
   PSTH.(stm_name)=mean(psth,1);
  SPK_trg.(stm_name)=[trg_data trg_label2];% for all cells
  %subplot(1,2,trig)
  
%   line([(SPK_trg(:,1)),(SPK_trg(:,1))]',[SPK_trg(:,2),1+SPK_trg(:,2)]','Marker', '.', 'MarkerSize', 2, ...
%       'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
%   line([0,0],[0,max(trg_label2)],'Color','r','MarkerSize', 4)
%   title([stm_name,' Spike raster ',fname(1:end-4)])
%   view([0,-90])
end
 

end
PSTH_all=[PSTH.('Flash') PSTH.('Chirp') bar_PSTH PSTH.('BG')  ];
%PSTH_all=[PSTH.('Chirp') bar_PSTH PSTH.('BG')  ];

SPK_all=SPK_trg;
% figure
% hold on
% 
% [all,t_all]=plot_visual_stimuli(1/bwid)
% plot(t_all,all*max(PSTH_all),':')
% %ylim([-2,2])
% xlim([-2,length(all)*bwid])
% 
% plot(t_all(1:length(PSTH_all)),PSTH_all);
%ylbl(cellfun(@isempty,ylbl))=[];


% visualize@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% ylbl = cellstr(ylbl(~cellfun('isempty', ylbl')));
% 
% ylabel ('Number of Channel','FontSize',22)
% xlabel('Time S','FontSize',22)
% xlim([-1,stim_dur])
% set(gca,'ytick',ytk,'yticklabel',ylbl,'FontSize',10)
% view([0,-90])
% set(hfig, 'PaperSize',[21 85])
% 
% folderName='Figures';
% fn = fullfile(fpath,folderName);
% if ~ exist(fn, 'dir' )
% mkdir(fpath, 'Figures')
% end
% print(hfig,'-dpng',[fpath,'\Figures\',stm_name,'Baseline removed raster plot ',fname(1:end-4),'.png'],'-r600')
% 
% 
% 
% % save@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% folderName='Matlab_results';
% fn = fullfile(fpath,folderName);
% if ~ exist(fn, 'dir' )
% mkdir(fpath, 'Matlab_results')
% end
% %save([fpath,'\Matlab_results\',stm_name,' ',fname(1:end-4),'.mat'],'Data')
% save([fpath,'\Matlab_results\',stm_name,'.mat'],'Data')

