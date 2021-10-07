jnk=3
fpath='D:\Hame2\Data\2018_06_13\2018_06_13_wl\'
fname='2018_06_13_wl.mat'
RawDat=load([fpath fname]); %load the workspace
vnames=sort(fieldnames(RawDat)); %get unit names
for i1=1:length(vnames)-jnk
    Data(i1).name=vnames(i1+jnk);
    eval(['Data(i1).spks=RawDat.',Data(i1).name{:},';']);
    
end

%%
%plot(diff(A2a))
stm=RawDat.A2a;
trl_nmber=4
stim_length=70;
Fs=50e3;

stm=stm(1:stim_length*trl_nmber)
stm(diff(stm)<4 )=[];

stim_dur=5;
between_trial=find(diff(stm)>2*stim_dur);
trial_times=[[1;between_trial+1], [between_trial(1);between_trial(2:end);between_trial(end)+stim_length]];
%%
% current_stim= tuningcurve(1:3:end,1);
% voltage_stim= VoltsUpDown(1:2:end,1);

load('D:\Hame2\current stimulus\colts_up_down.mat')
load('D:\Hame2\current stimulus\tuning_vurve.mat')
stimulus=voltage_stim;
stimulus=tuningcurve(1:3:end,1);

[v_val,v_idx]=sort(stimulus,'descend');
[v_u,i_u]=unique(v_val);
v_u=sort(v_u,'descend')

trials_number=[3,4];% number of current stimulus trials in recording files
trials_number=[1,2];% number of volatage stimulus trials in recording files
%%
stm=stm([trial_times(trials_number(1),1):trial_times(trials_number(end),2)]);
nstm=reshape(stm,[],2);
snstm=nstm(v_idx,:);
nstm(end+1,:)=nstm(end,:)+5;

trls=2;
%vol_times=stm(trial_times(voltage_trials,:));
for i2=1:length(Data)-1
        i2;
ch_name= char(Data(i2).name());
ch_name=ch_name(3:end);
   
data=Data(i2).spks;
data(data>=nstm(end,1) & data<nstm(1,2))=[];
data(data>=nstm(end,2) | data<nstm(1,1))=[];

PSTH=zeros(size(snstm,2)*10,length(v_u)*stim_dur*Fs);
t=length(stm)*stim_dur;
%ndata=data-stm(1);

    cnt=0;
    diffs=diff(stm);lbls=[];spk=[];lbls2=[];
    for stim=1:7
    stim_idx=stimulus==v_u(stim);
    sstm=nstm(stim_idx,:);
    cnt2=0;
    for trl=1:length(trials_number)
%         PSTH(trl,ndata snstm(trg,trl))=1;
% for stim=v_idx%length(unique(stimulus))
%             spk2=[spk2; data(data>=snstm(trg,trl) & data<snstm(trg,trl)+5)-snstm(1,trl)];


    for trg=1:10%length(unique(v_val))
        cnt=cnt+1;
        cnt2=cnt2+1;
        spk=[spk; data(data>=sstm(trg,trl) & data<sstm(trg,trl)+5)-sstm(trg,trl)];
        ndata=(data(data>=sstm(trg,trl) & data<sstm(trg,trl)+5)-sstm(trg,trl))*Fs;
        times=(stim*5-5)*Fs;
        trg_lbl(times+1)=1;
PSTH(cnt2,times+ceil(ndata))=1;
        lbl=ones(length(data(data>=sstm(trg,trl) & data<sstm(trg,trl)+5)),1)*cnt;
        lbl2=ones(length(data(data>=sstm(trg,trl) & data<sstm(trg,trl)+5)),1)*trl;

        lbls=[lbls;lbl];
        lbls2=[lbls2;lbl2];
        

    end
    end
    end
    SPKS=[spk lbls lbls2];
    a=mean(PSTH,1);
    %plot(mean(PSTH,1))
%B = smoothdata(mean(PSTH,1)','gaussian',Fs/10);
    h = fspecial('gaussian', Fs/10, .1);
    %fr=conv(mean(PSTH,1)',ones(1,Fs/10))/max(conv(mean(PSTH,1)',ones(1,Fs/10)));
    %plot(fr)
     y = normpdf([1:Fs/5],Fs/10,Fs/20);
         fr=conv(sum(PSTH,1)',y);
         fr=fr/max(fr);
         trg_lbl(length(fr))=0;
         
         fr(:,2)=trg_lbl;
         


  
    
%subplot(411)
ytic=[20:20:70*trls];
fig=figure;
subplot(211)
    line([SPKS(:,1)' ;SPKS(:,1)'],[SPKS(:,2)' ;SPKS(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
line([zeros(7,1)';5*ones(7,1)'],[ytic;ytic],'Marker', '.', 'MarkerSize', 2, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
xlim([0,5.04])
title('Spike raster for each stimulus 40 msec')
xlabel('Time sec')
ylabel('stimulation number ')
yticks(ytic)
yticklabels({'v = 100','v = 300','v = 500','v = 1000','v = 1500','v = 2000','v = 2500'})
yticklabels({'i = 1','i = 2.5','i = 5','i = 6.25','i =7.5','i = 8.7','i = 10'})

% psth
clear psth
% for trl= 1:size(snstm,2)
% psth(:,trl)=histcounts(data,nstm(:,trl))/mean(diff(nstm(:,trl)));
% spike_number=sum(sum(histcounts(data,nstm(:,trl))));
% end
% psth=psth(v_idx,:)
% clear spk_nmbr
% for trgs=1:7
%     tc=[1:10;11:20;21:30;31:40;41:50;51:60;61:70];
%     spk_nmbr(trgs)=mean(mean(psth(tc(trgs,:),:)));
% end
subplot(212)
plot(fr,'LineWidth',1)
xticks(linspace(Fs*2.5,length(fr),7))
xticklabels({'x = 100','x = 300','x = 500','x = 1000','x = 1500','x = 2000','x = 2500'})
yticklabels({'i = 1','i = 2.5','i = 5','i = 6.25','i =7.5','i = 8.7','i = 10'})

%plot(spk_nmbr)
    title(['psth estimates' ,' Channel  ',ch_name(4:end)]);
    legend('Response','Stimulus')
    
      folderName='Figures';
    fn = fullfile(fpath,folderName);
    if ~ exist(fn, 'dir' )
        mkdir(fpath, 'Figures')
    end
    print(fig,'-dpng',[fpath,'\Figures\','Electrical_current ',fname(1:end-4),'  ',ch_name,'2.png'],'-r300')


end





