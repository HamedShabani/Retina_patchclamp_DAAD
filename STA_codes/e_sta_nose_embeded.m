function sta=e_sta_nose_embeded(fpath,fname,jnk,nt,lock_out,trg_period)
%% this codes plots the sta of noise embeded sinusoids.
%nt is the filter lengeh
%jnk is 3
%;lock_out is .010 (10 ms)


RawDat=load([fpath fname]); %load the workspace
vnames=sort(fieldnames(RawDat)); %get unit names
for i1=1:length(vnames)-jnk
    Data(i1).name=vnames(i1+jnk);
    eval(['Data(i1).spks=RawDat.',Data(i1).name{:},';']);
    
end

stm=RawDat.A2a;
%stm(102324:104993)=[]

stm_missed=[ find(diff(stm)>.045 & diff(stm)<.1)'];
stm1=[];
%stm(1:280)=[];
for i = 1: length(stm)% add missed triggers
    
    if sum(i==stm_missed)
        stm1=[stm1;stm(i);stm(i)+.04];
    else
        stm1=[stm1;stm(i)];
    end
end
stm=stm1;stm1=[];
%stm=stm(87501:120177);after drug
%stm=stm(1:40000);% before drug
if ~isempty(trg_period)
    stm=stm(trg_period);
end

stm_missed1=[ find(diff(stm1)>.04 & diff(stm1)<.1)'];
trl_idx=[find(diff(stm)>3)'];
stm_idx=[ find(diff(stm)<3)'];
%stm(15000:17676)=[];

trls_times1=[[stm(1);stm(find(diff(stm)>.3)+1)] ,[stm(diff(stm)>.3);stm(end)]];
betweem_trls_times1=[stm(diff(stm)>.3) ,stm(find(diff(stm)>.3)+1)];
fs=25;
%stm(102324:104993)=[]
%%
load('D:\Hame2\Scripts\enoise_embeded\sin_noise_embeded.mat')
%load('D:\Hame2\Data\Sydney_data\191022\egwn_1_Freq =25_Mean =400_contrast =140.mat')
for i2=1:length(Data)-1
    i2
    ch_name= char(Data(i2).name());
    ch_name=ch_name(3:end);
    
    data=Data(i2).spks;
    data(data<stm(1) | data>[stm(end)+1/fs])=[];% rejecte data before and after stimulations
    for i =1:length(betweem_trls_times1) % reject data between trials
        data(data>=betweem_trls_times1(i,1)& data<=betweem_trls_times1(i,2))=[];
        %trl_data(i,:)= data(data>trls_times1(i,1) & data< trls_times1(i,2));
    end
    
       
        figure
    subplot(131)
dif2=diff(data)
histogram(dif2,[0:.005:.2])
xlabel('Time sec')
ylabel('number of Spikes')
title('ISI Histogram all data   ') 
    
    %lock_out=0;
    for st=1: length(stm)-1% reject direct spikes
        
        data(data>=stm(st) & data<stm(st)+lock_out)=[];
        %data(data>stm(st)+lock_out & data<=stm(st+1))=[];

    end

subplot(132)
dif2=diff(data)
histogram(dif2,[0:.005:.2])
xlabel('Time sec')
ylabel('number of Spikes')
title('ISI Histogram for Indirect reponse   ') 


% Burst Correction Derived from Nima
    diff_spiketimes = diff(data); % diff_spiketimes
     Burst_End = find(diff_spiketimes >.1); % p.ISI_gap is what we set based on the scatter plot of pre and post ISI
     new_spike_times = [];
    for iter = 1:length(Burst_End) - 1
     new_spike_times_holder = repmat(data(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
     new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
    end
    data=new_spike_times;
subplot(133)
dif2=diff(data)
dif2(dif2==0)=[];
[a]=histogram(dif2,[0:.005:.2])
xlabel('Time sec')
ylabel('number of Spikes')
title('ISI Histogram for Indirect reponse with Burst Correction   ')    
    
    lbls=[];spk=[];
    diffs=diff(stm);
    for trg=1:length(stm)-1% plot single stimulus
        spk=[spk; data(data>=stm(trg) & (data<stm(trg)+diffs(trg)))-stm(trg)];
        lbl=ones(length(data(data>=stm(trg) & data<stm(trg)+diffs(trg))),1)*trg;
        lbls=[lbls;lbl];
    end
    SPKS=[spk lbls];
    fig=figure('Units', 'Normalized', 'Position', [.1 .1 .8 .8])
    if length(SPKS)>1
        subplot(411)
        line([SPKS(:,1)' ;SPKS(:,1)'],[SPKS(:,2)' ;SPKS(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
        xlim([0,0.04])
        title('Spike raster for each stimulus 40 msec')
        xlabel('Time sec')
        ylabel('stimulation number ')
        
        
        stm_times1=stm(diff(stm)<3);
        [spk_counts]=histcounts(data,stm);
        
        trl_idx=[0;find(diff(stm)>3)];
        spk_counts_trls=[];trgrs=[];
        lbls_trl=[];spk_trl=[];
        sin_period=[475:500 ,599:625 ,736:750, 862:875, 987:1000, 1112:1124, 1243:1250, 1368:1375 ...
            1491:1500, 1616:1625, 1745:1750, 1870:1875, 1994:2000, 2119:2125, 2246:2250, 2372:2375];
        
%         sin_periods=[475,500 ,599,625 ,736,750, 862,875, 987,1000, 1112,1124, 1243,1250, 1368,1375 ...
%             1491,1500, 1616,1625, 1745,1750, 1870,1875, 1994,2000, 2119,2125, 2246,2250, 2372,2375];
%         for ii=1:2:length(sin_periods)/2
%         data(data>=(stm(sin_periods(ii))) & (data<stm(sin_periods(ii+1))))=[];
%         end
        
        for trl=1:length(trl_idx)  % plot full trial
            [spk_counts_trls]=[spk_counts_trls; histcounts(data,stm(trl_idx(trl)+1:trl_idx(trl)+2500))];
            trgrs=[trgrs stm(trl_idx(trl)+1:trl_idx(trl)+2500-1)];
            
            spk_trl=[spk_trl; data(data>=(stm(trl_idx(trl)+1)) & (data<stm(trl_idx(trl)+2500)))-stm(trl_idx(trl)+1)];
            lbl=ones(length(data(data>=(stm(trl_idx(trl)+1)) & (data<stm(trl_idx(trl)+2500)))),1)*trl;
            lbls_trl=[lbls_trl;lbl];
            
        end
        if length(spk_trl)>1
            
            SPKS_trl=[spk_trl lbls_trl];
            subplot(412)
            
            line([SPKS_trl(:,1)' ;SPKS_trl(:,1)'],[SPKS_trl(:,2)' ;SPKS_trl(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
            title('Spike raster for each trial 100 sec')
            xlim([0,100.04])
            
            ylabel('trial number ')
        end
        spk_counts_trls(:,sin_period-1)=0;
        spk_counts=reshape(spk_counts_trls',length(trgrs(:)),1);
        %spk_counts=histcounts(data,stm(rex2(:,1)))
        
        Stimulus_all_trials=repmat(noiseembeded(1:end-1)',trl,1);
        if sum(spk_counts)>10
            [sta,stc,rawmu,rawcov] = simpleSTC(Stimulus_all_trials,spk_counts,nt);
            %[sta,stc,rawmu,rawcov] = simpleSTC(Stimulus_all_trials,spk_counts',nt);
            
            [u,S,V] = svd(stc); % Compute eigenvectors of STC matrix
            
            % Compute iSTAC estimator
            nfilts = 1;
            stim_fs=25;

            % number of temporal elements of filter
%            tvec = (-nt/2+1:nt/2)'; % vector of time indices (in units of stim frames)
            tvec = (-nt/2+1:nt/2)'*1/stim_fs-.5/stim_fs; % vector of time indices (in units of stim frames)

            %istacfilts = compiSTAC(sta, stc, rawmu, rawcov, nfilts);
            spk_nmbr=sum(spk_counts);
            subplot(414)
            plot( tvec, sta,'LineWidth',1)
            patch([tvec(1) tvec(nt/2) tvec(nt/2) tvec(1)],[min(sta) min(sta) max(sta)+10  max(sta)+10],'red')
            xlabel('time before spike (bins)'); ylabel('filter coeff');
            title(['filter estimates' ,' Channel  ',ch_name(4:end),'  spike number is ',num2str(spk_nmbr) ]);
            alpha .1
            %ylim([-1300,-300])

            subplot(413)
            npsth=mean(spk_counts_trls)/max(mean(spk_counts_trls));
            
            plot(npsth*5*std(noiseembeded)+mean(noiseembeded),'LineWidth',1)
            ylabel({'Average of PSTH of all trials(Scaled)';'Voltage of stimulus mv'})
            hold on
            title('Stimulus VS Response noiseembed')
            %plot((rexp54Frozen/max(abs(rexp54Frozen))+.4))
            plot(noiseembeded)
            alpha .1
            legend('Response','Stimulus')

            %plot( tvec, sta./norm(sta), tvec, u(:,end), tvec, istacfilts);
            %legend('STA', 'STC', 'iSTAC', 'location', 'northwest');
            
            %     subplot(122)
            % for st=1: length(stm)-1% label ttl spikes
            % %     trgs_spk=[trgs_spk; data(data>stm(st) & data<=stm(st+1))-stm(st)];
            % %     trgs_label=[trgs_label;st*ones(length(data(data>stm(st) & data<=stm(st+1))),1)];
            %     trgs_spk_cnt(st,:) = histcounts(data(data>stm(st) & data<=stm(st+1))-stm(st),...
            %     [0:1/(40*fs):1/fs]);
            % end
            % plot(linspace(0,1/fs,length(mean(trgs_spk_cnt,1))),mean(trgs_spk_cnt,1)*fs)
            %     xlabel('time sec'); ylabel('Rate');
            
            folderName='Figures';
            fn = fullfile(fpath,folderName);
            if ~ exist(fn, 'dir' )
                mkdir(fpath, 'Figures')
            end
            print(fig,'-dpng',[fpath,'\Figures\','Electrical_sta_noiseembed_sinusoid_removed  ',fname(1:end-4),'  ',ch_name,'  lt=' ,num2str(lock_out),'.png'],'-r300')
                        saveas(fig,[fpath,'\Figures\','Electrical_sta_noiseembed_sinusoid_removed  ',fname(1:end-4),'  ',ch_name,'  lt=' ,num2str(lock_out),'.fig'])

        end
    end
end
