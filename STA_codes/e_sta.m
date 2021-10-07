function Data=e_sta(fpath,fname,jnk,nt,lock_out,trg_period,stim_file,meafig,trls)
%% sta
set (0, 'DefaultTextInterpreter' , 'none' )



RawDat=load([fpath fname]); %load the workspace
vnames=sort(fieldnames(RawDat)); %get unit names
for i1=1:length(vnames)-jnk
    Data(i1).name=vnames(i1+jnk);
    eval(['Data(i1).spks=RawDat.',Data(i1).name{:},';']);
    
end

stm=RawDat.A2a;
stm_missed=[ find(diff(stm)>.045 & diff(stm)<.1)'];
stm1=[];
for i = 1: length(stm)% add missed triggers
    
    if sum(i==stm_missed)
        stm1=[stm1;stm(i);stm(i)+.04];
    else
        stm1=[stm1;stm(i)];
    end
end
stm=stm1;stm1=[];
%stm=stm(120178:160177); %after drug
%stm=stm(40001:87500);% before drug
%stm=stm(1:40000);% before drug
if ~isempty(trg_period)
    
stm=stm(trg_period);
end
stm_missed1=[ find(diff(stm1)>.04 & diff(stm1)<.1)'];
trl_idx=[find(diff(stm)>1)'];
stm_idx=[ find(diff(stm)<1)'];

trls_times1=[[stm(1);stm(find(diff(stm)>.3)+1)] ,[stm(diff(stm)>.3);stm(end)]];
betweem_trls_times1=[stm(diff(stm)>.3) ,stm(find(diff(stm)>.3)+1)];
fs=25;


%%
%load('rexp54Frozen.mat')
load(stim_file)

%load('D:\Hame2\Data\Sydney_data\Stimili\Current\egwn_1_Freq =25_Mean =2_contrast =0.7.mat')
%load('D:\Hame2\Data\Sydney_data\Stimili\Current\egwn_1_Freq =25_Mean =4_contrast =1.4.mat')

%load('D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =200_contrast =70')
%load('D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =400_contrast =140')
%load('D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =1000_contrast =350')
rexp54Frozen=r';
stim=repmat(rexp54Frozen,[1],length(trl_idx)+1)';
[v,idx]=sort(stim,'descend');
for i2=1:length(Data)-1
    i2
    ch_name= char(Data(i2).name());
    ch_name=ch_name(3:end);
    data=Data(i2).spks;
    Data(i2).name=[fname(1:end-4),'_',ch_name];

%     figure
%     subplot(131)
% dif2=diff(data)
% histogram(dif2,[0:.005:.2])
% xlabel('Time sec')
% ylabel('number of Spikes')
% title('ISI Histogram all data   ')

    
    data(data<stm(1) | data>[stm(end)+1/fs])=[];% rejecte data before and after stimulations
    for i =1:length(betweem_trls_times1) % reject data between trials
        data(data>=betweem_trls_times1(i,1)& data<=betweem_trls_times1(i,2))=[];
        %trl_data(i,:)= data(data>trls_times1(i,1) & data< trls_times1(i,2));
    end
    %lock_out=0;
    for st=1: length(stm)% reject direct spikes
        
        data(data>=stm(st) & data<stm(st)+lock_out)=[];
    end
    
% subplot(132)
% dif2=diff(data)
% histogram(dif2,[0:.005:.2])
% xlabel('Time sec')
% ylabel('number of Spikes')
% title('ISI Histogram for Indirect reponse   ')

% Burst Correction Derived from Nima
    diff_spiketimes = diff(data); % diff_spiketimes
     Burst_End = find(diff_spiketimes >.1); % p.ISI_gap is what we set based on the scatter plot of pre and post ISI
     new_spike_times = [];
    for iter = 1:length(Burst_End) - 1
     new_spike_times_holder = repmat(data(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
     new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
    end
 %   data=new_spike_times;
% subplot(133)
% dif2=diff(data)
% dif2(dif2==0)=[];
% [a]=histogram(dif2,[0:.005:.2])
% xlabel('Time sec')
% ylabel('number of Spikes')
% title('ISI Histogram for Indirect reponse with Burst Correction   ')
    
    lbls=[];spk=[];lbl=[];SPKS=[];
    diffs=diff(stm);
    diffs=.03985;
    for trg=1:length(stm)
        %spk=[spk; data(data>=stm(trg) & (data<stm(trg)+diffs(trg)))-stm(trg)];
        spk=[spk; data(data>=stm(idx(trg)) & (data<stm(idx(trg))+diffs))-stm(idx(trg))];
        lbl=ones(length(data(data>=stm(trg) & data<stm(trg)+diffs)),1)*trg;
        %lbl=ones(length(data(data>=stm(trg) & data<stm(trg)+diffs(trg))),1)*trg;
        
        lbls=[lbls;lbl];
    end
    SPKS=[spk lbls];
        
%% plot psth ans raster of visusal stimuli
%        trls=2
%        meafig='HT_2019_07_16_lw';
%trls=2;
       %meafig=load(meafig)%'HT_2019_12_22wr';
       stim_dur=[4,32,12];
       %trls=2;
       %meafig='HT_2018_06_13_wl';
       bwid=.02
       
       
       [PSTH_all, SPK_all]=single_cell_visualisze(fpath,fname,jnk,stim_dur,bwid,trls,meafig,i2);
       [all,t_all]=plot_visual_stimuli(1/bwid);
       fig=figure('Units', 'Normalized', 'Position', [.1 .1 .8 .8])
subplot(5,3,[1:3])% psth
       plot(t_all,all*max(PSTH_all),':')
       %ylim([-2,2])
       hold on
       plot(t_all(1:length(PSTH_all)),PSTH_all);
       xlim([-2,length(all)*bwid])
       %ylim([-max(PSTH_all)-10,max(PSTH_all)+10])
       ylabel('FR (spike/sec)')
       title(['Visual stimulus versus response','   ',fname(1:end-4),'  Channel  ',ch_name(4:end)])
       %legend('Stimulus','PSTH')
       patch([0 2.03 2.03 0],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'y','EdgeColor','none')
       patch([2.09 4 4 2.09],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
       patch([7.14 10 10 7.14],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'y','EdgeColor','none')
       patch([10 13 13 10],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
       patch([38 42 42 38],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'r','EdgeColor','none')

       patch([42.1 45.1 45.1 42.1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'g','EdgeColor','none')
       patch([45.1 48.82-1 48.82-1 45.1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
       patch([48.82-1 51.93-1 51.93-1 48.82-1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'b','EdgeColor','none')
       patch([51.95-1 54.98-1 54.98-1 51.95-1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')



%        patch([0 2 2 0],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'y','EdgeColor','none')
%        patch([2 4 4 2],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
%        patch([7 10 10 7],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'y','EdgeColor','none')
%        patch([10 13 13 10],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
% 
% %       patch([3 6 6 3],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'y','EdgeColor','none')
% %       patch([6 9 9 6],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
%        patch([9 11 11 9],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
%        patch([33 38 38 33],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'r','EdgeColor','none')
% 
%        patch([38 41 41 38],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'g','EdgeColor','none')
%        patch([41 44 44 41],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
%        patch([44 47 47 44],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'b','EdgeColor','none')
%        patch([47 50 50 47],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')


alpha .05
xlim([0,50])
       hold off
      % ylim([0,max(PSTH_all)+max(PSTH_all)/4])

       v_stimuli(1).names='Flash'
       v_stimuli(2).names='Chirp'
       v_stimuli(3).names='BG'
       c=0;
       for plt=13:15
           subplot(5,3,plt)% raster
           c=c+1;
           stname=v_stimuli(c).names;
           
           SPK_trg=SPK_all.(stname);
           line([(SPK_trg(:,1)),(SPK_trg(:,1))]',[SPK_trg(:,2),1+SPK_trg(:,2)]','Marker', '.', 'MarkerSize', 2, ...
               'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
           line([0,0],[0,max(SPK_trg(:,2))],'Color','r','MarkerSize', 4)
           title([stname,' Spike raster ',fname(1:end-4)])
           %view([0,-90])
           xlim([-1,stim_dur(c)])
       end
       clear c
%% Plot raster to triggers 40msec         
        
      if length(SPKS)>5
      
subplot(5,3,[11:12])
        line([SPKS(:,1)' ;SPKS(:,1)'],[SPKS(:,2)' ;SPKS(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
        xlim([0,0.04])
        title('Spike raster for each 1 ms pulse sorted based on magnitude (40 ms)')
        xlabel('Time sec')
        ylabel('Stimulation Voltage mv ')
        yticks([0: trg/2:trg])
        yticklabels({num2str(round(v(1))),num2str(round(v(ceil(trg/2)))),num2str(round(v(trg)))})

        
        stm_times1=stm(diff(stm)<1);
        [spk_counts]=histcounts(data,stm);
        
        trl_idx=[0;find(diff(stm)>1)];
        spk_counts_trls=[];trgrs=[];
        lbls_trl=[];spk_trl=[];lbl=[];SPKS_trl=[];
        for trl=1:length(trl_idx)-1 %is applied only to dataset 2018_1_11WR
            [spk_counts_trls]=[spk_counts_trls; histcounts(data,stm(trl_idx(trl)+1:trl_idx(trl)+2500))];
            trgrs=[trgrs stm(trl_idx(trl)+1:trl_idx(trl)+2500-1)];
            
            spk_trl=[spk_trl; data(data>=stm(trl_idx(trl)+1) & (data<stm(trl_idx(trl)+2500)))-stm(trl_idx(trl)+1)];
            lbl=ones(length(data(data>=stm(trl_idx(trl)+1) & (data<stm(trl_idx(trl)+2500)))),1)*trl;
            lbls_trl=[lbls_trl;lbl];
            
        end
        %if length(spk_trl)>1
            
            SPKS_trl=[spk_trl lbls_trl];
%% PLot raters of 100 sec trials
subplot(5,3,[7:9])
            
            line([SPKS_trl(:,1)' ;SPKS_trl(:,1)'],[SPKS_trl(:,2)' ;SPKS_trl(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
            title('Spike raster for each trial of noise stimulus (100 sec)')
            xlim([0,100.04])
            ylabel('Trial number ')
            xlabel('Time s')
            %set(gca,'FontSize',8)
       % end
%% Compute STA
        spk_counts=reshape(spk_counts_trls',length(trgrs(:)),1);
        %spk_counts=histcounts(data,stm(rex2(:,1)))
        
        stim_fs=25;
        Stimulus_all_trials=repmat(rexp54Frozen(1:end-1)',trl,1);
        if sum(spk_counts)>1
            [sta,stc,~,rawcov] = simpleSTC_hamed(Stimulus_all_trials,spk_counts,nt);
            %[sta,stc,rawmu,rawcov] = simpleSTC(Stimulus_all_trials,spk_counts',nt);
            
            %[u,~,~] = svd(stc); % Compute eigenvectors of STC matrix
            
            % Compute iSTAC estimator
            
            % number of temporal elements of filter
            tvec = (-nt/2+1:nt/2)'*1/stim_fs-.5/stim_fs; % vector of time indices (in units of stim frames)
            
            %istacfilts = compiSTAC(sta, stc, rawmu, rawcov, nfilts);
            spk_nmbr=sum(spk_counts);
subplot(5,3,10)
            plot( tvec, sta,'LineWidth',2)
            hold on
            patch([-2 0 0 -2],[-1200 -1200 0  0],'black')
            %patch([0 2 2 0],[-1200 -1200 0  0],'black')

            xlabel('time before spike (sec)'); ylabel('E-STA (mv)');
            title([' Channel  ',ch_name(4:end),'  ',num2str(spk_nmbr) , '  Spikes _NO_BC' ]);
            line([tvec(1),tvec(end)],[mean(rexp54Frozen) ,[mean(rexp54Frozen)]],'Color','k','LineStyle','--')
            line([tvec(nt/2)+.5/stim_fs tvec(nt/2)+.5/stim_fs], [-1200 0],'Color','k','LineStyle','--')

            %set(gca,'FontSize',12)
            ylim([min(sta)+min(sta)/100,max(sta)-max(sta)/100])
            xlim([-1,.5])
            alpha .05
            %grid on
            Data(i2).time_sta=tvec;
            Data(i2).e_sta=sta;
        end
%% Plot Electrical stimulus and response            
subplot(5,3,[4:6])
            
            npsth=mean(spk_counts_trls)/max(mean(spk_counts_trls));
            t_stim=linspace(0,1/stim_fs*length(npsth),length(npsth))
            
            plot(t_stim,rexp54Frozen(1:end-1),':')
            alpha .1
            hold on
            
            plot(t_stim,npsth*5*std(rexp54Frozen)+mean(rexp54Frozen),'LineWidth',1)
            ylabel({'Voltage of stimulus mv'})
            yyaxis right
            ylabel('Average of PSTHs(Scaled)','Color','k')
            set( gca, 'YTick', [] )

            hold off
            title('Electrical stimulus VS response')
            legend('Stimulus','Response')
end
            %plot((rexp54Frozen/max(abs(rexp54Frozen))+.4))
            

            %set(gca,'FontSize',12)

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
    print(fig,'-dpng',[fpath,'\Figures\','Electrical_sta ',fname(1:end-4),'  ',ch_name,' lt= ', num2str(lock_out),' NO_BurtsCorrection 5 plots for paper  .png'],'-r300')
                           %saveas(fig,[fpath,'\Figures\','Electrical_sta_',fname(1:end-4),'  ',ch_name,'  lt=' ,num2str(lock_out),' Bernstein.fig'])

        
    
   %close all

close all
   end
 save([fpath,'\Matlab_results\','e_STA_NO_BurtsCorrection ',fname(1:end-4),'.mat'],'Data')
