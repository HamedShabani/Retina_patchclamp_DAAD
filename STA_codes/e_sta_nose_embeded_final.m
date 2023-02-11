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
    
       
        %figure
    %subplot(131)
dif2=diff(data)
%histogram(dif2,[0:.005:.2])
%xlabel('Time sec')
%ylabel('number of Spikes')
%title('ISI Histogram all data   ') 
    
    %lock_out=0;
    for st=1: length(stm)-1% reject direct spikes
        
        data(data>=stm(st) & data<stm(st)+lock_out)=[];
        %data(data>stm(st)+lock_out & data<=stm(st+1))=[];

    end

%subplot(132)
dif2=diff(data)
%histogram(dif2,[0:.005:.2])
%xlabel('Time sec')
%ylabel('number of Spikes')
%title('ISI Histogram for Indirect reponse   ') 


% Burst Correction Derived from Nima
    diff_spiketimes = diff(data); % diff_spiketimes
     Burst_End = find(diff_spiketimes >.1); % p.ISI_gap is what we set based on the scatter plot of pre and post ISI
     new_spike_times = [];
    for iter = 1:length(Burst_End) - 1
     new_spike_times_holder = repmat(data(Burst_End(iter) + 1), Burst_End(iter + 1) - Burst_End(iter), 1);
     new_spike_times = vertcat(new_spike_times, new_spike_times_holder);
    end
    %data=new_spike_times;
%subplot(133)
dif2=diff(data)
dif2(dif2==0)=[];
%[a]=histogram(dif2,[0:.005:.2])
%xlabel('Time sec')
%ylabel('number of Spikes')
%title('ISI Histogram for Indirect reponse without Burst Correction   ')    
    
    lbls=[];spk=[];
    diffs=diff(stm);
    for trg=1:length(stm)-1% plot single stimulus
        spk=[spk; data(data>=stm(trg) & (data<stm(trg)+diffs(trg)))-stm(trg)];
        lbl=ones(length(data(data>=stm(trg) & data<stm(trg)+diffs(trg))),1)*trg;
        lbls=[lbls;lbl];
    end
    SPKS=[spk lbls];
    %fig=figure('Units', 'Normalized', 'Position', [.1 .1 .8 .8])
    if length(SPKS)>1
        %subplot(411)
        %line([SPKS(:,1)' ;SPKS(:,1)'],[SPKS(:,2)' ;SPKS(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
        %    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
        %xlim([0,0.04])
       % title('Spike raster for each stimulus 40 msec')
        %xlabel('Time sec')
        %ylabel('stimulation number ')
        
        
        stm_times1=stm(diff(stm)<3);
        [spk_counts]=histcounts(data,stm);
        
        trl_idx=[0;find(diff(stm)>3)];
        spk_counts_trls=[];trgrs=[];
        lbls_trl=[];spk_trl=[];
        %sin_period=[475:500 ,599:625 ,736:750, 862:875, 987:1000, 1112:1124, 1243:1250, 1368:1375 ...
        %   1491:1500, 1616:1625, 1745:1750, 1870:1875, 1994:2000, 2119:2125, 2246:2250, 2372:2375];
        
         sin_period=[475:501 ,599:626 ,736:751, 862:876, 987:1001, 1112:1125, 1243:1251, 1368:1376 ...
            1491:1501, 1616:1626, 1745:1751, 1870:1876, 1994:2001, 2119:2126, 2246:2251, 2372:2376];
               
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
            %subplot(412)
            
            %line([SPKS_trl(:,1)' ;SPKS_trl(:,1)'],[SPKS_trl(:,2)' ;SPKS_trl(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
            %    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
            %title('Spike raster for each trial 100 sec')
            %xlim([0,100.04])
            
            %ylabel('trial number ')
        end
        spk_counts_trls(:,sin_period-1)=0;
        spk_counts=reshape(spk_counts_trls',length(trgrs(:)),1);
        %spk_counts=histcounts(data,stm(rex2(:,1)))
        
        Stimulus_all_trials=repmat(noiseembeded(1:end-1)',trl,1);
        if sum(spk_counts)>10
            [sta,stc,rawmu,rawcov] = simpleSTC_hamed(Stimulus_all_trials,spk_counts,nt);
            %[sta,stc,rawmu,rawcov] = simpleSTC(Stimulus_all_trials,spk_counts',nt);
            
            [u,S,V] = svd(stc); % Compute eigenvectors of STC matrix
            
            % Compute iSTAC estimator
            nfilts = 1;
            stim_fs=25;

            % number of temporal elements of filter
%            tvec = (-nt/2+1:nt/2)'; % vector of time indices (in units of stim frames)
            %tvec = (-nt/2+1:nt/2)'*1/stim_fs-.5/stim_fs; % vector of time indices (in units of stim frames)
            tvec = (-nt/2+1:nt/2)'*1/stim_fs-.5/stim_fs
            %istacfilts = compiSTAC(sta, stc, rawmu, rawcov, nfilts);
            spk_nmbr=sum(spk_counts);
    %        subplot(414)

            %ylim([-1300,-300])
            %%
            fig=figure('units','normalized','outerposition',[0.05 0.05 .9 .9]);
            cutperiod=[1:2500]
            SPKS_trl_Sample=SPKS_trl;
            SPKS_trl_Sample(:,1)=SPKS_trl_Sample(:,1)/(40/1000);
            t_hist=linspace(0,20,2000)
            
            t_stim = stm(trl_idx(1)+1:trl_idx(1)+2500);
            t_stim=t_stim(end-length(cutperiod)+1:end);
            t_stim=t_stim-t_stim(1);


            x_lim=[18,38]
            pois=[0.100    0.9500    0.7750    0.8150]
            for iplot= 1:4
            axes('Position',[pois(1) pois(2)-iplot*.2 pois(3) .2]);

            cutperiodt=cutperiod*40/1000;
            
            %plot(npsth(cutperiod)*5+SPKS_trl(end,end)+1,'LineWidth',1)
            maskidx=SPKS_trl(:,1)>=cutperiodt(1) & SPKS_trl(:,1)<cutperiodt(end);            
            plot(t_stim(2:end),histcounts([SPKS_trl(maskidx,1)-cutperiodt(1)],t_stim),'LineWidth',1)
            hold on
            box off
            
            %plot((rexp54Frozen/max(abs(rexp54Frozen))+.4))
            plot(t_stim,[(noiseembeded(cutperiod)- mean(noiseembeded(cutperiod)))/800]*3+SPKS_trl(end,end)+5,'LineWidth',1,'Color','r')
            alpha .1 
            
            if iplot==1
                title('Stimulus VS Response noise embedded')
            end
            
            hold on
            line([SPKS_trl(maskidx,1)'-cutperiodt(1) ;SPKS_trl(maskidx,1)'-cutperiodt(1)],[SPKS_trl(maskidx,2)' ;SPKS_trl(maskidx,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
            xlim(x_lim)
      

            set(gca,'ytick',[],'FontSize',14)
            if iplot<4
                box off;
                set(gca,'xtick',[])
                set(gca,'ytick',[])
            else 
                xticks([x_lim(1):5:x_lim(2)])
                xticklabels([x_lim(1):5:x_lim(2)]-x_lim(1)) 
            end
            %cutperiod=cutperiod+500;           %xlim([0,500])

            x_lim=x_lim+20;
            end

            axes('Position',[.45 0.03 .152 .1]);

            plot( tvec, sta,'LineWidth',2)
            %patch([tvec(1) tvec(nt/2) tvec(nt/2) tvec(1)],[min(sta) min(sta) max(sta)+10  max(sta)+10],'red')
            
            

            
            hold on
            patch([-2 0 0 -2],[-1200 -1200 0  0],'black')
            %patch([0 2 2 0],[-1200 -1200 0  0],'black')

            line([tvec(1),tvec(end)],[mean(noiseembeded) ,[mean(noiseembeded)]],'Color','k','LineStyle','--')
            line([tvec(nt/2)+.5/stim_fs tvec(nt/2)+.5/stim_fs], [-1200 0],'Color','k','LineStyle','--')

            %set(gca,'FontSize',12)
            ylim([min(sta)+min(sta)/100,max(sta)-max(sta)/100])
            xlim([-1,.5])
            alpha .05
   
            
            
            
            %xlabel('time before spike (bins)'); 
            ylabel('filter coeff');
            title(['filter estimates' ,' Channel  ',ch_name(4:end),'  spike number is ',num2str(spk_nmbr) ]);
            %alpha .1
            
            
            %%%% plot psth ans raster of visusal stimuli
%        trls=2
%        meafig='HT_2019_07_16_lw';
%trls=2;
       %meafig=load(meafig)%'HT_2019_12_22wr';
       stim_dur=[4,12];
       %trls=2;
       %meafig='HT_2018_06_13_wl';
       bwid=.02
       
       
       [PSTH_all, SPK_all]=single_cell_visualisze_nombar_nochirp_nocolor(fpath,fname,jnk,stim_dur,bwid,i2);
       [all,t_all]=plot_visual_stimuli(1/bwid);
%        fig=figure('Units', 'Normalized', 'Position', [.1 .1 .8 .8])
% subplot(5,3,[1:3])% psth
%        plot(t_all,all*max(PSTH_all),':')
%        %ylim([-2,2])
%        hold on
%        
%        plot(t_all(1:length(PSTH_all)),PSTH_all);
%        xlim([-2,length(all)*bwid])
%        ylim([-max(PSTH_all)-10,max(PSTH_all)+10])
%        ylabel('FR')
%        title(['Visual stimulus vs reponse','   ',fname(1:end-4),'  Channel  ',ch_name(4:end)])
%        %legend('Stimulus','PSTH')
%        patch([0 2.03 2.03 0],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'y','EdgeColor','none')
%        patch([2.09 4 4 2.09],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
%        patch([7.14 10.18 10.18 7.14],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'y','EdgeColor','none')
%        patch([10.2 13.26 13.26 10.2],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
% 
%        patch([42.1 45.1 45.1 42.1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'g','EdgeColor','none')
%        patch([45.1 48.82-1 48.82-1 45.1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
%        patch([48.82-1 51.93-1 51.93-1 48.82-1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'b','EdgeColor','none')
%        patch([51.95-1 54.98-1 54.98-1 51.95-1],[-max(PSTH_all) -max(PSTH_all) max(PSTH_all) max(PSTH_all)],'k','EdgeColor','none')
% 
%        alpha .05
%        hold off
      % ylim([0,max(PSTH_all)+max(PSTH_all)/4])

       v_stimuli(1).names='Flash'
       v_stimuli(2).names='BG'
%       v_stimuli(3).names='BG'
       c=0;
       %for plt=[13]
          % subplot(5,3,plt)% raster
           axes('Position',[.15 0.03 .152 .1]);

           c=c+1;
           stname=v_stimuli(c).names;
           
           SPK_trg=SPK_all.(stname);
           if length(SPK_trg)>0
               line([(SPK_trg(:,1)),(SPK_trg(:,1))]',[SPK_trg(:,2),1+SPK_trg(:,2)]','Marker', '.', 'MarkerSize', 2, ...
                   'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
               line([0,0],[0,max(SPK_trg(:,2))],'Color','r','MarkerSize', 4)
               title([stname,' Spike raster ',fname(1:end-4)])
               view([0,-90])
               xlim([-1,stim_dur(c)])
         %  end
       end
       clear c
            
            
            
            folderName='Figures';
            fn = fullfile(fpath,folderName);
            if ~ exist(fn, 'dir' )
                mkdir(fpath, 'Figures')
            end
            print(fig,'-dpng',[fpath,'\Figures\','Electrical_sta_noiseembed_sinusoid_detialed_removed_visual  ',fname(1:end-4),'  ',ch_name,'  lt=' ,num2str(lock_out),'.png'],'-r300')
           % saveas(fig,[fpath,'\Figures\','Electrical_sta_noiseembed_sinusoid_detialed_removed  ',fname(1:end-4),'  ',ch_name,'  lt=' ,num2str(lock_out),'.fig'])
%%

%%
%             fig=figure('units','normalized','outerposition',[0.05 0.05 .9 .9]);
%             pois=[0.100    0.9900    0.7750    0.8150]
% 
%             psth=histcounts([SPKS_trl(maskidx,1)-cutperiodt(1)],t_stim)';
%             psth=(psth-min(psth))/max(psth-min(psth));
%             
%             filter=flip(sta(1:25)-mean(sta));
%             %filter=(filter-min(filter))/max(filter);
%             
%             x1 = sameconv(psth,filter);             
% 
%             stim=(noiseembeded'-min(noiseembeded))/max(noiseembeded'-min(noiseembeded));
%             
%                       hold on
% 
%             
%             axes('Position',[pois(1) pois(2)-1*.2 pois(3) .2]);
%             plot(stim)
%             
%            
%             x2 = sameconv(stim,filter);
%             axes('Position',[pois(1) pois(2)-2*.2 pois(3) .2]);            
%             plot(x2)
%            
%             
%             axes('Position',[pois(1) pois(2)-3*.2 pois(3) .2]);
%             plot(x1);            hold on;
%             plot(psth*20)
%             
%             slen=length(stim);
%             r = 10*max(x2,0);   % Half-wave rectified nonlinearity
%             % Generate Poisson spike response ---------------
%             RefreshRate = 100;  % Stim refresh rate (Hz)
%             dtbin = .01;          % binsize for Poisson spike generation
%             rbig = repmat(r'/RefreshRate*dtbin,1./dtbin,1); % make Poisson spike train
%             sp = sum(rand(size(rbig))<rbig)';
%             % Hamed add
%             spks=rand(size(rbig))<rbig;
%             tVec = 0:dtbin:slen-dtbin;
%             axes('Position',[pois(1) pois(2)-4.5*.2 pois(3) .3]);
%             title('Raster Plot of generated data (stim*filter)')
%             plotRaster(spks,tVec)
%             xlim([0,slen*dtbin])
% 
%             %edges=[1:10:max(diff(find(spks(:)==1)))/10];
%             %isi_hist=histcounts(diff(find(spks==1)),edges)/sum(histcounts(diff(find(spks==1)),edges));
%             %subplot(223)
%             %bar(isi_hist);title('ISI Histogram of poisson spikes')
% 
%   
%            print(fig,'-dpng',[fpath,'\Figures\','Spike generator model stimulus_filter  ',fname(1:end-4),'  ',ch_name,'  lt=' ,num2str(lock_out),'.png'],'-r300')

        end
    end
    close all
end
