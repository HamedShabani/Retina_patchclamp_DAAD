function sta=e_sta_sps(fpath,fname,fnamesps,jnk,nt,lock_out,trg_period,stim_file,just_direct)
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
    %ch_name=ch_name(3:end);
    data=Data(i2).spks;

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
    if lock_out~=0;
        for st=1: length(stm)% reject direct spikes
            if just_direct~=1
                data(data>=stm(st) & data<(stm(st)+lock_out))=[];
            else
                data(data<(stm(st)+1/fs) & data>(stm(st)+lock_out))=[];
            end


        end
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
    %data=new_spike_times;
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
        spk=[spk; data(data>=stm(idx(trg)) & (data<stm(idx(trg))+5*diffs))-stm(idx(trg))];
        lbl=ones(length(data(data>=stm(idx(trg)) & (data<stm(idx(trg))+5*diffs))),1)*trg;
        %lbl=ones(length(data(data>=stm(trg) & data<stm(trg)+diffs(trg))),1)*trg;
        
        lbls=[lbls;lbl];
    end
    SPKS=[spk lbls];
    spk_counts_trls=[];trgrs=[];
    %if length(SPKS)>0
        
       bwid=.02
       fig=figure('Units', 'Normalized', 'Position', [.1 .1 .8 .8])

%% Plot raster to triggers 40msec         
        
        
subplot(3,4,[1,2,5,6,9,10])
        line([SPKS(:,1)' ;SPKS(:,1)'],[SPKS(:,2)' ;SPKS(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
        xlim([0,0.2])
        title('Spike raster for each stimulus 40 msec')
        xlabel('Time sec')
        ylabel('stimulation Amplitude  mv/uamp ')
        yticks([0: trg/5:trg])
        yticklabels({num2str((v(1))),num2str((v(ceil(trg/5)))),...
            num2str((v(2*ceil(trg/5)))),num2str((v(3*ceil(trg/5)))),...
            num2str((v(4*ceil(trg/5)))),num2str((v(trg)))})

        
        stm_times1=stm(diff(stm)<1);
        [spk_counts]=histcounts(data,stm);
        
        trl_idx=[0;find(diff(stm)>1)];
        
        lbls_trl=[];spk_trl=[];lbl=[];SPKS_trl=[];
        for trl=1:length(trl_idx)%-1 %is applied only to dataset 2018_1_11WR
            [spk_counts_trls]=[spk_counts_trls; histcounts(data,stm(trl_idx(trl)+1:trl_idx(trl)+2500))];
            trgrs=[trgrs stm(trl_idx(trl)+1:trl_idx(trl)+2500-1)];
            
            spk_trl=[spk_trl; data(data>=stm(trl_idx(trl)+1) & (data<stm(trl_idx(trl)+2500)))-stm(trl_idx(trl)+1)];
            lbl=ones(length(data(data>=stm(trl_idx(trl)+1) & (data<stm(trl_idx(trl)+2500)))),1)*trl;
            lbls_trl=[lbls_trl;lbl];
            
        end
        %if length(spk_trl)>1
            
            SPKS_trl=[spk_trl lbls_trl];
            %%
subplot(3,4,[7,8,11,12])
line([SPKS(:,1)' ;SPKS(:,1)'],[SPKS(:,2)' ;SPKS(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
        xlim([0,0.005])
        title('Spike raster for each stimulus 40 msec (zoomed on first 5 ms)')
        xlabel('Time sec')
        ylabel('stimulation Amplitude  mv/uamp ')
        yticks([0: trg/5:trg])
        yticklabels({num2str(round(v(1))),num2str(round(v(ceil(trg/5)))),...
            num2str(round(v(2*ceil(trg/5)))),num2str(round(v(3*ceil(trg/5)))),...
            num2str(round(v(4*ceil(trg/5)))),num2str(round(v(trg)))})
%% PLot raters of 100 sec trials
% subplot(3,4,[7,8,11,12])
%             
%             line([SPKS_trl(:,1)' ;SPKS_trl(:,1)'],[SPKS_trl(:,2)' ;SPKS_trl(:,2)'+1],'Marker', '.', 'MarkerSize', 2, ...
%                 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color','k')
%             title('Spike raster for each trial 100 sec')
%             xlim([0,100.04])
%             ylabel('trial number ')
%             xlabel('Time sec')
%             set(gca,'FontSize',12)
      % end
%% Compute STA
        spk_counts=reshape(spk_counts_trls',length(trgrs(:)),1);
        %spk_counts=histcounts(data,stm(rex2(:,1)))
        
        stim_fs=25;
        Stimulus_all_trials=repmat(rexp54Frozen(1:end-1)',trl,1);
        if sum(spk_counts)>10
            [sta,stc,~,rawcov] = simpleSTC(Stimulus_all_trials,spk_counts,nt);

            lengthX = length(sta);
            x = 1:lengthX;
            samplingRateIncrease = 10;
            newXSamplePoints = linspace(1, lengthX, lengthX * samplingRateIncrease);
            usta = spline(x, sta, newXSamplePoints);
            
            
            nusta=usta-min(usta);
            thr=mean(rexp54Frozen)%-min(usta);%+2*std(sta(26:end))

            [valmin,locmin]=min(usta((1:260)));% find negative peak
            ir=1
            while (usta(locmin+ir))<(thr+(valmin-thr)/2)
                ir=ir+1
            end
            
            il=1
            while (usta(locmin-il))<(thr+(valmin-thr)/2)
                il=il+1
            end   
            
            if usta(locmin-il)> usta(locmin+ir)
                yshiftn=ir;
            else 
                yshiftn=-il
            end
            
            [valmax,locmax]=max(usta((1:260)));% find negative peak

            irp=1
            while (usta(locmax+irp))>(thr+(valmax-thr)/2)
                irp=irp+1
            end
            
            ilp=1
            while (usta(locmax-ilp))>(thr+(valmax-thr)/2)
                ilp=ilp+1
            end 
            
            if usta(locmin-ilp)> usta(locmin+irp)
                yshiftp=-ilp;
            else 
                yshiftp=irp
            end
            
            
            latency_p= (250-locmax)*4;% 250th sample is in the middle (eaual to time zero)
            latency_n= (250-locmin)*4;% 250th sample is in the middle (eaual to time zero)
            halfmaxwidth=(il+ir)*4;

%             [halfw_r_val,loc_halfw_r]=min(distancer) ;% find the right side of half width
%             [halfw_l_val,loc_halfw_l]=min(distancel) ;% find the left side of half width
            

            
         

            
            
            tvec = (-nt/2+1:nt/2)'*1/stim_fs-.5/stim_fs; % vector of time indices (in units of stim frames)



            %[sta,stc,rawmu,rawcov] = simpleSTC(Stimulus_all_trials,spk_counts',nt);
            
            %[u,~,~] = svd(stc); % Compute eigenvectors of STC matrix
            
            % Compute iSTAC estimator
            
            % number of temporal elements of filter
 %           tvec = (-nt/2+1:nt/2)'*1/stim_fs; % vector of time indices (in units of stim frames)
 %           [pks,locs,w,p]= findpeaks (sta,tvec, 'Annotate' , 'extents', 'MinPeakHeight',thr )
            [pksn,locsn,wn,pn]= findpeaks (abs(sta),tvec, 'Annotate' , 'extents', 'MinPeakHeight',max(-sta) )
            T = table (locsn, pksn, wn, 'VariableNames' , { 'Location' , 'PeakHeight' , 'Width' });
% 
%             [peakmaxp, peaklocsp]=max(pks)
%             [peakmaxn, peaklocsn]=max(pksn)
%             %istacfilts = compiSTAC(sta, stc, rawmu, rawcov, nfilts);
            spk_nmbr=sum(spk_counts);
subplot(3,4,[3,4])
            plot( tvec, sta,'LineWidth',2)
%             findpeaks (sta,tvec, 'Annotate' , 'extents' ,'MinPeakHeight',thr )
%             text(locs(peaklocsp),pks(peaklocsp),num2str(w(peaklocsp)))
            text(.004*(locmin+ir-220),thr+(min(usta)-thr)/2,['FWHM(-)= ',num2str((il+ir)*4),' ms'])
            %text(.004*(locmin+ir-220),nusta(locmin+il)+min(usta)-15,['FWHM= ',num2str(round(wn*1000)),' ms(findpeaks)'])
            text(.004*(locmax+irp-220),thr+(max(usta)-thr)/2,['FWHM(+)= ',num2str((ilp+irp)*4),' ms'])

            text(.004*(locmin+ir-220),thr+(min(usta)-thr)/2+(min(usta)-thr)/5,['Latency(-)= ',num2str(latency_n),' ms'])
            text(.004*(locmax+irp-220),thr+(max(usta)-thr)/2+(max(usta)-thr)/5,['Latency(+)= ',num2str(latency_p),' ms'])

            
            
            line([.004*(locmin-il-250),.004*(locmin+ir-250)],[usta(locmin+yshiftn),usta(locmin+yshiftn)],'color','r','LineWidth',2)
            line([.004*(locmax-ilp-250),.004*(locmax+irp-250)],[usta(locmax+yshiftp),usta(locmax+yshiftp)],'color','r','LineWidth',2)

            hold on

            %patch([tvec(1) tvec(nt/2) tvec(nt/2) tvec(1)],[min(sta) min(sta) max(sta)  max(sta)],'red')
            patch([-1 0 0 -1],[min(sta) min(sta) max(sta)  max(sta)],'red')
            xlabel('time before spike (sec)'); ylabel('filter coeff');
            title(['filter estimates ' ,' Channel  ',ch_name,'  spike number is ',num2str(spk_nmbr) ]);
            line([tvec(1),tvec(end)],[mean(rexp54Frozen) ,[mean(rexp54Frozen)]],'Color','k','LineStyle','--')
            %set(gca,'FontSize',12)
            %ylim([min(sta)-25,max(sta)+25])
            alpha .05
        end
        
%% SPS


%       fig
%     h = get(fh,'Children');
%     newh = copyobj(h,2)
% for i=1:2
%     
% posnewh = get(newh(i),'Position');
% possub  = get(ax(i),'Position');
% set(newh(i),'Position',[posnewh(1) possub(2) posnewh(3) possub(4)])
%     end
% end
        %% Plot Electrical stimulus and response
% subplot(5,3,[4:6])
%             
%             npsth=mean(spk_counts_trls)/max(mean(spk_counts_trls));
%             t_stim=linspace(0,1/stim_fs*length(npsth),length(npsth))
%             
%             plot(t_stim,rexp54Frozen(1:end-1),':')
%             alpha .1
%             hold on
%             
%             plot(t_stim,npsth*5*std(rexp54Frozen)+mean(rexp54Frozen),'LineWidth',1)
%             ylabel({'Average of PSTH of all trials(Scaled)';'Voltage of stimulus mv'})
%             hold off
%             title('Electrical Stimulus VS Response')
%             legend('Stimulus','Response')

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
    
    if just_direct
    	typeofplot='just direct'
    else
        typeofplot='idirect'
    end
        
    print(fig,'-dpng',[fpath,'\Figures\','Electrical_sta ',fname(1:end-4),'  ',ch_name,' lt= ', num2str(lock_out),' ',typeofplot,' long plot .png'],'-r300')
%    saveas(fig,[fpath,'\Figures\','Electrical_sta_',fname(1:end-4),'  ',ch_name,'  lt=' ,num2str(lock_out),' No Burst Correction .fig'])
%just direct response
        
    end
   %close all
   end
