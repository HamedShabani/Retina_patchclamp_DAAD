function [fh,hst,PSTH]=rPSTH(spikes,trgs,xlimm,binw,xtk,xlbl,ytk,ylbl,fname,filtw)
% [fh,hst,PSTH]=rPSTH(spikes,trgs,xlim,binw,xtk,xlbl,ytk,ylbl,fname,filtw)
% Generalized Rastergram with associated PSTH
% [fh,hst,PSTH]=rPSTH(spikes,trgs,[0 4.1],.01,[0,203.1,406.2],{'ON','OFF',''},ytk,ylbl,fname,5)
% Use binw=.002; filtw=2 for anything less than 500 msec total.
% fh=figure handle
% hst=raw histograms for each stimulus
% PSTH=the resultant peristimulus time histogram (smoothed if filtw given)
% spikes = list of spike times for a single cell - in seconds
% trgs = nx2 list of trigger times (actually start and end times for the
%   rPSTH calculation relative to each true trigger time) - in seconds
%   Will be plotted from top of rastergram to bottom in the order given.
% xlim = limits for the x axes of rPSTH - in seconds
% binw = binning width for the PSTH - in seconds
% xtk = x tick mark list for the PSTH - in seconds
% xlbl = labels for the x tick marks
% ytk = list of stimulus block indices
%   NOTE:  Should be provided in order of bottom to top of rastergram
% ylbl = names for the stimulus blocks
% fname = file name (for the title of the rPSTH)
% filtw = (optional) width (SD) of the PSTH smoothing filter in multiples
%   of the bin width.
%Called by:  ArchAna1.m ArchAna4.m FlashAna.m
%  A stripped down version is used as a subfunction in EpiThreshAna10.m
%Calls:  suptitle.m
%Created by Daniel L. Rathbun 2012????
%REVISIONS
% 201209??  by Daniel L. Rathbun  Fixed the rastergram xlim to accept offsets.
% 20140807  by Daniel L. Rathbun  Removed the +1 from the ylim of the
% rastergram to get rid of the extra row.
xlimm2=[0,.1]
fh=figure('Units', 'Normalized', 'Position', [.1 .1 .8 .8])
 hold on; %initialize the figure.
% plot the rastergram for each stimulus
SPKS=[];
for i1=1:size(trgs,1)
    spks=spikes(spikes>=trgs(i1,1)& spikes<trgs(i1,2))-trgs(i1,1);
    subplot(4,2,[1,3,5]), line([spks';spks'],[ones(1,length(spks))*size(trgs,1)-i1; ones(1,length(spks))*size(trgs,1)-i1+1],'color','k');
    SPKS=[SPKS;spks [ones(length(spks),1)*size(trgs,1)-i1]];
end

set(gca,'ytick',ytk,'yticklabel',ylbl,'xtick',xtk,'xticklabel',xlbl,'xlim',xlimm,'ylim',[0 size(trgs,1)],'FontSize',13); %clean up axis
xlabel('time (ms)')
%set(gca,'ytick',ytk,'yticklabel',ylbl,'xtick',xtk,'xticklabel',xlbl,'ylim',[0 size(trgs,1)],'xlim',xlimm); %clean up axis
title('Long plot')
ylabel('mv/uamp');

%% hamed comes!

%     lbls=[];spk=[];lbl=[];SPKS=[];
%     for i1=1:size(trgs,1)
%         %spk=[spk; data(data>=stm(trg) & (data<stm(trg)+diffs(trg)))-stm(trg)];
%         spk=[spk; spikes(spikes>=trgs(i1,1)& spikes<trgs(i1,2))-trgs(i1,1)];
%         lbl=[lbl; ones(length(spikes(spikes>=trgs(i1,1)& spikes<trgs(i1,2))),1)*i1];
%         %lbl=ones(length(data(data>=stm(trg) & data<stm(trg)+diffs(trg))),1)*trg;
%         
%         %lbls=[lbls;lbl];
%     end
%     SPKS=[spk lbl];
%     
%     figure
%         line([SPKS(:,1)' ;SPKS(:,1)'],[SPKS(:,2)' ;SPKS(:,2)'+1],'Color','k')
% 





delete spks
for i1=1:size(trgs,1)
    spks=spikes(spikes>=trgs(i1,1)& spikes<trgs(i1,2))-trgs(i1,1);
    subplot(4,2,[2,4,6]), line([spks';spks'],[ones(1,length(spks))*size(trgs,1)-i1; ones(1,length(spks))*size(trgs,1)-i1+1],'color','k');
end
set(gca,'ytick',ytk,'yticklabel',ylbl,'xtick',xtk,'xticklabel',xlbl,'ylim',[0 size(trgs,1)],'FontSize',13); %clean up axis
title('Zoomed plot(100ms)')
xlabel('time (ms)')
xlim(xlimm2)

















% create a histogram for each stimulus
for i1=1:size(trgs,1)
    hst(:,i1)=histc(spikes,trgs(i1,1):binw:trgs(i1,2));
end
PSTH=mean(hst,2)/binw; %convert the histogram to Hz
subplot(4,2,7), plot(xlimm(1)+.5*binw:binw:xlimm(2),PSTH(1:end-1),'color','k'); ax1=gca; %plot the PSTH.
set(gca,'xlim',xlimm,'FontSize',13); ylabel('Hz'); %clean up axis
%xlabel('time (ms)');
%if requested, add a filtered PSTH

if nargin > 11
    SmFilt=pdf('norm',-filtw*3:filtw*3,0,filtw); %create the smoothing filter
    PSTH=conv([PSTH(end-filtw*3:end-1);PSTH(1:end-1);PSTH(1:filtw*3)],SmFilt,'valid'); %Assuming it's cyclic, we can pad the ends. [fix this for non-cyclic data]
    ax2=axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none','Ycolor','r'); %initialize 2nd axis
    line(linspace(xlimm(1),xlimm(2),length(PSTH)),PSTH,'color','r','parent',ax2); %plot smoothed PSTH
    ylim=get(ax2,'ylim'); set(ax2,'ylim',[0 ylim(2)]); %clean 2nd Y axis
end
%
PSTH=mean(hst,2)/binw; %convert the histogram to Hz
subplot(4,2,8), plot(linspace(xlimm(1),xlimm(2),length(PSTH)),PSTH,'color','k'); ax1=gca; %plot the PSTH.
set(gca,'xlim',xlimm2,'FontSize',13); ylabel('Hz'); %clean up axis
%xlabel('time (ms)');
%if requested, add a filtered PSTH
if nargin > 11
    SmFilt=pdf('norm',-filtw*3:filtw*3,0,filtw); %create the smoothing filter
   % PSTH=conv([PSTH(end-filtw*3:end-1);PSTH(1:end-1);PSTH(1:filtw*3)],SmFilt,'valid'); %Assuming it's cyclic, we can pad the ends. [fix this for non-cyclic data]
    ax2=axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none','Ycolor','r'); %initialize 2nd axis
    line(linspace(xlimm(1),xlimm(2),length(PSTH)),PSTH,'color','r','parent',ax2); %plot smoothed PSTH
    ylim=get(ax2,'ylim'); set(ax2,'ylim',[0 ylim(2)]); %clean 2nd Y axis
end
set(gca,'xlim',xlimm2,'FontSize',13); %clean joint x axis

suptitle(fname); %set(gcf,'position',[-1279 206 1280 636]); %figure title
%fhh = ancestor(fh, 'figure');
%fhh=get(gcf)
%fh=[];
