function [fh,hst,PSTH]=plot_current_sps(datafile,recs,name,thr,SampRate)
% This function loads and plots the sopike raster of patch clamp data with current stimulation
% 
% datafile: location of the file
%recs: name of the recording session
%Name: specify the name of stimulation paradigm(Current or Voltage)
%thr: threshold for spike detection
% SampRate: Sampling rate

CurrSeq=[2:2:48]; %long 1 msec pulses
reps=5;
FullCurrSeq=sort(repmat(CurrSeq,[1,reps]));
load(datafile)

eval(['traces=',recs{1},';']); %rename raw data for generalized code.
VoltTrace=reshape(traces(:,1,:),[size(traces,1)*size(traces,3),1]);
TTLtrace=reshape(traces(:,2,:),[size(traces,1)*size(traces,3),1]);
t=linspace(0,length(VoltTrace)/SampRate,length(VoltTrace));



ttls=find(TTLtrace>2); %get list of all time points during the ttl (stimulus on).
trgs=ttls(find(diff(ttls)>(SampRate/50))+1); %onset timestamps for most TTLs
trgs=sort([trgs; ttls(1)]); %include first TTL timestamp.
trgs=trgs/SampRate; %convert to units = seconds.
d = designfilt('lowpassfir', ...
    'PassbandFrequency',0.1,'StopbandFrequency',0.2, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple')
trs=filtfilt(d,VoltTrace);
%trs=VoltTrace;

%APs=find(VoltTrace>-40); %get list of all time points during the Action Potentials.
[pks1, locs1] = findpeaks(trs,  'MinPeakHeight' , thr,'MinPeakDistance',SampRate*.005);
spks=locs1;
shrtpnts=spks(diff(spks)<(SampRate*.003));
spks(diff(spks)<(SampRate*.003))=[];

    figure; 
    plot([VoltTrace,TTLtrace]);hold on
    plot (spks,trs(spks),'*') 
    plot( shrtpnts,VoltTrace(shrtpnts),'o')
    title(['Raw Traces ',name]);
    xlabel('time sec'); ylabel('mV');
%spks=APs(find(diff(APs)>10)+1); %onset timestamps for the APs.
%spks=sort([spks; APs(1)]); %include first spike timestamp.
spks=spks/SampRate; %convert to units = seconds.

[B,SortInd]=sort(VoltTrace(1:length(trgs)));
trgs=trgs(SortInd);

reps=ceil(length(trgs)/7);
RepMod=mod(length(trgs),7);
if RepMod==0
    ytk=[reps:reps:length(trgs)];
else
ytk=[reps:reps:reps*RepMod, reps*RepMod+(reps-1):reps-1:length(trgs)];
end
xlim=[-.05 2]; binw=.001; xtk=[0:.05:.2, .4:.2:2]; xlbl=[-50:50:200, 400:200:2000]; 
filtw=5; ylbl=num2str(CurrSeq(1:7)');
trgs=[trgs+xlim(1),trgs+xlim(2)]; %add beginning and end times for response binning and plotting

[fh,hst,PSTH]=rPSTH(spks,trgs,xlim,binw,xtk,xlbl,ytk,ylbl,name,filtw);  %
xlabel('time (ms)'); %set(gcf,'position',[297 1379 441 346]);
saveas(gcf,name);
saveas(gcf,[name,'.jpg']);