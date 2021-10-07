function [fh,hst,PSTH]=plot_voltage_sps(datafile,recs,name,thr,SampRate)
% This function reads the data with voltage stimulation and plots the spike
% rasteres
VoltSeq=[-100, -300, -500 -1000, -1500, -2000, -2500, -2500, -2000, -1500, -1000, -500, -300, -100,...
    -100, -300, -500 -1000, -1500, -2000, -2500, -2500, -2000, -1500, -1000, -500, -300, -100,...
    -100, -300, -500 -1000, -1500, -2000, -2500, -2500, -2000, -1500, -1000, -500, -300, -100,...
    -100, -300, -500 -1000, -1500, -2000, -2500, -2500, -2000, -1500, -1000, -500, -300, -100,...
    -100, -300, -500 -1000, -1500, -2000, -2500, -2500, -2000, -1500, -1000, -500, -300, -100];

 
 load(datafile)

 
eval(['traces=',recs{1},';']); %rename raw data for generalized code.
VoltTrace=reshape(traces(:,1,:),[size(traces,1)*size(traces,3),1]);
TTLtrace=reshape(traces(:,2,:),[size(traces,1)*size(traces,3),1]);
t=linspace(0,length(VoltTrace)/SampRate,length(VoltTrace));

%%TOGGLE to plot raw data
trs=VoltTrace;


ttls=find(TTLtrace>2); %get list of all time points during the ttl (stimulus on).
trgs=ttls(find(diff(ttls)>200)+1); %onset timestamps for most TTLs
trgs=sort([trgs; ttls(1)]); %include first TTL timestamp.
trgs=trgs/SampRate; %convert to units = seconds.

%[pks1ttls, locs1ttls] = findpeaks (TTLtrace,  'MinPeakHeight' , 2);
%trgs=locs1ttls/SampRate;

%APs=find(VoltTrace>-40); %get list of all time points during the Action Potentials.
[pks1, locs1] = findpeaks (VoltTrace,  'MinPeakHeight' , thr,'MinPeakDistance',SampRate*.002);
spks=locs1;


%figure; plot(t,[VoltTrace,TTLtrace]);hold on
shrtpnts=spks(diff(spks)<(SampRate*.002));
spks(diff(spks)<(SampRate*.002))=[];

    figure; 
    plot([VoltTrace,TTLtrace]);hold on
    plot (spks,trs(spks),'*') 
    plot( shrtpnts,VoltTrace(shrtpnts),'o')
    title(['Raw Traces ',name]);
    xlabel('time sec'); ylabel('mV');spks=spks/SampRate; %convert to units = seconds.

[B,SortInd]=sort(VoltSeq(1:length(trgs)));
trgs=trgs(SortInd);

reps=ceil(length(trgs)/7);
RepMod=mod(length(trgs),7);
if RepMod==0
    ytk=[reps:reps:length(trgs)];
else
ytk=[reps:reps:reps*RepMod, reps*RepMod+(reps-1):reps-1:length(trgs)];
end
xlim=[-.05 2]; binw=.001; xtk=[0:.05:.2, .4:.2:2]; xlbl=[-50:50:200, 400:200:2000]; 
filtw=5; ylbl=num2str(VoltSeq(1:7)');
trgs=[trgs+xlim(1),trgs+xlim(2)]; %add beginning and end times for response binning and plotting

[fh,hst,PSTH]=rPSTH(spks,trgs,xlim,binw,xtk,xlbl,ytk,ylbl,name,filtw);  %
xlabel('time (ms)'); %set(gcf,'position',[297 1379 441 346]);
saveas(gcf,name);
saveas(gcf,[name,'.jpg']);
end