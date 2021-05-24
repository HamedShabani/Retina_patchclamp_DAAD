

%note:  make sure that the thresholds for TTL and spike detection are valid
%for each recording.
%Sampled at 50kHz. 20 samples during TTL.  Must have been a post time of
%1000 microsec after the 1000 microsec voltage stimulus pulse.
fs=50000
load('D:\Hame2\Data\Sydney_data\191022\191022 data files.mat');
date='2019_10_22'
for ii=2:2
    switch ii
        case 1
            recs={'OFFT_191022_NoiseEmbeddedTrial1_BD_400mV'};
            name='OFFT_191022_BD_400mV_new';
        case 2
            recs={'OFFT_191022_NoiseEmbeddedTrial1_BD_800mV'};
            name='OFFT_191022_BD_800mV_new';
        case 3
            recs={'OFFT_191022_NoiseEmbeddedTrial1_BD_1000mv'};
            name='OFFT_191022_BD_1000mV_new';     
        case 4
            recs={'OFFT_191022_NoiseEmbeddedTrial1_BD_2uA'};
            name='OFFT_191022_BD_2uA_new'; 
        case 5
            recs={'OFFT_191022_NoiseEmbeddedTrial1_BD_4uA'};
            name='OFFT_191022_BD_4uA_new';
            

    end
%%
trgs1=[]; spks1=[];
reps=length(recs);
for i1=1:reps
    eval(['traces=',recs{i1},';']); %rename raw data for generalized code.
    
%%TOGGLE to plot raw data
%     figure; plot(traces); %take a gander at the raw patch recording (col 1) and ttl signal (col 2).
%     title(['Raw Traces ',name,'__',num2str(i1)]);
%     xlabel('samples'); ylabel('mV');
%     saveas(gcf,['Raw Traces ',name,'__',num2str(i1)]);
%     saveas(gcf,['Raw Traces ',name,'__',num2str(i1),'.jpg']);
    
    ttls=find(traces(:,2)>.5); %get list of all time points during the ttl (stimulus on).
    trgs=ttls(find(diff(ttls)>200)+1); %onset timestamps for most TTLs
    skipd=ttls(find(diff(ttls)>2000)+1)-2000; %timestamps for 8 missing TTLs (where the voltage pulse had and sctual amplitude of 0V, so no TTL was delivered)
    trgs=sort([trgs;  ttls(1)]); %Concatenate all 2500 TTL onsets.
    trgs=trgs/fs; %convert to units = seconds.
    trgs1=[trgs1;trgs+1000*(i1-1)]; %add 1000 seconds to each subsequent repetition

    %APs=find(traces(:,1)>-40); %get list of all time points during the Action Potentials.
    trs=traces(:,1);
    [pks1, locs1] = findpeaks (traces(:,1),  'MinPeakHeight' , -40);
    spks=locs1;
    shrtpnts=spks(diff(spks)<50.5);
    spks(diff(spks)<50.5)=[];
    
    t=linspace(0,length(traces(:,1))/fs,length(traces(:,1)));
    figure; plot([traces(:,1),traces(:,2)]);hold on
    plot (spks,trs(spks),'*') 
    plot( shrtpnts,trs(shrtpnts),'o')
    title(['Raw Traces ',name]);
    xlabel('time sec'); ylabel('mV');
    %if isempty(APs); else %in case there are no spikes
    %spks=APs(find(diff(APs)>10)+1); %onset timestamps for the APs.
    %spks=sort([spks; APs(1)]); %include first spike timestamp.
    spks=spks/fs; %convert to units = seconds.
    spks1=[spks1;spks+1000*(i1-1)]; %add 1000 seconds to each subsequent repetition   
    %end

end
end
load('2019_10_09.mat')
A2a=trgs1;
adch_12a=spks1;
clearvars -except A1a A2a A3a adch_12a trgss name 
eval([name,'=adch_12a'])
clear adch_12a name

