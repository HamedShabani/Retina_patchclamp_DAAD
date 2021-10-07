

%note:  make sure that the thresholds for TTL and spike detection are valid
%for each recording.
%Sampled at 50kHz. 20 samples during TTL.  Must have been a post time of
%1000 microsec after the 1000 microsec voltage stimulus pulse.

load('D:\Hame2\Data\Sydney_data\191015\esta\191015 Data Files.mat');
fs=50000
for ii=2:2
    switch ii
        case 1
            recs={'ONS_191015_NoiseEmbeddedTrial1_AD','ONS_191015_NoiseEmbeddedTrial2_AD','ONS_191015_NoiseEmbeddedTrial3_AD'...
                'ONS_191015_NoiseEmbeddedTrial4_AD','ONS_191015_NoiseEmbeddedTrial5_AD'};
            name='ONS_191015_AD';
        case 2
            recs={'ONS_191015_NoiseEmbeddedTrial1_BD','ONS_191015_NoiseEmbeddedTrial2_BD','ONS_191015_NoiseEmbeddedTrial3_BD'...
                'ONS_191015_NoiseEmbeddedTrial4_BD','ONS_191015_NoiseEmbeddedTrial5_BD'};
            name='ONS_191015_BD';

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

