%% 1000
FolderName='D:\Hame2\Data\Sydney_data\191022\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =1000_contrast =350.mat'
%stim_dur=4;% stimulus duration sec 
%sta=e_sta_noise(FolderName,'2019_10_09.mat',jnk,50,.01,[]);
%sta=e_sta_nose_embeded(FolderName,'2019_10_22_1000_new.mat',jnk,50,0.005,[]);
%sta=e_sta(FolderName,'2019_10_22_1000_new.mat',jnk,50,0.0,[]);
sta=e_sta_sps(FolderName,'2019_10_22_1000_new.mat',[],jnk,50,0.01,[],stimfile)
%% 400
close all

FolderName='D:\Hame2\Data\Sydney_data\191022\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =400_contrast =140.mat'

sta=e_sta_sps(FolderName,'2019_10_22_400_new.mat',[],jnk,50,0,[],stimfile)
sta=e_sta_sps(FolderName,'2019_10_22_400_new.mat',[],jnk,50,.010,[],stimfile)

%% 800
FolderName='D:\Hame2\Data\Sydney_data\191022\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =800_contrast =240.mat'
sta=e_sta_sps(FolderName,'2019_10_22_800_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'2019_10_22_800_new.mat',[],jnk,50,0,[],stimfile)
%% 2u
FolderName='D:\Hame2\Data\Sydney_data\191022\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Current\egwn_1_Freq =25_Mean =2_contrast =0.7.mat'
sta=e_sta_sps(FolderName,'2019_10_22_2u_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'2019_10_22_2u_new.mat',[],jnk,50,0,[],stimfile)
%% 4u
FolderName='D:\Hame2\Data\Sydney_data\191022\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Current\egwn_1_Freq =25_Mean =4_contrast =1.4.mat'
sta=e_sta_sps(FolderName,'2019_10_22_4u_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'2019_10_22_4u_new.mat',[],jnk,50,0,[],stimfile)
%% 800 embedded
FolderName='D:\Hame2\Data\Sydney_data\191009\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\sin_noise_embeded_800.mat'
sta=e_sta_sps(FolderName,'2019_10_09_Cell1__191009__BD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'2019_10_09_Cell1__191009__BD_new.mat',[],jnk,50,0,[],stimfile)

%% 200 embedded
FolderName='D:\Hame2\Data\Sydney_data\191024\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =200_contrast =70.mat'
sta=e_sta_sps(FolderName,'2019_10_24_OFFS_191024_200mv_BD2_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'2019_10_24_OFFS_191024_200mv_BD2_new.mat',[],jnk,50,0,[],stimfile)
%% 400 embedded
FolderName='D:\Hame2\Data\Sydney_data\191024\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Current\egwn_1_Freq =25_Mean =4_contrast =1.4.mat'
sta=e_sta_sps(FolderName,'2019_10_24_ONS_191024_4uA_BD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'2019_10_24_ONS_191024_4uA_BD_new.mat',[],jnk,50,0,[],stimfile)
%%

%% 800 embedded fs =10k

close all
FolderName='D:\Hame2\Data\Sydney_data\190821\sta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='sin_noise_embeded_800'
stimfile='D:\Hame2\Data\Sydney_data\sin_noise_embeded_800.mat'

sta=e_sta_sps(FolderName,'OFFT1_190821_NoiseEmbedd_BD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT1_190821_NoiseEmbedd_BD_new.mat',[],jnk,50,0,[],stimfile)

sta=e_sta_sps(FolderName,'OFFT1_190821_NoiseEmbedd_AD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT1_190821_NoiseEmbedd_AD_new.mat',[],jnk,50,0,[],stimfile)


sta=e_sta_sps(FolderName,'OFFT2_190821_NoiseEmbedd_BD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT2_190821_NoiseEmbedd_BD_new.mat',[],jnk,50,0,[],stimfile)


%% 800 embedded fs =10k

close all
FolderName='D:\Hame2\Data\Sydney_data\190819\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stimfile='sin_noise_embeded_800'
stimfile='D:\Hame2\Data\Sydney_data\sin_noise_embeded_800.mat'

sta=e_sta_sps(FolderName,'OFFT_190819_NoiseEmbedd_BD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_190819_NoiseEmbedd_BD_new.mat',[],jnk,50,0,[],stimfile)

sta=e_sta_sps(FolderName,'OFFT_190819_NoiseEmbedd_AD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_190819_NoiseEmbedd_AD_new.mat',[],jnk,50,0,[],stimfile)


sta=e_sta_sps(FolderName,'ONS_190819_NoiseEmbedd_BD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'ONS_190819_NoiseEmbedd_BD_new.mat',[],jnk,50,0,[],stimfile)

sta=e_sta_sps(FolderName,'ONS_190819_NoiseEmbedd_AD_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'ONS_190819_NoiseEmbedd_AD_new.mat',[],jnk,50,0,[],stimfile)

%% 800 embedded fs =50k

close all
FolderName='D:\Hame2\Data\Sydney_data\191022\Second_cell\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;

stimfile='sin_noise_embeded_800'
stimfile='D:\Hame2\Data\Sydney_data\sin_noise_embeded_800.mat'
stimfile='D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =800_contrast =240.mat'
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_800mv_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_800mv_new.mat',[],jnk,50,0,[],stimfile)


stimfile='egwn_1_Freq =25_Mean =400_contrast =140.mat'
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_400mv_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_400mv_new.mat',[],jnk,50,0,[],stimfile)



stimfile='egwn_1_Freq =25_Mean =200_contrast =70.mat'
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_200mv_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_200mv_new.mat',[],jnk,50,0,[],stimfile)


stimfile='egwn_1_Freq =25_Mean =1000_contrast =350.mat'
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_1000mv_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_1000mv_new.mat',[],jnk,50,0,[],stimfile)


stimfile='egwn_1_Freq =25_Mean =4_contrast =1.4.mat'
sta=e_sta_sps(FolderName,'OFFT_191022__second_cell_4uA_new.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_191022__second_cell_4uA_new.mat',[],jnk,50,0,[],stimfile)


stimfile='D:\Hame2\Data\Sydney_data\Stimili\Voltage\egwn_1_Freq =25_Mean =800_contrast =240.mat'
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_800mv_relocatedstim.mat',[],jnk,50,.010,[],stimfile)
sta=e_sta_sps(FolderName,'OFFT_191022_Second_cell_BD_800mv_relocatedstim.mat',[],jnk,50,0,[],stimfile)



%%
FolderName='C:\Users\abcd1234\Desktop\Hamed\RGC\Sydney_data\Mohit Collab\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stim_dur=4;% stimulus duration sec 
%sta=e_sta_noise(FolderName,'2019_10_09.mat',jnk,50,.01,[]);
sta=e_sta_nose_embeded(FolderName,'2019_08_19.mat',jnk,50,0.005,[]);

%%

FolderName='C:\Users\abcd1234\Desktop\Hamed\RGC\Sydney_data\Mohit Collab\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stim_dur=4;% stimulus duration sec 
%sta=e_sta_noise(FolderName,'2019_10_09.mat',jnk,50,.01,[]);
%sta=e_sta_nose_embeded(FolderName,'2019_07_21_wl_secondmouce.mat',jnk,50,0.00,[]);
sta=e_sta_nose_embeded(FolderName,'2019_07_21_wl_secondmouce.mat',jnk,50,0.00,[1:40000]);
%%

FolderName='D:\Hame2\Data\Sydney_data\191022\esta\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stim_dur=4;% stimulus duration sec 
%sta=e_sta_noise(FolderName,'2019_10_09.mat',jnk,50,.01,[]);
%sta=e_sta_nose_embeded(FolderName,'2019_07_21_wl_secondmouce.mat',jnk,50,0.00,[]);
%sta=e_sta_nose_embeded(FolderName,'2019_10_22_400.mat',jnk,50,0.005,[]);
sta=e_sta(FolderName,'2019_10_22_400.mat',jnk,50,0.0,[],stim_file);
sta=e_sta(FolderName,'2019_10_22_800.mat',jnk,50,0.01,[]);
sta=e_sta(FolderName,'2019_10_22_1000.mat',jnk,50,0.0,[]);

%%
FolderName='D:\Hame2\Data\Sydney_data\191022\current\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stim_dur=4;% stimulus duration sec 
sta=e_sta(FolderName,'2019_10_22_4u.mat',jnk,50,0.0,[]);
%%
FolderName='D:\Hame2\Data\Sydney_data\191022\Second_cell\';cd(FolderName)
bwid=.1;%binwidth sec
jnk=3;
stim_dur=4;% stimulus duration sec 
sta=e_sta(FolderName,'2019_10_22_4u.mat',jnk,50,0.0,[]);
