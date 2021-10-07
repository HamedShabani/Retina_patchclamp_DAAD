
FolderName='D:\Hame2\Data\2020_07_23\2020_07_23_l1\'
fname='2020_07_23_l1_cathodic'
trgs_name='Heval_Cathodic'
set (0, 'DefaultTextInterpreter' , 'none' )

jnk=3;
Fs=50e3
stim_dur=5
fpath=FolderName
RawDat=load([fpath fname]); %load the workspace
vnames=sort(fieldnames(RawDat)); %get unit names
Ythix={'v = -111','i = -250','i = -428','i = -666','i =-1000','i = 666','i = 428','i = 250','i = 111'};
Xthix={'i = 1','i = 2.5','i = 5','i = 6.25','i =7.5','i = 8.7','i = 10'};
for i1=1:length(vnames)-jnk
    Data(i1).name=vnames(i1+jnk);
    eval(['Data(i1).spks=RawDat.',Data(i1).name{:},';']);
    
end
stm=RawDat.A2a;
close all
for i2=1:length(Data)-1
ch_name= char(Data(i2).name());
ch_name=ch_name(3:end);
   
data=Data(i2).spks;

trls=length(stm)/9
stim_order=[1:9];
spk=[];
lbl=[]
k=0
A2a=stm
for ord=1:9
    for kk= 1:trls
        k=(k+1);
        trl=(kk-1)*9;
    spk=[spk;data(data>A2a(trl+stim_order(ord))& data<[A2a(trl+stim_order(ord))+5])-A2a(trl+stim_order(ord))];
    tmp=data(data>A2a(trl+stim_order(ord))& data<[A2a(trl+stim_order(ord))+5])-A2a(trl+stim_order(ord));
    lbl=[lbl;k*ones(length(tmp),1)];
    end
    
end
fig=figure
line([spk'; spk'],[lbl';lbl'+1],'Color','k','MarkerSize',4)
hold on
line([zeros(1,9); ones(1,9)*5],[[trls+1:trls:10*trls];[trls+1:trls:10*trls]],'Color','k','MarkerSize',4)
ytic=[trls/2:trls:10*trls];

xlim([-.1,2])
ylim([0,trls*9])
yticks(ytic)

yticklabels(Ythix)
%xticklabels(Xthix)
%plot(spk_nmbr)
title(['Spike Raster' ,' Channel  ',ch_name(4:end),'   ',fname]);
%legend('Response','Stimulus')
xlabel('Time s')

folderName='Figures';
fn = fullfile(fpath,folderName);
    if ~ exist(fn, 'dir' )
        mkdir(fpath, 'Figures')
    end
print(fig,'-dpng',[fpath,'\Figures\',trgs_name,' stimulus ',fname(1:end-4),'  ',ch_name,'2.png'],'-r300')


end


    
