clear
trgs_name_all={'Heval_Anodic','Heval_Cathodic'}
trgs_ch={'RawDat.A2a_a','RawDat.A2a_c'}

FolderName='D:\Hame2\Data\2020_08_04\2020_08_l1\'
fname='2020_08_l1_all'
fpath=FolderName
RawDat=load([fpath fname]); %load the workspace
set (0, 'DefaultTextInterpreter' , 'none' )
Fs=50e3
stim_dur=5
        vnames=sort(fieldnames(RawDat)); %get unit names
        
        %Xthix={'i = 1','i = 2.5','i = 5','i = 6.25','i =7.5','i = 8.7','i = 10'};
        jnk=6
        for i1=1:length(vnames)-jnk
            Data(i1).name=vnames(i1+jnk);
            eval(['Data(i1).spks=RawDat.',Data(i1).name{:},';']);
            
        end
%stm=RawDat.A2a_a;
for i2 =1:length(Data)-1
fig = figure('units','normalized','outerposition',[0.1 0.1 .5 .5],'visible','off');

    for ii=1:2
        stm=eval(trgs_ch{ii});
        trgs_name=trgs_name_all{ii};
        
        if contains(trgs_name, 'Anodic')
            load('D:\Hame2\Scripts\Heval_stim\Anodic_Heval_10_rep.mat')
%            Ythix={'v = 111/-1000(10%)','v = 250/-1000(20%)','v = 428/-1000(30%)','v = 666/-1000(40%)','v =1000/-1000(50%)','v =1000/-666(60%)','v =1000/-428(70%)','v = -1000/250(80%)','v = -1000/111(90%)'};
            Ythix={'v = 1000/-111(10%)','v = 1000/-250(20%)','v = 1000/-428(30%)','v = 1000/-666(40%)','v =1000/-1000(50%)','v =666/-1000(60%)','v =428/-1000(70%)','v = 250/-1000(80%)','v = 111/-1000(90%)'};

            jnk=6;
            
        elseif contains(trgs_name, 'Cathodic')
            Ythix={'v = -1000/111(10%)','v = -1000/250(20%)','v = -1000/428(30%)','v = -1000/666(40%)','v =-1000/1000(50%)','v =-666/1000(60%)','v =-428/1000(70%)','v = -250/1000(80%)','v = -111/1000(90%)'};
%            Ythix={'v = -111/1000(10%)','v = -1000/250(20%)','v = -1000/428(30%)','v = -1000/666(40%)','v =-1000/1000(50%)','v =1000/-666(60%)','v =1000/-428(70%)','v = 1000/-250(80%)','v = 1000/-111(90%)'};

            load('D:\Hame2\Scripts\Heval_stim\Cathodic_Heval_10_rep.mat')
            
            jnk=6;
        end
        
        [val ,idx]=sort(DD(:,[3])./DD(:,[5]));
        DDD=DD(idx,:);

        %[val ,idx]=sort(abs(DD(:,[3,5])));
        
        ch_name= char(Data(i2).name());
        ch_name=ch_name(3:end);
        
        data=Data(i2).spks;
        
        trls=length(stm)/9;
        stim_order=idx;
        spk=[];
        lbl=[];
        k=0;
        A2a=stm;
        ord=0;
        for ordd=1:9
            for kk= 1:trls
                ord=ord+1;
                
                k=(k+1);
                trl=0;
                spk=[spk;data(data>A2a(trl+stim_order(ord))& data<[A2a(trl+stim_order(ord))+5])-A2a(trl+stim_order(ord))];
                tmp=data(data>A2a(trl+stim_order(ord))& data<[A2a(trl+stim_order(ord))+5])-A2a(trl+stim_order(ord));
                lbl=[lbl;k*ones(length(tmp),1)];
            end
            
        end
        subplot(1,2,ii)
%        set(fig, 'PaperSize',[21 85])
        line([spk'; spk'],[lbl';lbl'+1],'Color','k','MarkerSize',5)
        hold on
        line([zeros(1,9); ones(1,9)*5],[[trls+1:trls:10*trls];[trls+1:trls:10*trls]],'Color','k','MarkerSize',4)
        ytic=[trls/2:trls:10*trls];
        
        xlim([-.1,.5])
        ylim([0,trls*9])
        yticks(ytic)
        
        yticklabels(Ythix)
 set(gca,'Fontsize',11);        %xticklabels(Xthix)
        %plot(spk_nmbr)
        title(['Channel  ',ch_name(4:end),'   ',fname,'  ',trgs_name(6:end)]);
        %legend('Response','Stimulus')
        xlabel('Time s')
        

        
    end
    
            folderName='Figures';
        fn = fullfile(fpath,folderName);
        if ~ exist(fn, 'dir' )
            mkdir(fpath, 'Figures')
        end
        print(fig,'-dpng',[fpath,'\Figures\',' Heval stimulus ',fname(1:end-4),'  ',ch_name,' dual.png'],'-r300')
        close all
    
end