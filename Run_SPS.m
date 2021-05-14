clear
set (0, 'DefaultTextInterpreter' , 'none' )
datafile='D:\Hame2\Data\Sydney_data\191009\esta\191009 Data Files.mat'
SampRate=10000;% sample rate
thr=-40% threshold for spike detection
close all

for ii=[3:4]
    
    switch ii
        case 1
            recs={'Cell1_191009_SPSVoltage_BD'};
            name='Cell1_191009_BD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)

        case 2
            recs={'Cell1_191009_SPSVoltage_AD'};
            name='Cell1_191009_AD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)

           
       case 3
           SampRate=50000;

            recs={'Cell1_191009_SPSCurrent1ms_AD'};
            name='Cell1_191009_AD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
           
        case 4
            SampRate=50000;

            recs={'Cell1_191009_SPSCurrent1ms_BD'};
            name='Cell1_191009_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        
    end




end

 %%
 close all
clear
 set (0, 'DefaultTextInterpreter' , 'none' )

datafile=('D:\Hame2\Data\Sydney_data\191011\esta\191011Experiments.mat')
SampRate=10000;
thr=-40
 for ii=3:3
%%TOGGLE: To run either the short or long pulses

    switch ii
        case 1
            SampRate=50000;
            recs={'OFFS_191011_SPSCurrent1ms_BD'};
            name='OFFS_191011_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        case 2
            SampRate=50000;
            recs={'OFFT_191011_SPSCurrent1ms_BD'};
            name='OFFT_191011_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        case 3
            SampRate=50000;
            recs={'ONOFF_191011_SPSCurrent1ms_BD'};
            name='ONOFF_191011_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
            
            
            
        case 4
            recs={'OFFS_191011_SPSVoltage_BD'};
            name='OFFS_191011_BD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        case 5
            recs={'OFFT_191011_SPSVoltage_BD'};
            name='OFFT_191011_BD';
 
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        case 6
            recs={'ONOFF_191011_SPSVoltage_BD'};
            name='ONOFF_191011_BD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)            
       
    end

 end
 

 %% 20190821
 close all
clear
 set (0, 'DefaultTextInterpreter' , 'none' )

datafile=('D:\Hame2\Data\Sydney_data\190821\sta\190821 Data Files.mat')
SampRate=10000;
thr=-40
 for ii=8:8
%%TOGGLE: To run either the short or long pulses

    switch ii
        case 1
            recs={'OFFT1_190821_SPSCurrent100us_BD'};
            name='OFFT1_190821_BD100us';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        case 2
            thr=-25
            recs={'OFFT1_190821_SPSCurrent100us_AD'};
            name='OFFT1_190821_BD100us';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
            
        case 3
            recs={'OFFT1_190821_SPSCurrent1ms_BD'};
            name='OFFT1_190821_BD1ms';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        case 4
            recs={'OFFT1_190821_SPSCurrent1ms_AD'};
            name='OFFT1_190821_AD1ms';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)   
            
        case 5
            recs={'OFFT2_190821_SPSCurrent100us_BD'};
            name='OFFT2_190821_BD100us';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        case 6
            recs={'OFFT2_190821_SPSCurrent100us_AD'};
            name='OFFT2_190821_AD100us';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        case 7
            thr=-40
            recs={'OFFT2_190821_SPSCurrent1ms_BD'};
            name='OFFT2_190821_BD1ms';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
        case 8
            recs={'OFFT2_190821_SPSCurrent1ms_AD'};
            name='OFFT2_190821_AD1ms';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)   
                        
            
            
            
            
            
        case 9
            recs={'OFFT1_190821_SPSVoltage_BD'};
            name='OFFT1_190821_BD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        case 10
            recs={'OFFT1_190821_SPSVoltage_AD'};
            name='OFFT1_190821_AD';
 
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        case 11
            recs={'OFFT2_190821_SPSVoltage_BD'};
            name='OFFT2_190821_BD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)            
       
       case 12
           thr=-50
            recs={'OFFT2_190821_SPSVoltage_AD'};
            name='OFFT2_190821_AD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)     
            
            
    end

 end

  %% 20190819
close all
clear
set (0, 'DefaultTextInterpreter' , 'none' )

datafile=('D:\Hame2\Data\Sydney_data\190819\esta\OFFT_190819_SPSVoltage_AD.mat')
SampRate=10000;
thr=-40
 for ii=1:1
%%TOGGLE: To run either the short or long pulses
    switch ii
        case 1
            recs={'OFFT_190819_SPSVoltage_AD'};
            name='OFFT_190819_SPS_AD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)

        case 2
            recs={'OFFT_190819_SPSCurrent100us_BD'};
            name='OFFT_190819_SPS100us_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)

    end
 end
 
 %%
  %% 20190822
close all
clear
set (0, 'DefaultTextInterpreter' , 'none' )

datafile=('D:\Hame2\Data\Sydney_data\191022\Second_cell\esta\191022 data files cell 2.mat')
SampRate=10000;
thr=-30
 for ii=2:2
%%TOGGLE: To run either the short or long pulses

    switch ii
        case 1
            recs={'OFFT_191022_SPSVoltage_BDrelocatedstim'};
            name='OFFT_191022_SPS_BDrelocatedstim';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        
       
        case 2
            SampRate=50000;
            recs={'OFFT_191022_SPSCurrent1ms_BD'};
            name='OFFT_191022_SPS1ms_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)

    end
 end
 
   %% 20191022
   
close all
clear
set (0, 'DefaultTextInterpreter' , 'none' )

datafile=('D:\Hame2\Data\Sydney_data\191022\191022 data files.mat')
SampRate=10000;
thr=-30
 for ii=1:2
%%TOGGLE: To run either the short or long pulses

    switch ii
        case 1
            recs={'OFFT_191022_SPSVoltage_BD'};
            name='OFFT_191022_SPS_BD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        
       
        case 2
            SampRate=50000;
            recs={'OFFT_191022_SPSCurrent1ms_BD'};
            name='OFFT_191022_SPS1ms_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)

    end
 end
 
  
   %% 20191015
   
   close all
clear
set (0, 'DefaultTextInterpreter' , 'none' )

datafile=('D:\Hame2\Data\Sydney_data\191015\esta\191015 data files.mat')
SampRate=10000;
thr=-30
 for ii=1:4
%%TOGGLE: To run either the short or long pulses

    switch ii
        case 2
            recs={'ONS_191015_SPSVoltage_AD'};
            name='ONS_191015_SPS_AD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        case 1
            recs={'ONS_191015_SPSVoltage_BD'};
            name='ONS_191015_SPS_BD';
            Name=['Voltage  ',name];
            plot_voltage_sps(datafile,recs,Name,thr,SampRate)
        
       
        case 3
            SampRate=50000;
            recs={'ONS_191015_SPSCurrent1ms_BD'};
            name='ONS_191015_SPS1ms_BD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)
            
        case 4
            SampRate=50000;
            recs={'ONS_191015_SPSCurrent1ms_AD'};
            name='ONS_191015_SPS1ms_AD';
            Name=['Current  ',name];
            plot_current_sps(datafile,recs,Name,thr,SampRate)    

    end
 end
 
  
