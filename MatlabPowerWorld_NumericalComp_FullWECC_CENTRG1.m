% Jacob Eisenbarth, May 2019
% This script will be used to run multiple different PSLF cases within
% Matlab. First, the script writes a PSLF .dyd file from Matlab using a
% series of functions. Next,it will start an EPCL code which will create a
% channel file. Finally, the channel file will be read in then
% stored as a structure which can be saved as a .mat file. Included into this
% code is the addition of changing individual parameters to MINIMIZE THE COST
% Function then continuing the simulation in PSLF.

%% Initialize Matlab
clc,close all,clear all, format longG,format compact

tic
%% Load Event Data
% data_event=udread("C:\upslf19\MyPSLF\16ls1a_DropGen_PALOVRD1.chf",[]);
% data_event.Data([243,17573],:)=[];  %Get rid of multiple points at same time

%% Establish a connection with PowerWorld / SimAuto
disp('>> Connecting to PowerWorld Simulator / SimAuto...')
SimAuto = actxserver('pwrworld.SimulatorAuto');
disp('Connection established')

%% Data to be Entered into dyd File
%Model genrou
genrou(1)=5.8; %Tpdo
% genrou(2)=0.015; %Tppdo
genrou(2)=0.016; %Tppdo
genrou(3)=0.6; %Tpqo
genrou(4)=0.03; %Tppqo
genrou(5)=3.25; %H
genrou(6)=0; %D
genrou(7)=2.05; %Ld
genrou(8)=1.95; %Lq
genrou(9)=0.42; %Lpd
genrou(10)=0.65; %Lpq
genrou(11)=0.24; %Lppd
genrou(12)=0.12; %Ll
genrou(13)=0.125; %S1
genrou(14)=0.33; %S12
genrou(15)=0.0019; %Ra
genrou(16)=0; %Rcomp
genrou(17)=0.063; %Xcomp

index_genrou=[1:5,7:15,17];     %Index for numerical parameters to edit

%Model exac8b
exac8b(1)=0.02; %Tr
exac8b(2)=200; %Kvp
exac8b(3)=0; %Kvi
exac8b(4)=60; %Kvd
exac8b(5)=0.02; %Tvd
exac8b(6)=999; %Vimax
exac8b(7)=0.02; %Ta
exac8b(8)=10; %Vrmax
exac8b(9)=-10; %Vrmin
exac8b(10)=1; %Ke
exac8b(11)=1.5; %Te
exac8b(12)=0.15; %Kc
exac8b(13)=0.45; %Kd
exac8b(14)=6.5; %E1
exac8b(15)=0.3; %S(E1)
exac8b(16)=9; %E2
exac8b(17)=3; %S(E2)
exac8b(18)=0; %limflg

index_exac8b=[1:2,4:17];     %Index for numerical parameters to edit

%Model pss2a
pss2a(1)=2; %J1
pss2a(2)=0; %K1
pss2a(3)=3; %J2
pss2a(4)=0; %K2
pss2a(5)=15; %Tw1
pss2a(6)=5; %Tw2
pss2a(7)=15; %Tw3
pss2a(8)=1; %Tw4
pss2a(9)=1; %T6
pss2a(10)=5; %T7
pss2a(11)=0.3077; %Ks2
pss2a(12)=1; %Ks3
pss2a(13)=1; %Ks4
pss2a(14)=0.4; %T8
pss2a(15)=1; %T9
pss2a(16)=1; %n
pss2a(17)=1; %m
pss2a(18)=2; %Ks1
pss2a(19)=0.4; %T1
pss2a(20)=0.2; %T2
pss2a(21)=0.4; %T3
pss2a(22)=0.2; %T4
pss2a(23)=0.05; %Vstmax
pss2a(24)=-0.05; %Vstmin
pss2a(25)=1; %a
pss2a(26)=0.4; %Ta
pss2a(27)=0.2; %Tb

index_pss2a=[5:15,18:27];     %Index for numerical parameters to edit

genrou_original=genrou;
exac8b_original=exac8b;
pss2a_original=pss2a;


index=struct('genrou',index_genrou,'exac8b',index_exac8b,'pss2a',index_pss2a);

clear index_genrou index_exac8b index_oel1 index_pss2a index_hygovr

% list=[ones(length(index.genrou),1),[1:length(index.genrou)]',ones(length(index.genrou),1)*1;...
%     2*ones(length(index.exac8b),1),[1:length(index.exac8b)]',ones(length(index.exac8b),1)*1;...
%     3*ones(length(index.pss2a),1),[1:length(index.pss2a)]',ones(length(index.pss2a),1)*1;];

list=[1,1,1;
    1,5,1;
    1,13,1;
    2,2,1;
    2,5,1;
    2,10,1;
    3,13,1;
    3,14,1;]
    
    
% list=1;
for k=1:size(list,1)
    numericalsims=100;
    for x=1:numericalsims
        %% Write Data Aux File for Centralia Event PlayIn
        datacsv=load(['D:\Users\JEisenbarth\Desktop\PowerWorld Files\FULL WECC and CENTG1 PlayIn\CJ_EventData.Mat']);
        
        texpected=[0:.004:20];
        datacsv.Data=TrimEventData(datacsv.Data,texpected);
        ndx1=PWFind(datacsv,'Bus ',' 47741 ','V pu');
        ndx2=PWFind(datacsv,'Bus ',' 47741 ','V angle No shift');
        ndx3=PWFind(datacsv,'Bus ',' 47741 ','Frequency in PU');
        
        t1=datacsv.Data(:,1);
        v1=datacsv.Data(:,ndx1);
        vang1=datacsv.Data(:,ndx2);
        finitial=datacsv.Data(1,ndx3);     %PU
        forig=CalcFfromVang(datacsv.Data(:,ndx2),datacsv.Data(:,1),finitial)';
        
        %%ADD Filtered Measurement Noise to V and F
        vmeasnoise=.001*randn(length(v1),1);
        vangmeasnoise=.25*randn(length(vang1),1);
        
        [b,a]=butter(3,.05);
        
        vmeasnoise=filter(b,a,vmeasnoise);
        vangmeasnoise=filter(b,a,vangmeasnoise);
        
        v1=v1+vmeasnoise;
        vang1=vang1+vangmeasnoise;
        
        f1=CalcFfromVang(vang1,datacsv.Data(:,1),finitial)';

        filename=['C:\MatlabPowerWorld_NumericalComp_FullWECC_CENTRG1\PlayInData.aux'];
        
        WritePlayInAux(filename,t1,v1,f1)
        data_PlayInAux(k,x).t1=t1;
        data_PlayInAux(k,x).v1=v1;
        data_PlayInAux(k,x).f1=f1;
        data_PlayInAux(k,x).vang1=vang1;
        
        
%         figure
%         subplot(3,1,1)
%         plot(t1,(v1+vmeasnoise),t1,v1)
%         legend('Added Measurement Noise','Original')
%         title('CENTR G1 Event 3: Voltage')
%         grid
%         xlim([0,20])
%         
%         subplot(3,1,2)
%         plot(t1,vang1,t1,(vang1-vangmeasnoise))
%         legend('Added Measurement Noise','Original')
%         title('CENTR G1 Event 3: Frequency')
%         grid
%         xlim([0,20])
%         
%         subplot(3,1,3)
%         plot(t1,(f1),t1,forig)
% 
%         legend('Added Measurement Noise','Original')
%         title('CENTR G1 Event 3: Frequency')
%         grid
%         xlim([0,20])

        
        %% Setup to Run to Minimize Cost Function
        %Setup Column Vector of Parameter to Adjust
        %         theta_indicies=[ones(length(index.genrou),1),[1:length(index.genrou)]',ones(length(index.genrou),1)*1;...
        %             2*ones(length(index.exac8b),1),[1:length(index.exac8b)]',ones(length(index.exac8b),1)*1;...
        %             3*ones(length(index.pss2a),1),[1:length(index.pss2a)]',ones(length(index.pss2a),1)*1;];
        
        %         theta_indicies=[[1,1,1];[1,2,1];[1,3,1];[1,4,1];[1,5,1];[1,6,1];[1,7,1];[1,9,1];[1,10,1];[1,11,1];[1,12,1];[1,13,1];[1,14,1];[1,15,1];[1,16,1];[2,1,1];[2,3,1];[2,6,1];[2,7,1];[2,8,1];[2,9,1];[4,3,1];[4,4,1];[4,7,1];[4,8,1];[4,9,1];[4,10,1];[4,11,1];[4,12,1];[4,13,1];[4,14,1];[4,15,1];[4,16,1];[4,17,1];[4,18,1];[4,19,1];[4,20,1];[4,24,1];[5,1,1];[5,2,1];[5,3,1];[5,4,1];[5,5,1];[5,6,1];[5,12,1];[5,15,1];[5,16,1];[5,19,1];[5,20,1];[5,21,1];[5,22,1];[5,24,1];[5,25,1];[5,26,1];[5,27,1]];
        
        theta_indicies=list(k,:);
        
        %             theta_indicies=[1,4,1;1,5,1];    %1st column is model,2nd column is numerical parameter,3rd column is what residual vector to use 1=P 2=Q 3=P&Q
        %Ex. [1,5,1]->model=genrou, parameter=H, P for
        % residual calculations.
        %Ex. [2,2,2]->model=exac8b, parameter=Kr, Q for
        % residual calculations.
        
        %Setup theta Vectors
        theta=zeros(size(theta_indicies,1),1);
        for b=1:length(theta)
            if theta_indicies(b,1)==1
                theta(b)= genrou_original(index.genrou(theta_indicies(b,2)));
            elseif theta_indicies(b,1)==2
                theta(b)= exac8b_original(index.exac8b(theta_indicies(b,2)));
            elseif theta_indicies(b,1)==3
                theta(b)= pss2a_original(index.pss2a(theta_indicies(b,2)));
            end
        end
        
        percentnominal=abs(.01*theta);
        
        %% Run Minimizing Cost Function
        %dyd and chf file names
        filenamedyd='C:\MatlabPowerWorld_NumericalComp_FullWECC_CENTRG1\CENTRG1_PlayIn.dyd';
        PQ_Flag=2;
        
    
        %% Run Simulation w Original theta in model
        [data_orig(k,x)] = PowerWorld_WriteDYD_Run(filenamedyd,genrou_original,exac8b_original,pss2a_original,SimAuto);

        
        
        
        %                         opts=optimoptions(@lsqnonlin,'TolFun',1e-12,'Display','iter','Diagnostics','off','Tolx',1e-12,'MaxFunEvals',50000,'SpecifyObjectiveGradient',true);
        opts=optimoptions(@lsqnonlin,'TolFun',1e-12,'Display','iter','Diagnostics','off','Tolx',1e-12,'MaxFunEvals',50000,'DiffMinChange',percentnominal,'SpecifyObjectiveGradient',true);
        %         opts=optimoptions(@lsqnonlin,'TolFun',1e-12,'Display','iter','Diagnostics','off','Tolx',1e-12,'MaxFunEvals',0,'MaxIterations',0,'SpecifyObjectiveGradient',true);
        
        %         residual = @(theta) residual_PowerWorld(theta,theta_indicies,index,datacsv,filenamedyd,genrou_original,exac8b_original,pss2a_original,filenamechf,PQ_Flag,SimAuto);
        residual = @(theta) residual_Jacobian_PowerWorld(theta,theta_indicies,index,datacsv,filenamedyd,genrou_original,exac8b_original,pss2a_original,PQ_Flag,SimAuto,percentnominal);
        %         [final_theta(k,x),resnorm(k,x),residual,exitflag,output(k,x),lambda,Jacobian] = lsqnonlin(residual,theta,[],[],opts);
        if (list(k,:)==[1,2,1])
            ub=Inf(1);
            lb=.016;
        elseif (list(k,:)==[3,3,1])
            lb=.016;
            ub=Inf(1);
            
%         elseif (list(k,:)==[1,11,1])
%             lb=-1*Inf(1);
%             ub=.192;
        else
            lb=-1*Inf(1);
            ub=Inf(1);
        end
        %             lb=(-1* Inf(length(theta),1));
        %         lb([1:4,16,19,21,25,32:36,45,47,52])=1/60;
        %         [final_theta(k,x),resnorm(k,x),residual,exitflag,output(k,x),lambda,Jacobian] = lsqnonlin(residual,theta,[],[],opts);
        [final_theta(k,x),resnorm(k,x),residual,exitflag,output(k,x),lambda,Jacobian] = lsqnonlin(residual,theta,lb,[],opts);
        
        
        %         %%Load channel file from simulation.
        %         data_modified(k,x)=udread(filenamechf,[]);
        
        x
              
        filenamemat=['D:\Users\JEisenbarth\Desktop\PowerWorld Files\CENTRG1 Parameter Testing\Single Parameter Test\FinalTheta_CENTRG1_PowerWorld_WithVFNoiseSingleParam.mat']
        %% Put thetas into model
        for m=1:size(theta_indicies,1)
            if theta_indicies(m,1)==1
                genrou(index.genrou(theta_indicies(m,2)))=final_theta(k,m);
            elseif theta_indicies(m,1)==2
                exac8b(index.exac8b(theta_indicies(m,2)))=final_theta(k,m);
            elseif theta_indicies(m,1)==3
                pss2a(index.pss2a(theta_indicies(m,2)))=final_theta(k,m);
            end
        end
        %% Run Simulation w final theta in model
        [data(k,x)] = PowerWorld_WriteDYD_Run(filenamedyd,genrou,exac8b,pss2a,SimAuto);
        %% Set Model as Original
        genrou=genrou_original;
        exac8b=exac8b_original;
        pss2a=pss2a_original;
        
%         save(filenamemat,'final_theta','resnorm','output','list','data','data_PlayInAux','data_orig')
                
        %         filenamemoddata=['D:\Users\JEisenbarth\Desktop\Full WECC Filtered Noise Parameter Testing\Full Numberical Comp for Noise and No Noise\FinalThetaNoNoise_ModifiedData.mat']
        %         filenamemoddata=['D:\Users\JEisenbarth\Desktop\Full WECC Filtered Noise Parameter Testing\Full Numberical Comp for Noise and No Noise\FinalThetaNoise',num2str(x),'_ModifiedData.mat']
        %         filenamemoddata=['D:\Users\JEisenbarth\Desktop\Full WECC Filtered Noise Parameter Testing\Sensitivity Test for 109 cases\singleparam_',num2str(k),'_ModifiedData.mat']
        
%                 save(filenamemoddata,'data_modified')
        %         clear final_theta
    toc
    end
    
    
end

% %% Plots of Real Power and Reactive Power
% %Plot Real Power
% % ndxP=jfind(data_original,'41742:COULEE22:500:1 :imetr   :pbr');
% % ndxP_event=jfind(data_event,'41742:COULEE22:500:1 :imetr   :pbr');
% Legend_label={};
% figure
% subplot(2,1,1)
% hold on
%
% jplot(data_event,'41742:COULEE22:500:1 :imetr   :pbr ')
% % for k=1:size(data_modified,2)
% for k=1:1
%     jplot_playin(data_modified(k),'41742:COULEE22:500:1 :imetr   :pbr ')
%     Legend_label{k}=['Sens. Test ',num2str(k),': H=',num2str(final_theta(k))];
% end
% hold off
% xlabel('Time')
% ylabel('Real Power MW')
% grid
% legend(['Event: No Noise',Legend_label],'Location','Best','NumColumns',2,'FontSize',6)
% title('Real Power Comparison')
% xlim([0,90])
%
% Legend_label={};
% subplot(2,1,2)
% hold on
%
% jplot(data_event,'41742:COULEE22:500:1 :imetr   :qbr ')
% % for k=1:size(data_modified,2)
% for k=1:1
%     jplot_playin(data_modified(k),'41742:COULEE22:500:1 :imetr   :qbr ')
%     Legend_label{k}=['Sens. Test ',num2str(k),': H=',num2str(final_theta(k))];
% end
% hold off
% xlabel('Time')
% ylabel('Real Power MVAR')
% grid
% legend(['Event: No Noise',Legend_label],'Location','Best','NumColumns',2,'FontSize',6)
% title('Reactive Power Comparison')
% xlim([0,90])