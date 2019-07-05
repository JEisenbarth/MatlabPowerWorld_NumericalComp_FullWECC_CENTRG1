function [residual] = residual_PowerWorld(theta,theta_indicies,index,data_event,filenamedyd,genrou,exac8b,pss2a,PQ_Flag,SimAuto)
%residual_PowerWorld This function calculates the P, Q, or PQ residual based on a given theta. 
%   PQ_Flag=0 then use P Residual
%   PQ_Flag=1 then use Q Residual
%   PQ_Flag=2 then use PQ Residual

%% Put thetas into model.
for m=1:length(theta)
    if theta_indicies(m,1)==1
        genrou(index.genrou(theta_indicies(m,2)))=theta(m);
    elseif theta_indicies(m,1)==2
        exac8b(index.exac8b(theta_indicies(m,2)))=theta(m);
    elseif theta_indicies(m,1)==3
        pss2a(index.pss2a(theta_indicies(m,2)))=theta(m);

    end           
end
theta
%% Write dyd and run simulation
[data]=PowerWorld_WriteDYD_Run(filenamedyd,genrou,exac8b,pss2a,SimAuto);

%% Find indicies for residual calculations then calc residual.
ndxP=PWFind(data,'Branch ',' 47741 47740 1 ','MW To');
ndxP_event=PWFind(data_event,'Branch ',' 47741 47740 1 ','MW To');
ndxQ=PWFind(data,'Branch ',' 47741 47740 1 ','Mvar To');
ndxQ_event=PWFind(data_event,'Branch ',' 47741 47740 1 ','Mvar To');

%% Calculate residual based on PQ_Flag.
if(PQ_Flag==0)
    residual=data_event.Data(:,ndxP_event)-data.Data(:,ndxP);
elseif(PQ_Flag==1)
    residual=data_event.Data(:,ndxQ_event)-data.Data(:,ndxQ);
elseif(PQ_Flag==2)
    residual=[data_event.Data(:,ndxP_event);data_event.Data(:,ndxQ_event)]-[data.Data(:,ndxP);data.Data(:,ndxQ);];
end


% %% Plot Check
% figure
% subplot(3,1,1)
% hold on
% plot(data_event.Data(:,1),data_event.Data(:,ndxP_event),'LineWidth',1,'DisplayName','Event')
% plot(data.Data(:,1),data.Data(:,ndxP),'LineWidth',1,'DisplayName','PlayIn')
% hold off
% title('P Plot')
% legend();
% 
% subplot(3,1,2)
% hold on 
% plot(data_event.Data(:,1),data_event.Data(:,ndxQ_event),'LineWidth',1,'DisplayName','Event')
% plot(data.Data(:,1),data.Data(:,ndxQ),'LineWidth',1,'DisplayName','PlayIn')
% hold off
% title('Q Plot')
% legend();
% 
% subplot(3,1,3)
% hold on
% plot(residual,'DisplayName',['KS1=',num2str(theta)])
% hold off
% title('Residual')
% legend();
% 
% ndxV=PWFind(data,'Bus ',' 47741 ','V pu');
% ndxV_event=PWFind(data_event,'Bus ',' 47741 ','V pu');
% figure
% subplot(2,1,1)
% hold on
% plot(data_event.Data(:,1),data_event.Data(:,ndxV_event),'LineWidth',1,'DisplayName','Event')
% plot(data.Data(:,1),data.Data(:,ndxV),'LineWidth',1,'DisplayName','PlayIn')
% hold off
% title('V Plot')
% legend();
% 
% ndxVang=PWFind(data,'Bus ',' 47741 ','V angle No shift');
% ndxVang_event=PWFind(data_event,'Bus ','V angle No shift');
% subplot(2,1,2)
% hold on
% plot(data_event.Data(:,1),data_event.Data(:,ndxV_event),'LineWidth',1,'DisplayName','Event')
% plot(data.Data(:,1),data.Data(:,ndxV),'LineWidth',1,'DisplayName','PlayIn')
% hold off
% title('Vang Plot')
% legend();


end