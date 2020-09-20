%% all figures of the paper are made using this script, loading all the results from the analysis
load('all_results_1may2020.mat')
addpath(genpath('/home/np451/Desktop/dataset alignment'));


%% polarisation vs shear stress


figure(10)
ms=8 ; lw=1
x= shear_tap;
y_90= POL(ind_tap); ey= emPOL(ind_tap);
errorbar(x,y_90,ey,'bd','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','b');hold on;
%---exponential fit---%
ft = fittype('1- exp(-x/xc)','independent','x');
fo = fitoptions('Method','NonLinearLeastSquares','Weights',1./ey.^2,'StartPoint',[0.1]);
[fit_out,errfit] = fit(x',y_90,ft,fo);
xc = fit_out.xc;
exc = confint(fit_out); exc= abs(exc(2)-exc(1))/2;

plot([0:0.01:1],fit_out([0:0.01:1]),'b','LineWidth',2)
set(gca,'XScale', 'lin', 'YScale', 'lin');
xlabel('Shear stress $\tau$ [dyne/c$m^2$]','Interpreter','latex') ;ylabel('Alignment $\Phi$ ','Interpreter','latex');
 ylim([-0.5,1.4]); xlim([0.00,0.9])
 x0=0;y0=0;width=600;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
legend({'DIV$>$15','fit DIV$>$15'},'Interpreter','latex','Location','se')

title(strcat('fit x0=',num2str(xc),'+',num2str(exc)));


ms=8 ; lw=1;
figure(11)
x= shear_tap_div;
y_90= POL(ind_tap_div); ey= emPOL(ind_tap_div);
errorbar(x,y_90,ey,'rd','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','r');hold on;
%---exponential fit---%
ft = fittype('1- exp(-x/xc)','independent','x');
fo = fitoptions('Method','NonLinearLeastSquares','Weights',1./ey.^2,'StartPoint',[0.1]);
[fit_out,errfit] = fit(x',y_90,ft,fo);
xc_div = fit_out.xc;
exc_div = confint(fit_out); exc_div= abs(exc_div(2)-exc_div(1))/2;

plot([0:0.01:1],fit_out([0:0.01:1]),'r','LineWidth',2)
set(gca,'XScale', 'lin', 'YScale', 'lin');
xlabel('Shear stress $\tau$ [dyne/c$m^2$]','Interpreter','latex') ;ylabel('Alignment $\Phi$ ','Interpreter','latex');
 ylim([-0.5,1.4]); xlim([0.00,0.9])
 x0=0;y0=0;width=600;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
set(gca,'XScale', 'lin', 'YScale', 'lin');
legend({'DIV$=$8','fit DIV$=$8'},'Interpreter','latex','Location','se')
title(strcat('fit x0=',num2str(xc),'+',num2str(exc),' DIV=8',' x0=',num2str(xc_div),'+',num2str(exc_div)));

%% polarisazion vs shear stress ---- probability distribution

figure();ms=10;lw=2;
ind_plot =  M>0  & Flowtime>1 & Div>14

x=Shear_ind(ind_plot);
shears=unique(x)
cc=1;std_theta=[];
for s=1:numel(shears)
%subplot(4,4,s)
%hold on;
figure()
y= Theta(Shear_ind==shears(s) & ind_plot);
std_theta(cc)=std(y);
y_90 = angle(exp(i*y)*exp(-i*pi/2))
polarhistogram(y_90,'Normalization','pdf');
set(gca,'ThetaZeroLocation','top','ThetaDir','counterclockwise');
x0=0;y0=0;width=400;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);
title(strcat('Shear',num2str(shears(s))));
name= strcat('Shear',num2str(shears(s)));
%saveas(gcf,strcat('/u/homes/np451/Dropbox/alignment with the flow/figure_alignment_shear/',name,'_2.pdf'))
%close('all')
cc=cc+1;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% here I report the code to gather all the analysis results.
% the result is the file all_result loaded at the beginning 
%this part to be runned it need to have all the analysis,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  plot for polarisation vs shear stress with all the data 27.10.19
addpath(genpath('/home/np451/Documents'));
clear all;
close all;
subdirOct={...
    '11.10.19/control/FL','11.10.19/5rpm/FL','11.10.19/5rpm_2/FL'...
    '18.10.19/8rpm/FL','18.10.19/8rpm_2/FL','18.10.19/control/FL','22.10.19/control/FL'...
    '25.10.19/8rpm/FL','25.10.19/8rpm_2/FL','25.10.19/control/FL',...
     '22.10.19/40rpm/FL','22.10.19/40rpm_2/FL','22.10.19/control/FL'...   %%%% data with 40 rpm for one day
}; 

holesOct={...
        'none',[1:8000],'none',...
    'none','none','none','none',...
    'none','none','none',...
   'none','none','none'...
    };

subdirOct = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/Octorber2019/PIV/',subdirOct);
    
subdirJuly={...   
   '16.7.18/0/FL'...
   '29.7.18/0/FL'...
   '16.7.18/5A/FL', '16.7.18/5B/FL', '16.7.18/5C/FL'... %%% straight
  }; 
subdirJuly = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/July2018/PIV/',subdirJuly);
 
subdirJune={...   
   '13.6.18/0/FL'...
   '18.6.18/0/FL'...
...   '18.6.18/100/FL',
   '18.6.18/200/FL'...  %% straight
   }; 
subdirJune = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/June2018/PIV/',subdirJune);

holesPast={...
    'none','none','none','none',...
    'none','none','none','none',...   
};

subdirMarch={...   
'22.3.19/beads/40rpm/FL','22.3.19/beads/40rpm2/FL','22.3.19/beads/control/FL'...
'29.3.19/beads/40rpm/FL','29.3.19/beads/40rpm2/FL','29.3.19/beads/control/FL'...
'5.4.19/beads/40rpm2/FL','5.4.19/beads/40rpm/FL'...,'5.4.19/beads/control/FL'...
'12.4.19/beads/40rpm/FL','12.4.19/beads/40rpm2/FL'...,'12.4.19/beads/control/FL'... 
}; 

subdirMarch = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/PIV/',subdirMarch);

holesMarch={...
    'none','none','none'...  
    [1:2000,4900:10000],'none','none'...% [1:500,500:10000],'none','none'...
    'none','none'...,'none'...  
    'none','none'...  
};

%holesPast={...
%    [1:2000,4900:10000],'none','none'...% [1:500,500:10000],'none','none'...
%    'none',[3865:4075,6415:7195],...
%    'none','none','none','none',...
%    'none','none','none'... %%% added for 22.3.19
%    'none','none','none'... %%% for 29.3.19
% 'none','none','none','none' %%% for 5.4.19
%    };



subdir=[subdirOct,subdirMarch,subdirJuly,subdirJune];

holes=[holesOct,holesPast,holesMarch];

shear_ind=[];div=[];
pol=[];pp=[]; f=[];rho=[]; posx=[]; uv=[]; mv=[]; m=[]; theta=[];
good=[];subdir_ind= []; fov_ind = []; 
flow=[]; flowtime=[]; av_flow_baseline=[];
av_flowtime=[];av_rho=[];
npo=[];pol_sign_tot=[];
shear_um=1700;
cc=1;
px2mu = 0.28;
for jj=1:numel(subdir)
    cd(subdir{jj}); 
    if exist('Res3.mat'); load('Res3.mat');
    else load('Res2.mat');
    end
    pol_sign=1;   %%% sign that accounts for the direction of the flow respective to camera 

    allfolders= regexp(Res.insert_name,'/','split'); 
    if   isempty(allfolders{end})
        allfolders=allfolders(1:end-1);
    end
    flow_string=allfolders{end-1};
    if strcmp(allfolders{end-2},'beads'); time_string=allfolders{end-3};
        else time_string=allfolders{end-2};
    end

   
     N_fov=numel(Res.fov);
%     if jj==1
%     N_fov=numel(Res.fov)-1;
%     end
    
   for fov=1:N_fov;   %starting the cycle over all the field of view

%        if fov==1
%          cd(Res.fov(fov).file_mat(1:end- numel(Res.fov(fov).filename)));
%          fl_var=load(Res.fov(fov).filename); fps=fl_var.mo.FrameRate;
%        end

    clc;
    disp(fov)
       
   % this is to fix a problem with the validation method    
   flowtime_temp=3; %%%%% days of flow

   original_path_PIV =Res.fov(fov).file_mat(1:end- numel(Res.fov(fov).filename));
   pattern= 'D:/np451/';
   if strfind(original_path_PIV,pattern)==1;
       path_PIV = strcat('/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/',original_path_PIV(numel(pattern):end));
       cd(path_PIV);
   else cd(original_path_PIV);
   end

   fl_var=load(Res.fov(fov).filename); fps=fl_var.mo.FrameRate;

    ulim= [-50,50];
    vlim= [-50,50];
 

%    X=fl_var.X;Y=fl_var.Y;U=fl_var.U;V=fl_var.V;
%    [U1,V1]= PIV_Validation(X,Y,U,V,ulim,vlim);
%    [x,y,u1,v1] =  PIV_ChangeFormat(X,Y,U1,V1);   

%----- filter data with no validation --------------------
    u=fl_var.u;v=fl_var.v; x=fl_var.x;y_90=fl_var.y;
    ind_out = u<ulim(1) | u>ulim(2) | v<vlim(1) | v>vlim(2);
    u(ind_out)=nan;v(ind_out)=nan;
    u1=u;v1=v;
    u1m=nanmean(u1,3);     
    v1m=nanmean(v1,3);
%    [vortex,info]=IDvortex(x(1,:),y(:,1),u1m,v1m);
%    nvort=cat(1,nvort,numel(info.Loc))
%---------data with validation-------------------------------   
%   u1m=Res.fov(fov).u1m;
%   v1m=Res.fov(fov).v1m;


%-----------calculate variables ----------------------------------

   M1m = sqrt(u1m.^2+ v1m.^2);
   
   
%   ind = Res.fov(fov).ind & M1m(:)'*0.28*fps >2;  %%% only the good one selected from the field of view through frequency filter
   ind = Res.fov(fov).ind;  %%% only the good one selected from the field of view through frequency filter
   nv= v1m(ind)./M1m(ind);   % normalised vector in the direction of the external flow.
   nu= u1m(ind)./M1m(ind);   % normalised vector in the direction perpendicular the external flow.
   IF= Res.fov(fov).freq.IF;
%%%%% do a new plot with the not validated arrows

%==========================================================================   
% ----------- make plot of the field of view and save it----------------
% 
%     ind2 = ind% & M1m(:)'*0.28*fps >2;
%     figure(5);
% %    imshow(ss);hold on;
%     imagesc(IF); colormap(gray);colorbar();
%     hold on; axis image
%     quiver(x(ind2),y(ind2),u1m(ind2),v1m(ind2),'r','LineWidth',2);
% %    quiverc(x(ind),y(ind),u1m(ind),v1m(ind));hold on;
% 
%     xlabel('X [$\mu$m]','Interpreter','latex');
%     ylabel('Y[$\mu$m]','Interpreter','latex');
%     title(strcat('Not validated X',num2str(Res.fov(fov).posx),' Y',num2str(Res.fov(fov).posy),...
%         ' Pol=',num2str(nanmean(nv),2),' good= ',num2str(Res.fov(fov).good)),'Interpreter','latex');
%     
%     x0=0;y0=0;width=1000;height=500;
%     set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
%     saveas(gcf,strcat(Res.fov(fov).filename(1:end-4),'no validated_PIV_std.png'));
%     close(5)   
%    
%=========================================================================   



% here depending on the experiment we assign a flow, there is some mess becasue we want to use controls from other directories    

      if strcmp(time_string,'18.6.18')  & strcmp(flow_string,'0')
          shear_temp= shear_calculate('tapered',40,shear_um*(floor(1*(shear_um-1)/shear_um)+1)); Res.flow= 0; div_temp=21; seg =1;
 
      elseif strcmp(time_string,'13.6.18')  & strcmp(flow_string,'0')
          shear_temp= shear_calculate('tapered',40,shear_um*(floor(2*(shear_um-1)/shear_um)+1)); Res.flow= 0; div_temp=15; seg =2;

      elseif strcmp(time_string,'16.7.18')  & strcmp(flow_string,'0')
         shear_temp= shear_calculate('tapered',40,shear_um*(floor(3*(shear_um-1)/shear_um)+1)); Res.flow= 0;  div_temp=15; seg =3;
 
      elseif strcmp(time_string,'29.7.18')  & strcmp(flow_string,'0')
        shear_temp= shear_calculate('tapered',40,shear_um*(floor(4*(shear_um-1)/shear_um)+1)); Res.flow= -10;   div_temp=15; seg =4;
 
      elseif strcmp(flow_string,'100') | strcmp(flow_string,'200') 
         shear_temp= shear_calculate('straight',10);  Res.flow= -1; pol_sign=-1;  div_temp=15; seg = -1;

      elseif strcmp(flow_string,'5A') | strcmp(flow_string,'5B') | strcmp(flow_string,'5C')
% this one the name is maybe mistaken      shear_temp= shear_calculate('straight',5);  Res.flow= -1;  % flow = -1 for straight channels
         shear_temp= shear_calculate('straight',10);  Res.flow= -1;  div_temp=15; seg = -1;% flow = -1 for straight channels

         
%%%%% experiment october         
         
      elseif strcmp(time_string,'11.10.19')
          shear_temp= shear_calculate('tapered',8,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 8; div_temp=8; seg= floor(Res.fov(fov).posx/shear_um);  %%% maybe resFLow=5
      
      elseif strcmp(time_string,'11.10.19') & strcmp(flow_string,'control')
          shear_temp= shear_calculate('tapered',8,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 0; div_temp=8; seg= floor(Res.fov(fov).posx/shear_um); %%% maybe resFLow=5

      elseif strcmp(time_string,'18.10.19')
          shear_temp= shear_calculate('tapered',8,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 8; div_temp=15; seg= floor(Res.fov(fov).posx/shear_um);
      elseif strcmp(time_string,'18.10.19') & strcmp(flow_string,'control')
          shear_temp= shear_calculate('tapered',8,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 0; div_temp=15; seg= floor(Res.fov(fov).posx/shear_um);
 
      elseif strcmp(time_string,'25.10.19')
          shear_temp= shear_calculate('tapered',8,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 8;  div_temp=21; seg= floor(Res.fov(fov).posx/shear_um);
   elseif strcmp(time_string,'25.10.19') & strcmp(flow_string,'control')
          shear_temp= shear_calculate('tapered',8,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 0;  div_temp=21; seg= floor(Res.fov(fov).posx/shear_um);

      elseif strcmp(time_string,'22.10.19') & ~strcmp(flow_string,'control')  %%% this is the experiment where flow was only one day
          shear_temp= shear_calculate('tapered',40,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 40; flowtime_temp = 1 ;  div_temp=15; seg= floor(Res.fov(fov).posx/shear_um);
      elseif strcmp(time_string,'22.10.19') & strcmp(flow_string,'control')  %%% this is the experiment where flow was only one day
          shear_temp= shear_calculate('tapered',40,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 0; flowtime_temp = 3 ;  div_temp=28; seg= floor(Res.fov(fov).posx/shear_um);
     

%%%% experiment March      
      
      elseif strcmp(time_string,'22.3.19') 
            shear_temp= shear_calculate('straight',40);  Res.flow= 40;  div_temp=8;  seg=-1  
      elseif strcmp(time_string,'22.3.19') & strcmp(flow_string,'control')
            shear_temp= shear_calculate('straight',40);  Res.flow= 0;  div_temp=8;    seg =-1

            
      elseif strcmp(time_string,'29.3.19') & strcmp(flow_string,'control')
         shear_temp= shear_calculate('tapered',40,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 0;  div_temp=15; seg= floor(Res.fov(fov).posx/shear_um);
    elseif strcmp(time_string,'29.3.19') 
         shear_temp= shear_calculate('tapered',40,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 40;  div_temp=15; seg= floor(Res.fov(fov).posx/shear_um);

         
      elseif strcmp(time_string,'12.4.19')
         shear_temp= shear_calculate('tapered',40,shear_um*(floor(Res.fov(fov).posx/shear_um))); Res.flow= 40;  div_temp=21; seg= floor(Res.fov(fov).posx/shear_um);

      elseif strcmp(time_string,'5.4.19') 
         shear_temp= shear_calculate('straight',40);  Res.flow= 40;  div_temp=28; seg =-1;%% here Res.flow should -1
      elseif strcmp(time_string,'5.4.19') & strcmp(flow_string,'control')
         shear_temp= shear_calculate('straight',40);  Res.flow= 0;  div_temp=28; seg =-1;%% here Res.flow should -1
 

         
       end
      
      
      if strcmp(flow_string,'control') | strcmp(flow_string,'0'); Res.flow= 0; end
         
      
      subdir_ind= cat(1,subdir_ind, jj*ones(size(nv(:))));
      div=cat(1,div,div_temp*ones(size(nv(:))));
      fov_ind= cat(1,fov_ind , fov*ones(size(nv(:))));
      shear_ind=cat(1,shear_ind, shear_temp*ones(size(nv(:))));
      
      %shear_ind=cat(1,shear_ind,shear_temp);
      posx_temp= Res.fov(fov).posx;          
      dummy=posx_temp*ones(size(nv(:)));
      posx=cat(1,posx,dummy(:));

      flow=cat(1,flow,Res.flow*ones(size(nv(:)))); %%%% the flow in rpm , need to distinguish control from exp,  is -1 for straight channel

      flowtime=cat(1,flowtime,flowtime_temp*ones(size(nv(:))));  %%%% time of applied flow
     
   
      if  isa(holes{jj},'char'); if strcmp(holes{jj},'none');good=cat(1,good,1*ones(size(nv(:))));end
      elseif any(Res.fov(fov).posx == holes{jj}); good=cat(1,good,0*ones(size(nv(:)))); 
      else good=cat(1,good,1*ones(size(nv(:)))); 
      end


% in general I always applied on the second experiment the flow in the other direction    

    if strcmp(flow_string,'8rpm_2') | strcmp(flow_string,'5rpm_2') | strcmp(flow_string,'40rpm2') | strcmp(flow_string,'40rpm_2');
        pol_sign=-1;
    end
    pol_sign_tot(cc)=pol_sign; 
      

      mv_temp=sqrt(v1m(ind).^2);
      m_temp = sqrt(v1m(ind).^2 + u1m(ind).^2);
      
      m_temp = m_temp*fps*px2mu;   %%% velocity in um/s
      mv_temp=mv_temp*fps*px2mu;  %%% velocity in um/s
      
      mv=cat(1,mv, mv_temp(:));     
      m=cat(1,m, m_temp(:));
      
         
   %%%%%% velocity along polarisation axis
       uv_temp=v1m(ind); %%%%% this is the velocity and you need to ind them
       uv_temp=pol_sign*uv_temp*fps*px2mu;  %%% velocity in um/s
       uv=cat(1,uv, uv_temp(:));
      
    %%%%% polarisation and frequency 
      pp= cat(1,pp,pol_sign*nv(:)); 
      theta= cat(1,theta,atan2(pol_sign*nv(:),nu(:)));
      F=Res.fov(fov).freq.F32(ind);
      f=cat(1,f,F(:));

    %%%%% for nematic order parameter  
%       s_nematic=  
    
     %%%%% in order to find the local density 
      
      fXX=Res.fov(fov).x(ind);
      fYY=Res.fov(fov).y(ind);
      F32= Res.fov(fov).freq.F32(ind);   
    d_lim=70;
    near= find_nearcilia(fXX(:),fYY(:),d_lim);
%    near_pic= zeros(size(Res.fov(fov).freq.F32));
%    near_pic(ind)=near;
    [near2,df] = find_nearcilia_df(fXX(:),fYY(:),F32,d_lim);
    rho= cat(1,rho, near(:)./(pi*(d_lim/32)^2) );
      
      
%    time_temp=ones([numel(fXX),1])*time_temp;
%   time=cat(1,time,time_temp);     
%      pol= cat(1,pol,pol_sign*(Res.fov(fov).Polx));
%      pol(cc)= pol_sign*(Res.fov(fov).Polx);
    average_density=  numel(nv(:))/(size(u,1)*size(u,2));
    if average_density <0.01 ; good(end)=0; end

%--------- load average value over the field of view -----------%      
      if good(end)==1 & flowtime(end)==3;  %%% only the one with flow for 3 days and good
         av_shear(cc) = shear_ind(end);
         av_posx(cc)  = posx(end);
         av_seg(cc) = seg;
         av_pol(cc)=    nanmedian(pol_sign*nv(:));
         av_flow(cc)=   flow(end);
         av_rho(cc)=        nanmedian(rho);
         av_flow_baseline(cc)= nanmedian(M1m(~ind));
         av_pol_rho(cc)= (nanmedian(pol_sign*nv(near<median(near(:))))...
             - nanmedian(pol_sign*nv(near>median(near(:)))))/av_pol(cc);
         av_div(cc)=div_temp;
         av_n(cc) = numel(nv(:))/(size(u,1)*size(u,2)) ;
         cc=cc+1;
         
      end


   end
end 
    
    good=logical(good);
    Pol=pp(good);
    Rho=rho(good);
    Posx= posx(good);
    Div= div(good);
    Flowtime= flowtime(good);
    Uv=uv(good);
    Mv=mv(good);
    M=m(good);
%    Npo=npo(good);
%    Time=time(good);
    Flow= flow(good);
    Shear_ind=shear_ind(good);
    
    Subdir_ind= subdir_ind(good);
    Fov_ind= fov_ind(good);
    Theta=theta(good);

    

%--------- calculate average value over all the field of views------------

    [G,idshear,idflow,iddiv]=findgroups(av_shear,av_flow,(av_div>14)*15);
    POL=accumarray(G',av_pol',[],@median); 
    emPOL=accumarray(G',av_pol',[],@std)./sqrt(accumarray(G',av_pol',[],@numel)) ;
    ePOL=accumarray(G',av_pol',[],@std) ;
    RHO_pol = accumarray(G',av_pol_rho',[],@median);
    emRHO_pol=accumarray(G',av_pol_rho',[],@std)%./sqrt(accumarray(G',av_pol_rho',[],@numel)) ;
    RHO = accumarray(G',av_rho',[],@median);
    Nfov=   accumarray(G',av_n',[],@median)
    

%-------------------------------------------------------------------------
%    [G,idshear,idflow]=findgroups(Shear_ind,Flow);
%    POL=accumarray(G,Pol,[],@median); 
%    ePOL=accumarray(G,Pol,[],@std)./sqrt(accumarray(G,Pol,[],@numel)) ;

%    NPO=accumarray(G,Npo,[],@median); 
%    eNPO=accumarray(G,Npo,[],@std)./sqrt(accumarray(G,Npo,[],@numel)) ;

    
  ind_c= (idflow==0) & (iddiv>14);
  ind_tap = (idflow>1) & (iddiv>14); 
  ind_str = (idflow==-1) & (iddiv>14)%(idflow~=40) & (idflow~=0); 
  ind_tap_div= (idflow>1) & (iddiv<14)
  ms=10;lw=2;
  
    shear_tap= idshear(ind_tap);
    shear_c= idshear(ind_c);
    shear_str= idshear(ind_str);
    shear_tap_div= idshear(ind_tap_div);

save('all_results_1may2020.mat')



