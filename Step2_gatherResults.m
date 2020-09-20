

clear all
% load the essential functions store in the directory Matlab_essentials
addpath(genpath('/home/np451/Documents'));

ppm= 0.14*2;  %%% 20X

path_dir = '/home/np451/Desktop/dataset alignment';
cd(path_dir); 


% folders with the analysis
subdir={...   
    'DAY15-highshear/FL/','DAY15-control/FL/'
}; 
subdir=strcat('PIV/',subdir)


posrange={...
    'all','all'
    }

% folders with the BF images
BFdir={... 
    'DAY15-highshear/FL/','DAY15-control/FL/'
};

%%
for jj=1:numel(subdir)
    
start=0;
cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');

clear xx;clear yy;clear uu;clear vv;clear uu1;clear vv1;
clear XX;clear YY;clear UU;clear VV;clear UU1;clear VV1;

clear Res;clear FF;


% all the results from an insert (or cell culture is stored in Res)
Res.insert_name= strcat(path_dir,'/',subdir{jj});
Res.posrange=posrange{jj};

%%
for ii=1:size(d,1)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    %% start loop for over the fov in the directory
    filename= d(ii).name;

    cd(path_dir);cd((subdir{jj}));
    fl_var= load(filename);
    u= fl_var.u;v= fl_var.v;u1= fl_var.u1;v1= fl_var.v1; x= fl_var.x;y= fl_var.y;
    um=nanmean(u,3); vm=nanmean(v,3); u1m=nanmean(u1,3); v1m=nanmean(v1,3);
    fl_var.mo;

    % get the position of the stage from the name of the file to find the
    % corresponding BF image
    
    indxL=strfind(filename,'_X')+2;
    posx=[];cc=0; while ~isempty(str2num(filename(indxL+cc))) |filename(indxL+cc)=='-' ; 
    posx=strcat(posx,filename(indxL+cc));cc=cc+1;end;posx=str2num(posx);
    indyL=strfind(filename,'Y')+1;
    posy=[];cc=0; while  ~isempty(str2num(filename(indyL+cc)))| filename(indyL+cc)=='-'; 
    posy=strcat(posy,filename(indyL+cc));cc=cc+1;end;cc=0;posy=str2num(posy);

    if  isa(Res.posrange,'char'); if Res.posrange=='all';Res.fov(ii).good=1;end
    elseif posx>=Res.posrange(1) & posx<=Res.posrange(2); Res.fov(ii).good=1;
    else  Res.fov(ii).good=0;   
    end
    Res.fov(ii).posx=posx;Res.fov(ii).posy=posy; Res.fov(ii).filename=filename;
    Res.fov(ii).file_mat= strcat(path_dir,'/',subdir{jj},'/',filename);

    
    
    
    %% load the std from BF videos and calculate mask and frequency map
  

    cd(path_dir);cd(BFdir{jj});
    d_BF=dir(strcat('*X',num2str(posx),'*Y',num2str(posy),'*.movie') );
    mo=moviereader(d_BF(1).name);
    f_lim=[16,35];prob_cilia=0.8;box_size=4;area_min=10;
    
    % first it find the position of the cilia based on the fact that they move within acertain range of frequency 
    % find the beating frequency of these cilia and make a map of it.
    
    [F_temp,s,BW_temp] = find_cilia(mo,f_lim,prob_cilia,box_size);
    [F4,BW] = remove_debris(F_temp,area_min,f_lim);
    ss=imresize(BW,box_size,'nearest');    
    IF= imresize(F4,box_size,'nearest');
    
    if size(ss,2)<fl_var.mo.width; dummy=zeros([fl_var.mo.height,fl_var.mo.width]); 
        dummy_fs=zeros([fl_var.mo.height,fl_var.mo.width,size(fs,3)],'uint8');
        dummy_fs(:,end-size(ss,2)+1:end,:)=fs; fs=dummy_fs;
        dummy(:,end-size(ss,2)+1:end)= ss; ss=logical(dummy);
    end; 
    
    
    % average the beating frequency inside the box of the PIV (from IF to F32)
    
    box_ss=[];F32=[];bs=16; n_good_px= 0.2;
    for b=1:numel(x);
        y_box=floor(y(b)-bs);
        x_box=floor(x(b)-bs);
        temp_ss =ss(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
        temp_f32=IF(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
        box_ss(b)=mean(temp_ss(:));
        if box_ss(b)>n_good_px
            F32(b)=nanmedian(temp_f32(:));
        else F32(b)=nan;
        end
    end
    ind= box_ss>n_good_px;
    
    cd(path_dir);cd((subdir{jj}));
    F32=reshape(F32,size(x));
    IF(isnan(IF))=0;
    
     Res.fov(ii).freq.F4=F4; Res.fov(ii).freq.BW =BW; Res.fov(ii).freq.IF=IF;
     Res.fov(ii).freq.F_temp=F_temp; Res.fov(ii).freq.BW_temp =BW_temp;
     Res.fov(ii).freq.area_min=area_min;Res.fov(ii).freq.f_lim=f_lim;
     Res.fov(ii).freq.F32=F32;
     
    %% find polarisation for the field of view
    
    M= sqrt(u1m(ind).^2 +v1m(ind).^2);Res.fov(ii).Mm=median(M);
    nu= u1m(ind)./M; nv= v1m(ind)./M;
    pol=sqrt(nanmean(nu)^2+ nanmean(nv)^2); polx= nanmean(nv); 
    Res.fov(ii).Pol= pol;   Res.fov(ii).x=x;   Res.fov(ii).y=y; Res.fov(ii).u1m=u1m; Res.fov(ii).v1m=v1m;
    Res.fov(ii).ind=ind;    Res.fov(ii).nu=nu; Res.fov(ii).nv=nv; Res.fov(ii).Polx= polx;

    %% make plot of the field of view with both frequency and beating direction    
    figure(5);
%    imshow(ss);hold on;
    imagesc(IF); colormap(gray);colorbar();
    hold on; axis image
    quiver(x(ind),y(ind),u1m(ind),v1m(ind),'r','LineWidth',2);
    xlabel('X [$\mu$m]','Interpreter','latex');
    ylabel('Y[$\mu$m]','Interpreter','latex');
    title(strcat('X',num2str(posx),' Y',num2str(posy),' good= ',num2str(Res.fov(ii).good),' Pol=',num2str(polx,2)),'Interpreter','latex');
    
    x0=0;y0=0;width=1000;height=500;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    saveas(gcf,strcat(filename,'_PIV_std.png'));
    close(5)

  
    %% find orientation correlation fov (not used in the paper)

    bin_res= 32; 
    [idr,idth,cc,ecc,n_cc] = corr_orientation_theta_func(x(ind),y(ind),nu,nv,bin_res);
    Res.fov(ii).idr =idr; Res.fov(ii).idth =idth; Res.fov(ii).cc=cc; Res.fov(ii).ecc=ecc; Res.fov(ii).n_cc=n_cc;
    
    %%%%% fit correlation fucntion with exponential to find tipic lenght 
    if numel(idr)>6
    
    ft = fittype('a*exp(-x/b)+c','independent','x');
    fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,idr(2),-1],...
    'Upper',[1,Inf,1],...
    'StartPoint',[0.9, idr(2)*10, 0]);
 
    fit_bins= idr; fit_cc=cc;max_cc=floor(numel(cc)/2);
    fit_out = fit(fit_bins(1:max_cc),fit_cc(1:max_cc),ft,fo);
    Res.fov(ii).fit_out=fit_out;
   
    figure(4)
     plot(fit_out,idr,cc,'ko');hold on;
    else
     figure(4)
     plot(idr,cc,'ko');hold on;  
     Res.fov(ii).fit_out=nan;
    end   
     legend({'freq_corr','orientation_corr'});
     xlabel('pixel [$\mu$m]','Interpreter','latex');
     ylabel(' correlation','Interpreter','latex');
     title('Orientation correlation function','Interpreter','latex');
     
    x0=0;y0=0;width=1000;height=500;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    saveas(gcf,strcat(filename,'_correlations.png'));
    close(4)

        
    %% store the information on this fov 

    xx(start+1:start+numel(x(ind))) = x(ind) +posy/ppm;
    yy(start+1:start+numel(x(ind))) = y(ind) -posx/ppm;

    
    uu(start+1:start+numel(x(ind))) = um(ind);
    vv(start+1:start+numel(x(ind))) = vm(ind);
    uu1(start+1:start+numel(x(ind))) = u1m(ind);
    vv1(start+1:start+numel(x(ind))) = v1m(ind);
    FF(start+1:start+numel(x(ind))) = F32(ind);  %%%%% load all the frequencies

%%%%%%%copiato
       
    start= start+ numel(x(ind));
end
Res.xx=xx;Res.yy=yy;Res.uu1=uu1;Res.vv1=vv1; Res.FF=FF;

XX= xx; YY= yy; UU= uu; VV= vv; UU1= uu1; VV1= vv1;


% find average properties of all the culture

M= sqrt(UU1.^2 +VV1.^2);nu= UU1./M;nv= VV1./M;
Res.nu=nu;Res.nv=nv;Res.M=M;
Mm=median(M);
pol=sqrt(nanmean(nu)^2+ nanmean(nv)^2);
polx= nanmean(nu); 
UM= nanmean(nu);
VM= nanmean(nv);
Res.Pol=pol;
Res.Polx=polx;


figure(2);
subplot(2,2,1)
quiver(XX,YY,UU1,VV1);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern with validation');
axis image

subplot(2,2,2)
quiver(XX,YY,nu,nv,.6);
xlabel('X [um]');
ylabel(' Y[um]');
title('Flow pattern validation and normalisation');
axis image


subplot(2,2,4)
rose(angle(complex(nu,nv)),18); 
title(strcat('Orientation distribution;  mean pol=',num2str(polx))); 

x0=0;y0=0;width=1000;height=1000;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',15);


    saveas(gcf,'flow_pattern.png');
    close all;
   

    save('Res3.mat','Res');

    cd(path_dir);
end


