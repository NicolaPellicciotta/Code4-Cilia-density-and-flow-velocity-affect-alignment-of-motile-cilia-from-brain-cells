%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---First step:analyse FL images with PIVlab------%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the essential functions store in the directory Matlab_essentials
addpath(genpath('/home/np451/Desktop/dataset alignment'));


%% Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=256;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=64;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=3;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=128;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=64;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=32;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower

%%% Standard image preprocessing settings
p = cell(8,1);
%Parameter                       %Setting           %Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=0;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denaoise filter, 0 = disable
p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size

%% directories to analise

long={...
    'DAY15-highshear/FL/','DAY15-control/FL/'
}

   

for insert=long; % cycle on each experiment directory, where insert is a single culture 

%% load files
insert=insert{1};
%directory with the raw data
data_dir = strcat('/home/np451/Desktop/dataset alignment/',insert);
%analysis folder
a_folder=strcat('/home/np451/Desktop/dataset alignment/PIV/',insert);
mkdir(a_folder); 

cd(data_dir) 
suffix='.movie';

direc = dir('*FL*.movie');
N_files= size(direc,1);

%%

for i=1:N_files % for on each field of view
    cd(data_dir)
    exp_name = direc(i).name;
    exp_name=exp_name(1:end-6);


        mo=moviereader(direc(i).name);
        movie_path=strcat(data_dir,mo.Filename)
        frame_stack=mo.read();
        frame_stack=(frame_stack)-movmin(frame_stack,40,3);
        max_int= double(max(frame_stack(:)));
        
        
        if max_int>256        
        frame_stack= uint8((double(frame_stack)/max_int)*255);
        else frame_stack= uint8(frame_stack);
        end
%        filter high frequency spatial noise
        for kk=1:size(frame_stack,3); 
            frame_stack(:,:,kk)= wiener2(frame_stack(:,:,kk),[5,5]);
        end
        
        
        %%  PIV

        [X,Y,U,V] = PIV_GetData(frame_stack,s,p);
        [x,y,u,v] =  PIV_ChangeFormat(X,Y,U,V);

        close('all');
        
        cd(a_folder);

        %% PostProcessing, removing no sense velocities

        ulim= [-50,50];
        vlim= [-50,50];

        [U1,V1]= PIV_Validation(X,Y,U,V,ulim,vlim);
        [x,y,u1,v1] =  PIV_ChangeFormat(X,Y,U1,V1);


        %%%%%% making a standard deviation figure
        ss= std(double(frame_stack),[],3);
        figure(1);
        subplot(2,2,1)
        title(strcat(exp_name,'frame 1'));
        imagesc(ss);
      
        %%%%%% look at the data before processing
        subplot(2,2,2)
        title(strcat(exp_name,' PreProcessing'));
        quiver(x,-y,nanmean(u,3),-nanmean(v,3),3);
  
        %%%% look at results after processing 
        subplot(2,2,3)
        title(strcat(exp_name,' PostProcessing'));
        quiver(x,-y,mean(u1,3),-mean(v1,3),3);
        fig=figure(1);
        saveas(fig,strcat(exp_name,'_PostPro.png'));
        close(1); 
  
        clear frame_stack; 
        % save results of the analysis
        save(strcat(exp_name,'.mat'));
        
        
end

end


