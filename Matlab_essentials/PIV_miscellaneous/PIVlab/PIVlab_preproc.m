function out = PIVlab_preproc (in,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize)
if size(in,3)>1
    in(:,:,2:end)=[];
end
%this function preprocesses the images
if numel(roirect)>0
    x=roirect(1);
    y=roirect(2);
    width=roirect(3);
    height=roirect(4);
else
    x=1;
    y=1;
    width=size(in,2)-1;
    height=size(in,1)-1;
end
%roi (x,y,width,height)
in_roi=in(y:y+height,x:x+width);

if intenscap == 1
    %Intensity Capping: a simple method to improve cross-correlation PIV results
    %Uri Shavit � Ryan J. Lowe � Jonah V. Steinbuck
    n = 2; 
    up_lim_im_1 = median(double(in_roi(:))) + n*std2(in_roi); % upper limit for image 1
    brightspots_im_1 = find(in_roi > up_lim_im_1); % bright spots in image 1
    capped_im_1 = in_roi; capped_im_1(brightspots_im_1) = up_lim_im_1; % capped image 1
    in_roi=capped_im_1;
end
if clahe == 1
    numberoftiles1=round(size(in_roi,1)/clahesize);
    numberoftiles2=round(size(in_roi,2)/clahesize);
    if numberoftiles1 < 2
    numberoftiles1=2;
    end
    if numberoftiles2 < 2
    numberoftiles2=2;
    end
    in_roi=adapthisteq(in_roi, 'NumTiles',[numberoftiles1 numberoftiles2], 'ClipLimit', 0.01, 'NBins', 256, 'Range', 'full', 'Distribution', 'uniform');
end

if highp == 1
    h = fspecial('gaussian',highpsize,highpsize);
    in_roi=double(in_roi-(imfilter(in_roi,h,'replicate')));
    in_roi=in_roi/max(max(in_roi))*255;
end

if wienerwurst == 1
    in_roi=wiener2(in_roi,[wienerwurstsize wienerwurstsize]);
end

out=in;
out(y:y+height,x:x+width)=in_roi;
% out=uint8(out); %this is bad as only works properly with images that started off as uint8
out = gray2ind(mat2gray(out),256);