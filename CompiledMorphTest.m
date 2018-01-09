%% Compiled commands use to achieve morphing

% Figure & Initialize
do_trig=0;
h = figure(2); clf;
whitebg(h,[0 0 0]);

img = (imread('project2_testimg.png'));

% Control points
p1 = [1 1; 257 1; 1 257; 257 257; 129 129];
p2(1) = {[1 1; 257 1; 1 257; 257 257; 129 33]};
p2(2) = {[1 1; 257 1; 1 257; 257 257; 33 129]};
p2(3) = {[1 1; 257 1; 1 257; 257 257; 129 223]};
p2(4) = {[1 1; 257 1; 1 257; 257 257; 223 129]};
p2(5) = {cell2mat(p2(1))};

%VIDEO STUFF
if do_trig
    fname = 'MorphingResult_Eval.avi';
else
    fname = 'MorphingResult_Eval.avi';
end

try
    % VideoWriter based video creation
    h_avi = VideoWriter(fname, 'Uncompressed AVI');
    h_avi.FrameRate = 10000;
    h_avi.open();
catch
    % Fallback deprecated avifile based video creation
    h_avi = avifile(fname,'fps',10000);
end


for j=1:4
    im_source = img;
    im_target = img;
    
    %[selectedMovingPoints,selectedFixedPoints] = cpselect(im_source,im_target,'Wait',true)
    source_pts=p1;
    target_pts=cell2mat(p2(j));

    %Initialize intermediate output
    im_intermed = zeros(size(im_source));

    %triangulation for mean shape
    MeanShape = (0.5*source_pts)+(0.5*target_pts);
    TRI = delaunay(MeanShape(:,1),MeanShape(:,2));

    % number of triangles
    TriangleNum = size(TRI,1);  
    
%     %See figures to be morphed
%     figure; imagesc(im_source)
%     figure; imagesc(im_target)

    %ADD LOOP FOR WARP_FRAC AND DISSOLVE_FRAC
    for warp_frac=linspace(0,1,60)
        dissolve_frac=0;

        % new intermediate shape according to warp_frac
        InterPts = ((1-warp_frac).*source_pts)+(warp_frac.*target_pts); 

        % find coordinates in images source and target
        CordInSource = zeros(3,3,TriangleNum);
        CordInTarget = zeros(3,3,TriangleNum);
        CordInIntermed = zeros(3,3,TriangleNum);

        for i =1:TriangleNum
          for j=1:3
            CordInSource(:,j,i) = [source_pts(TRI(i,j),:)'; 1]; %coordinates of each triangle vertex in source image
            CordInTarget(:,j,i) = [target_pts(TRI(i,j),:)'; 1]; %coordinates of each triangle vertex in target image
            CordInIntermed(:,j,i) = [InterPts(TRI(i,j),:)'; 1]; %coordinates of each triangle vertex in intermediate image
          end
        end

        % create a grid for the morphed image
        [x,y] = meshgrid(1:size(im_intermed,2),1:size(im_intermed,1));
        XYpix=[x(:) y(:)];

        % for each element of the grid of the morphed image, find  which triangle it falls in
        TM = tsearchn([InterPts(:,1) InterPts(:,2)],TRI,XYpix);

        % YOUR CODE STARTS HERE

        PixInTRIAll=[];
        SourceWarpAll=[];
        TargetWarpAll=[];

        for i=1:TriangleNum
                PixInTRI= XYpix(find(TM==i),:); %finds all pixels in triangle
                PixInTRIAll=vertcat(PixInTRIAll,PixInTRI); %save subscripts of pixels in each triangle
                IndPixInTRI=(horzcat(PixInTRI, ones(size(PixInTRI,1),1)))'; %indices of all pixels in triangle
                bary=CordInIntermed(:,:,i)\IndPixInTRI; %calculate barycentric coord for each point in triangle
                SourceWarp=round(CordInSource(:,:,i)*bary);
                SourceWarpAll=horzcat(SourceWarpAll,SourceWarp);
                TargetWarp=round(CordInTarget(:,:,i)*bary);
                TargetWarpAll=horzcat(TargetWarpAll,TargetWarp);
        end

        SourceWarpAll(find(SourceWarpAll > size(im_source,2)))=[size(im_source,2)];
        %SourceWarpAll(find(SourceWarpAll(1,:) > size(im_source,1)),:)=[size(im_source,1)];
        SourceWarpAll(SourceWarpAll == 0)=[1];

        TargetWarpAll(find(TargetWarpAll > size(im_target,2)))=[size(im_target,2)];
        %SourceWarpAll(find(SourceWarpAll(1,:) > size(im_source,1)),:)=[size(im_source,1)];
        TargetWarpAll(TargetWarpAll == 0)=[1];

        IntermedInd=(sub2ind([size(im_intermed,1) size(im_intermed,2)],PixInTRIAll(:,2),PixInTRIAll(:,1)))';
        SourceWarpInd=sub2ind([size(im_source,1) size(im_source,2)],SourceWarpAll(2,:),SourceWarpAll(1,:));
        TargetWarpInd=sub2ind([size(im_target,1) size(im_target,2)],TargetWarpAll(2,:),TargetWarpAll(1,:)); 

        %repeat for target
        rs = im_source(:,:,1);
        gs = im_source(:,:,2);
        bs = im_source(:,:,3);

        im_intermed_rs=zeros(size(im_intermed,1),size(im_intermed,2));
        im_intermed_gs=zeros(size(im_intermed,1),size(im_intermed,2));
        im_intermed_bs=zeros(size(im_intermed,1),size(im_intermed,2));

        im_intermed_rs(IntermedInd)=rs(SourceWarpInd);
        im_intermed_gs(IntermedInd)=gs(SourceWarpInd);
        im_intermed_bs(IntermedInd)=bs(SourceWarpInd);

        rt = im_target(:,:,1);
        gt = im_target(:,:,2);
        bt = im_target(:,:,3);

        im_intermed_rt=zeros(size(im_intermed,1),size(im_intermed,2));
        im_intermed_gt=zeros(size(im_intermed,1),size(im_intermed,2));
        im_intermed_bt=zeros(size(im_intermed,1),size(im_intermed,2));

        im_intermed_rt(IntermedInd)=rt(TargetWarpInd);
        im_intermed_gt(IntermedInd)=gt(TargetWarpInd);
        im_intermed_bt(IntermedInd)=bt(TargetWarpInd);

        SourceWarpInter=uint8(cat(3,im_intermed_rs,im_intermed_gs,im_intermed_bs));
        %figure; imagesc(SourceWarpInter);

        TargetWarpInter=uint8(cat(3,im_intermed_rt,im_intermed_gt,im_intermed_bt));
        %figure; imagesc(TargetWarpInter);

        %imwrite(SourceWarpColor,'source1.jpg')
        morphed_im=(1-dissolve_frac).*SourceWarpInter+ dissolve_frac.* TargetWarpInter;
        imagesc(morphed_im);
            axis image; axis off;drawnow;
            try
                % VideoWriter based video creation
                h_avi.writeVideo(getframe(gcf));
            catch
                % Fallback deprecated avifile based video creation
                h_avi = addframe(h_avi, getframe(gcf));
            end
    end
end

try
    % VideoWriter based video creation
    h_avi.close();
catch
    % Fallback deprecated avifile based video creation
    h_avi = close(h_avi);
end
clear h_avi;