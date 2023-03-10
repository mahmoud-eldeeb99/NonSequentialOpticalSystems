classdef NonSequentialOpticalModel
   properties

   end
   methods
       

%%%%%%%%% <<<<<<<<<<< simple raytracing funcations >>>>>>>>>>>.



       % free space
function raysOut = freeSpace(z,raysIn,L)

M = [1 L;
     0 1];
 raysOut = M*raysIn;
end




% create the rays 
function [rays rayColors]=createRays(z)

numX = 4; % number of spatial points
maxSize = 0.1; % Image size in mm
X = -maxSize/2:maxSize/numX:maxSize/2;
numU = 4; % number of rays from each spatial point
maxAngle = pi/72; % maximum divergence of each group of rays
U = -maxAngle/2:maxAngle/numU:maxAngle/2;
c = [1.0 0.0 0.0;
     0.5 0.5 0.0;
     0.0 1.0 0.0;
     0.0 0.5 0.5;
     0.0 0.0 1.0;
     0.5 0.5 0.5];
rayCount = 1;
for i = 1:numX+1
    for j = 1:numU+1
        rayColors(:,rayCount) = c(i,:); % Rays from each point get the same color
        rays(1,rayCount)=X(i);
        rays(2,rayCount)=U(j);
        rayCount = rayCount + 1;
    end
end


fprintf("Created %d rays of %d angles from %d positions\n", rayCount-1, numU, numX)
end





% draw rays
function [status] = drawRays(z,raysIn,raysOut,LIn,LOut,rayColors)

status = 0;  %Starting drawing
size(raysIn,2)
for i = 1:size(raysIn,2)
    X = [LIn LOut];
    Y = [raysIn(1,i) raysOut(1,i)];
    line(X,Y,'LineStyle','-','Color',rayColors(:,i))
end
status = 1; % Drawing completed
end




% thin lense
function raysOut = thinLens(z,raysIn,f)

M = [1      0;
    -(1/f) 1];
 
 raysOut = M*raysIn;
end







%%%%%%%%>>>>>>>>>>>>> aberations funcations >>>>>>>>>>>>>>>>>>>>>
% thin lense
function raysOut = Ab_thinLens(z,raysIn,f1,f2)

M1 = [1      0;
    (1/f1) 1];

M2 = [1      0;
    -(1/f2) 1];

len=length(raysIn)/2

for i= 1:len+1
    i
    raysOut(:,i)=M1*raysIn(:,i)
    raysOut(:,i+1)=M2*raysIn(:,i+1)

%     raysOut(:,4)=M2*raysIn(:,4)
%     raysOut(:,3)=M1*raysIn(:,3)
end

end




% 
% create the rays 
function [rays rayColors]=Ab_createRays(z)

numX = 0; % number of spatial points
maxSize = 0.1; % Image size in mm
X = -maxSize/2:maxSize/numX:maxSize/2;
numU = 3; % number of rays from each spatial point
maxAngle = pi/72; % maximum divergence of each group of rays
U = -maxAngle/2:maxAngle/numU:maxAngle/2;
c = [1.0 0.0 0.0;
     0.5 0.5 0.0;
     0.0 1.0 0.0;
     0.0 0.5 0.5;
     0.0 0.0 1.0;
     0.5 0.5 0.5];
rayCount = 1;
for i = 1:numX+1
    for j = 1:numU+1
        rayColors(:,rayCount) = c(i,:); % Rays from each point get the same color
        rays(1,rayCount)=X(i);
        rays(2,rayCount)=U(j);
        rayCount = rayCount + 1;
    end
end
end
      
      


%%%%% >>>>>>>>>>>> 4F Imaging system >>>>>>>>>>>>>>
%>> filters 
% H_Single
function H_Single=H_SingleSlit(z,dim,center)
H_Single = zeros(dim,dim);
v_single_slit_width   = 100;   % width pixels
H_SingleSlit_hight  = 10;    %  hight pixels

H_Single((center-H_SingleSlit_hight):(center+H_SingleSlit_hight),(center-v_single_slit_width):(center+v_single_slit_width)) =ones(2*H_SingleSlit_hight+1,2*v_single_slit_width+1);
end



function horizontal_double=horizontal_double_slit(z,dim,center)
horizontal_double = zeros(dim,dim);
horizontal_double_width   = 100;   
horizontal_double_hight  = 20;    
horizontal_double_gap     = 50;   
horizontal_double(((center-horizontal_double_gap/2)-horizontal_double_hight):((center-horizontal_double_gap/2)+horizontal_double_hight),...
    (center-horizontal_double_width):(center+horizontal_double_width)) = ...
                                          ones(2*horizontal_double_hight+1,2*horizontal_double_width+1);
                                      
horizontal_double(((center+horizontal_double_gap/2)-horizontal_double_hight):((center+horizontal_double_gap/2)+horizontal_double_hight),...
    (center-horizontal_double_width):(center+horizontal_double_width)) = ...
                                          ones(2*horizontal_double_hight+1,2*horizontal_double_width+1);
    
end




%
function vertical_single=vertical_single_slit(z,dim,center)
vertical_single = zeros(dim,dim);
v_single_slit_width   = 20;   % pixels
vertical_single_hight  = 100;    % pixels
vertical_single((center-vertical_single_hight):(center+vertical_single_hight),(center-v_single_slit_width):(center+v_single_slit_width)) = ...
                                          ones(2*vertical_single_hight+1,2*v_single_slit_width+1);
 
end



%
function vertical_double=vertical_double_slit(z,dim,center)
vertical_double = zeros(dim,dim);
v_single_slit_width   = 20;   % pixels
vertical_single_slit_hight  = 100;    % pixels
v_spacing     = 60;   % pixels
vertical_double((center-vertical_single_slit_hight):(center+vertical_single_slit_hight),...
    ((center-v_spacing/2)-v_single_slit_width):((center-v_spacing/2)+v_single_slit_width)) = ...
                                          ones(2*vertical_single_slit_hight+1,2*v_single_slit_width+1);
                                      
vertical_double((center-vertical_single_slit_hight):(center+vertical_single_slit_hight),...
    ((center+v_spacing/2)-v_single_slit_width):((center+v_spacing/2)+v_single_slit_width)) = ...
                                          ones(2*vertical_single_slit_hight+1,2*v_single_slit_width+1);

    
end


%
function pinhole=pinhole_filter(z,dim,f)
Xfreq = ((-dim/2):(dim/2-1))*f;
Yfreq = -Xfreq;
[FX,FY] = meshgrid(Xfreq,Yfreq);  %2-D arrays hold fx location and fy location of all points
freq_rad = sqrt(FX.^2 + FY.^2);
maxfreq = (dim/2-1)*f;

cutoff_freq1 = 0.1*maxfreq;
pinhole = double(freq_rad <= cutoff_freq1); % pinhole
        
    
end





      
      
   end
end