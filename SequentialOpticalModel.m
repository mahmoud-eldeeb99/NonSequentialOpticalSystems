classdef SequentialOpticalModel < handle
    properties
        x_left = 0.0;
        x_right = 1.0;
        y_bottom = -0.5;
        y_top = 0.5;
        y_rays_min = 0;
        y_rays_max = 0;

        ray_color_policy = 'red';       % affects colors or rays, can be '%color_name%', 'for-each-ray', 'wavelength-dependent', 'random'
        autoscale_enabled = 1;      % boolean, can be 0 (false) or 1 (true)
        autoplot_enabled = 1;        % set 0 for manual plotting

        default_k = 2*pi/650*10^(9);

        x_current = 0.0;
        x_prev = 0.0;
        n_visible_rays = 1;
        rays = [0; 0; 1; 0; 2*pi/650*10^(9)];       % [y, angle, amplitude, phase, wave vector]. Note: wave amplitude = amplitude * exp(- i * (wave vector) * (delta x))
        rays_prev = [0; 0; 1; 0; 2*pi/650*10^(9)];
        intensity;
    end

    methods

%%
% Constructors

        function obj = SequentialOpticalModel
            axis([obj.x_left, obj.x_right, obj.y_bottom, obj.y_top]);
        end

%%
% Rays creating

        function status = createRays_new(obj, rays_new)
            obj.setRays(rays_new);
            obj.n_visible_rays = size(rays_new, 2);
            status = 0;
        end


        function status = createRaysFromTemplate(obj, template, n_y = 3, n_angle = 5, max_half_angle = pi/180)
            y_third = (obj.y_top - obj.y_bottom) / 3;
            switch template
                case 'single'
                    obj.setRays([obj.y_bottom + y_third * 3 / 2; 0; 1; 0; obj.default_k]);
                    status = 0;
                    obj.n_visible_rays = 1;
                case 'flat'
                    obj.setRays([linspace(obj.y_bottom + y_third, obj.y_top - y_third, n_y);
                        zeros(1, n_y);
                        ones(1, n_y);
                        zeros(1, n_y);
                        obj.default_k * ones(1, n_y)]);
                    obj.n_visible_rays = n_y;
                    status = 0;
                case 'spherical'
                    if fix(n_angle/2) ~= 0
                        if rem(n_angle, 2)
                            angles = max_half_angle / fix(n_angle/2) * [-fix(n_angle/2):fix(n_angle/2)];
                        else
                            angles = linspace(-max_half_angle, max_half_angle, n_angle);
                        end
                    else
                        angles = 0;
                    end
                    obj.setRays([obj.y_bottom + y_third * 3 / 2 * ones(1, n_angle);
                        angles;
                        ones(1, n_angle);
                        zeros(1, n_angle);
                        obj.default_k * ones(1, n_angle)]);
                    status = 0;
                    obj.n_visible_rays = n_angle;
                case 'combined'
                    if fix(n_angle/2) ~= 0
                        if rem(n_angle, 2)
                            angles = max_half_angle / fix(n_angle/2) * [-fix(n_angle/2):fix(n_angle/2)];
                        else
                            angles = linspace(-max_half_angle, max_half_angle, n_angle);
                        end
                    else
                        angles = 0;
                    end
                    obj.setRays([repelem(linspace(obj.y_bottom + y_third, obj.y_top - y_third, n_y), 1, n_angle);
                        repmat(angles, 1, n_y);
                        ones(1, n_y * n_angle);
                        zeros(1, n_y * n_angle);
                        obj.default_k * ones(1, n_y * n_angle)]);
                    status = 0;
                    obj.n_visible_rays = n_y * n_angle;
                otherwise
                    status = 1;
            end
        end


        function status = createHiddenRays(obj, n_new_hidden_rays_per_existed_ray = 2, policy = 'random', max_divergency_angle = atan((obj.y_top - obj.y_bottom) / (obj.x_right - obj.x_current)))
            n_half = fix(n_new_hidden_rays_per_existed_ray/2);
            n_old_rays = size(obj.rays, 2);
            new_rays = repmat(obj.rays, 1, 2*n_half);
            switch policy
                case 'uniform'
                    angles_addition = repmat(linspace(0, 1, n_half), n_old_rays, 1);;
                    angles_addition = max_divergency_angle / 2 * [-angles_addition, angles_addition];
                case 'random_forward'
                    angles_addition = exp(-rand(n_old_rays, n_half).^2*5);
                    angles_addition = max_divergency_angle / 2 * [-angles_addition, angles_addition];
                case 'random'
                    angles_addition = rand(n_old_rays, n_half);
                    angles_addition = max_divergency_angle / 2 * [-angles_addition, angles_addition];
                otherwise
                    status = 1;

                    status = 0;
            end
            amplitude_inv_norm = 1 ./ sqrt(1 + sum(cos(angles_addition).^2, 2));
            amplitude_inv_norm = repmat(transpose(amplitude_inv_norm), 1, 2*   n_half + 1);
            new_rays(2,:) = new_rays(2,:) + transpose(angles_addition(:));
            new_rays(3,:) = new_rays(3,:) .* cos(transpose(angles_addition(:)));
            obj.rays = [obj.rays, new_rays];
            obj.rays(3,:) = obj.rays(3,:) .* amplitude_inv_norm;
        end


%%
% Propagators through media and optical elements

        function status = freeSpace_new(obj, L = obj.x_right - obj.x_current)
            propagator_first_order = [1, L, 0, 0, 0;
                0, 1, 0, 0, 0;
                0, 0, 1, 0, 0;
                0, 0, 0, 1, 0;
                0, 0, 0, 0, 1];

            phase_change = obj.rays(5, :) .* L .* (1 + obj.rays(2, :).^2 / 2);      % this this provides whole Fraunhofer's difraction (first nonzero order)
            phase_change = rem(phase_change, 2*pi);       % renormalization
            calculated_propagator_second_order = [zeros(3, size(obj.rays, 2)); phase_change; zeros(1, size(obj.rays, 2))];        % good news: all other terms of second order with respect to y and angle are zero

            obj.setRays(propagator_first_order * obj.rays + calculated_propagator_second_order);
            obj.x_prev = obj.x_current;
            obj.x_current = obj.x_current + L;
            if obj.autoplot_enabled
                obj.drawRays_new;
            end
            status = 0;
        end


        function status = thinLens_new(obj, f)
            propagator = [1, 0, 0, 0, 0;
                -1/f, 1, 0, 0, 0;
                0, 0, 1, 0, 0;
                0, 0, 0, 1, 0;
                0, 0, 0, 0, 1];

            obj.setRays(propagator * obj.rays);
            if obj.autoplot_enabled
                obj.drawLens(f);
            end
            status = 0;
        end


%%
% Helpers

        function status = setBorders(obj, borders);
            if nargin > 0
                obj.x_left = borders(1);
                obj.x_right = borders(2);
                obj.y_bottom = borders(3);
                obj.y_top = borders(4);
                axis(borders);
                status = 0;
            else
                status = 1;
            endif
        end


        function status = autoScale(obj, absolute_scaling)
            if obj.autoscale_enabled
                y_min = min(obj.rays(1, 1:obj.n_visible_rays));
                y_max = max(obj.rays(1, 1:obj.n_visible_rays));
                if absolute_scaling
                    y_min = min(y_min, obj.y_rays_min);
                    y_max = max(y_max, obj.y_rays_max);
                end
                y_bottom = y_min - (y_max - y_min)/2;
                y_top = y_max + (y_max - y_min)/2;
                if y_top ~= y_bottom
                    obj.setBorders([obj.x_left, obj.x_right, y_bottom, y_top]);
                endif
            end
        end

        function status = start(obj)
            obj.x_current = obj.x_left;
            obj.y_rays_min = min(obj.rays(1, 1:obj.n_visible_rays));
            obj.y_rays_max = max(obj.rays(1, 1:obj.n_visible_rays));
        end


        function status = setRays(obj, new_rays)
            obj.rays_prev = obj.rays;
            obj.rays = new_rays;
            obj.y_rays_min = min([obj.y_rays_min, new_rays(1, 1:obj.n_visible_rays)]);
            obj.y_rays_max = max([obj.y_rays_max, new_rays(1, 1:obj.n_visible_rays)]);
        end


        function status = calcIntensity(obj, window = 1, resolution = max(19, min(151, size(obj.rays, 2)/100)))
            obj.intensity = zeros(1,resolution);
            delta = (obj.y_top - obj.y_bottom) / resolution;
            for j = 1:resolution
                y = obj.y_bottom + (j - 0.5) * delta;
                mask = (y - window * delta / 2 < obj.rays(1,:)) & (obj.rays(1,:) <= y + window * delta / 2);
                obj.intensity(j) = sum(mask .* obj.rays(3,:) .* exp(i * obj.rays(4,:)));
                %obj.intensity(:,j) = sum(([1; 1] * (mask .* obj.rays(3,:))) .* [cos(obj.rays(4,:)); sin(obj.rays(4,:))], 2);
            end
            %obj.intensity = obj.intensity(:,1).^2 + obj.intensity(:,2).^2;
            obj.intensity = abs(obj.intensity).^2;
            max_intensity = max(obj.intensity);
            if max_intensity ~= 0
                obj.intensity = 1 / max_intensity * obj.intensity;
            end
            if obj.autoplot_enabled
                obj.drawIntensity(resolution);
            end
            status = 0;
        end
%%
% Drawing

        function drawRays_new(obj, x_a, x_b, rays_a, rays_b)
            n_rays = obj.n_visible_rays;
            X = [obj.x_prev * ones(1, n_rays);
                obj.x_current * ones(1, n_rays)];
            Y = [obj.rays_prev(1, 1:n_rays);
                obj.rays(1, 1:n_rays)];
            switch obj.ray_color_policy
                case 'for-each-ray'
                    line_color = [linspace(1, 0, n_rays - fix(n_rays/3)), linspace(0, 0, fix(n_rays/3));
                        linspace(0, 0, n_rays - fix(n_rays/3)), linspace(0, 1, fix(n_rays/3));
                        linspace(0, 1, n_rays - fix(n_rays/3)), linspace(1, 0, fix(n_rays/3))];
                    for i = 1:n_rays
                        line(X(:,i), Y(:,i), 'Color', line_color(:,i));
                    end
                case 'random'
                    line_color = [linspace(1, 0, n_rays - fix(n_rays/3)), linspace(0, 0, fix(n_rays/3));
                        linspace(0, 0, n_rays - fix(n_rays/3)), linspace(0, 1, fix(n_rays/3));
                        linspace(0, 1, n_rays - fix(n_rays/3)), linspace(1, 0, fix(n_rays/3))];
                    line_color = line_color(:, randperm(n_rays));
                    for i = 1:n_rays
                        line(X(:,i), Y(:,i), 'Color', line_color(:,i));
                    end
                case 'wavelength-dependent'
                        ;
            otherwise
                line(X, Y, 'Color', obj.ray_color_policy);
            end
        end


        function drawLens(obj, f)
            obj.autoScale(0);
            y_delta = 0.9 * (obj.y_top - obj.y_bottom)/4;
            line([obj.x_current; obj.x_current], [obj.y_bottom + y_delta; obj.y_top - y_delta], 'Color', [0, 0, 0], 'LineWidth', 2);
            text(obj.x_current, obj.y_bottom + y_delta*2/3, ['f = ', num2str(f)], 'FontSize', 14);
            obj.autoScale(1);
        end


        function drawIntensity(obj, resolution)
            line([obj.x_right, obj.x_right], [obj.y_bottom, obj.y_top], 'Color', [0,0,0]);
            axis([obj.x_left, obj.x_right + 1, obj.y_bottom, obj.y_top]);
            X = obj.x_right + obj.intensity;
            Y = obj.y_bottom + (obj.y_top - obj.y_bottom) / resolution * ([1:resolution] - 0.5);
            plot(X, Y);
        end















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Old functions


function raysOut = freeSpace(z, raysIn, L)

  PropagatorM = [1 L;
    0 1];

  raysOut = PropagatorM * raysIn;

end


function raysOut = thinLens(z,raysIn,f)

M = [1      0;
    -(1/f) 1];

 raysOut = M*raysIn;
end


function raysOut = thickLens(z,raysIn,p1,p2,d,n)

M = [1-p2*(d/n)      p1+p2-p1*p2*(d/n);
    -d/n             1-p1*(d/n)];

 raysOut = M*raysIn;
end





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









function [rays, rayColors]=createRays(z)

numX = 2; % number of spatial points
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

fprintf('Created %d rays of %d angles from %d positions\n', rayCount-1, numU, numX)
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
    raysOut(:,1)=M1*raysIn(:,1)
    raysOut(:,2)=M1*raysIn(:,2)

    raysOut(:,4)=M2*raysIn(:,4)
    raysOut(:,3)=M2*raysIn(:,3)
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
