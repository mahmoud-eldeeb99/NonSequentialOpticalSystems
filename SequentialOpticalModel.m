classdef SequentialOpticalModel < handle
    properties (Constant)
        ULTRAVIOLET_BORDER_HARD = 2 * pi / 400 * 10^9;
        ULTRAVIOLET_BORDER_SOFT = 2 * pi / 450 * 10^9;
        INFRARED_BORDER_SOFT = 2 * pi / 650 * 10^9;
        INFRARED_BORDER_HARD = 2 * pi / 800 * 10^9;
    end
    properties
        x_left = 0.0;
        x_right = 1.0;
        y_bottom = -0.5;
        y_top = 0.5;
        y_rays_min = 0;
        y_rays_max = 0;

        ray_color_policy = 'red';       % affects colors or rays, can be '%color_name%', 'for-each-ray', 'wavelength-dependent', 'natural', 'random'
        autoscale_enabled = 1;      % boolean, can be 0 (false) or 1 (true)
        autoplot_enabled = 1;        % set 0 for manual plotting

        default_k = 2 * pi / 650 * 10^9;

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


        function status = createRaysFromTemplate(obj, template, n_y, n_angle, max_half_angle)
            if nargin < 5
                max_half_angle = (obj.y_top - obj.y_bottom) / (obj.x_right - obj.x_left) / 3;
                if nargin < 4
                    n_angle = 5;
                    if nargin < 3
                        n_y = 3;
                        if nargin < 2
                            template = 'single';
                        end
                    end
                end
            end

            y_third = (obj.y_top - obj.y_bottom) / 3;
            switch template
                case 'single'
                    obj.setRays([obj.y_bottom + y_third * 3 / 2; 0; 1; 0; obj.default_k]);
                    status = 0;
                    obj.n_visible_rays = 1;
                case 'flat'
                    obj.setRays([obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y);
                        zeros(1, n_y);
                        ones(1, n_y);
                        zeros(1, n_y);
                        obj.default_k * ones(1, n_y)]);
                    obj.n_visible_rays = n_y;
                    status = 0;
                case 'spherical'
                    obj.setRays([obj.y_bottom + y_third * 3 / 2 * ones(1, n_angle);
                        obj.centralizedDistribution(-max_half_angle, max_half_angle, n_angle);
                        ones(1, n_angle);
                        zeros(1, n_angle);
                        obj.default_k * ones(1, n_angle)]);
                    status = 0;
                    obj.n_visible_rays = n_angle;
                case 'combined'
                    obj.setRays([repelem(obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y), 1, n_angle);
                        repmat(obj.centralizedDistribution(-max_half_angle, max_half_angle, n_angle), 1, n_y);
                        ones(1, n_y * n_angle);
                        zeros(1, n_y * n_angle);
                        obj.default_k * ones(1, n_y * n_angle)]);
                    status = 0;
                    obj.n_visible_rays = n_y * n_angle;
                case 'rainbow'
                    obj.setRays([obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y);
                        zeros(1, n_y);
                        ones(1, n_y);
                        zeros(1, n_y);
                        obj.centralizedDistribution(obj.INFRARED_BORDER_SOFT, obj.ULTRAVIOLET_BORDER_SOFT, n_y)]);
                    obj.n_visible_rays = n_y;
                    status = 0;
                case 'visible'
                    obj.setRays([obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y);
                        zeros(1, n_y);
                        ones(1, n_y);
                        zeros(1, n_y);
                        obj.centralizedDistribution(obj.INFRARED_BORDER_HARD, obj.ULTRAVIOLET_BORDER_HARD, n_y)]);
                    obj.n_visible_rays = n_y;
                    status = 0;
                otherwise
                    status = 1;
            end
        end


        function status = createHiddenRays(obj, n_new_hidden_rays_per_existed_ray, policy, max_divergency_angle)
            if nargin < 4
                max_divergency_angle = atan((obj.y_top - obj.y_bottom) / (obj.x_right - obj.x_current));
                if nargin < 4
                    policy = 'random';
                    if nargin < 2
                        n_new_hidden_rays_per_existed_ray = 2;
                    end
                end
            end

            n_half = fix(n_new_hidden_rays_per_existed_ray/2);
            n_old_rays = size(obj.rays, 2);
            new_rays = repmat(obj.rays, 1, 2*n_half);
            switch policy
                case 'uniform'
                    angles_addition = repmat(linspace(0, 1, n_half), n_old_rays, 1);
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

        function status = freeSpace_new(obj, L)
            if nargin < 2
                L = obj.x_right - obj.x_current;
            end

            propagator_first_order = [1, L, 0, 0, 0;
                0, 1, 0, 0, 0;
                0, 0, 1, 0, 0;
                0, 0, 0, 1, 0;
                0, 0, 0, 0, 1];

            phase_change = obj.rays(5, :) .* L .* (1 + obj.rays(2, :).^2 / 2);      % this provides whole Fraunhofer's difraction (first nonzero order)
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


        function status = thinLens_new(obj, f, half_diameter, refraction_coeff)
            if nargin < 4
                refraction_coeff = 0;
                if nargin < 3
                    half_diameter = 1.05 * max(max(obj.rays(1,:)), -min(obj.rays(1,:)));
                end
            end
            if half_diameter <= 0
                max_rays_y = max(max(obj.rays(1,:)), -min(obj.rays(1,:)));
                if max_rays_y == 0
                    max_rays_y = (obj.y_top - obj.y_bottom) / 6;
                end
                half_diameter = 1.05 * max_rays_y;
            end
            propagator = [1, 0, 0, 0, 0;
                -1/f, 1, 0, 0, 0;
                0, 0, 1, 0, 0;
                0, 0, 0, 1, 0;
                0, 0, 0, 0, 1];

            affected_rays = abs(obj.rays(1,:)) <= half_diameter;
            unaffected_rays = ~affected_rays;
            sample_array = [1: size(obj.rays, 2)];
            affected_rays =  sample_array(affected_rays);
            unaffected_rays = sample_array(unaffected_rays);

            if size(refraction_coeff, 2) > 1
                chrom_aberration_coeff = zeros(1, size(obj.rays, 2));
                for iOrder = 1 : (size(refraction_coeff, 2) - 1)
                    chrom_aberration_coeff = chrom_aberration_coeff + refraction_coeff(iOrder + 1) * (2 * pi)^(-2*iOrder) * (obj.rays(5, affected_rays)).^(2*iOrder);
                end
                chrom_aberration_angle_change = - 1 / f / (refraction_coeff(1) - 1) * chrom_aberration_coeff .* obj.rays(1, affected_rays);
            else
                chrom_aberration_angle_change = zeros(1, size(affected_rays, 2));
            end
            calculated_propagator_second_order = [zeros(1, size(affected_rays, 2));
                chrom_aberration_angle_change;
                zeros(3, size(affected_rays, 2))];

            new_rays = zeros(5, size(obj.rays, 2));
            new_rays(:, affected_rays) = propagator * obj.rays(:, affected_rays) + calculated_propagator_second_order;
            new_rays(:, unaffected_rays) = obj.rays(:, unaffected_rays);
            obj.setRays(new_rays);

            if obj.autoplot_enabled
                obj.drawLens(f, half_diameter);
            end
            status = 0;
        end


        function status = obstacle(obj, y_bottom, y_top)
            rays_to_keep = (obj.rays(1,:) < y_bottom) | (y_top < obj.rays(1,:));
            obj.deleteRays(rays_to_keep);

            if obj.autoplot_enabled
                obj.drawWall(y_bottom, y_top);
            end
            status = 0;
        end


%%
% Helpers

        function distribution = centralizedDistribution(obj, left, right, n)
            if n > 1
                if rem(n, 2)
                    n_half = fix(n / 2);
                    distribution = (left + right) / 2 + (right - left) / (2 * n_half) * [-n_half : n_half];
                else
                    distribution = linspace(left, right, n);
                end
            else
                distribution = (right + left) / 2;
            end
        end


        function status = setBorders(obj, borders)
            if nargin > 0
                obj.x_left = borders(1);
                obj.x_right = borders(2);
                obj.y_bottom = borders(3);
                obj.y_top = borders(4);
                axis(borders);
                status = 0;
            else
                status = 1;
            end
        end


        function status = autoScale(obj, absolute_scaling)
            if obj.autoscale_enabled
                y_min = min(obj.rays(1, 1:obj.n_visible_rays));
                y_max = max(obj.rays(1, 1:obj.n_visible_rays));
                if absolute_scaling
                    y_min = min(y_min, obj.y_rays_min);
                    y_max = max(y_max, obj.y_rays_max);
                end
                y_bottom = y_min - 0.2 * (y_max - y_min);
                y_top = y_max + 0.2 * (y_max - y_min);
                if y_top ~= y_bottom
                    obj.setBorders([obj.x_left, obj.x_right, y_bottom, y_top]);
                end
            end
        end

        function status = start(obj)
            obj.x_current = obj.x_left;
            obj.y_rays_min = min(obj.y_bottom + 1/7 * (obj.y_top - obj.y_bottom), min(obj.rays(1, 1:obj.n_visible_rays)));
            obj.y_rays_max = max(obj.y_top - 1/7 * (obj.y_top - obj.y_bottom), max(obj.rays(1, 1:obj.n_visible_rays)));
        end


        function status = setRays(obj, new_rays)
            obj.rays_prev = obj.rays;
            obj.rays = new_rays;
            obj.y_rays_min = min([obj.y_rays_min, new_rays(1, 1:obj.n_visible_rays)]);
            obj.y_rays_max = max([obj.y_rays_max, new_rays(1, 1:obj.n_visible_rays)]);
        end


        function status = deleteRays(obj, rays_to_keep)
            obj.n_visible_rays = sum(rays_to_keep(1:obj.n_visible_rays));
            new_rays = obj.rays(:, rays_to_keep);
            obj.setRays(new_rays);
        end


        function status = calcIntensity(obj, resolution, window)
            if nargin < 3
                window = 1;
                if nargin < 2
                    resolution = max(19, min(151, fix(size(obj.rays, 2)/(2*10^4))));
                end
            end
            obj.intensity = zeros(1,resolution);
            delta = (obj.y_top - obj.y_bottom) / resolution;
            amplitude = exp(i * obj.rays(4,:)) .* obj.rays(3,:);
            rays_y = obj.rays(1,:);
            Y = obj.y_bottom + ([1:resolution] - 0.5) * delta;
            for j = 1:resolution
                y = obj.y_bottom + (j - 0.5) * delta;
                mask = (y - window * delta / 2 < rays_y) & (rays_y <= y + window * delta / 2);
                obj.intensity(j) = abs(sum(amplitude(mask)))^2;
            end
            max_intensity = max(obj.intensity);
            if max_intensity ~= 0
                obj.intensity = 1 / max_intensity * obj.intensity;
            end
            if obj.autoplot_enabled
                obj.drawIntensity(resolution);
            end
            status = 0;
        end


        function color_returned = colorFromWavevector(obj, k, invisible_color)
            switch invisible_color
                case 'black'
                    r0 = 0;
                    g0 = 0;
                    b0 = 0;
                case 'white'
                    r0 = 1;
                    g0 = 1;
                    b0 = 1;
                otherwise
                    r0 = 0;
                    g0 = 0;
                    b0 = 0;
            end
            red = obj.INFRARED_BORDER_SOFT;
            green = (obj.INFRARED_BORDER_SOFT + obj.ULTRAVIOLET_BORDER_SOFT) / 2;
            blue = obj.ULTRAVIOLET_BORDER_SOFT;
            if (k <= obj.INFRARED_BORDER_HARD) || (obj.ULTRAVIOLET_BORDER_HARD <= k)
                r = r0;
                g = g0;
                b = b0;
            elseif (obj.INFRARED_BORDER_HARD < k) && (k <= red)
                r = (k - obj.INFRARED_BORDER_HARD) / (red - obj.INFRARED_BORDER_HARD) + r0 * (red - k) / (red - obj.INFRARED_BORDER_HARD);
                g = g0 * (red - k) / (red - obj.INFRARED_BORDER_HARD);
                b = b0 * (red - k) / (red - obj.INFRARED_BORDER_HARD);
            elseif (red < k) && (k <= green)
                r = (green - k) / (green - red);
                g = (k - red) / (green - red);
                b = 0;
            elseif (green < k) && (k <= blue)
                r = 0;
                g = (blue - k) / (blue - green);
                b = (k - green) / (blue - green);
            elseif (blue < k) && (k < obj.ULTRAVIOLET_BORDER_HARD)
                r = cos(((acos(r0)/pi + 0.5) * (k - blue) / (obj.ULTRAVIOLET_BORDER_HARD - blue) - 0.5) * pi);
                g = g0 * (k - blue) / (obj.ULTRAVIOLET_BORDER_HARD - blue);
                b = (obj.ULTRAVIOLET_BORDER_HARD - k) / (obj.ULTRAVIOLET_BORDER_HARD - blue) + g0 * (k - blue) / (obj.ULTRAVIOLET_BORDER_HARD - blue);
            end
            color_returned = [r; g; b];
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
                case 'natural'
                    for i = 1:n_rays
                        line(X(:,i), Y(:,i), 'Color', obj.colorFromWavevector(obj.rays(5, i), 'white'));
                    end
                case 'wavelength-dependent'
                    for i = 1:n_rays
                        line(X(:,i), Y(:,i), 'Color', obj.colorFromWavevector(obj.rays(5, i), 'black'));
                    end
            otherwise
                line(X, Y, 'Color', obj.ray_color_policy);
            end
            obj.autoScale(1);
        end


        function drawLens(obj, f, half_diameter)
            y_center = (obj.y_top + obj.y_bottom)/2;
            line([obj.x_current; obj.x_current], [y_center - half_diameter; y_center + half_diameter], 'Color', [0, 0, 0], 'LineWidth', 2);
            text(obj.x_current, y_center - half_diameter - 0.05 * (obj.y_top - obj.y_bottom), ['f = ', num2str(f)], 'FontSize', 14);
            obj.autoScale(1);
        end


        function drawWall(obj, y_min, y_max)
            line([obj.x_current; obj.x_current], [y_min; y_max], 'Color', [0, 0, 0], 'LineWidth', 2);
        end


        function drawIntensity(obj, resolution)
            line([obj.x_right, obj.x_right], [obj.y_bottom, obj.y_top], 'Color', [0,0,0]);
            axis([obj.x_left, obj.x_right * 1.1, obj.y_bottom, obj.y_top]);
            X = obj.x_right + obj.intensity * 0.1 * obj.x_right;
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
