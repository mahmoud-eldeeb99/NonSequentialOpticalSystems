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
        cluster_number;
        intensity;
        intensity_meas_point;
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
            obj.cluster_number = [1 : size(rays_new, 2)];
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
                    obj.cluster_number = 1;
                    status = 0;
                case 'flat'
                    obj.setRays([obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y);
                        zeros(1, n_y);
                        ones(1, n_y);
                        zeros(1, n_y);
                        obj.default_k * ones(1, n_y)]);
                    obj.n_visible_rays = n_y;
                    obj.cluster_number = ones(1, n_y);
                    status = 0;
                case 'spherical'
                    obj.setRays([obj.y_bottom + y_third * 3 / 2 * ones(1, n_angle);
                        obj.centralizedDistribution(-max_half_angle, max_half_angle, n_angle);
                        ones(1, n_angle);
                        zeros(1, n_angle);
                        obj.default_k * ones(1, n_angle)]);
                    obj.n_visible_rays = n_angle;
                    obj.cluster_number = ones(1, n_angle);
                    status = 0;
                case 'combined'
                    obj.setRays([repelem(obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y), 1, n_angle);
                        repmat(obj.centralizedDistribution(-max_half_angle, max_half_angle, n_angle), 1, n_y);
                        ones(1, n_y * n_angle);
                        zeros(1, n_y * n_angle);
                        obj.default_k * ones(1, n_y * n_angle)]);
                    obj.n_visible_rays = n_y * n_angle;
                    obj.cluster_number = repelem([1 : n_y], 1, n_angle);
                    status = 0;
                case 'rainbow'
                    obj.setRays([obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y);
                        zeros(1, n_y);
                        ones(1, n_y);
                        zeros(1, n_y);
                        obj.centralizedDistribution(obj.INFRARED_BORDER_SOFT, obj.ULTRAVIOLET_BORDER_SOFT, n_y)]);
                    obj.n_visible_rays = n_y;
                    obj.cluster_number = [1 : n_y];
                    status = 0;
                case 'visible'
                    obj.setRays([obj.centralizedDistribution(obj.y_bottom + y_third, obj.y_top - y_third, n_y);
                        zeros(1, n_y);
                        ones(1, n_y);
                        zeros(1, n_y);
                        obj.centralizedDistribution(obj.INFRARED_BORDER_HARD, obj.ULTRAVIOLET_BORDER_HARD, n_y)]);
                    obj.n_visible_rays = n_y;
                    obj.cluster_number = [1 : n_y];
                    status = 0;
                otherwise
                    status = 1;
            end
        end


        function status = createHiddenRays(obj, n_new_hidden_rays_per_existed_ray, policy, max_divergency_angle)
            if nargin < 4
                max_divergency_angle = pi / 4;
                %atan((obj.y_top - obj.y_bottom) / (obj.x_right - obj.x_current));
                if nargin < 3
                    policy = 'random';
                    if nargin < 2
                        n_new_hidden_rays_per_existed_ray = 2;
                    end
                end
            end

            n_half = fix(n_new_hidden_rays_per_existed_ray/2);
            n_old_rays = size(obj.rays, 2);
            old_rays = obj.rays;
            delta_x = - 1 / 2 * old_rays(1, :) .* old_rays(2, :);
            old_rays = obj.freeSpacePropagate(old_rays, delta_x);
            delta_x_phi = - old_rays(4, :) / old_rays(5, :) * (1 - old_rays(2, :).^2 / 2);
            old_rays = obj.freeSpacePropagate(old_rays, delta_x_phi);
            delta_x = delta_x + delta_x_phi;
            new_rays = repmat(old_rays, 1, 2*n_half);
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
            amplitude_inv_norm = 1 ./ sqrt(1 + sum(((1 + cos(angles_addition)) / 2).^2, 2));
            %amplitude_inv_norm = 1 ./ sqrt(1 + sum(cos(angles_addition).^2, 2));
            amplitude_inv_norm = repmat(transpose(amplitude_inv_norm), 1, 2*   n_half + 1);
            new_rays(2,:) = new_rays(2,:) + transpose(angles_addition(:));
            new_rays(3,:) = new_rays(3,:) .* (1 + cos(transpose(angles_addition(:)))) / 2;
            %new_rays(3,:) = new_rays(3,:) .* cos(transpose(angles_addition(:)));
            mask = [0: 2*n_half - 1] * n_old_rays;
            for i_ray = 1 : n_old_rays
                %for i_sub_ray = 1 : 2*n_half
                %    new_rays(:, i_ray * i_sub_ray) = obj.freeSpacePropagate(new_rays(:, i_ray * i_sub_ray), -delta_x(i_ray));
                %end
                new_rays(:, i_ray + mask) = obj.freeSpacePropagate(new_rays(:, i_ray + mask), -delta_x(i_ray));
            end
            obj.rays = [obj.rays, new_rays];
            obj.rays(3,:) = obj.rays(3,:) .* amplitude_inv_norm;
            obj.cluster_number = repmat([1 : n_old_rays], 1, 2*n_half + 1);
        end


%%
% Propagators through media and optical elements

        function rays_output = freeSpacePropagate(obj, rays_input, L)
            if size(L, 2) == 1
                propagator_first_order = [1, L, 0, 0, 0;
                    0, 1, 0, 0, 0;
                    0, 0, 1, 0, 0;
                    0, 0, 0, 1, 0;
                    0, 0, 0, 0, 1];

                first_order_term = propagator_first_order * rays_input;
            else
                first_order_term(1, :) = rays_input(1, :) + L .* rays_input(2, :);
                first_order_term(2:5, :) = rays_input(2:5, :);
            end

            phase_change = rays_input(5, :) .* L .* (1 + rays_input(2, :).^2 / 2);      % this provides angle-dependent phase change
            phase_change = rem(phase_change, 2*pi);       % renormalization
            second_order_term = [zeros(3, size(rays_input, 2)); phase_change; zeros(1, size(rays_input, 2))];        % good news: all other terms of second order with respect to y and angle are zero

            rays_output = first_order_term + second_order_term;
        end


        function status = freeSpace_new(obj, L)
            if nargin < 2
                L = obj.x_right - obj.x_current;
            end
            obj.setRays(obj.freeSpacePropagate(obj.rays, L));
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
                for i_order = 1 : (size(refraction_coeff, 2) - 1)
                    chrom_aberration_coeff = chrom_aberration_coeff + refraction_coeff(i_order + 1) * (2 * pi)^(-2*i_order) * (obj.rays(5, affected_rays)).^(2*i_order);
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
            rays_to_keep = (obj.rays(1,:) <= y_bottom) | (y_top <= obj.rays(1,:));
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


        function status = deleteRays(obj, rays_to_keep_log)
            obj.n_visible_rays = sum(rays_to_keep_log(1:obj.n_visible_rays));
            new_rays = obj.rays(:, rays_to_keep_log);
            obj.cluster_number = obj.cluster_number(rays_to_keep_log);
            obj.setRays(new_rays);
        end


        function status = calcIntensity(obj, resolution, method, max_rays_in_batch, coherence_length)
            status = 1;
            if nargin < 5
                coherence_length = 2*pi / min(obj.rays(5,:));
                if nargin < 4
                    max_rays_in_batch = 100;
                    if nargin < 3
                        method = 'cluster';
                        if nargin < 2
                            resolution = max(19, min(151, fix(size(obj.rays, 2)/(2*10^4))));
                        end
                    end
                end
            end
            obj.intensity = zeros(1,resolution);
            obj.intensity_meas_point = linspace(obj.y_bottom, obj.y_top, resolution);

            rays_y = obj.rays(1,:);
            rays_angles = obj.rays(2,:);
            rays_amp = obj.rays(3,:);
            numeric_array = [1 : size(rays_y, 2)];
            numeric_array_aux = [1 : resolution];

            switch method
                case 'cluster'
                    y_meas_points = linspace(obj.y_bottom, obj.y_top, resolution);
                    window = (obj.y_top - obj.y_bottom) / resolution;
                    amplitude = zeros(1,resolution);
                    k = obj.rays(5, 1);

                    for i_cluster = 1 : max(obj.cluster_number)
                        mask_cluster_log = (obj.cluster_number == i_cluster);
                        mask_rays_numeric = numeric_array(mask_cluster_log);
                        %mask_cluster_numeric = numeric_array(mask_cluster_log);
                        %mask_rays_log = (obj.y_bottom <= rays_y(mask_cluster_numeric)) & (rays_y(mask_cluster_numeric) <= obj.y_top);
                        %mask_rays_numeric = mask_cluster_numeric(mask_rays_log);
                        [y_min, i_min] = min(rays_y(mask_rays_numeric));
                        [y_max, i_max] = max(rays_y(mask_rays_numeric));
                        angle_max = rays_angles(mask_rays_numeric(i_max));
                        angle_min = rays_angles(mask_rays_numeric(i_min));
                        mask_points_log = (y_min <= y_meas_points) & (y_meas_points <= y_max);
                        mask_points_numeric = numeric_array_aux(mask_points_log);
                        if (size(mask_points_numeric, 2) == 0) && (obj.y_bottom < y_min) && (y_max < obj.y_top)
                            resolution = resolution + 1;
                            numeric_array_aux(end + 1) = resolution;
                            y_meas_points(end + 1) = (y_min + y_max) / 2;
                            amplitude(end + 1) = 0;
                            [y_meas_points, index] = sort(y_meas_points);
                            amplitude = amplitude(index);
                            mask_points_log = (y_min <= y_meas_points) & (y_meas_points <= y_max);
                            mask_points_numeric = numeric_array_aux(mask_points_log);
                        end

                        r_cluster = cos((angle_max + angle_min) / 2) * (y_max - y_min) / (angle_max - angle_min);
                        if isnan(r_cluster)
                            amplitude(mask_points_numeric) = sum(rays_amp(mask_rays_numeric) .* exp(i * obj.rays(4, mask_rays_numeric)));
                        elseif r_cluster == 0
                            amplitude(mask_points_numeric) = sum(rays_amp(mask_rays_numeric) .* exp(i * obj.rays(4, mask_rays_numeric)));
                        else
                            if (cos(angle_max) - cos(angle_min)) * r_cluster > 0
                                i_ref = i_min;
                            else
                                i_ref = i_max;
                            end
                            y_ref = rays_y(mask_rays_numeric(i_ref));
                            angle_ref = obj.rays(2, mask_rays_numeric(i_ref));
                            phase_ref = obj.rays(4, mask_rays_numeric(i_ref));
                            del_y = y_meas_points(mask_points_numeric) - y_ref;
                            del_r = angle_ref * del_y + 1/2 * del_y.^2 / r_cluster;
                            phases = rem(k * del_r, 2*pi) + phase_ref;
                            normale = cos(angle_ref) ./ (1 + del_r ./ r_cluster);
                            amplitude_init = [];
                            for i_y = 1 : size(y_meas_points(mask_points_numeric), 2)
                                diff_y = rays_y(mask_rays_numeric) - y_meas_points(mask_points_numeric(i_y));
                                [y_above, i_above] = min(diff_y(diff_y >= 0));
                                [y_below, i_below] = max(diff_y(diff_y <= 0));
                                linear_coeff = (rays_amp(mask_rays_numeric(i_above)) - rays_amp(mask_rays_numeric(i_below))) / (y_above - y_below);
                                if isnan(linear_coeff) || isinf(linear_coeff)
                                    amplitude_init(end + 1) = rays_amp(mask_rays_numeric(i_below));
                                else
                                    amplitude_init(end + 1) = rays_amp(mask_rays_numeric(i_below)) - y_below * linear_coeff;
                                end
                            end
                            amplitude(mask_points_numeric) = amplitude(mask_points_numeric) + amplitude_init .* sqrt(1 ./ (y_max - y_min) ./ (1 + del_r ./ r_cluster)) .* exp(i * phases);
                        end
                    end
                    intensity_y = y_meas_points;
                    intensity_value = abs(amplitude).^2;

                otherwise
                    amplitude = exp(i * obj.rays(4,:)) .* obj.rays(3,:);

                    power_y = (obj.y_top + obj.y_bottom) / 2;
                    window_half = (obj.y_top - obj.y_bottom) / 2;
                    mask_log = ((power_y(1) - window_half(1) < rays_y) & (rays_y <= power_y(1) + window_half(1)));
                    mask_list = {numeric_array(mask_log)};
                    n_rays_inside = size(rays_y(mask_list{1}), 2);
                    [max_rays_inside, index_max_elem] = max(n_rays_inside);

                    while max_rays_inside > max_rays_in_batch
                        [max_rays_inside, index_max_elem] = max(n_rays_inside);
                        while (window_half(index_max_elem) < coherence_length) && (max_rays_inside > 0)
                            n_rays_inside(index_max_elem) = 0;
                            [max_rays_inside, index_max_elem] = max(n_rays_inside);
                        end

                        window_half(index_max_elem) = window_half(index_max_elem) / 2;
                        window_half(end + 1) = window_half(index_max_elem);
                        power_y(end + 1) = power_y(index_max_elem) + window_half(index_max_elem);
                        power_y(index_max_elem) = power_y(index_max_elem) - window_half(index_max_elem);

                        mask_left_log = (power_y(index_max_elem) - window_half(index_max_elem) < rays_y(mask_list{index_max_elem})) & (rays_y(mask_list{index_max_elem}) <= power_y(index_max_elem) + window_half(index_max_elem));
                        mask_right_log = ~mask_left_log;
                        mask_left_numeric = mask_list{index_max_elem}(mask_left_log);
                        mask_right_numeric = mask_list{index_max_elem}(mask_right_log);
                        mask_list{index_max_elem} = mask_left_numeric;
                        mask_list{end + 1} = mask_right_numeric;

                        n_rays_inside(index_max_elem) = size(rays_y(mask_left_numeric), 2);
                        n_rays_inside(end + 1) = size(rays_y(mask_right_numeric), 2);
                    end

                    for j = 1 : size(power_y, 2)
                        power_value(j) = abs(sum(amplitude(mask_list{j})))^2;
                    end

                    delta_half = (obj.y_top - obj.y_bottom) / resolution / 2;
                    intensity_y = obj.y_bottom + 2 * delta_half * ([1:resolution] - 0.5);
                    for j = 1:resolution
                        mask = (intensity_y(j) - delta_half < power_y) & (power_y <= intensity_y(j) + delta_half);
                        intensity_value(j) = sum(power_value(mask));
                    end
            end

            max_intensity = max(intensity_value);
            if max_intensity ~= 0
                intensity_value = 1 / max_intensity * intensity_value;
                obj.intensity = intensity_value;
                obj.intensity_meas_point = intensity_y;
            end

            if obj.autoplot_enabled
                obj.drawIntensity(obj.intensity, obj.intensity_meas_point);
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
            %n_rays = size(obj.rays, 2);
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
            x_center = obj.x_current;
            width = half_diameter / 4;
            width_curve = width / 3;
            width_inner = width - 2 * width_curve;
            n_pieces = 16;
            r = (half_diameter^2 + width_curve^2) / (2 * width_curve);
            max_angle = asin(half_diameter / r);
            angle = linspace(-max_angle, max_angle, n_pieces);
            if f >= 0
                arc_y_left = y_center + r * sin(angle);
                arc_x_left = x_center + r * (cos(max_angle) - cos(angle)) - width_inner/2;
                arc_y_right = y_center - r * sin(angle);
                arc_x_right = x_center - r * (cos(max_angle) - cos(angle)) + width_inner/2;
            else
                arc_y_left = y_center + r * sin(angle);
                arc_x_left = x_center - r * (1 - cos(angle)) - width_inner/2;
                arc_y_right = y_center - r * sin(angle);
                arc_x_right = x_center + r * (1 - cos(angle)) + width_inner/2;
            end
            arc_y_left = [arc_y_left(1), arc_y_left];
            arc_x_left = [x_center, arc_x_left];
            arc_y_right(end + 1) = arc_y_right(end);
            arc_x_right(end + 1) = x_center;

            arc_y = [arc_y_left, arc_y_right];
            arc_x = [arc_x_left, arc_x_right];
            area(arc_x, arc_y, 'FaceColor', 'blue','FaceAlpha', 0.5, 'EdgeColor', 'black', 'EdgeAlpha', 1, 'LineWidth', 0.5);
            text(obj.x_current, y_center - half_diameter - 0.05 * (obj.y_top - obj.y_bottom), ['f = ', num2str(f)], 'FontSize', 14);
            obj.autoScale(1);
        end


        function drawWall(obj, y_min, y_max)
            line([obj.x_current; obj.x_current], [y_min; y_max], 'Color', [0, 0, 0], 'LineWidth', 2);
        end


        function drawIntensity(obj, intensity, Y)
            line([obj.x_right, obj.x_right], [obj.y_bottom, obj.y_top], 'Color', [0,0,0]);
            axis([obj.x_left, obj.x_right * 1.1, obj.y_bottom, obj.y_top]);
            X = obj.x_right + intensity * 0.1 * obj.x_right;

            plot(X, Y);
        end















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Old functions



function raysOut = thickLens(z,raysIn,p1,p2,d,n)

M = [1-p2*(d/n)      p1+p2-p1*p2*(d/n);
    -d/n             1-p1*(d/n)];

 raysOut = M*raysIn;
end


%lense surface power
function out= LensSurfacePower(z,n1,n2,R)
out=(n2-n1)/R;
end



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
