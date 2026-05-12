%**************************************************************************
%*****************************PM-ABL v2.0**********************************
%**************************************************************************
% 2-D Langevin model to calculate ensemble of particle trajectories
% from U and V winds, variances, and relaxation timescales. Particles are
% binned to convert particle positions into perturbed aerosol concentrations.
%--------------------------------------------------------------------------
% This particular routine is designed to run the particle model for every case
% with median cloud fraction > 50% in the NEPac LES library.
% 2D fields of solar insolation (SOLIN) and net shortwave at TOA (SWNT) are
% used to calculate initial reflectance values to be perturbed by the ship
% aerosols from the particle model. Twomey forcing is calculated for 
% all cases and all injection
% timing scenarios (default: at sunrise). Testing of the injection duration
% can be done using the 
% variable injection_duration. All tests assume random ship motion within 
% the assigned grid and boundary conditions are doubly periodic.

% Inversion heights are read in from ZINV for the aerosol concentration 
% calculation.
% Boundary-layer averaged accumulation mode aerosol is used in the Twomey
% calculation (as Nd_background).
%--------------------------------------------------------------------------
% This version calculates repeated diurnal cycles until equilibrium 
% aerosol conditions are reached (Day 1 and Day 2 indices are 
% changed using diurnal_indices_day1 and diurnal_indices_subsequent).
%--------------------------------------------------------------------------

% C_0 parameter changes as a function of precipitation. Case is considered
% precipitating if the surface rain rate is greater than 0.01 mm/day at any
% point
prec_thresh = 0.01; % prec. threshold in mm/day
C_0_prec = 0.15; % C_0 for precipitating cases (McMichael et al., 2024)
C_0_noprec = 0.5; % C_0 for non-precipitating cases (McMichael et al., 2024)

% Particles decay to represent wet/dry deposition processes on an e-folding timescale
% set by tau_decay (default = 2 days). This can be conceptualized as the
% accumulated probability of a removal event occurring over a particle batch lifetime.
% tau_decay = [0.25, 0.5, 1, 2, 4] --> num_days = [1, 2, 4, 8, 16]
tau_decay = 0.5; % days

% Array of all LES cases with median cloud fraction > 50%
case_name = {'01_20180627', '02_20200830', '02_20210819', '02_20210629', ...
             '03_20180617', '03_20180623', '03_20190729', '03_20200821', ...
             '04_20180731', '04_20180804', '04_20190624', '04_20200804', ...
             '05_20210806', '06_20180704', '06_20190821', '06_20200609', ...
             '06_20210823'};

num_les_cases = length(case_name);

% Ship and domain parameters
num_ships = 10; % number of ships
part_num = 1000; % number of particles injected at each time step
x_range = 200; % km
y_range = 200; % km

dt = 600; % time step in seconds
injection_rate = 10^16; % particles per second
rho_mg = 1.15*10^6; % mg/m^3
total_deleted_from_first_group = 0.0;

% Time parameters
time_s = double((time - time(1)) * 3600.0 * 24.0);
dt_les = round(time_s(2) - time_s(1));
start_t = floor(min(time_s));

% Plotting parameters
save_interval = 1200; % in seconds [this matches the 2-D data]
day_s = 86400; % seconds in a day
num_days = 2; % number of diurnal cycles

% Create time arrays
time_les = transpose(start_t:dt_les:num_days*day_s);
new_time = transpose(start_t:dt:num_days*day_s);
new_time_hr = new_time / 3600.0;

% Define injection cases relative to dawn
injection_offsets = [0]; % Hours relative to dawn: -4, 0, +4
num_cases = length(injection_offsets);
injection_duration = 12 * 3600; % injection duration in seconds (12 hours)

% Grid parameters for particle density calculation
grid_size = x_range; % 1 x 1 km grid for now
x_edges = linspace(0, x_range * 1000, grid_size + 1);
y_edges = linspace(0, y_range * 1000, grid_size + 1);

x_range_m = x_range * 1000; % Convert to meters
y_range_m = y_range * 1000; % Convert to meters

block_size = 1; % 1-km coarse graining
x_blocks = floor(x_range / block_size);
y_blocks = floor(y_range / block_size);

% Initialize arrays for all cases
save_size = ((day_s * num_days) / save_interval) + 1;

domain_twomey_all = zeros(length(time_les), num_cases, num_les_cases);
particle_conc_2D_all = zeros(x_blocks, y_blocks, save_size, num_cases, num_les_cases);

density_gif_initialized = false;
prob_gif_initialized = false;
density_filename = 'particle_density.gif';
prob_filename = 'prob_conc.gif';

% Loop over all LES cases in case_name
for les_idx = 1:num_les_cases

    current_case = case_name{les_idx};
    disp(['Processing LES case: ', current_case]);

    % Load statistics from NEP_'case_name'.mat
    stats_file = ['NEP_', current_case, '.mat'];
    stats_data = load(stats_file);
    U_bl = stats_data.U_bl;
    U2_bl = stats_data.U2_bl;
    V_bl = stats_data.V_bl;
    V2_bl = stats_data.V2_bl;
    TLu_bl = stats_data.TLu_bl;
    TLv_bl = stats_data.TLv_bl;
    PREC = stats_data.PREC;
    ZINV = stats_data.ZINV;
    NAc_bl = stats_data.NAc_bl;
    RADSWDN = stats_data.RADSWDN;
    time = stats_data.time;

    % Load SOLIN and SWNT for Twomey calculations from NEP_'case_name'_2D.nc
    nc_file = ['NEP_', current_case, '_2D.nc'];
    SOLIN = ncread(nc_file, 'SOLIN'); % [x, y, time]
    SWNT = ncread(nc_file, 'SWNT'); % [x, y, time]

    % Make RADSWDN a 1D array
    radswdn_1d = RADSWDN(1,:);

    % Extract the first 24 hours (indices 1:72) and repeat for num_days
    diurnal_indices_day1 = 1:73; % Select 24 hour period (including t=0)
    diurnal_indices_subsequent = 2:73; % Subsequent day indices

    % Extract the second 24 hours (indices 72:144) and repeat for num_days
    %diurnal_indices_day1 = 72:144; % Select 24 hour period (including t=0)
    %diurnal_indices_subsequent = 73:144; % Subsequent day indices

    diurnal_day1 = length(diurnal_indices_day1);
    diurnal_ind = length(diurnal_indices_subsequent);

    tot_elements = diurnal_day1 + (num_days-1)*diurnal_ind;

    % Preallocate repeated arrays
    U_bl_repeated = zeros(tot_elements, 1);
    U2_bl_repeated = zeros(tot_elements, 1);
    V_bl_repeated = zeros(tot_elements, 1);
    V2_bl_repeated = zeros(tot_elements, 1);
    TLu_bl_repeated = zeros(tot_elements, 1);
    TLv_bl_repeated = zeros(tot_elements, 1);
    PREC_repeated = zeros(tot_elements, 1);
    ZINV_repeated = zeros(tot_elements, 1);
    NAc_bl_repeated = zeros(tot_elements, 1);
    RADSWDN_repeated = zeros(tot_elements, 1);

    % Preallocate 2D arrays
    SOLIN_repeated = zeros(size(SOLIN, 1), size(SOLIN, 2), tot_elements);
    SWNT_repeated = zeros(size(SWNT, 1), size(SWNT, 2), tot_elements);

    % Preallocate
    SOLIN_twomey = zeros(size(SOLIN, 1), size(SOLIN, 2), tot_elements, num_les_cases);
    SWNT_twomey = zeros(size(SWNT, 1), size(SWNT, 2), tot_elements, num_les_cases);
   
    % Initialize 2D arrays
    SOLIN_day1 = SOLIN(:,:,1:diurnal_day1);
    SOLIN_diurnal = SOLIN(:,:,1:diurnal_ind);
    SWNT_day1 = SWNT(:,:,1:diurnal_day1);
    SWNT_diurnal = SWNT(:,:,1:diurnal_ind);

    % Initialize 2D arrays
    % SOLIN_day1 = SOLIN(:,:,72:144);
    % SOLIN_diurnal = SOLIN(:,:,73:144);
    % SWNT_day1 = SWNT(:,:,72:144);
    % SWNT_diurnal = SWNT(:,:,73:144);

    % Initialize repeated arrays 
    U_bl_repeated(1:diurnal_day1) = U_bl(diurnal_indices_day1)';
    U2_bl_repeated(1:diurnal_day1) = U2_bl(diurnal_indices_day1)';
    V_bl_repeated(1:diurnal_day1) = V_bl(diurnal_indices_day1)';
    V2_bl_repeated(1:diurnal_day1) = V2_bl(diurnal_indices_day1)';
    TLu_bl_repeated(1:diurnal_day1) = TLu_bl(diurnal_indices_day1)';
    TLv_bl_repeated(1:diurnal_day1) = TLv_bl(diurnal_indices_day1)';
    PREC_repeated(1:diurnal_day1) = PREC(diurnal_indices_day1)';
    ZINV_repeated(1:diurnal_day1) = ZINV(diurnal_indices_day1)';
    NAc_bl_repeated(1:diurnal_day1) = NAc_bl(diurnal_indices_day1)';
    RADSWDN_repeated(1:diurnal_day1) = radswdn_1d(diurnal_indices_day1)';
    
    %Initialize 2D arrays
    SOLIN_repeated(:,:,1:diurnal_day1) = SOLIN_day1;
    SWNT_repeated(:,:,1:diurnal_day1) = SWNT_day1;

    % Initialize start index for subsequent days
    start_idx = diurnal_day1 + 1; 

    % Loop for subsequent days
    for d = 2:num_days

        end_idx = start_idx + diurnal_ind - 1; 

        U_bl_repeated(start_idx:end_idx) = U_bl(diurnal_indices_subsequent)';
        U2_bl_repeated(start_idx:end_idx) = U2_bl(diurnal_indices_subsequent)';
        V_bl_repeated(start_idx:end_idx) = V_bl(diurnal_indices_subsequent)';
        V2_bl_repeated(start_idx:end_idx) = V2_bl(diurnal_indices_subsequent)';
        TLu_bl_repeated(start_idx:end_idx) = TLu_bl(diurnal_indices_subsequent)';
        TLv_bl_repeated(start_idx:end_idx) = TLv_bl(diurnal_indices_subsequent)';
        PREC_repeated(start_idx:end_idx) = PREC(diurnal_indices_subsequent)';
        ZINV_repeated(start_idx:end_idx) = ZINV(diurnal_indices_subsequent)';
        NAc_bl_repeated(start_idx:end_idx) = NAc_bl(diurnal_indices_subsequent)';
        RADSWDN_repeated(start_idx:end_idx) = radswdn_1d(diurnal_indices_subsequent)';

        SOLIN_repeated(:,:,start_idx:end_idx) = SOLIN_diurnal;
        SWNT_repeated(:,:,start_idx:end_idx) = SWNT_diurnal;

        start_idx = end_idx + 1; % Update start index for next day
    end

    SOLIN_twomey(:,:,:,les_idx) = SOLIN_repeated;
    SWNT_twomey(:,:,:,les_idx) = SWNT_repeated;

    % Interpolate fields
    MEAN_U_interp = interp1(time_les, U_bl_repeated, new_time);
    %MEAN_U_interp = zeros(size(new_time));
    MEAN_U2_interp = interp1(time_les, U2_bl_repeated, new_time);
    MEAN_V_interp = interp1(time_les, V_bl_repeated, new_time);
    %MEAN_V_interp = zeros(size(new_time));
    MEAN_V2_interp = interp1(time_les, V2_bl_repeated, new_time);
    TS_u_interp = interp1(time_les, TLu_bl_repeated, new_time).*3600.0;
    TS_v_interp = interp1(time_les, TLv_bl_repeated, new_time).*3600.0;
    Prec_interp = interp1(time_les, PREC_repeated, new_time);
    Bl_depth = interp1(time_les, ZINV_repeated, new_time).*1000.0;
    NAc_interp = interp1(time_les, NAc_bl_repeated, new_time);
    RADSW_interp = interp1(time_les, RADSWDN_repeated, new_time);

    % Determine nighttime periods for Twomey calculation
    night_ind = (RADSWDN_repeated(:) == 0);
    dusk_Twomey = find(diff([0 night_ind']) == 1);
    dawn_Twomey = find(diff([0 night_ind']) == -1);

    if (RADSWDN_repeated(end) == 0)
        dawn_Twomey = [dawn_Twomey, length(RADSWDN_repeated)];
    end

    % Determine nighttime periods for injection schedule
    night_ind = (RADSW_interp(:) == 0);
    dusk = find(diff([0 night_ind']) == 1);
    dawn = find(diff([0 night_ind']) == -1);

    if (RADSW_interp(end) == 0)
        dawn = [dawn, length(RADSW_interp)];
    end

    % Smooth diurnal cycle to avoid abrupt changes from one cycle to
    % another
    smooth_window = 20;
    MEAN_U_interp = movmean(MEAN_U_interp, smooth_window);
    MEAN_U2_interp = movmean(MEAN_U2_interp, smooth_window);
    MEAN_V_interp = movmean(MEAN_V_interp, smooth_window);
    MEAN_V2_interp = movmean(MEAN_V2_interp, smooth_window);
    TS_u_interp = movmean(TS_u_interp, smooth_window);
    TS_v_interp = movmean(TS_v_interp, smooth_window);
    Prec_interp = movmean(Prec_interp, smooth_window);
    Bl_depth = movmean(Bl_depth, smooth_window);
    NAc_interp = movmean(NAc_interp, smooth_window);

    % Compute injection periods
    num_dawns = length(dawn); % Number of dawns detected
    injection_periods = cell(num_cases, num_dawns);
    for case_idx = 1:num_cases
        offset = injection_offsets(case_idx) * 3600;
        for d = 1:num_dawns
            start_time = new_time(dawn(d)) + offset;
            end_time = start_time + injection_duration;
            if end_time > new_time(end)
                end_time = end_t;
            end
            injection_periods{case_idx, d} = [start_time, end_time];
        end
    end

    % Loop over injection cases
    for case_idx = 1:num_cases
        % Initialize arrays
        part_pos_x = [];
        part_pos_y = [];
        part_vel_u = [];
        part_vel_v = [];
        part_release_time = [];

        % Initialize particle position width
        std_dev_x = 1.0 * 1000.0;
        std_dev_y = 1.0 * 1000.0;

        % Initialize ship positions
        initial_positions_x = x_range * 1000 * rand(num_ships, 1);
        initial_positions_y = y_range * 1000 * rand(num_ships, 1);

        % Random ship motions
        %10 m/s ship velocity
        ship_motion_y_options = [-10, 0, 10];
        ship_motion_x_options = [-10, 0, 10];
        %5 m/s ship velocity
        %ship_motion_y_options = [-5, 0, 5];
        %ship_motion_x_options = [-5, 0, 5];
        %3 m/s ship velocity
        %ship_motion_y_options = [-3, 0, 3];
        %ship_motion_x_options = [-3, 0, 3];
        %1 m/s ship velocity
        %ship_motion_y_options = [-1, 0, 1];
        %ship_motion_x_options = [-1, 0, 1];
        %0 m/s ship velocity
        %ship_motion_y_options = [-0.01, 0, 0.01];
        %ship_motion_x_options = [-0.01, 0, 0.01];


        [ship_motion_x_grid, ship_motion_y_grid] = meshgrid(ship_motion_x_options, ship_motion_x_options);
        valid_ship_motion_indices = find(~(ship_motion_x_grid == 0 & ship_motion_y_grid == 0));
        num_valid_motions = length(valid_ship_motion_indices);
        random_indices = randi(num_valid_motions, num_ships, 1);
        selected_ship_motion_indices = valid_ship_motion_indices(random_indices);
        ship_motion_x = ship_motion_x_grid(selected_ship_motion_indices);
        ship_motion_y = ship_motion_y_grid(selected_ship_motion_indices);

        % Initialize particle velocities
        std_dev_vel = 0.0;
        mean_u_vel = MEAN_U_interp(1);
        mean_v_vel = MEAN_V_interp(1);

        % Preallocate particle concentration
        particle_conc_2D = zeros(x_blocks, y_blocks, save_size);
        num_step = 1;

        % Main time loop
        for n = 1:length(new_time)
            current_time = new_time(n);

            % Check if current time is within any injection period for this case
            is_injecting = false;
            for d = 1:num_dawns
                period = injection_periods{case_idx, d};
                if (current_time >= period(1) && current_time <= period(2))
                    is_injecting = true;
                    break;
                end
            end

            if is_injecting
                % Inject new particles at the ship's current position for
                % each ship
                for m = 1:num_ships
                    % Update ship position based on its motion
                    initial_positions_x(m) = initial_positions_x(m) + ship_motion_x(m) * dt;
                    initial_positions_y(m) = initial_positions_y(m) + ship_motion_y(m) * dt;

                    % Check if the ship has left the domain
                    if initial_positions_x(m) < 0 || initial_positions_x(m) > x_range * 1000 || ...
                            initial_positions_y(m) < 0 || initial_positions_y(m) > y_range * 1000
                        
                        % Regenerate the ship's position and motion
                        initial_positions_x(m) = x_range * 1000 * rand();
                        initial_positions_y(m) = y_range * 1000 * rand();

                        random_index = randi(num_valid_motions, 1);
                        selected_ship_motion_index = valid_ship_motion_indices(random_index);
                        ship_motion_x(m) = ship_motion_x_grid(selected_ship_motion_index);
                        ship_motion_y(m) = ship_motion_y_grid(selected_ship_motion_index);
                    end

                    part_pos_init_x = initial_positions_x(m) + std_dev_x * randn(part_num, 1);
                    part_pos_init_y = initial_positions_y(m) + std_dev_y * randn(part_num, 1);
                    part_vel_u_init = mean_u_vel + std_dev_vel * randn(part_num, 1);
                    part_vel_v_init = mean_v_vel + std_dev_vel * randn(part_num, 1);

                    part_pos_x = [part_pos_x; part_pos_init_x];
                    part_pos_y = [part_pos_y; part_pos_init_y];
                    part_vel_u = [part_vel_u; part_vel_u_init];
                    part_vel_v = [part_vel_v; part_vel_v_init];
                    part_release_time = [part_release_time; repmat(current_time / 3600, part_num, 1)];
                end
            end

            % Update particle positions and velocities for all particles
            for i = 1:length(part_pos_x)

                % Check for precipitation to determine C_0 value
                if Prec_interp(n) > prec_thresh
                    TLu = TS_u_interp(n) * (1 / (0.75 * C_0_prec));
                    TLv = TS_v_interp(n) * (1 / (0.75 * C_0_prec));
                else
                    TLu = TS_u_interp(n) * (1 / (0.75 * C_0_noprec));
                    TLv = TS_v_interp(n) * (1 / (0.75 * C_0_noprec));
                end

                % Zonal velocities and positions
                part_vel_u(i) = part_vel_u(i) + (MEAN_U_interp(n) - part_vel_u(i)) * (dt / TLu) ...
                    + sqrt((2 * MEAN_U2_interp(n) * dt) / TLu) * randn();
                part_pos_x(i) = part_pos_x(i) + part_vel_u(i) * dt;

                % Apply periodic boundary conditions in x direction
                if part_pos_x(i) < 0
                    part_pos_x(i) = part_pos_x(i) + x_range * 1000;
                elseif part_pos_x(i) > x_range * 1000
                    part_pos_x(i) = part_pos_x(i) - x_range * 1000;
                end

                % Meridional velocities and positions
                part_vel_v(i) = part_vel_v(i) + (MEAN_V_interp(n) - part_vel_v(i)) * (dt / TLv) ...
                    + sqrt((2 * MEAN_V2_interp(n) * dt) / TLv) * randn();
                part_pos_y(i) = part_pos_y(i) + part_vel_v(i) * dt;

                % Apply periodic boundary conditions in y direction
                if part_pos_y(i) < 0
                    part_pos_y(i) = part_pos_y(i) + y_range * 1000;
                elseif part_pos_y(i) > y_range * 1000
                    part_pos_y(i) = part_pos_y(i) - y_range * 1000;
                end
            end

            % Track remaining particles after accounting for deposition
            % processes
            remaining_particles = true(1, length(part_release_time));
            unique_release_times = unique(part_release_time);

            for rt = length(unique_release_times):-1:1
                release_time = unique_release_times(rt);
                indices = find(part_release_time == release_time);
                particle_age = current_time - release_time * 3600.0;
                decay_prob = 1 - exp(-particle_age / (tau_decay * day_s));
                remaining_fraction = 1 - decay_prob;

                if ~isempty(indices)
                    num_particles_remaining = round(remaining_fraction * (num_ships * part_num));
                    num_particles_to_delete = length(indices) - num_particles_remaining;

                    if num_particles_to_delete > 0
                        particles_to_delete = indices(randperm(length(indices), num_particles_to_delete));
                        remaining_particles(particles_to_delete) = false;
                    end
                end
            end

            % Delete the particles after marking
            part_pos_x(~remaining_particles) = [];
            part_pos_y(~remaining_particles) = [];
            part_vel_u(~remaining_particles) = [];
            part_vel_v(~remaining_particles) = [];
            part_release_time(~remaining_particles) = [];

            % Print remaining particles after deletion
            disp(['Case ', num2str(case_idx), ' - Number of particles after deletion: ', num2str(length(part_release_time))]);

            % Calculate and save particle concentration every save_interval
            if mod(current_time, save_interval) == 0

                % figure('visible', 'off');
                % Calculate particle density/concentration
                particle_density = histcounts2(part_pos_x, part_pos_y, x_edges, y_edges);
                dx = x_edges(2) - x_edges(1);
                dy = y_edges(2) - y_edges(1);
                cell_area = dx * dy;
                count_per_particle = (injection_rate * dt) / part_num;
                part_dens = (particle_density * count_per_particle) / cell_area;
                particle_concentration = part_dens ./ (Bl_depth(n) * rho_mg);

                % Coarse grain the particle_concentration matrix
                coarse_particle_concentration = zeros(x_blocks, y_blocks);
                for i = 1:x_blocks
                    for j = 1:y_blocks
                        x_start = (i-1)*block_size + 1;
                        x_end = i*block_size;
                        y_start = (j-1)*block_size + 1;
                        y_end = j*block_size;
                        block = particle_concentration(x_start:x_end, y_start:y_end);
                        coarse_particle_concentration(i, j) = mean(block(:));
                    end
                end

                % Save 2D particle concentration
                particle_conc_2D(:, :, num_step) = coarse_particle_concentration;
                particle_conc_2D_all(:, :, num_step, case_idx, les_idx) = coarse_particle_concentration;
                num_step = num_step + 1;

                % cmap = colormap('jet');
                % cmap(1, :) = [0 0 0];
                % 
                % imagesc('XData', x_edges / 1000, 'YData', y_edges / 1000, 'CData', coarse_particle_concentration');
                % set(gca, 'FontSize', 10, 'FontName', 'Georgia');
                % colorbar;
                % colormap(cmap);
                % caxis([0 100]);
                % xlim([0 x_range])
                % ylim([0 y_range])
                % xlabel('x [km]','FontName', 'Georgia');
                % ylabel('y [km]','FontName', 'Georgia');
                % title(['Particle Concentration (#/mg) at Time = ', num2str(current_time / 86400), ' Day']);
                % grid on;
                % 
                % % Capture the density plot as an image
                % frame = getframe(gcf);
                % im = frame2im(frame);
                % [imind, cm] = rgb2ind(im, 256);
                % 
                % % Write to the density GIF file
                % if ~density_gif_initialized
                %     imwrite(imind, cm, density_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
                %     density_gif_initialized = true;
                % else
                %     imwrite(imind, cm, density_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
                % end
                % 
                % close(gcf);

            end
        end
    end

    %----------------------------------------------------------------------------
    % Compute Twomey forcing
    %----------------------------------------------------------------------------
    delta_ac_all = zeros(200, 200, length(time_les), num_cases);
    rad_forcing_all = zeros(200, 200, length(time_les), num_cases);

    for case_idx = 1:num_cases

        for t = dawn_Twomey(num_days):dusk_Twomey(num_days+1)
            
            % Find the corresponding particle concentration slice
            part_conc_2D_t = particle_conc_2D_all(:, :, t, case_idx, les_idx);

            % Compute reflectance from 2-D LES data
            SWUP = SWNT_repeated(:, :, t) - SOLIN_repeated(:, :, t);
            REFLECT = -SWUP ./ SOLIN_repeated(:, :, t);

            % Tile and clip reflectance to match particle concentration domain
            tiled_REFLECT = repmat(REFLECT, 4, 4);
            clipped_REFLECT = tiled_REFLECT(1:2000, 1:2000);
            clipped_REFLECT(clipped_REFLECT < 0.07) = 0;

            tiled_SOLIN = repmat(SOLIN_repeated(:, :, t), 4, 4);
            clipped_SOLIN = tiled_SOLIN(1:2000, 1:2000);

            % Coarse-grain reflectance and SOLIN
            coarse_REFLECT = zeros(200, 200);
            coarse_SOLIN = zeros(200, 200);
            block_size_les = 10;

            for i = 1:200
                for j = 1:200
                    block = clipped_REFLECT((i-1)*block_size_les + 1:i*block_size_les, ...
                        (j-1)*block_size_les + 1:j*block_size_les);
                    coarse_REFLECT(i, j) = mean(block(:));

                    block_solin = clipped_SOLIN((i-1)*block_size_les + 1:i*block_size_les, ...
                        (j-1)*block_size_les + 1:j*block_size_les);
                    coarse_SOLIN(i, j) = mean(block_solin(:));
                end
            end

            % Compute delta_ac for the current time step
            Nd_background = NAc_bl_repeated(t); %account for varying background aerosol
            phi = 1.0; % Atmospheric correction factor
            delta_ac = zeros(200, 200);

            %Twomey calculation
            for i = 1:200
                for j = 1:200
                    if part_conc_2D_t(i, j) ~= 0.0
                        delta_ac(i, j) = (coarse_REFLECT(i, j)*(1-coarse_REFLECT(i, j))*...
                            ((((part_conc_2D_t(i, j)+Nd_background)/Nd_background)^(1/3))-1)) / ...
                            (1 + coarse_REFLECT(i, j)*((((part_conc_2D_t(i, j)+Nd_background)/Nd_background)^(1/3))-1));
                    else
                        delta_ac(i, j) = 0.0;
                    end
                end
            end
            delta_ac_all(:, :, t, case_idx) = delta_ac;

            % Compute radiative forcing
            rad_forcing = coarse_SOLIN .* delta_ac;
            rad_forcing_all(:, :, t, case_idx) = -rad_forcing;

            % Domain-average Twomey forcing
            domain_twomey_all(t, case_idx, les_idx) = mean(rad_forcing(:));
        end
    end

end

mean_Twomey = mean(fillmissing(domain_twomey_all(end-71:end, :, :), 'constant', 0), [1 3]);