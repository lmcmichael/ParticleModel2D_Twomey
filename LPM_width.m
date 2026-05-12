%**************************************************************************
%*****************************PM-ABL v2.0**********************************
%**************************************************************************
%**************************************************************************
% 2-D Langevin model to calculate ensemble of particle trajectories
% from U and V winds, variances, and relaxation timescales
%--------------------------------------------------------------------------
% This routine calculates the east-west and north-south plume widths of
% initial cloud of injected particles

% C_0 parameter changes as a function of precipitation. Case is considered
% precipitating if the surface rain rate is greater than 0.01 mm/day

% no particle decay accounted for in plume width estimations

% mean winds are set to 0. 
%--------------------------------------------------------------------------
%region and date of NEPac LES simulation
%.mat file created with NEP_simulations.m
    NEPac_array = {'NEP_01_20180627','NEP_01_20180815','NEP_01_20190801','NEP_01_20190803','NEP_01_20200607',...
                   'NEP_01_20200822','NEP_01_20210601','NEP_01_20210807','NEP_01_20210808'...
                   'NEP_02_20180808','NEP_02_20190605','NEP_02_20190725','NEP_02_20190829',...
                   'NEP_02_20200614','NEP_02_20200830','NEP_02_20210629','NEP_02_20210819',...
                   'NEP_02_20210828','NEP_03_20180617','NEP_03_20180621','NEP_03_20180623',...
                   'NEP_03_20190729','NEP_03_20200729','NEP_03_20200821','NEP_03_20210608',...
                   'NEP_03_20210807','NEP_03_20210815','NEP_04_20180622','NEP_04_20180731',...
                   'NEP_04_20180804','NEP_04_20180822','NEP_04_20190624','NEP_04_20190808',...
                   'NEP_04_20190814','NEP_04_20200804','NEP_04_20210828','NEP_05_20180601','NEP_05_20180702',...
                   'NEP_05_20180813','NEP_05_20190715','NEP_05_20190720',...
                   'NEP_05_20200627','NEP_05_20200730','NEP_05_20210806','NEP_05_20210816',...
                   'NEP_06_20180704','NEP_06_20180708','NEP_06_20190619','NEP_06_20190803',...
                   'NEP_06_20190821','NEP_06_20200609','NEP_06_20200819','NEP_06_20210704',...
                   'NEP_06_20210823'};

num_cases = length(NEPac_array);
max_time_steps = 144; 

zonal_width = NaN(max_time_steps, num_cases);
meridional_width = NaN(max_time_steps, num_cases);
average_width = NaN(max_time_steps, num_cases);
average_widths = NaN(max_time_steps, num_cases);
spreading_rates = NaN(max_time_steps - 1, num_cases);
spread_rates_smooth = NaN(max_time_steps - 1,num_cases);
spread_rates_smooth_low = NaN(max_time_steps - 1,num_cases);
spread_rates_smooth_mid = NaN(max_time_steps - 1,num_cases);
median_CF_values = NaN(num_cases, 1);
median_LWP_values = NaN(num_cases, 1);
median_NAc_values = NaN(num_cases, 1);
median_ZINV_values = NaN(num_cases, 1);

for case_idx = 1:num_cases
    case_name = NEPac_array{case_idx};
    filename = strcat(case_name, '.mat');
    
    load(filename);
   
    % Run Particle Model
    % Parameters
    ens_num = 1; %number of ships
    part_num = 5000; %number of particles
    x_range = 100000; % km
    y_range = 100000; % km
    %particle decay timescale
    tau_decay = 2; %days
    %C_0 for precipitating cases
    C_0_prec = 0.15;
    %C_0 for non-precipitating cases
    C_0_noprec = 0.5;
    dt = 1200; % time step of particle model in seconds
    injection_rate = 10^16; %particles per second
    rho_mg = 1.15*10^6; %g/m3
    total_deleted_from_first_group = 0.0;
    % Define the injection schedule
    injection_duration = 900; % hours --> seconds
    total_cycle_duration = 48 * 3600; % hours --> seconds

    % Time parameters
    day_s = 86400; %number of seconds in day
    num_days = 2; %number of diurnal cycles to be simulated
    time_s = double((time - time(1)) * 3600.0 * 24.0);
    dt_les = round(time_s(2) - time_s(1));
    start_t = floor(min(time_s));
    end_t = time_s(end);
    %create num_days long array of dt_les
    time_les = transpose(0:dt_les:end_t);
    new_time = transpose(start_t:dt:end_t);
    new_time_hr = new_time / 3600.0;
    num_step = 1;

    %cloud fraction bins
    %<0.25, 0.25-0.75, >0.75
    % Interpolate
    low_thresh = 0.50;
    high_thresh = 1.0;
    CF_interp = interp1(time_les, CF, new_time);
    median_CF_values(case_idx) = median(CF_interp);

    low_CF = find(median_CF_values < low_thresh);
    mid_CF = find(median_CF_values >= low_thresh & median_CF_values <= high_thresh);
    high_CF = find(median_CF_values > high_thresh);

    %LWP bins (g/m2)
    %<10, 10-50, >50
    % Interpolate
    low_thresh = 10;
    high_thresh = 50;
    LWP_interp = interp1(time_les, LWP, new_time);
    median_LWP_values(case_idx) = median(LWP_interp);

    low_LWP = find(median_LWP_values < low_thresh);
    mid_LWP = find(median_LWP_values >= low_thresh & median_LWP_values <= high_thresh);
    high_LWP = find(median_LWP_values > high_thresh);

    %NAc bins (#/mg)
    %<20, 20-100, >100
    % Interpolate
    low_thresh = 20;
    high_thresh = 100;
    NAc_interp = interp1(time_les, NAc_bl, new_time);
    median_NAc_values(case_idx) = median(NAc_interp);

    low_NAc = find(median_NAc_values < low_thresh);
    mid_NAc = find(median_NAc_values >= low_thresh & median_NAc_values <= high_thresh);
    high_NAc = find(median_NAc_values > high_thresh);

    %Boundary layer depth (km)
    %<0.8, 0.8-1.2, >1.2
    % Interpolate
    low_thresh = 0.8;
    high_thresh = 1.2;
    ZINV_interp = interp1(time_les, ZINV, new_time);
    median_ZINV_values(case_idx) = median(ZINV_interp);

    low_ZINV = find(median_ZINV_values < low_thresh);
    mid_ZINV = find(median_ZINV_values >= low_thresh & median_ZINV_values <= high_thresh);
    high_ZINV = find(median_ZINV_values > high_thresh);

    % arrays to store plume widths
    zonal_width = zeros(length(new_time),1);
    meridional_width = zeros(length(new_time),1);
    average_width = zeros(length(new_time),1);

    % Variables to store standard deviation 
    ship_std_x = zeros(ens_num, 1); 
    ship_std_y = zeros(ens_num, 1); 

    % Interpolate U_bl, U2_bl, V_bl, V2_bl, TLu_bl, and TLv_bl
    MEAN_U_interp = interp1(time_les, U_bl, new_time);
    MEAN_U_interp(:) = 0;
    MEAN_U2_interp = interp1(time_les, U2_bl, new_time);
    MEAN_V_interp = interp1(time_les, V_bl, new_time);
    MEAN_V_interp(:) = 0;
    MEAN_V2_interp = interp1(time_les, V2_bl, new_time);
    TS_u_interp = interp1(time_les, TLu_bl, new_time).*3600.0; %convert to s
    TS_v_interp = interp1(time_les, TLv_bl, new_time).*3600.0; %convert to s

    % Interpolate precipitation array
    prec_thresh = 0.01; %prec. threshold in mm/day
    Prec_interp = interp1(time_les, PREC, new_time);

    % Interpolate inversion height for particle density calculation
    Bl_depth = interp1(time_les, ZINV, new_time).*1000.0; %height in meters

    % Initialize arrays to store particle positions and velocities
    part_pos_x = [];
    part_pos_y = [];
    part_vel_u = [];
    part_vel_v = [];
    part_release_time = [];

    % Initialize particle positions based on first LES datapoint (t = 2.0 hr)
    std_dev_x = 1.0 * 100.0; % initial x-plume width
    std_dev_y = 1.0 * 100.0; % initial y-plume width

    % Initialize ship at center of domain
    initial_positions_x = (x_range * 1000)/2.0 ; % random x positions in meters
    initial_positions_y = (y_range * 1000)/2.0; % random y positions in meters

    % Random ship motions
    ship_motion_y_options = [0, 0, 0];
    ship_motion_x_options = [0, 0, 0];
    [ship_motion_x_grid, ship_motion_y_grid] = meshgrid(ship_motion_x_options, ship_motion_y_options);
    valid_ship_motion_indices = find(~(ship_motion_x_grid == 0 & ship_motion_y_grid == 0));

    % Select random indices for ship motions
    ship_motion_x = 0;
    ship_motion_y = 0;

    % Initialize particle velocities
    std_dev_vel = 0.0;
    mean_u_vel = MEAN_U_interp(1);
    mean_v_vel = MEAN_V_interp(1);

    % Define grid parameters for density calculation
    grid_size = x_range ; % 1 x 1 km grid
    x_edges = linspace(0, x_range * 1000, grid_size + 1);
    y_edges = linspace(0, y_range * 1000, grid_size + 1);

    x_range_m = x_range * 1000; % Convert to meters
    y_range_m = y_range * 1000; % Convert to meters

    % Main time loop to calculate each ensemble particle trajectory
    for n = 1:length(new_time)
        current_time = new_time(n);

        % Determine the current time within the total cycle duration
        time_in_cycle = mod(current_time, total_cycle_duration);

        if time_in_cycle <= injection_duration
            % Inject new particles at the ship's current position for each ensemble member
            for m = 1:ens_num
                part_pos_init_x = initial_positions_x(m) + std_dev_x * randn(part_num, 1); % particle x position in meters
                part_pos_init_y = initial_positions_y(m) + std_dev_y * randn(part_num, 1); % particle y position in meters
                part_vel_u_init = mean_u_vel + std_dev_vel * randn(part_num, 1);
                part_vel_v_init = mean_v_vel + std_dev_vel * randn(part_num, 1);

                part_pos_x = [part_pos_x; part_pos_init_x];
                part_pos_y = [part_pos_y; part_pos_init_y];
                part_vel_u = [part_vel_u; part_vel_u_init];
                part_vel_v = [part_vel_v; part_vel_v_init];
                part_release_time = [part_release_time; repmat(current_time / 3600, part_num, 1)]; % release time in hours
            end
        end

        % Update particle positions and velocities for all particles
        for i = 1:length(part_pos_x)

            %check for precipitation to determine C_0 value
            if Prec_interp(n) > prec_thresh
                % Calculate relaxation timescale from k/eps from LES
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
                part_pos_x(i) = part_pos_x(i) + x_range * 1000; % wrap to the right side of the domain
            elseif part_pos_x(i) > x_range * 1000
                part_pos_x(i) = part_pos_x(i) - x_range * 1000; % wrap to the left side of the domain
            end

            % Meridional velocities and positions
            part_vel_v(i) = part_vel_v(i) + (MEAN_V_interp(n) - part_vel_v(i)) * (dt / TLv) ...
                + sqrt((2 * MEAN_V2_interp(n) * dt) / TLv) * randn();

            part_pos_y(i) = part_pos_y(i) + part_vel_v(i) * dt;

            % Apply periodic boundary conditions in y direction
            if part_pos_y(i) < 0
                part_pos_y(i) = part_pos_y(i) + y_range * 1000; % wrap to the top side of the domain
            elseif part_pos_y(i) > y_range * 1000
                part_pos_y(i) = part_pos_y(i) - y_range * 1000; % wrap to the bottom side of the domain
            end

        end

        % Print remaining particles after deletion
        disp(['Number of particles after deletion: ', num2str(length(part_release_time))]);

        first_injection_indices = find(part_release_time == 0); % Indices of particles that were initially injected

        % Calculate widths
        for m = 1:ens_num
            ship_pos_x = part_pos_x(first_injection_indices);
            ship_pos_y = part_pos_y(first_injection_indices);

            delta_x = ship_pos_x - initial_positions_x(m);
            delta_y = ship_pos_y - initial_positions_y(m);

            ship_std_x_test = 2.0*std(delta_x);
            ship_std_y_test = 2.0*std(delta_y);
        end

        zonal_width(n) = ship_std_x_test;
        meridional_width(n) = ship_std_y_test;
        average_width(n) = (ship_std_x_test + ship_std_y_test) ./ 2.0;

    end
    
    % Compute spreading rates 
    average_widths(:, case_idx) = average_width;
    spreading_rates(:, case_idx) = diff(average_widths(:, case_idx)) ./ dt * 3600 / 1000; % km/hr

end

%determine nighttime periods
night_ind = (RADSWDN(1,:)==0);
dusk = find(diff([0 night_ind]) == 1);
dawn = find(diff([0 night_ind]) == -1);

if(RADSWDN(1,end)==0) 
    dawn = [dawn, length(time_new)];
end

% Plot plume widths
figure;
colors = lines(num_cases);
for case_idx = 1:num_cases
    plot(new_time_hr, average_widths(:,case_idx) ./ 1000, 'Color', colors(case_idx,:), 'LineWidth', 2);
    hold on;
end
xlim([0 48]);
ylim([0 250]);
current_ylim = ylim;
for i=1:length(dusk)
    patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
        [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
end
xlabel('Time [hours]');
ylabel('Plume Width [km]');
grid on;
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
pbaspect([2.0 1 1]);


%legend(NEPac_array, 'NumColumns', 1, 'Location','eastoutside','Interpreter','none');


% Plot plume widths
figure;
colors = lines(length(low_CF));
for i=1:length(low_CF)
plot(new_time_hr, average_widths(:,low_CF(i)) ./ 1000, 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
end
xlim([0 48]);
ylim([0 250]);
current_ylim = ylim;
for i=1:length(dusk)
    patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
        [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
end
xlabel('Time [hours]');
ylabel('Plume Width [km]');
title('Median Cloud Fraction below 50%','FontName', 'Georgia');
grid on;
num_low_CF = numel(low_CF);
text(5, max(ylim)-20, ['# of cases: ', num2str(num_low_CF)], 'FontSize', 18, 'FontName', 'Georgia');
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
pbaspect([2.0 1 1]);

%legend(NEPac_array{low_CF}, 'NumColumns', 1, 'Location','eastoutside','Interpreter','none');

figure;
colors = lines(length(mid_CF));
for i=1:length(mid_CF)
plot(new_time_hr, average_widths(:,mid_CF(i)) ./ 1000, 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
end
xlim([0 48]);
ylim([0 140]);
current_ylim = ylim;
for i=1:length(dusk)
    patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
        [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
end
pbaspect([1.0 1 1]);
xlabel('Time [hours]');
ylabel('Plume Width [km]');
title('Median Cloud Fraction above 50%','FontName', 'Georgia');
grid on;
num_mid_CF = numel(mid_CF);
text(5, max(ylim)-20, ['# of cases: ', num2str(num_mid_CF)], 'FontSize', 18, 'FontName', 'Georgia');
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
pbaspect([2.0 1 1]);

%legend(NEPac_array{mid_CF}, 'NumColumns', 1, 'Location','eastoutside','Interpreter','none');

% Plot plume widths
% figure;
% colors = lines(length(high_CF));
% for i=1:length(high_CF)
% plot(new_time_hr, average_widths(:,high_CF(i)) ./ 1000, 'Color', colors(i,:), 'LineWidth', 2);
%     hold on;
% end
% xlim([0 48]);
% ylim([0 180]);
% current_ylim = ylim;
% for i=1:length(dusk)
%     patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%         [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% end
% xlabel('Time (hours)');
% ylabel('Plume Width (km)');
% title('Median Cloud Fraction above 75%','FontName', 'Georgia');
% grid on;
% num_high_CF = numel(high_CF);
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_high_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% legend(NEPac_array{high_CF}, 'NumColumns', 1, 'Location','eastoutside','Interpreter','none');

%plot spreading rate
for case_idx = 1:num_cases
    spread_rates_smooth(:,case_idx) = smooth(spreading_rates(:,case_idx),10, 'moving');
end

figure;
colors = lines(num_cases);
for case_idx = 1:num_cases
    plot(new_time_hr(1:end-1), spread_rates_smooth(:, case_idx), 'Color', colors(case_idx, :), 'LineWidth', 2);
    hold on;
end
xlim([0 48]);
ylim([0 7]);
current_ylim = ylim;
for i=1:length(dusk)
    patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
        [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
end
grid on;
xlabel('Time [hours]');
ylabel('Spreading Rate [km/hr]');
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
pbaspect([2.0 1 1]);
%legend(NEPac_array, 'NumColumns', 1, 'Location','eastoutside','Interpreter', 'none');

%plot spreading rate
for i = 1:length(mid_CF)
    spread_rates_smooth_mid(:,i) = smooth(spreading_rates(:,mid_CF(i)),10, 'moving');
end

figure;
colors = lines(length(mid_CF));
for i = 1:length(mid_CF)
    plot(new_time_hr(1:end-1), spread_rates_smooth_mid(:, i), 'Color', colors(i, :), 'LineWidth', 2);
    hold on;
end
xlim([0 48]);
ylim([0 6]);
current_ylim = ylim;
for i=1:length(dusk)
    patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
        [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
end
grid on;
pbaspect([2.0 1 1]);
num_mid_CF = numel(mid_CF);
text(35, max(ylim)-20, ['# of cases: ', num2str(num_mid_CF)], 'FontSize', 14, 'FontName', 'Georgia');
xlabel('Time [hours]');
ylabel('Spreading Rate [km/hr]');
title('Median Cloud Fraction above 50%','FontName', 'Georgia','LineWidth',2);
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
pbaspect([2.0 1 1]);
%legend(NEPac_array{mid_CF}, 'NumColumns', 1, 'Location','eastoutside','Interpreter', 'none');

% plot spreading rate
for i = 1:length(low_CF)
    spread_rates_smooth_low(:,i) = smooth(spreading_rates(:,low_CF(i)),10, 'moving');
end

figure;
colors = lines(length(low_CF));
for i = 1:length(low_CF)
    plot(new_time_hr(1:end-1), spread_rates_smooth_low(:, i), 'Color', colors(i, :), 'LineWidth', 2);
    hold on;
end
xlim([0 48]);
ylim([0 7]);
current_ylim = ylim;
for i=1:length(dusk)
    patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
        [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
end
grid on;
pbaspect([1.0 1 1]);
num_low_CF = numel(low_CF);
text(35, max(ylim)-20, ['# of cases: ', num2str(num_low_CF)], 'FontSize', 14, 'FontName', 'Georgia');
xlabel('Time [hours]');
ylabel('Spreading Rate [km/hr]');
title('Median Cloud Fraction below 50%','FontName', 'Georgia');
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
pbaspect([2.0 1 1]);
%legend(NEPac_array{low_CF}, 'NumColumns', 1, 'Location','eastoutside','Interpreter', 'none');

% plot median plume width and shade 25th-75th percentile region
median_widths = median(average_widths(:,mid_CF), 2)./1000.0;
P25_widths = prctile(average_widths(:,mid_CF), 25, 2)./1000.0;
P75_widths = prctile(average_widths(:,mid_CF), 75, 2)./1000.0;

x_patch = [new_time_hr, fliplr(new_time_hr)];
y_patch = [P25_widths', fliplr(P75_widths')];

initial_width = 2; % Initial width in km
const_rate = 1.5; % Spreading rate in km/hr
total_time = 48; % Total time in hours

% Time vector from 0 to total_time in hours
time1 = 1:total_time; % From 0 to 48 hours

% Calculate the width as a function of time
width_const = initial_width + const_rate * time1;

median_widths_s = smooth(median_widths, 10);
P25_widths_s = smooth(P25_widths, 10);
P75_widths_s = smooth(P75_widths, 10);

%day1 = mean(median_widths_s(1:72));
 
figure;
for i = 1:length(dusk)
plot(new_time_hr, median_widths_s, 'black', 'LineWidth', 3);
hold on;
%plot(new_time_hr, P25_widths, 'r', 'LineWidth', 2,'LineStyle','--');
%plot(new_time_hr, P75_widths, 'r', 'LineWidth', 2,'LineStyle','--');
plot(time1, width_const,'r','LineWidth', 3,'LineStyle','-.')
fill([new_time_hr; flipud(new_time_hr)], [P25_widths_s; flipud(P75_widths_s)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
xlabel('Time [hours]','FontName', 'Georgia');
ylabel('Plume Width [km]','FontName', 'Georgia');
xlim([0 48]);
ylim([0.0 100]);
% Set tick label font style
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
title('Median Cloud Fraction above 50%','FontName', 'Georgia');
current_ylim = ylim;
patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
    [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
pbaspect([2.0 1 1]);
grid on;
legend('median width', '1.5 km/hr','Location','southeast');
end

%plot median plume width and shade 25th-75th percentile region
median_widths = median(average_widths(:,low_CF), 2)./1000.0;
P25_widths = prctile(average_widths(:,low_CF), 25, 2)./1000.0;
P75_widths = prctile(average_widths(:,low_CF), 75, 2)./1000.0;

x_patch = [new_time_hr, fliplr(new_time_hr)];
y_patch = [P25_widths', fliplr(P75_widths')];

initial_width = 2; % Initial width in km
const_rate = 1.5; % Spreading rate in km/hr
total_time = 48; % Total time in hours

% Time vector from 0 to total_time in hours
time1 = 1:total_time; % From 0 to 48 hours

% Calculate the width as a function of time
width_const = initial_width + const_rate * time1;

median_widths_s = smooth(median_widths, 10);
P25_widths_s = smooth(P25_widths, 10);
P75_widths_s = smooth(P75_widths, 10);

figure;
for i = 1:length(dusk)
plot(new_time_hr, median_widths_s, 'black', 'LineWidth', 3);
hold on;
%plot(new_time_hr, P25_widths, 'r', 'LineWidth', 2,'LineStyle','--');
%plot(new_time_hr, P75_widths, 'r', 'LineWidth', 2,'LineStyle','--');
plot(time1, width_const,'r','LineWidth', 3,'LineStyle','-.')
fill([new_time_hr; flipud(new_time_hr)], [P25_widths_s; flipud(P75_widths_s)], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
xlabel('Time [hours]','FontName', 'Georgia');
ylabel('Plume Width [km]','FontName', 'Georgia');
xlim([0 48]);
ylim([0.0 100]);
title('Median Cloud Fraction below 50%','FontName', 'Georgia');
current_ylim = ylim;
patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
    [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
pbaspect([2.0 1 1]);
grid on;
legend('median width', '1.5 km/hr','Location','southeast');
end
% 
% %plot median spreading rate and shade 25th-75th percentile region
% median_sr = median(spread_rates_smooth, 2);
% P25_widths = prctile(spread_rates_smooth, 25, 2);
% P75_widths = prctile(spread_rates_smooth, 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr, 'black', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 

median_sr_high = median(spread_rates_smooth(:,mid_CF), 2);
P25_widths = prctile(spread_rates_smooth_mid(:,:), 25, 2);
P75_widths = prctile(spread_rates_smooth_mid(:,:), 75, 2);

median_sr_high = median(spread_rates_smooth(:,mid_CF), 2);
mean_day1_sr = mean(median_sr_high(1:72));
mean_day1_sr_ind = mean(spread_rates_smooth(1:72,mid_CF),1);
mean_day1_sr = mean(mean_day1_sr_ind);
mean_day2_sr_ind = mean(spread_rates_smooth(72:143,mid_CF),1);

delta_sr = spread_rates_smooth(1:72,mid_CF) - spread_rates_smooth(72:143,mid_CF);
mean_day2_sr = mean(mean_day2_sr_ind);

delta_spreadingrate = abs(mean_day1_sr_ind - mean_day2_sr_ind);

mean_day2_sr = mean(median_sr_high(72:143));

x_patch = [new_time_hr, fliplr(new_time_hr)];
y_patch = [P25_widths', fliplr(P75_widths')];

num_mid_CF = numel(mid_CF);

median_sr_high_s = smooth(median_sr_high,10);
P25_widths_s = smooth(P25_widths,10);
P75_widths_s = smooth(P75_widths,10);

figure;
for i = 1:length(dusk)
plot(new_time_hr(1:end-1), median_sr_high_s, 'black', 'LineWidth', 3);
hold on;
%plot(new_time_hr(1:end-1), P25_widths, 'b', 'LineWidth', 2,'LineStyle','--');
%plot(new_time_hr(1:end-1), P75_widths, 'b', 'LineWidth', 2,'LineStyle','--');
fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths_s; flipud(P75_widths_s)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
xlabel('Time [hours]','FontName', 'Georgia');
ylabel('Spreading Rate [km/hr]','FontName', 'Georgia');
title('Median Cloud Fraction above 50%','FontName', 'Georgia');
xlim([0 48]);
ylim([0.0 4]);
set(gca, 'FontSize', 14, 'FontName', 'Georgia', 'LineWidth', 2);
current_ylim = ylim;
%text(35, max(ylim)-1, ['# of cases: ', num2str(num_mid_CF)], 'FontSize', 14, 'FontName', 'Georgia');
patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
    [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
pbaspect([2.0 1 1]);
grid on;
end

median_sr_high = median(spread_rates_smooth(:,low_CF), 2);
P25_widths = prctile(spread_rates_smooth_low(:,:), 25, 2);
P75_widths = prctile(spread_rates_smooth_low(:,:), 75, 2);

x_patch = [new_time_hr, fliplr(new_time_hr)];
y_patch = [P25_widths', fliplr(P75_widths')];

median_sr_high_s = smooth(median_sr_high,10);
P25_widths_s = smooth(P25_widths,10);
P75_widths_s = smooth(P75_widths,10);

num_low_CF = numel(low_CF);

figure;
for i = 1:length(dusk)
plot(new_time_hr(1:end-1), median_sr_high_s, 'black', 'LineWidth', 3);
hold on;
%plot(new_time_hr(1:end-1), P25_widths, 'g', 'LineWidth', 2,'LineStyle','--');
%plot(new_time_hr(1:end-1), P75_widths, 'g', 'LineWidth', 2,'LineStyle','--');
fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths_s; flipud(P75_widths_s)], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
xlabel('Time [hours]','FontName', 'Georgia');
ylabel('Spreading Rate [km/hr]','FontName', 'Georgia');
title('Median Cloud Fraction below 50%','FontName', 'Georgia');
xlim([0 48]);
ylim([0.0 5]);
set(gca, 'FontSize', 14, 'FontName', 'Georgia','LineWidth',2);
current_ylim = ylim;
%text(35, max(ylim)-1, ['# of cases: ', num2str(num_low_CF)], 'FontSize', 14, 'FontName', 'Georgia');
patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
    [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
pbaspect([2.0 1 1]);
grid on;
end
% 
% median_sr_mid = median(spread_rates_smooth(:,mid_CF), 2);
% P25_widths = prctile(spread_rates_smooth(:,mid_CF), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,mid_CF), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_mid_CF = numel(mid_CF);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_mid, 'r', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median Cloud Fraction between 25% and 75%','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_mid_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% median_sr_low = median(spread_rates_smooth(:,low_CF), 2);
% P25_widths = prctile(spread_rates_smooth(:,low_CF), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,low_CF), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_low_CF = numel(low_CF);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low, 'g', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median Cloud Fraction below 25%','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_low_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low, 'g', 'LineWidth', 3);
% hold on;
% plot(new_time_hr(1:end-1), median_sr_mid, 'r', 'LineWidth', 3);
% plot(new_time_hr(1:end-1), median_sr_high, 'b', 'LineWidth', 3);
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% legend('Low CF','Mid CF','High CF', 'Location','northeast');
% end
% 
% %plot median spreading rate and shade 25th-75th percentile region
% %plot LWP-binned results
% median_sr_high_lwp = median(spread_rates_smooth(:,high_LWP), 2);
% P25_widths = prctile(spread_rates_smooth(:,high_LWP), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,high_LWP), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_high_CF = numel(high_LWP);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_high_lwp, 'b', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'b', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'b', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median LWP above 50 g/m^2','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_high_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% median_sr_mid_lwp = median(spread_rates_smooth(:,mid_LWP), 2);
% P25_widths = prctile(spread_rates_smooth(:,mid_LWP), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,mid_LWP), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_mid_CF = numel(mid_LWP);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_mid_lwp, 'r', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median LWP between 10 and 50 g/m^2','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_mid_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% median_sr_low_lwp = median(spread_rates_smooth(:,low_LWP), 2);
% P25_widths = prctile(spread_rates_smooth(:,low_LWP), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,low_LWP), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_low_CF = numel(low_LWP);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low_lwp, 'g', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median LWP below 10 g/m^2','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_low_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low_lwp, 'g', 'LineWidth', 3);
% hold on;
% plot(new_time_hr(1:end-1), median_sr_mid_lwp, 'r', 'LineWidth', 3);
% plot(new_time_hr(1:end-1), median_sr_high_lwp, 'b', 'LineWidth', 3);
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% legend('Low LWP','Mid LWP','High LWP', 'Location','northeast');
% end
% 
% %plot median spreading rate and shade 25th-75th percentile region
% %plot NAc-binned results
% median_sr_high_nac = median(spread_rates_smooth(:,high_NAc), 2);
% P25_widths = prctile(spread_rates_smooth(:,high_NAc), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,high_NAc), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_high_CF = numel(high_NAc);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_high_nac, 'b', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'b', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'b', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median NAc above 80 #/mg','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_high_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% median_sr_mid_nac = median(spread_rates_smooth(:,mid_NAc), 2);
% P25_widths = prctile(spread_rates_smooth(:,mid_NAc), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,mid_NAc), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_mid_CF = numel(mid_NAc);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_mid_nac, 'r', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median NAc between 20 and 80 #/mg','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_mid_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% median_sr_low_nac = median(spread_rates_smooth(:,low_NAc), 2);
% P25_widths = prctile(spread_rates_smooth(:,low_NAc), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,low_NAc), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_low_CF = numel(low_NAc);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low_nac, 'g', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median NAc below 20 #/mg','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_low_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low_nac, 'g', 'LineWidth', 3);
% hold on;
% plot(new_time_hr(1:end-1), median_sr_mid_nac, 'r', 'LineWidth', 3);
% plot(new_time_hr(1:end-1), median_sr_high_nac, 'b', 'LineWidth', 3);
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% legend('Low NAc','Mid NAc','High NAc', 'Location','northeast');
% end
% 
% %plot median spreading rate and shade 25th-75th percentile region
% %plot NAc-binned results
% median_sr_high_zinv = median(spread_rates_smooth(:,high_ZINV), 2);
% P25_widths = prctile(spread_rates_smooth(:,high_ZINV), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,high_ZINV), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_high_CF = numel(high_ZINV);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_high_zinv, 'b', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'b', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'b', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median BLD above 1.2 km','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_high_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% median_sr_mid_zinv = median(spread_rates_smooth(:,mid_ZINV), 2);
% P25_widths = prctile(spread_rates_smooth(:,mid_ZINV), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,mid_ZINV), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_mid_CF = numel(mid_ZINV);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_mid_zinv, 'r', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'r', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median BLD between 0.8 and 1.2 km','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_mid_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% median_sr_low_zinv = median(spread_rates_smooth(:,low_ZINV), 2);
% P25_widths = prctile(spread_rates_smooth(:,low_ZINV), 25, 2);
% P75_widths = prctile(spread_rates_smooth(:,low_ZINV), 75, 2);
% 
% x_patch = [new_time_hr, fliplr(new_time_hr)];
% y_patch = [P25_widths', fliplr(P75_widths')];
% 
% num_low_CF = numel(low_ZINV);
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low_zinv, 'g', 'LineWidth', 2);
% hold on;
% plot(new_time_hr(1:end-1), P25_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% plot(new_time_hr(1:end-1), P75_widths, 'g', 'LineWidth', 2,'LineStyle','--');
% fill([new_time_hr(1:end-1); flipud(new_time_hr(1:end-1))], [P25_widths; flipud(P75_widths)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Confidence interval fill
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% title('Median BLD below 0.8 km','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% text(35, max(ylim)-1, ['# of cases: ', num2str(num_low_CF)], 'FontSize', 14, 'FontName', 'Georgia');
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% end
% 
% figure;
% for i = 1:length(dusk)
% plot(new_time_hr(1:end-1), median_sr_low_zinv, 'g', 'LineWidth', 3);
% hold on;
% plot(new_time_hr(1:end-1), median_sr_mid_zinv, 'r', 'LineWidth', 3);
% plot(new_time_hr(1:end-1), median_sr_high_zinv, 'b', 'LineWidth', 3);
% xlabel('Time (hours)','FontName', 'Georgia');
% ylabel('Spreading Rate (km)','FontName', 'Georgia');
% xlim([0 48]);
% ylim([0.0 5]);
% set(gca, 'FontSize', 12, 'FontName', 'Georgia');
% current_ylim = ylim;
% patch([time_new(dusk(i)) time_new(dusk(i)) time_new(dawn(i)) time_new(dawn(i)) time_new(dusk(i))],[current_ylim(1) current_ylim(2) current_ylim(2) current_ylim(1) current_ylim(1)], ...
%     [0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor', 'none');
% pbaspect([1.0 1 1]);
% grid on;
% legend('Low BLD','Mid BLD','High BLD', 'Location','northeast');
% end