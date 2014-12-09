clear
close all

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction  [m]
dy = 0.1e-3;        % grid point spacing in the y direction  [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% testing
c0 = 1500;
rho0 = 1.2;

% define the properties of the propagation medium
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 0.75;

% define the ratio between the barrier and background sound speed and density
barrier_scale = 20;

% create the time array using the barrier sound speed
t_end = 40e-6;                % [s]
CFL = 0.5;                    % Courant–Friedrichs–Lewy number
kgrid.t_array = makeTime(kgrid, c0*barrier_scale, CFL, t_end);

% create a mask of a barrier with a slit
slit_thickness = 2;                     % [grid points]
slit_width = 20;                        % [grid points]
slit_x_pos = Nx - Nx/4;                 % [grid points]
slit_offset = Ny/2 - slit_width/2 - 1;  % [grid points]
slit_mask = zeros(Nx, Ny);
slit_mask(slit_x_pos:slit_x_pos + slit_thickness, 1:1 + slit_offset) = 1;
slit_mask(slit_x_pos:slit_x_pos + slit_thickness, end - slit_offset:end) = 1;

% assign the slit to the properties of the propagation medium
medium.sound_speed = c0*ones(Nx, Ny);
medium.density = rho0*ones(Nx, Ny);
medium.sound_speed(slit_mask == 1) = barrier_scale*c0;
medium.density(slit_mask == 1) = barrier_scale*rho0;

% define a single source point
source.p_mask = zeros(Nx, Ny);
source.p_mask(end, Ny*3/4) = 1;

% define a time varying sinusoidal source
source_freq = 0.25e6;   % [Hz]
source_mag = 100;         % [Pa]
source.p = source_mag*sin(2*pi*source_freq*kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

% define four sensor points centered about source.p0
sensor_radius = 40; % [grid points]
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2 + sensor_radius, Ny/2) = 1;
sensor.mask(Nx/2 - sensor_radius, Ny/2) = 1;
sensor.mask(Nx/2, Ny/2 + sensor_radius) = 1;
sensor.mask(Nx/2, Ny/2 - sensor_radius) = 1;

% set the acoustic variables that are recorded
sensor.record = {'p', 'u'};

% set the input options
input_args = {'PMLInside', false, 'PlotPML', false, ...
    'DisplayMask', slit_mask, 'DataCast', 'single'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

plot(kgrid.t_array, sensor_data.p);