%% parameters
K = 864; %m/day     10;  %Transmissivity = 1m/day x 10m = 10m2/day
S = 0.002;%storativity (S = Ss*b + Sy) 
b = 0.25; %meter %aquifer thickness




%% load pumping data

%pumping data in L/week

pump = load('pump_data.csv');
well_rates = [50 100 200 0];
%columns:
%landunit_id  well_type  hours   wellDepth
for i=1:size(pump,1);
    Qwell(i) = -pump(i,3).*well_rates(pump(i,2)); %L/day
    Qwell(i) = Qwell(i)/1000; %m3/day
end
Qwell =Qwell*1e4; %scaling for "effect" (not real)

pump_cell = pump(:,1);

%% load ground elevation data
g_elev = load('elevation.csv');

%smooth it out to look less sharply varying:
gg = reshape(g_elev,100,100);
g0 = ones(120)*mean(g_elev); 
g0(11:110,11:110) = gg;
g0(:,1:10) = 900;
gg = conv2(g0,ones(10)/100,'same');
gg = gg(11:110,11:110);
g_elev = gg(:);

%evaluate depth to well (pump setting) based on ground elevation and well
%depth:
well_elev = g_elev(pump(:,1)) - pump(:,4);



%% grid and time step setup: 
dx = 135;    %requires square cells (dx=dy)

%number of cells in each direction (not required to be square)
Nx = 100;
Ny = 100;

dt = 10;
tmin = 0;
tmax = 380;

err_tol = 5;  %used for auto scaling time to reduce errors if needed


%% Set boundary conditions for cells

%initialize matricies to store prescribed head and flux values
bH = zeros(Ny,Nx);
bq = zeros(Ny,Nx);
cells = ones(Ny,Nx);

% set no flow everywhere except at bottom with fixed head zone and top with
% fixed flux

cells(:,[1,Nx]) = 3; %no flow boundary
cells([1,Ny],:) = 3;

%bq(1,:) = 0.0001;              %flux boundaries for inflows to watershed
%bq(1:50,[1,Nx])=0.0001;

%bq(25:75,1) = -4.8e-3;%0.04;      %set inflow for flux boundary

cells(25:75,1) = 2;
bH(25:75,1) = 1500;              %fixed boundary
%bH = g_elev(25:75,1)*0.9;      %boundary that matches elevation

%take out the corners (not constrained)
cells([1,Ny],[1,Nx]) = 2;
bH([1,Ny],[1,Nx]) = 0;


%% set initial values for cells:
Ho =1500*ones(Nx*Ny,1);
%Ho = g_elev - 2;% mean(g_elev(25:75,1))*1 + (g_elev-mean(g_elev(:)))*0.95;;%- (g_elev-50)*0.2;  %use this one for an
%initial value that looks like the surface


%% Set source terms
R = g_elev.*(g_elev>1400)/1e6;


%% RUN MODEL

gwfd_homo_trans_v5;

%get rid of corners (not properly constrained):plot(H(1000,:))
H(H<1e-6 & H>-1e-6) = nan;


%% data vis

plot(H(1000,:));
pause; 

for i=1:10:length(t); 
    surf(flipud(reshape(H(:,i),100,100)));
    shading flat;
    title(i);
    axis([0 100 0 100 1230 1550]);
    colorbar;
    pause(0.1);
end