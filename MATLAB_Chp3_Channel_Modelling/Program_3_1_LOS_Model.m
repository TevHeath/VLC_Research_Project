% Gathered from Program 3.1 - LOS Channel Gain

%% Variables:
theta = 70; % semi-angle at half power for the light source
m = -log10(2)/log10(cosd(theta)); % Lambertian order of emission (defines how light intensity decreases with angle)
P_total = 20; % total transmitted power in watts
Adet = 1e-4; % physical area of the photodetector (PD) in square meters

%% Optics parameters:
Ts = 1; % optical system transmission coefficient (assumed ideal, hence 1)
index = 1.5; % refractive index of the concentrator (material-dependent)
FOV = 60*pi/180; % field of view of the receiver in radians
G_Con = (index^2) / sin(FOV); % gain of the optical concentrator, depends on the refractive index and FOV

%% Room Dimension:
lx = 10; ly = 10; lz = 10; % room dimensions in meters (length, width, height)
h = 5.15; % distance between the source and receiver plane (in meters)
Nx = lx * 20; Ny = ly * 20; % number of grid points in the receiver plane for X and Y dimensions

XT = 0; YT = 0; % coordinates of the transmitter (light source) at the center of the room
x = -lx/2 : lx/Nx : lx/2; % X-axis grid for the receiver plane
y = -ly/2 : ly/Ny : ly/2; % Y-axis grid for the receiver plane
[XR, YR] = meshgrid(x, y); % creates a meshgrid for the X and Y coordinates in the receiver plane

D1 = sqrt((XR - XT(1,1)).^2 + (YR - YT(1,1)).^2 + h^2); 
% calculates the Euclidean distance from the transmitter to each point on the receiver grid
cosphi_A1 = h ./ D1; % calculates the cosine of the angle of incidence for each point

%% Channel Gain Calculation:
H_A1 = (m + 1) * Adet .* cosphi_A1.^(m + 1) ./ (2 * pi .* D1.^2); 
% calculates the channel DC gain for line-of-sight (LOS) propagation from the light source to each point

P_rec = P_total .* H_A1 .* Ts .* G_Con; % received optical power at each point on the receiver plane
P_rec_dBm = 10 * log10(P_rec); % converts the received power to dBm for easier visualization

%% Visualization:
meshc(x, y, P_rec_dBm); % generates a 3D mesh plot with contours for the received power distribution
xlabel('X (m)'); % labels the X-axis
ylabel('Y (m)'); % labels the Y-axis
zlabel('Received power (dBm)'); % labels the Z-axis (received power in dBm)
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(P_rec_dBm)) max(max(P_rec_dBm))]); 
% sets the axis limits based on the room dimensions and the range of received power values
