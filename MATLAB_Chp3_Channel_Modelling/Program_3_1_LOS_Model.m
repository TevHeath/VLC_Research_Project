%Gathered from Program 3.1 - LOS Channel Gain
%Variables:
theta = 70; % semi-angle at half power
m = -log10(2)/log10(cosd(theta)); %Lambertian order of emission
P_total = 20;
Adet=1e-4; %detector physical area of a PD

%% Optics parameters
Ts = 1;
index = 1.5;
FOV = 60*pi/180;
G_Con = (index^2)/sin(FOV); %gain in an optical concentrator

%% Room Dimension
lx = 5; ly = 5; lz = 3; %Room dimension in meter
h = 2.15; %the distance berween the soruce and receiver plane
Nx = lx*20; Ny = ly*20; %number of grid in the receiver plane

XT = 0; YT = 0;
x=-lx/2:lx/Nx:lx/2;
y=-ly/2:ly/Ny:ly/2;
[XR,YR]=meshgrid(x,y);

D1 = sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2);
% distance verctor from source 1
cosphi_A1=h./D1; % angle vector

%%
H_A1=(m+1)*Adet.*cosphi_A1.^(m+1)./(2*pi.*D1.^2);
% channel DC gain for source 1
P_rec=P_total.*H_A1.*Ts.*G_Con; % received power from source 1;
P_rec_dBm=10*log10(P_rec);

meshc(x,y,P_rec_dBm);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Received power (dBm)');
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(P_rec_dBm)) max(max(P_rec_dBm))]);
