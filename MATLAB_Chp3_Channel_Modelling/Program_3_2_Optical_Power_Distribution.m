P_total=1; % Total transmitted power

rho=0.8; % Reflection coefficient of the walls

% Room dimensions (lx, ly, lz) in meters
lx=5; 
ly=5; 
lz=2.15; % Height of the room

% Number of grid points along x, y, and z axes
Nx=lx*30; 
Ny=ly*30; 
Nz=round(lz*30);

% Calculation of grid area for reflections
dA=lz*ly/(Ny*Nz); % Area element of the surface

FOV=60; % Field of view (FOV) of the receiver in degrees

% Create grid points along x, y, and z axes
x=linspace(-lx/2,lx/2,Nx); % X-axis grid points
y=linspace(-ly/2,ly/2,Ny); % Y-axis grid points
z=linspace(-lz/2,lz/2,Nz); % Z-axis grid points

% Create a 3D mesh grid for receiver positions on the lower surface (ZR = -lz/2)
[XR,YR,ZR] = meshgrid(x,y,-lz/2); 

% Transmitter position at the center of the room's ceiling
TP1=[0 0 lz/2]; % Transmitter position (x=0, y=0, z=lz/2)

%%%%%%%%%%%%%%% Calculation for wall 1 (North face) %%%%%%%%%%%%%%%%%%

% Loop over x and y grid points for the receiver plane
for ii=1:Nx
    for jj=1:Ny
        RP=[x(ii) y(jj) -lz/2]; % Receiver position vector
        h1(ii,jj)=0; % Initialize reflection from North face (wall 1)

        % Loop over y and z grid points for wall 1
        for kk=1:Ny
            for ll=1:Nz
                WP1=[-lx/2 y(kk) z(ll)]; % Point of incidence on wall 1 (North face)

                % Distance from the transmitter to the point of incidence (WP1) on wall 1
                D1=sqrt(dot(TP1-WP1,TP1-WP1)); 
                
                % Cosine of the angle between the transmitter and the wall point (vertical)
                cos_phi= abs(WP1(3)-TP1(3))/D1; 

                % Cosine of the angle for the reflection from the wall (horizontal)
                cos_alpha = abs(TP1(1)-WP1(1))/D1; 

                % Distance from the point of incidence (WP1) to the receiver (RP)
                D2=sqrt(dot(WP1-RP,WP1-RP)); 

                % Cosine of the angle between the wall point and receiver (horizontal)
                cos_beta=abs(WP1(1)-RP(1))/D2;

                % Cosine of the angle between the wall point and receiver (vertical)
                cos_psi=abs(WP1(3)-RP(3))/D2;

                % Check if the reflection falls within the receiver's field of view (FOV)
                if abs(acosd(cos_psi))<=FOV 
                    % Update the channel gain (h1) for the reflection from wall 1
                    h1(ii,jj)=h1(ii,jj)+(m+1)*Adet*rho*dA*...
                    cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%% Repeat the above calculation for remaining walls %%%%%%%%%%%%%%%%%
% The code to calculate the channel gain from other walls (h2, h3, and h4) should follow a similar logic.
% Wall 2 (South face), Wall 3 (East face), and Wall 4 (West face) calculations should be added here.
% For each wall, adjust the position vectors (WP1) and the corresponding angles.

% Calculate total received power from all reflections (h1, h2, h3, h4)
P_rec_A1=(h1)*P_total.*Ts.*G_Con; % Total received power including reflections from all walls
