C=3e8*1e-9; % Speed of light in ns (converting from meters per second to nanoseconds)
theta=70; % Semi-angle at half power (in degrees)
m=-log10(2)/log10(cosd(theta)); % Lambertian order of emission, calculated using the angle
P_total=1; % Total normalized transmitted power
Adet=1e-4; % Detector physical area of a photodiode (PD)
rho=0.8; % Reflection coefficient of the walls
Ts=1; % Gain of an optical filter
index=1.5; % Refractive index of a lens at a PD
FOV=60; % Field of view (FOV) of the receiver in degrees
G_Con=(index^2)/(sind(FOV).^2); % Gain of the optical concentrator

%% Room dimensions (lx, ly, lz) in meters
lx=5; 
ly=5; 
lz=3-0.85; % Height of the room minus the receiver height
Nx=lx*30; Ny=ly*30; Nz=round(lz*30); % Number of grid points in x, y, and z directions

dA=lz*ly/(Ny*Nz); % Calculation of grid area (used for reflections)
x=-lx/2:lx/Nx:lx/2; % X-axis grid points
y=-ly/2:ly/Ny:ly/2; % Y-axis grid points
z=-lz/2:lz/Nz:lz/2; % Z-axis grid points

%% First transmitter setup
TP1=[0 0 lz/2]; % Position of the first transmitter in the room
TPV=[0 0 -1]; % Transmitter position vector
RPV=[0 0 1]; % Receiver position vector
WPV1=[1 0 0]; % Wall position vector for wall 1
WPV2=[0 1 0]; % Wall position vector for wall 2
WPV3=[-1 0 0]; % Wall position vector for wall 3
WPV4=[0 -1 0]; % Wall position vector for wall 4

delta_t=1/2; % Time resolution in ns, used in the form of 1/2^m

% Loop over x and y grid points to calculate channel gain
for ii=1:Nx+1
    for jj=1:Ny+1
        RP=[x(ii) y(jj) -lz/2]; % Receiver position in the grid
        t_vector=0:25/delta_t; % Time vector in ns
        h_vector=zeros(1,length(t_vector)); % Initialize channel gain vector
        
        % LOS (Line-of-Sight) channel gain calculation
        D1=sqrt(dot(TP1-RP,TP1-RP)); % Distance between transmitter and receiver
        cosphi= lz/D1; % Cosine of the angle between the transmitter and receiver
        tau0=D1/C; % Time delay of the LOS path
        index=find(round(tau0/delta_t)==t_vector); % Find the index in the time vector for tau0
        if abs(acosd(cosphi))<=FOV % Check if the angle is within the FOV
            h_vector(index)=h_vector(index)+(m+1)*Adet.*cosphi.^(m+1)./(2*pi.*D1.^2); % Update channel gain
        end
        
        % Reflection from the first wall (loop over y and z grid points)
        count=1; % Counter for wall reflections
        for kk=1:Ny+1
            for ll=1:Nz+1
                WP1=[-lx/2 y(kk) z(ll)]; % Position on the first wall
                D1=sqrt(dot(TP1-WP1,TP1-WP1)); % Distance from the transmitter to the wall point
                cos_phi= abs(WP1(3)-TP1(3))/D1; % Cosine of the angle between transmitter and wall
                cos_alpha = abs(TP1(1)-WP1(1))/D1; % Cosine of the angle for reflection from the wall
                D2=sqrt(dot(WP1-RP,WP1-RP)); % Distance from the wall point to the receiver
                cos_beta=abs(WP1(1)-RP(1))/D2; % Cosine of the angle between the wall and receiver
                cos_psi=abs(WP1(3)-RP(3))/D2; % Cosine of the angle for the receiver
                tau1=(D1+D2)/C; % Time delay for the reflected path
                index=find(round(tau1/delta_t)==t_vector); % Find the index in the time vector for tau1
                if abs(acosd(cos_psi))<=FOV % Check if the angle is within the FOV
                    % Update channel gain for the reflected path
                    h_vector(index)=h_vector(index)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2); 
                end
                count=count+1; % Increment the counter for reflections
            end
        end
        
        % Similar calculations for reflections from remaining walls can be added here
    end
    
    % Time vector scaling
    t_vector=t_vector*delta_t; % Scale time vector by delta_t
    mean_delay(ii,jj)=sum((h_vector).^2.*t_vector)/sum(h_vector.^2); % Calculate mean delay
    Drms(ii,jj)= sqrt(sum((t_vector-mean_delay(ii,jj)).^2.*h_vector.^2)/sum(h_vector.^2)); % Calculate RMS delay spread
end

%% Plot RMS delay spread
surf(x,y, Drms); % 3D surface plot of RMS delay spread

% Optional: Uncomment to plot mean delay instead
% surf(x,y,mean_delay);

% Set plot axes limits
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(Drms)) max(max(Drms))]);
