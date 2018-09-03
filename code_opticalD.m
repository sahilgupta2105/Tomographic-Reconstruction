%%%%%CONSTANTS%%%%%%

%Uniform pressure field//////Units= atm
P= 1;

%Line strength for target frequency and temperature//// Units= cm-2/atm
S= 3.16e-02;

%Target center frequency///// Units= cm-1
mu_o= 4988.655;

%Uniform Temparature and molecular mass of gas
T=296;
M=44;

%broadening coefficients for air and gas////// Units= cm-1/atm
bc_air= 0.077;
bc_self= 0.104;

%%%%%%MESH%%%%%
%Mesh size

N=50;

n_ray=2*(N-1);

%% coordinates%%
x= linspace(0,1,N);
y= linspace(0,1,N);

%Plot the mesh grid
[xx,yy]=meshgrid(x,y);
figure, plot(xx,yy,'g');
hold on;
plot(yy,xx,'g');
set(gca,'Color',[0,0,0])
hold on;

%Defing array for concentration of species at cell center
conc = zeros(N-1,N-1);

dx= 1/N;

%coordinates of cell centers, lasers and detectors
x_cell = zeros(N-1,N-1);
y_cell = zeros(N-1,N-1);


%top hat profile defined centered at x=0.25, y=0.25 for radius= 0,125 units
for i= 1:1:N-1
    for j=1:1:N-1
        
        x_cell(i,j) = (x(j)+x(j+1))/2;% coordinates of cell center
        y_cell(i,j) = (y(i)+y(i+1))/2;
        
        %         distance= sqrt((x_cell(i,j)-0.25)^2 + (y_cell(i,j)-0.35)^2);
        %         %condition for concentration distribution
        %         if distance<=0.125
        %             conc(i,j)= 0.18*(0.5 + 0.5*erf(0.4*(0.25-distance))); % top-hat function
        %         end
        %         distance= sqrt((x_cell(i,j)-0.55)^2 + (y_cell(i,j)-0.75)^2);
        %         %condition for concentration distribution
        %         if distance<=0.075
        %             conc(i,j)= 0.18*(0.5 + 0.5*erf(0.4*(0.25-distance))); % top-hat function
        %         end
        distance=sqrt(x_cell(i,j)^2+y_cell(i,j)^2);
        if distance<0.5
            cosG=0.25*(1-cos(2*pi*(x_cell(i,j)+0.5)^0.8))*(1-cos(2*pi*(y_cell(i,j)+0.5)^(2/3)));
        else
            cosG=0;
        end
        cosGauss=1.09*(0.3*cosG + 0.8*(exp(-(9*(x_cell(i,j)-0.2)^2+(6*(y_cell(i,j)-0.1)^2)))+exp(-(8*(x_cell(i,j)-0.7)^2+(30*(y_cell(i,j)-0.7)^2)))));
        conc(i,j)= 0.18*cosGauss;
        
    end
end

fprintf('Concentration loaded\n');
%coordinates for lasers

prompt='No. of lasers:';
n_laser= input(prompt);
prompt1='Values of angles in degrees (Range:[-90,90]):';
angles=zeros(n_laser,1);
for i=1:1:n_laser
    angles(i,1)=input(prompt1);
end
x_laser=zeros(n_ray,n_laser);
y_laser=zeros(n_ray,n_laser);
x_det=zeros(n_ray,n_laser);
y_det=zeros(n_ray,n_laser);

for i=1:1:n_laser
    
    [var1,var2,var3,var4]=laser(angles(i),N);
    x_laser(:,i)=var1(:);
    y_laser(:,i)=var2(:);
    x_det(:,i)=var3(:);
    y_det(:,i)=var4(:);
    
end

fprintf('Coordinates of lasers loaded\n');
for i=1:1:n_laser
    plot(x_laser(:,i),y_laser(:,i),'r.');
end
for i=1:1:n_laser
    plot(x_det(:,i),y_det(:,i),'b.');
end
hold off;

%Surface plot of concentration distribution

figure,surface(x_cell,y_cell,conc);

%Line-shape function: VOIGT PROFILE%%%%

gamma_G= 3.581e-07*mu_o*sqrt(T/M); % HWHM for Doppler broadening line
C_L= zeros(N-1,N-1); % Weight of Lorentzian profile
C_G= zeros(N-1,N-1); % Weight of Gaussian profile
gamma_V= zeros(N-1,N-1); % HWHM Voigt Profile

for i= 1:1:N-1
    for j=1:1:N-1
        gamma_L= P*(bc_self*conc(i,j)+bc_air*(1-conc(i,j))); % HWHM pressure broadening
        d= (gamma_L-gamma_G)/(gamma_L+gamma_G);
        C_L(i,j)= 0.68188 + 0.61293*d - 0.18384*d*d - 0.11568*d*d*d;
        C_G(i,j)= 0.32460 - 0.61825*d + 0.17681*d*d + 0.12109*d*d*d;
        gamma_V(i,j)= 0.5346*gamma_L + sqrt(0.2166*gamma_L*gamma_L + gamma_G*gamma_G); % HWHM Voigt
    end
end

fprintf('Begin: calculation of ray length\n');

%%%calculating length of each ray inside each box

laser_len = zeros((N-1)*(N-1),n_ray,n_laser);

for k=1:1:n_laser
    for i=1:1:n_ray
        for j=1:1:(N-1)*(N-1)
            
            %coordinates of the cell
            ji=ceil(j/(N-1));
            jj=rem(j-1,N-1)+1;
            
            %definig the  rectangle
            xbox = [x(jj) x(jj) x(jj+1) x(jj+1) x(jj)];
            ybox = [y(ji) y(ji+1) y(ji+1) y(ji) y(ji)];
            
            % Laser 1: ith ray
            x_line_1 = [x_laser(i,k) x_det(i,k)];
            y_line_1 = [y_laser(i,k) y_det(i,k)];
            
            %Intersection with ray from Laser 1
            [xi_1, yi_1] = polyxpoly(x_line_1, y_line_1, xbox, ybox);
            
            %Length of ray inside the rectangle for LASER 1
            if isempty(xi_1)
                laser_len(j,i,k) = 0;
            else
                laser_len(j,i,k) = sqrt((xi_1(2)-xi_1(1))^2 + (yi_1(2)-yi_1(1))^2);
            end
        end
    end
end

fprintf('End: Calculation of ray length\n');


%Line funciton

VP = zeros(N-1,N-1); % Voigt profile

%Laser intensities

I = zeros(100,n_ray,n_laser);

fprintf('Begin: Wavenumber loop\n');

%%%%%%wave number iteration main loop

mu=linspace(4986.168,4991.1458,100);
I_o= zeros(1,100);
mT = zeros(N-1,N-1);% Modulation index
%modulation amplitude

a=0.2;


wn= mu_o;



% Voigt profile function
for i= 1:1:N-1
    for j=1:1:N-1
        VP(i,j)= (C_L(i,j)/pi)*((gamma_V(i,j)/((wn-mu_o)^2 + (gamma_V(i,j))^2))) + (C_G(i,j)/gamma_V(i,j))*sqrt(log(2)/pi)*exp((-log(2)*(wn-mu_o)^2)/(gamma_V(i,j))^2);
    end
end
%calculating optical depth%%%% Number of views= 3 %

% Optical Depth is calculated for each laser

OD = zeros(n_ray,n_laser);

%Optical depth measurement
for l=1:1:n_laser
    for i=1:1:n_ray
        for j=1:1:N-1
            for k=1:1:N-1
                OD(i,l)= OD(i,l)+(-1/pi)*(P*conc(j,k)*S*VP(j,k)*laser_len((j-1)*(N-1)+k,i,l)); %%%%%%LASER 1
            end
        end
    end
end

new_inverse
