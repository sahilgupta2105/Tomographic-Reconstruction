%%%%CONSTANTS%%%%%%

%Uniform pressure field//////Units= atm
P= str2num(get(handles.edit_pressure,'String'));

%Line strength for target frequency and temperature//// Units= cm-2/atm
S= str2num(get(handles.edit_line_strength,'String'));

%Target center frequency///// Units= cm-1
mu_o= str2num(get(handles.edit_freq,'String'));

%Uniform Temparature and molecular mass of gas
T=str2num(get(handles.edit_temperature,'String'));
M=str2num(get(handles.edit_mm,'String'));

%broadening coefficients for air and gas////// Units= cm-1/atm
bc_air= str2num(get(handles.edit_bc_air,'String'));
bc_self= str2num(get(handles.edit_bc_gas,'String'));

%%%%%%MESH%%%%%
%Mesh size
% prompt3='Mesh resolution: ';
% N=input(prompt3);
N=str2num(get(handles.edit_mesh,'String'));

set(handles.edit_notif,'String','Data Loaded...');
drawnow
done=0.1;
f=statusbar
f=statusbar('Working',f)
f=statusbar(0.1,f)
v=statusbar('on',f)
n_ray=2*(N-1);

%% coordinates%%
x= linspace(0,1,N);
y= linspace(0,1,N);

%Plot the mesh grid
[xx,yy]=meshgrid(x,y);


%Defing array for concentration of species at cell center
conc = zeros(N-1,N-1);

dx= 1/N;

%coordinates of cell centers, lasers and detectors
x_cell = zeros(N-1,N-1);
y_cell = zeros(N-1,N-1);

%conc. profile selection
conc_profile= get(handles.popup_conc,'Value');

switch conc_profile
    
    case 1
        set(handles.edit_notif,'String','Please select conc. profile...');
        drawnow
        
    case 2
        
        for i= 1:1:N-1
            for j=1:1:N-1
                x_cell(i,j) = (x(j)+x(j+1))/2;% coordinates of cell center
                y_cell(i,j) = (y(i)+y(i+1))/2;
                %%%%%%%%%%%%%%%%%%%%%%%%%Gauss-Phantom profile
                distance=sqrt(x_cell(i,j)^2+y_cell(i,j)^2);
                if distance<0.5
                    cosG=0.25*(1-cos(2*pi*(x_cell(i,j)+0.5)^0.8))*(1-cos(2*pi*(y_cell(i,j)+0.5)^(2/3)));
                else
                    cosG=0;
                end
                cosGauss=1.09*(0.3*cosG + 0.8*(exp(-(40*(x_cell(i,j)-0.3)^2+(15*(y_cell(i,j)-0.4)^2)))+exp(-(15*(x_cell(i,j)-0.7)^2+(60*(y_cell(i,j)-0.7)^2)))));
                conc(i,j)= 0.18*cosGauss;
            end
        end
        
     case 3
        
        for i= 1:1:N-1
            for j=1:1:N-1
                x_cell(i,j) = (x(j)+x(j+1))/2;% coordinates of cell center
                y_cell(i,j) = (y(i)+y(i+1))/2;

                        distance= sqrt((x_cell(i,j)-0.25)^2 + (y_cell(i,j)-0.35)^2);
                        %condition for concentration distribution
                        if distance<=0.125
                            conc(i,j)= 0.18*(0.5 + 0.5*erf(0.4*(0.25-distance))); % top-hat function
                        end
                        distance= sqrt((x_cell(i,j)-0.55)^2 + (y_cell(i,j)-0.75)^2);
                        %condition for concentration distribution
                        if distance<=0.075
                            conc(i,j)= 0.18*(0.5 + 0.5*erf(0.4*(0.25-distance))); % top-hat function
                        end
            end
        end
end

set(handles.edit_notif,'String','Concentration loaded...');
drawnow
done=0.2;
statusbar(done,f)
v=statusbar('on',f)
%coordinates for lasers
n_laser=str2num(get(handles.edit_lasers,'String'));
% prompt='No. of lasers:';
% n_laser= input(prompt);
% prompt1='Values of angles in degrees (Range:[-90,90]):';
% angles=zeros(n_laser,1);
% for i=1:1:n_laser
%     angles(i,1)=input(prompt1);
% end
eval(['dum1 = {' get(handles.edit_angles,'String') '}']);
angles=cell2mat(dum1);
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
axes(handles.fig_mesh)
plot(xx,yy,'g');
hold on;
plot(yy,xx,'g');
set(gca,'Color',[0,0,0])
hold on;

for i=1:1:n_laser
    plot(x_laser(:,i),y_laser(:,i),'r.');
   
end
for i=1:1:n_laser
    plot(x_det(:,i),y_det(:,i),'b.');
   
end
axis image;
drawnow
hold off;


set(handles.edit_notif,'String','Lasers loaded...');
drawnow



done=0.30;
statusbar(done,f)
v=statusbar('on',f)
%Surface plot of concentration distribution
axes(handles.fig_conc)
surface(x_cell,y_cell,conc);
axis image;
colorbar
drawnow

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

set(handles.edit_notif,'String','Ray length calculation started...');
drawnow
done=0.4
statusbar(done,f)
v=statusbar('on',f)
%%%calculating length of each ray inside each box
dt_status= 0.1/n_laser;
laser_len = zeros((N-1)*(N-1),n_ray,n_laser);

for k=1:1:n_laser
            done=done+dt_status;
            statusbar(done,f)
            v=statusbar('on',f)
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
done=0.50;
statusbar(done,f)
v=statusbar('on',f)
set(handles.edit_notif,'String','Ray calculation completed...');
drawnow

%Line funciton

VP = zeros(N-1,N-1); % Voigt profile

%Laser intensities

I = zeros(100,n_ray,n_laser);

set(handles.edit_notif,'String','Optical Depth calculation started...');
drawnow

%%%%%%wave number iteration main loop

delta_mu=0.0005*mu_o;
mu=linspace(mu_o-delta_mu,mu_o+delta_mu,100);
I_o= zeros(1,100);
mT = zeros(N-1,N-1);% Modulation index
%modulation amplitude

a=0.2;


wn= mu_o;
done=0.6
statusbar(done,f)
v=statusbar('on',f)
%  Voigt profile function
for i= 1:1:N-1
    for j=1:1:N-1
        mT(i,j)=a/gamma_V(i,j);
        fun = @(x)((C_L(i,j)/pi)*((gamma_V(i,j)./((mT(i,j)*cos(x)).^2 + (gamma_V(i,j))^2))) + (C_G(i,j)/gamma_V(i,j))*sqrt(log(2)/pi)*exp((-log(2)*(mT(i,j)*cos(x)).^2)./(gamma_V(i,j))^2))*cos(2*x);
        VP(i,j)= integral(fun,-pi,pi,'ArrayValued', true);
    end
end
%calculating optical depth%%%% Number of views= 3 %

done=0.7;
statusbar(done,f)
v=statusbar('on',f)

% Optical Depth is calculated for each laser

OD = zeros(n_ray,n_laser);

%Optical depth measurement
for l=1:1:n_laser
     done=done+dt_status;
            statusbar(done,f)
            v=statusbar('on',f)
    for i=1:1:n_ray           
        for j=1:1:N-1
            for k=1:1:N-1
                OD(i,l)= OD(i,l)+(-1/pi)*(P*conc(j,k)*S*VP(j,k)*laser_len((j-1)*(N-1)+k,i,l)); %%%%%%LASER 'l'
            end
        end
    end
end
done=0.9
statusbar(done,f)
v=statusbar('on',f)
set(handles.edit_notif,'String','TDLAS calculation completed...');
drawnow

%Intensity calculation for each laser

%The total intial intensity for each laser will be equally distributed
%among all the rays of the particular laser. In this case each laser has 2N
%rays

% I_o(:)= (1/(n_ray))*1*wn(:);

% for i=1:1:n_ray
%     I_1(ii,i)=I_o(ii)*exp(-OD_1(i,1));% LASER 1
%     I_2(ii,i)=I_o(ii)*exp(-OD_2(i,1));% LASER 2
%     I_3(ii,i)=I_o(ii)*exp(-OD_3(i,1));% LASER 3
%     I_4(ii,i)=I_o(ii)*exp(-OD_4(i,1));% LASER 4
% % end
% I_1= wn(50);
% I_1= I_1+
done=0.9;
statusbar(done,f)
v=statusbar('on',f)
axes(handles.fig_od)
for i=1:1:n_laser
plot(OD(:,i),'color',rand(1,3));
dum_label{i}= strcat('Laser ',num2str(i));
hold on;
end
hold off;
xlabel('Ray number');
ylabel('Optical Depth');
legend(dum_label);
drawnow
set(handles.edit_notif,'String','Dataset generation successful...');
drawnow

%%%%Weighted lengths inside each box
t_len=sum(laser_len,1);
laser_lenc=laser_len;
for j=1:1:n_laser
    for i=1:1:n_ray
        laser_lenc(:,i,j)=laser_lenc(:,i,j)/t_len(1,i,j);
    end
end

%%%Creating matrix 'W', containing the weighted lengths of rays inside
%%%boxes
l_mat=zeros(n_laser*n_ray,(N-1)*(N-1));
p=1;
for j=1:1:n_laser
    for i=1:1:n_ray
        l_mat(p,:)=laser_lenc(:,i,j);
        p=p+1;
    end
end


% % Intensity matrix
p_mat=zeros(n_ray*n_laser,1);
for i=1:1:n_laser
    for j=1:1:n_ray
        p_mat((i-1)*n_ray+j,1)=OD(j,i);% 'P' matrix
    end
end


done=1;
statusbar(done,f)
v=statusbar('on',f)
delete(statusbar)

save guioutput
