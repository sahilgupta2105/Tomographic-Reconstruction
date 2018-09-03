function [x_laser,y_laser,x_det,y_det]= laser(i,N)

%%%%%%function for calculating laser and detector locations
dx=1/N;
if i<=0
    m= tan(pi*i/180); %%%slope of laser line
    
    %%%calculating the extreme ends of line
    d= -3*dx*sqrt(1+m^2);
    e= 1- m + 3*dx*sqrt(1+m^2);
    
    x_laser=linspace((1-d)/(m+1/m),(1-d*m)/(1+m^2),2*(N-1)+2);
    y_laser=linspace(1+(d-1)/(1+m^2),(m+d)/(1+m^2),2*(N-1)+2);
    x_det=linspace((1-e)/(m+1/m),(1-e*m)/(1+m^2),2*(N-1)+2);
    y_det=linspace(1+(e-1)/(1+m^2),(m+e)/(1+m^2),2*(N-1)+2);
    
    x_laser(1)=[];
    y_laser(1)=[];
    x_det(1)=[];
    y_det(1)=[];
    x_laser(2*(N-1)+1)=[];
    y_laser(2*(N-1)+1)=[];
    x_det(2*(N-1)+1)=[];
    y_det(2*(N-1)+1)=[];
    
else
    m= tan(pi*-i/180); %%%slope of laser line
    
    %%%calculating the extreme ends of line
    d= -3*dx*sqrt(1+m^2);
    e= 1- m + 3*dx*sqrt(1+m^2);
    
    x_laser=linspace(1-(1-d)/(m+1/m),1-(1-d*m)/(1+m^2),2*(N-1)+2);
    y_laser=linspace(1+(d-1)/(1+m^2),(m+d)/(1+m^2),2*(N-1)+2);
    x_det=linspace(1-(1-e)/(m+1/m),1-(1-e*m)/(1+m^2),2*(N-1)+2);
    y_det=linspace(1+(e-1)/(1+m^2),(m+e)/(1+m^2),2*(N-1)+2);
    x_laser(1)=[];
    y_laser(1)=[];
    x_det(1)=[];
    y_det(1)=[];
    x_laser(2*(N-1)+1)=[];
    y_laser(2*(N-1)+1)=[];
    x_det(2*(N-1)+1)=[];
    y_det(2*(N-1)+1)=[];
    
    
end
