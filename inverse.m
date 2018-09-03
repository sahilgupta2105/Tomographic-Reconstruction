load('guioutput')

set(handles.edit_notif,'String','MAART: In process...');
drawnow
epsilon=str2num(get(handles.epsilon,'String'));
delta=str2num(get(handles.delta,'String'));
%MAART
%%% 'delta' is the Regularization parameter
A=l_mat;
b=p_mat;
[m,n] = size(A);
X = ones(n,1);
Xn=zeros(n,1);
beta=0.25;

for t=1:1:25000
    
    %%%Calculating lamda
    lamda= zeros(n,m);
    sum1=zeros(m,1);
    X_A=zeros(n,m);
    for i=1:1:m
        dum1=A(i,:)';
        sum1(i,1)= sum(X.*dum1);
        
        for j=1:1:n
            X_A(j,i)=(X(j)*A(i,j));
        end
    
    end
    I=find(sum1>0)';
        
    for i=I
        lamda(:,i)=(X_A(:,i))./sum1(i,1);
    end

    lamda=beta*lamda;
    %%%Calculating X*A_i
    dum=zeros(m,1);
    for k=1:1:m
        dum(k,1)=dot(X,A(k,:));
    end
    %%%Updating value of X
    Xn=X-lamda*(dum-b);
    %%%Non-Negative condition
    for l=1:1:n
        if Xn(l,1)<0
            Xn(l,1)=0.0;
        end
    end
    
    %%%Regularization of Xn
    
    for i=1:1:n
        nn=sqrt(n);
        yy=fix(i/nn);
        xx=rem(i,nn);
        if yy==0||yy==nn-1||xx==1||xx==0
        else
          mat=[Xn(i,1),Xn(i-nn-1,1),Xn(i-nn,1),Xn(i-nn+1,1),Xn(i-1,1),Xn(i+1,1),Xn(i+nn-1,1),Xn(i+nn,1),Xn(i+nn+1,1)];  
          Xn(i,1)=(1-delta)*mat(1)+delta*(sum(mat)-mat(1))/8;
        end
    end
    
    %%%Calculating the error
    err=norm(Xn(:,1)-X);
    ss=strcat('Error...',num2str(err));
    set(handles.edit_notif,'String',ss);
    drawnow
    X=Xn(:,1);
    
    if err<epsilon
        break;
    end
    
end 
%%%%%%%%%%%%%%%END of MAART%%%%%%%%%%%%%%%%%%%%%%%%%

x_conc= X;
set(handles.edit_notif,'String','Inverse calculation completed');
drawnow
conc_ans=zeros(N-1,N-1);
for j=1:1:(N-1)
    for  i=1:1:(N-1)
        conc_ans(i,j)=(-pi)*x_conc(i+(j-1)*(N-1))/(P*S*VP(i,j));
        if conc_ans(i,j)<0
            conc_ans(i,j)=0.0;
        end
    end
end

norm11=(conc_ans'-conc);
set(handles.edit_notif,'String','Calculation completed...');
drawnow
axes(handles.fig_inv)
surface(x_cell,y_cell,conc_ans');
colorbar
drawnow
guidata(hObject,handles);
axes(handles.fig_error)
surface(x_cell,y_cell,norm11);
colorbar

