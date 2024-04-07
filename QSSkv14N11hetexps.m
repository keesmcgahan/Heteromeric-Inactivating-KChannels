%[Taus,ACTdata,time,Kappa]= ActCurveKv14;
[Recoverext,timeext,yvals]=InacRecovExt;
[INACTdataext,timeext]=InActCurveExt;

 
function [Recover,time,yvals]=InacRecovExt
    prot=3;
    %t=[0:0.1:500];
    Vsteps=[.1,.5,1,5,10,20,50,100,500,1000,2000,5000,7000,10000,20000];
    %Vsteps=[.1];
     Recover=zeros(1,length(Vsteps));
     t=[0 5200];
     Vstep=0;
     [sol]=voltageclamp2(t,Vstep,prot);
 tic    
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
        
        
    
    endtime=6000+Vstep;
    [solext]=vclampextend(t,Vstep,prot,sol,endtime);
    
    time=solext.x;
    yvals=solext.y;
    
    Kappa=zeros(1,length(time));
    

    k1=.950; k2=.617; k3=.40; k4=.49; k5=.484; k6=.13; k7=.4;
    
    m1 = find(time>=20,1);
    m2 = find(time>=5200,1)-1;
    m3 = find(time>=5200+Vstep,1);
    V=-90;
    a1=k1*exp((V-(100*k3))/(100*k4));
    b1=k2*exp((V-(100*k3))/(100*k4)).*exp(-(V+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(V+(100*k5))/(100*k6)));
    k1f11=.12889*exp((V+45)/33.90877);
    k1b11=.12889*exp(-(V+45)/12.42101);
    Kappa(1,1:m1)=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);
    Kappa(1,m2:m3)=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);
    V=50;
    a1=k1*exp((V-(100*k3))/(100*k4));
    b1=k2*exp((V-(100*k3))/(100*k4)).*exp(-(V+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(V+(100*k5))/(100*k6)));
    k1f11=.12889*exp((V+45)/33.90877);
    k1b11=.12889*exp(-(V+45)/12.42101);
    Kappa(1,m1:m2)=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);
    Kappa(1,m3:end)=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);
    
    datapts=yvals(1,:)./Kappa;
     k = find(time>=5200,1);
     Recover(i)=max(datapts(1,k:end))/max(datapts(1,1:m1));
     figure(5)
     hold on
     plot(time,datapts(1,:))
 end
toc
  figure(6)
   plot(Vsteps,Recover,'o')
end   


function [output,time]=ActCurve
    prot=1;
    t=[0 500];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50,60,70,90];
    output=zeros(1,length(Vsteps));
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
  
    output(i)=S(end,1);
    figure(1)
    hold on
    plot(time,S(:,1))
  end
end  

function [INACTdata,time]=InActCurve
    prot=2;
    t=[0 6000];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    INACTdata=zeros(1,length(Vsteps));
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
    output(i)=S(end,1);
    k = find(time>=5000,1);
    INACTdata(i)=max(S(k:end,1));
    figure(3)
    hold on
    plot(time,S(:,1))
 end
 toc
  figure(4)
 plot(Vsteps,INACTdata,'o')
end  

function [INACTdata,time]=InActCurveExt
    prot=2;
    
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    %output=zeros(1,length(Vsteps));
    INACTdata=zeros(1,length(Vsteps));
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    t=[0 4900];
        
    [sol]=voltageclamp2(t,Vstep,prot);
    [maxsol]=voltageclamp2(t,50,prot);
    endtime=6000;
    [solext]=vclampextend(t,Vstep,prot,sol,endtime);
    
    time=solext.x;
    yvals=solext.y;
    
     Kappa=zeros(1,length(time));
    k1=.950; k2=.617; k3=.40; k4=.49; k5=.484; k6=.13; k7=.4;
    
    m1 = find(time>=20,1);
    m2 = find(time>=5000,1)-1;
    V=-90;
    a1=k1*exp((V-(100*k3))/(100*k4));
    b1=k2*exp((V-(100*k3))/(100*k4)).*exp(-(V+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(V+(100*k5))/(100*k6)));
    k1f11=.12889*exp((V+45)/33.90877);
    k1b11=.12889*exp(-(V+45)/12.42101);
    Kappa(1,1:m1)=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);
    a1=k1*exp((Vstep-(100*k3))/(100*k4));
    b1=k2*exp((Vstep-(100*k3))/(100*k4)).*exp(-(Vstep+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(Vstep+(100*k5))/(100*k6)));
    k1f11=.12889*exp((Vstep+45)/33.90877);
    k1b11=.12889*exp(-(Vstep+45)/12.42101);
    Kappa(1,m1:m2)=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);
    V=50;
    a1=k1*exp((V-(100*k3))/(100*k4));
    b1=k2*exp((V-(100*k3))/(100*k4)).*exp(-(V+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(V+(100*k5))/(100*k6)));
    k1f11=.12889*exp((V+45)/33.90877);
    k1b11=.12889*exp(-(V+45)/12.42101);
    Kappa(1,m2:end)=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);
    datapts=yvals(1,:)./Kappa;
    dataptsmaxV=maxsol.y(1,:)/Kappa(end);
     k = find(time>=5000,1);
     INACTdata(i)=max(datapts(1,k:end))/max(dataptsmaxV(1,:));
    figure(3)
    hold on
    plot(time,datapts)
 end
 toc
  figure(4)
 plot(Vsteps,INACTdata,'o')
end  

function [Recover,time]=InacRecovCurve
    prot=3;
    %t=[0:0.1:500];
    Vsteps=[.1,.5,1,5,10,20,50,100,500,1000,2000,5000,7000,10000,20000];
    %Vsteps=[.1];
     Recover=zeros(1,length(Vsteps));
     tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
        t=[0:.1: 6000+Vstep];
    [time,S]=voltageclamp(t,Vstep,prot);
    k = find(time>=5200,1);
    Recover(i)=max(S(k:end,1));
    figure(5)
    hold on
    plot(time,S(:,1))
 end
 toc
  figure(6)
   plot(Vsteps,Recover,'o')
end  

function [Taus,ACTdata,time,Kappa]=ActCurveKv14
    prot=4;
    t=[0 50];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    ACTdata=zeros(1,length(Vsteps));
    Taus=zeros(1,length(Vsteps));
    %Ninfs=zeros(1,length(Vsteps));  
    
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
    
    
    Kappa=zeros(1,length(time));
    k1=.950; k2=.617; k3=.40; k4=.49; k5=.484; k6=.13; k7=.4;
    output(i)=S(end,1);
    k = find(time>=20,1);
    j = find(time>=50,1)-1;
    V=-90;
    a1=k1*exp((V-(100*k3))/(100*k4));
    b1=k2*exp((V-(100*k3))/(100*k4)).*exp(-(V+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(V+(100*k5))/(100*k6)));
    Kappa(1,1:k)=(a1+b1).^4./a1.^4;
    
    a1=k1*exp((Vstep-(100*k3))/(100*k4));
    b1=k2*exp((Vstep-(100*k3))/(100*k4)).*exp(-(Vstep+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(Vstep+(100*k5))/(100*k6)));
    Kappa(1,k:end)=(a1+b1).^4./a1.^4;
    
    ACTdata(i)=max(S(k:j,1));
    
    figure(7)
    hold on
    plot(time,S(:,1)./Kappa')
    hold off
  

 end
 toc

 
%    figure(2)
%  plot(Vsteps,ACTdata,'o')
%  
end  



function [time,S]=voltageclamp(t,Vstep,prot)

S0 = [1,0]'; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8); 
[time,S]=ode15s(@(t,S)deRHS(t,S,Vstep,prot),t,S0,options);
S(end,:);

end

function [sol]=voltageclamp2(t,Vstep,prot)

S0 = [1,0]';
options = odeset('RelTol',1e-8,'AbsTol',1e-8); 
[sol]=ode15s(@(t,S)deRHS(t,S,Vstep,prot),t,S0,options);


end

function [solext]=vclampextend(t,Vstep,prot,sol,endtime)

solext=odextend(sol,@(t,S)deRHS(t,S,Vstep,prot),endtime);

end

function [V]=SSact(t,Vstep)

if t<300
    V=-100;
else
    V=Vstep;    
end

end



function [V]=SSinact(t,Vstep)

if t<20
    V=-90;
elseif t>=20 && t<=5000
    V=Vstep;
else
    V=50;    
end

end

function [V]=InacRecov(t,Vstep)

if t<20
    V=-90;
elseif t>=20 && t<5200
    V=50;
elseif t>=5200 && t<(Vstep+5200)
    V=-90;
else
    V=50;
end

end

function [V]=SSActkv14paper(t,Vstep)

    if t<20
        V=-90;
    elseif t>20 && t<50
        V=Vstep;
    elseif t>50 && t<80
        V=-50;
    else 
        V=-90;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%ODE%%%%%%%%%%%%%%%%%%%5

function s_prime=deRHS(t,S,Vstep,prot) 

    if prot==1
        V=SSact(t,Vstep);
    elseif prot==2
        V=SSinact(t,Vstep);
    elseif prot==3
        V=InacRecov(t,Vstep);
    elseif prot==4
        V=SSActkv14paper(t,Vstep);
    end

%%%Variables%%%
X=S(1);
I=S(2);

k1=.950;
k2=.617;
k3=.40;
k4=.49;
k5=.484;
k6=.13;
k7=.4;

a1=k1*exp((V-(100*k3))/(100*k4));
b1=k2*exp((V-(100*k3))/(100*k4)).*exp(-(V+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(V+(100*k5))/(100*k6)));
aI=14.24/1000;
bI=1.48/1000;

k1f11=.12889*exp((V+45)/33.90877);
k1b11=.12889*exp(-(V+45)/12.42101);
K=((b1+a1).^2.*(k1f11+k1b11).^2)./(k1f11.^2*a1.^2);


s_prime=[-2/4*aI*X/K+bI*I,...
   2/4*aI*X/K-bI*I]';  

end