% clear all
% clc
%[time,S]=voltageclamp(t);
%[Actend,time]=ActCurve;
%[Inactend]=InActCurve;
 %InacRecovCurve
   [output,ACTdata,time]= ActCurveKv14;
%   [INACTdata,timeIN]=InActCurve;
%[Recover,time]=InacRecovCurve;
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
    yvals=(solext.y);
     k = find(time>=5200,1);
     Recover(i)=(max(yvals(2,k:end))).^3.*(max(yvals(5,k:end)));
     figure(5)
     hold on
     plot(time,yvals(2,:).^3.*yvals(5,:))
 end
toc
  figure(6)
  hold on
   plot(Vsteps,Recover,'o')
end   



function [INACTdata,time]=InActCurveExt
    prot=2;
    
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    %output=zeros(1,length(Vsteps));
    INACTdata=zeros(1,length(Vsteps));
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    t=[0 5000];
        
    [sol]=voltageclamp2(t,Vstep,prot);
    endtime=6000;
    [solext]=vclampextend(t,Vstep,prot,sol,endtime);
    
    time=solext.x;
    yvals=(solext.y);
     k = find(time>=5000,1);
     INACTdata(i)=(max(yvals(2,k:end))).^3.*(max(yvals(5,k:end)));
    figure(3)
    hold on
    plot(time,yvals(2,:).^3.*yvals(5,:))
 end
 toc
  figure(4)
   hold on
 plot(Vsteps,INACTdata,'o')
end  


function [output,ACTdata,time]=ActCurveKv14
    prot=4;
    t=[0 50];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    ACTdata=zeros(1,length(Vsteps));
    Taus=zeros(1,length(Vsteps));
    Ninfs=zeros(1,length(Vsteps));  
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
    output(i)=S(end,2);
    k = find(time>=20,1);
    j = find(time>=50,1)-1;
    ACTdata(i)=(max(S(k:j,2))).^3.*(max(S(k:j,5))).^1;
    
%     fo = fitoptions('Method','NonlinearLeastSquares','Lower',[ACTdata(i).^(1/4),0],'Upper',[ACTdata(i).^(1/4),10],...
%                'StartPoint',[ACTdata(i).^(1/4),.1]);
%     ft = fittype('(a*(1-exp(-(x-20)/b)))^4','options',fo);
%     [curve4,fof1]=fit(time(k:j),(S(k:j,2)),ft);
%     Taus(i)=curve4.b;
%     Ninfs(i)=curve4.a;
    
    
 end
 toc
%  figure(1)
%    hold on
%     plot(Vsteps,Taus,'o')
 
   figure(2)
    hold on
 plot(Vsteps,ACTdata,'o')
 
end  



function [time,S]=voltageclamp(t,Vstep,prot)

S0 = [1,0,0,1,0]'; 
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[time,S]=ode15s(@(t,S)deRHS(t,S,Vstep,prot),t,S0,options);
S(end,:);

end

function [sol]=voltageclamp2(t,Vstep,prot)

S0 = [1,0,0,1,0]';
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
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

if t<5000
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



%%%for a 2 kv11: 2kv14Ctype het channel
pC=S(1);
pO=S(2);
pI=S(3);
pC11=S(4);
pO11=S(5);



%Best Fit Parms so far
%parms=[452.4278    1.3168    3.2838    7.5995    0.1221    0.5359];
parms=[429.42,  2.15,  8.56,   5.36,    .102,  .46];
m1=parms(1);
n1=parms(2);
c1=parms(3);
d1=parms(4);
m2=parms(5);
n2=parms(6);

%Rates Kv14
 a1=m1*exp(n1*V*.01)/1000;
 b1=c1*exp(-d1*V*.01)/1000;
 aI=m2/1000;
 bI=n2/1000;
 %Rates Kv11
 a2=.12889*exp((V+45)/33.90877);
 b2=.12889*exp(-(V+45)/12.42101);
 
%%2Kv1.1 and 2Kv1.4
 s_prime=[-a1*pC+b1*pO,...
          a1*pC-b1*pO+bI*pI-aI*pO,...
          aI*pO-bI*pI,...
          -a2*pC11+b2*pO11,...
          -b2*pO11+a2*pC11]';


end