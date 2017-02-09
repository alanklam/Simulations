%% 
% This Matlab script performs numerical simulation for a 1D model of neural firing described by the Hodgkin-Huxley equations 
% (https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model)
% The model includes also myelinated segment, formulated by compartmental model.

clear all
%%
% Tunable parameters of the model are defined here

%step size
dt=0.01;
t = 0:dt:100;
%number of compartments
k = 51;

%Array Declaration
Va = ones(k,length(t));
Vm = zeros(k,length(t));
Vt = Va + Vm;

HH = zeros(k,1);%Hodgkin Huxley term
It = zeros(k,k); %Cytoplasm node current operator
Itn = zeros(k,k); % Periaxonal node current operator
% Inj = zeros(k,length(t));%Injected Current in muA/cm2
Iax = zeros(k,length(t));%Current through axon membrane

n = ones(k, length(t));
m = ones(k, length(t));
h = zeros(k, length(t));
Kc = [1:k]; %Sodium channels per micro^2
Nac = [1:k]; %Potassium channels per micro^2
gK = [1:k]; %potassium conductance per micro^2
gNa = [1:k]; %sodium conductance per micro^2
x = [1:k]; %length of compartments
A = zeros(k,1); %1/Surface Area
NR = 1:k; %Nodes of Ranvier
LJ = NR - 1;%Left juxtapodal nodes
RJ = NR + 1;%Right juxtapodal nodes



%Constants

% Reversal potentials
Vl = 50.6; 
Vk = 72;   
Vna = -55; 
% Conductance/Resistance
gl = .3;%mS per cm^2
rmy = 0.104;%Myelin resistivity (kohm cm2)
rc = 2e-4;%Cytoplasm resistivity (Mohm cm)
rpa = 6e-4;%Periaxonal resistivity (Mohm cm)
gNach = 2;%e-8 mS=10pS per Na channel
gKch = 2;%e-8 mS per K channel
%Capacitance per micro^2
Ca = 1; % muF/cm^2 axon capacitance
Cm = 0.6;%0.001 muF/cm^2 myelin capacitance
%Node Concentrations (amount of channels present per micro^2)
KNR = 18;%18
NaNR = 60;%60
KJN = 500;%
NaJN = 25;%
KM = 10;%
NaM = 25;%
%Node Lengths
d = 1;% in 10 micro diameter
NRx = 20;%Node of Ranvier
JNx = 20;%Juxtapodal
Mx = 20;%Internode
Inj=10*exp(-((1:k)'/4).^2);

% Node Creation
for c = 1:k
    if ismember(c,NR)
        x(c) = NRx;
        gK(c) = gKch*KNR;
        gNa(c) = gNach*NaNR;
        A(c) = 1/(pi*d*x(c));
    elseif  ismember(c,LJ) || ismember(c,RJ)
        x(c) = JNx;
        gK(c) = gKch*KJN;
        gNa(c) = gNach*NaJN;
        A(c) = 1/(pi*d*x(c));
    else
        x(c) = JNx;
        gK(c) = gKch*KM;
        gNa(c) = gNach*NaM;
        A(c) = 1/(pi*d*x(c));
    
    end
end

 

%Creating Cytoplasm Current Matrix
It(1,1) = (-pi*d*d)/(4*rc*(x(1)/2+x(2)/2));
It(1,2) = (pi*d*d)/(4*rc*(x(1)/2+x(2)/2));
It(k,k) = (-pi*d*d)/(4*rc*(x(k-1)/2+x(k)/2));
It(k,k-1) = (pi*d*d)/(4*rc*(x(k-1)/2+x(k)/2));
for c=2:1:k-1
    It(c,c-1) =(pi*d^2)/(4*rc*(x(c-1)/2+x(c)/2));
    It(c,c) = -(pi*d^2)/(4*rc*(x(c-1)/2+x(c)/2))-(pi*d*d)/(4*rc*(x(c)/2+x(c+1)/2));
    It(c, c+1) = (pi*d^2)/(4*rc*(x(c)/2+x(c+1)/2));
end


%Creating Periaxonal Current Matrix

if ~ismember(1,NR)
    Itn(1,1) = -(pi*d*d)/(4*rpa*(x(1)/2+x(2)/2));
    Itn(1,2) = (pi*d*d)/(4*rpa*(x(1)/2+x(2)/2));
end

if ~ismember(k,NR);    
    Itn(k,k) = -(pi*d*d)/(4*rpa*(x(k-1)/2+x(k)/2));
    Itn(k,k-1) = (pi*d*d)/(4*rpa*(x(k-1)/2+x(k)/2));
end  
% for c=2:k-1
%     if  ismember(c,LJ)
%         Itn(c,c) =-(pi*d*d)/(4*rpa*(x(c-1)/2+x(c)/2));
%         Itn(c,c-1) = (pi*d*d)/(4*rpa*(x(c-1)/2+x(c)/2));
%     elseif  ismember(c,RJ)
%         Itn(c,c) = -(pi*d*d)/(4*rpa*(x(c)/2+x(c+1)/2));
%         Itn(c,c+1) = (pi*d*d)/(4*rpa*(x(c)/2+x(c+1)/2));
%     else
%         Itn(c,c) = -(pi*d*d)/(4*rpa*(x(c-1)/2+x(c)/2))-(pi*d*d)/(4*rpa*(x(c)/2+x(c+1)/2));
%         Itn(c,c-1) = (pi*d*d)/(4*rpa*(x(c-1)/2+x(c)/2));
%         Itn(c,c+1) = (pi*d*d)/(4*rpa*(x(c)/2+x(c+1)/2));
%     end
% end
for c=2:k-1
        Itn(c,c) = -(pi*d*d)/(4*rpa*(x(c-1)/2+x(c)/2))-(pi*d*d)/(4*rpa*(x(c)/2+x(c+1)/2));
        Itn(c,c-1) = (pi*d*d)/(4*rpa*(x(c-1)/2+x(c)/2));
        Itn(c,c+1) = (pi*d*d)/(4*rpa*(x(c)/2+x(c+1)/2));
end


% HH + myelin
for i = 1:(length(t)-1)
    Iax(:,i) = It*Vt(:,i);
      for j = 1:k
        n(j,i+1) = n(j,i) + dt*(An(Va(j,i))*(1-n(j,i)) - ...
            Bn(Va(j,i))*n(j,i));
        
        m(j,i+1) = m(j,i) + dt*(Am(Va(j,i))*(1-m(j,i)) - ...
            Bm(Va(j,i))*m(j,i));
        
        h(j,i+1) = h(j,i) + dt*(Ah(Va(j,i))*(1-h(j,i)) - ...
            Bh(Va(j,i))*h(j,i));
        HH(j) = Ion(n(j,i),m(j,i),h(j,i),Va(j,i),gl,gK(j),gNa(j),...
            Vl,Vk,Vna);
      end
      Va(:,i+1) = Va(:,i) + dt*(A.*Iax(:,i) - HH(:)+Inj)/Ca;
      Vm(:,i+1) = Vm(:,i) + dt*((1/Cm)*A.*(Iax(:,i)+Itn*Vm(:,i)) -...
          (1/rmy)*Vm(:,i)/Cm);
      Vm(NR,i+1)=0;
      Vt(:,i+1) = Va(:,i+1) + Vm(:,i+1);   
end


%function [I] = Ion(n,m,h,V,gl,gk,gna,Vl,Vk,Vna)
%I = gl*(V+Vl) + gk*(n^4)*(V + Vk) + gna*(m^3)*h*(V + Vna);
%end
% figure;
% mesh(Vt);
% figure
% mesh(Vm);

% plot(t,Va(1,:))
% 


%%
% Show action potential propagation along axon as movie

for i=1:50:length(t);
    plot(cumsum(x)'/1000,Vt(1:k,i)); ylim([-80,50]); hold on; 
    plot(cumsum(x)'/1000,Inj,'r'); hold off; %plot injection current
    xlabel('x (mm)')
    F(i)=getframe;
end;
movie(F,1,15)










