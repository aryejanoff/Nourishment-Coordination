%NOTE: turn on the first line below if you want to explore the nourishment efficiency
%for coordination; turn on the second line below if you want to explore the
%nourishment efficiencies for conservative non-coordination
load('future_erosion_sandcost_DATA.mat','gamma_vec','phi_vec','R2_coord','R3_coord');
% load('future_erosion_sandcost_DATA.mat','gamma_vec','phi_vec','R2_uncoord_cons','R3_uncoord_cons');
n2=length(gamma_vec);
n3=length(phi_vec);
N_efficiency_storage=NaN(n2,n3);

for iR2=1:n2
    N_efficiency_vector=NaN(1,n3);
    for iR3=1:n3
    %% Rotation Lengths by Coordination Scheme
    %NOTE: turn on the first two lines below if you want to explore the nourishment
    %efficiencies for coordination; turn on the third and fourth lines below if you
    %want to explore the nourishment efficiencies for conservative
    %non-coordination
        R2=R2_coord(iR2,iR3); %Coordinated Rotation Length Community 1
        R3=R3_coord(iR2,iR3); %Coordinated Rotation Length Community 2
%         R2=R2_uncoord_cons(iR2,iR3); %Uncoordinated Rotation Length Community 1
%         R3=R3_uncoord_cons(iR2,iR3); %Uncoordinated Rotation Length Community 2
        %% Input Physical Pmarameters %%
        lot1=30; lotj=30; lotm=30;
        lot_size=[lot1 lotj lotm];
        w_init=30;
        beta2=0.4; beta3=0.4; %beach width hedonic parameters
        PV2=1e6;
        alpha2=PV2/(w_init^beta2);
        PV3=2e6;
        alpha3=PV3/(w_init^beta3);
        s1=1500; sj=1500; sm=1500; 
        s=[s1,sj,sm]; %alongshore compartment lengths (m)
        rows_along2=s(2)/lot_size(2); rows_along3=s(3)/lot_size(3);
        rows_cross2=1; %# of cross-shore property rows 2
        rows_cross3=1; %# of cross-shore property rows 3
        comm_width2=rows_cross2*lot_size(1); comm_width3=rows_cross3*lot_size(1); %initial Community Width (m)
        psi=0;
        D=16; %depth of closure (m) NJ=16
        gamma=gamma_vec(iR2); %5; %erosion rate (m/yr) NJ=1
        d=600000; %alongshore flux coeff NJ=50000
        K=2000; %cross-shore flux coeff NJ=2000
        phi2=phi_vec(iR3); phi3=phi_vec(iR3);
        phi=[NaN phi2 phi3]; %sand cost ($/m^3) NJ=10
        c=1e6; %fixed nourishment cost ($) NJ=3e6
        rho=0.06; %discount rate
        xN2=50; xN3=50;
        xN=[NaN xN2 xN3]; %nourishment magnitude (m) NJ=50
        nu=0; %beach width decline beyond max threshold
        theta_eq=0.02; %equilibrium shoreface slope
        k2=0;
        k3=0;
        %% Computational Parameters %%
        tmax=50; dt=0.05; t=0:dt:tmax; n=length(t); 
        Smax=3; ds=1; S=1:ds:Smax; m=length(S);
        A2=alpha2*rho; A3=alpha3*rho; 
        theta=zeros(n,m); qL=zeros(n,m); qC=zeros(n,m); fS=zeros(n,m); fT=zeros(n,m); 
        xS=zeros(n,m); xT=zeros(n,m); xH=zeros(n,m); w=zeros(n,m);
        B=zeros(n,m); B_disc=zeros(n,m); TB=zeros(n,m);
        C=zeros(n,m); TC=zeros(n,m);
        nb=NaN(n,m); N_h=zeros(n,m);        
        volume=zeros(n,m); nE=zeros(n,m);
        V_nourish=zeros(1,n); V_loss=zeros(1,n); efficiency=zeros(1,n); 
        q_loss=zeros(1,n); q_along=zeros(1,n); q_cross=zeros(1,n); q_gamma=zeros(1,n);
        nb2=NaN(1); nb3=NaN(1); behavior=NaN(1); NB2=NaN(1); NB3=NaN(1); TNB=NaN(1);

        %% Initial Conditions %%
        xS(1,1)=comm_width2+w_init; xS(1,2)=comm_width2+w_init; xS(1,3)=comm_width3+w_init; 
        xT(1,:)=xS(1,:)+(D/theta_eq);
        xH(1,1)=comm_width2; xH(1,2)=comm_width2; xH(1,3)=comm_width3;
        N_h(1,1)=rows_cross2; N_h(1,2)=rows_cross2; N_h(1,3)=rows_cross3;
        w(1,:)=xS(1,:)-xH(1,:);

        %% Main Code %%
        for i=1:n-1
            for j=2:m-1                           
                %% Nourishment Initiation + Volume
                ff2=round(t(i)-k2*R2,4);
                if ff2==0
                    k2=k2+1; volume(i,2)=0.5*s(2)*D*xN(2); nE(i,2)=1;
                end
                ff3=round(t(i)-k3*R3,4);
                if ff3==0
                    k3=k3+1; volume(i,3)=0.5*s(3)*D*xN(3); nE(i,3)=1;
                end
                
                %% Fluxes (Along/Cross-shore) and Shoreface Dynamics
                qL(i,j)=d*((xS(i,j-1)-xS(i,j))/((s(j-1)+s(j))/2)); qL(i,1)=d*((xS(i,m)-xS(i,1))/((s(m)+s(1))/2)); qL(i,m)=d*((xS(i,m-1)-xS(i,m))/((s(m-1)+s(m))/2)); 
                theta(i,j)=D/(xT(i,j)-xS(i,j)); theta(i,1)=D/(xT(i,1)-xS(i,1)); theta(i,m)=D/(xT(i,m)-xS(i,m)); 
                qC(i,j)=K*(theta(i,j)-theta_eq); qC(i,1)=K*(theta(i,1)-theta_eq); qC(i,m)=K*(theta(i,m)-theta_eq);
                q_along(i)=((2*(qL(i,1)-qL(i,2)))/(s(2)+s(3))); q_cross(i)=((4*(qC(i,2)+qC(i,3)))/D); q_gamma(i)=gamma;
                q_loss(i)=q_along(i)+q_cross(i)+gamma; %
                V_loss(i+1)=V_loss(i)+dt*q_loss(i); %(dt/3)*(q_loss(1)+4*(sum(q_loss(2:2:end)))+2*sum(q_loss(2:1:end))+q_loss(end));
                efficiency(i)=V_nourish(i)./(V_nourish(i)+V_loss(i));
                
                %% ODE's ((2*d/(s^2))*(xS(i,m)-2*xS(i,1)+xS(i,2)))-(4*K*(theta(i,1)-theta_eq)/D)
                fT(i,j)=(4*K*(theta(i,j)-theta_eq)/D)-gamma; fT(i,1)=(4*K*(theta(i,1)-theta_eq)/D)-gamma; fT(i,m)=(4*K*(theta(i,m)-theta_eq)/D)-gamma;
                fS(i,j)=(2*d/s(j))*((xS(i,j-1)-xS(i,j))/((s(j-1)+s(j))/2)-(xS(i,j)-xS(i,j+1))/((s(j)+s(j+1))/2))-(4*K*(theta(i,j)-theta_eq)/D)-gamma; fS(i,1)=(2*d/s(1))*((xS(i,m)-xS(i,1))/((s(m)+s(1))/2)-(xS(i,1)-xS(i,2))/((s(1)+s(2))/2))-(4*K*(theta(i,1)-theta_eq)/D)-gamma; fS(i,m)=(2*d/s(m))*((xS(i,m-1)-xS(i,m))/((s(m-1)+s(m))/2)-(xS(i,m)-xS(i,1))/((s(m)+s(1))/2))-(4*K*(theta(i,m)-theta_eq)/D)-gamma; 

                %% Numerical Approximations

                xT(i+1,j)=xT(i,j)+dt*fT(i,j); xT(i+1,1)=xT(i,1)+dt*fT(i,1); xT(i+1,m)=xT(i,m)+dt*fT(i,m);
                xS(i+1,1)=xS(i,1)+dt*fS(i,1);  
                if xS(i,j)<=lot_size(1)
                    volume(i,j)=0; volume(i+1,j)=0; N_h(i,j)=0; xH(i,j)=0; N_h(i+1,j)=0; xH(i+1,j)=0;  
                else
                    xH(i+1,j)=xH(i,j);
                end
                if xS(i,m)<=lot_size(1)
                    volume(i,m)=0; volume(i+1,m)=0; N_h(i,m)=0; xH(i,m)=0; N_h(i+1,m)=0; xH(i+1,m)=0;  
                else
                    xH(i+1,m)=xH(i,m);
                end
                if volume(i,2)~=0
                    xS(i+1,2)=xS(i,2)+xN(2);
                elseif volume(i,2)==0
                    xS(i+1,2)=xS(i,2)+dt*fS(i,2);
                end
                if volume(i,3)~=0
                    xS(i+1,3)=xS(i,3)+xN(3); 
                elseif volume(i,3)==0
                    xS(i+1,3)=xS(i,3)+dt*fS(i,3);
                end
                if xS(i,j)<=xH(i,j)+5
                    xH(i+1,j)=xH(i,j)-lot_size(1); N_h(i+1,j)=N_h(i,j)-1;
                else
                    xH(i+1,j)=xH(i,j);
                end
                if xS(i,m)<=xH(i,m)+5
                    xH(i+1,m)=xH(i,m)-lot_size(1); N_h(i+1,m)=N_h(i,m)-1;
                else
                    xH(i+1,m)=xH(i,m);
                end
                w(i,j)=xS(i,j)-xH(i,j); w(i,1)=xS(i,1)-xH(i,1); w(i,m)=xS(i,m)-xH(i,m);
                V_nourish(i+1)=V_nourish(i)+(2*(volume(i,2)+volume(i,3))./((s(2)+s(3))*D));

                %% Housing Lines 
                N_h(i,j)=xH(i,j)/lot_size(1); N_h(i,1)=xH(i,1)/lot_size(1); N_h(i,m)=xH(i,m)/lot_size(1);
                
                %% Benefit 
                B(i,2)=rows_along2*A2*((w(i,2)).^beta2)*((N_h(i,2)).^psi)-(nu*(w(i,2).^2)); B(i,3)=rows_along3*A3*((w(i,3)).^beta3)*((N_h(i,3)).^psi)-(nu*(w(i,3).^2));
                B_disc(i,j)=B(i,j)*exp(-rho*t(i)); B_disc(i,m)=B(i,m)*exp(-rho*t(i));
                TB(i,j)=(dt/3)*(B_disc(1,j)+4*(sum(B_disc(2:2:end,j)))+2*sum(B_disc(2:1:end,j))+B_disc(end,j)); TB(i,m)=(dt/3)*(B_disc(1,m)+4*(sum(B_disc(2:2:end,m)))+2*sum(B_disc(2:1:end,m))+B_disc(end,m));

                %% Cost 
                C(i,j)=(c+volume(i,j)*phi(j))*exp(-rho*t(i)); C(i,m)=(c+volume(i,m)*phi(m))*exp(-rho*t(i));
                if volume(i,j)==0
                    C(i,j)=0;
                end
                if volume(i,m)==0
                    C(i,m)=0;
                end
                if N_h(i,j)==0
                    C(i,j)=0;
                end
                if N_h(i,m)==0
                    C(i,m)=0;
                end
                TC(i,j)=sum(C(:,j)); TC(i,m)=sum(C(:,m));

                %% marginal net benefit 
                nb(i,j)=TB(i,j)-TC(i,j); nb(i,m)=TB(i,m)-TC(i,m);

                %% Net Benefit
                NB2=nb(i,2); NB3=nb(i,3); TNB=NB2+NB3;
            end
        end
        
        %% Identify Behavior
        if (max(xS((tmax/dt-50):end,2))>comm_width2+w_init+xN(2) && max(xS((tmax/dt-50):end,3))>comm_width3+w_init+xN(3))
            behavior=0; %seaward growth
        elseif max(xS((tmax/dt-50):end,2))<=comm_width2+w_init+xN(2) && max(xS((tmax/dt-50):end,3))<=comm_width3+w_init+xN(3) && xH(end,2)-xH(1,2)==0 && xH(end,3)-xH(1,3)==0
            behavior=3; %hold the line
        elseif  (xH(end,2)-xH(1,2)<0 && xH(end,3)-xH(1,3)<0) && (isnan(R2)==0 || isnan(R3)==0)
            behavior=6; %slow retreat
        elseif xH(end,2)==0 && xH(end,3)==0 && isnan(R2)==1 && isnan(R3)==1 %xS(end,2)<=xH(1,2) && xS(end,3)<=xH(1,3)
            behavior=9; %retreat
        elseif (max(xS((tmax/dt-50):end,2))>comm_width2+w_init+xN(2) || max(xS((tmax/dt-50):end,3))>comm_width3+w_init+xN(3)) && ((max(xS((tmax/dt-50):end,2))<=comm_width2+w_init+xN(2) && xH(end,2)-xH(1,2)==0) || (max(xS((tmax/dt-50):end,3))<=comm_width3+w_init+xN(3) && xH(end,3)-xH(1,3)==0))
            behavior=1.5; %mixing: seaward growth/hold the line
        elseif ((max(xS((tmax/dt-50):end,2))<=comm_width2+w_init+xN(2) && xH(end,2)-xH(1,2)==0) || (max(xS((tmax/dt-50):end,3))<=comm_width3+w_init+xN(3) && xH(end,3)-xH(1,3)==0)) && (xH(end,2)-xH(1,2)<0 || xH(end,3)-xH(1,3)<0)
            behavior=4.5; %mixing: hold the line/slow retreat
        elseif (max(xS((tmax/dt-50):end,2))>comm_width2+w_init+xN(2) || max(xS((tmax/dt-50):end,3))>comm_width3+w_init+xN(3)) && (xH(end,2)-xH(1,2)<0 || xH(end,3)-xH(1,3)<0)
            behavior=10.5; %mixing: seaward growth/slow retreat
        end
              
        %% Nourishment Efficiency
        if efficiency(end-1)==0
            Eff=NaN;
        elseif efficiency(end-1)~=0
            Eff=efficiency(end-1);
        end

        %% Output Storage
        N_efficiency_vector(iR3)=Eff;
    end
    N_efficiency_storage(iR2,:)=N_efficiency_vector;
end

%% Figures
%NOTE: turn on figure 1 if you want to explore the nourishment efficiency
%for coordination; turn on figure 2 if you want to explore the nourishment
%efficiency for conservative non-coordination

figure (1) %figure 10d in paper
pcolor(phi_vec,gamma_vec,N_efficiency_storage)
shading flat
colormap jet
pbaspect([1 1 1])
set(gca,'FontSize',20)
xlabel('Background Erosion Rate')
ylabel({'Nourishment Efficiencies';'Sand Cost ($/m^3)'})
title('Coordination')
c1=colorbar;

% figure (2) %figure 10e in paper
% pcolor(phi_vec,gamma_vec,N_efficiency_storage)
% shading flat
% colormap jet
% pbaspect([1 1 1])
% set(gca,'FontSize',20)
% xlabel('Background Erosion Rate')
% ylabel('Sand Cost ($/m^3)')
% title('Non-coordination')
% c2=colorbar;

