function [Vratio_coord,Vratio_uncoord_cons,Vratio_uncoord,R2coord_star,R3coord_star,R2uncoord_star_cons,R3uncoord_star_cons,R2uncoord_star,R3uncoord_star,M_coord,M_uncoord_cons,M_uncoord,share2_coord,share3_coord,share2_uncoord_cons,share3_uncoord_cons,share2_uncoord,share3_uncoord,behavior_coord,behavior_uncoord_cons,behavior_uncoord,Vlost_coord,VLlost_uncoord_cons,VLlost_uncoord,V_equaleffort_coord]=maincode(d)


% Rotation Length Sensitivities
R_min=0.5;
R_max=50;
R2_vec=linspace(R_min,R_max,R_max/R_min); R3_vec=linspace(R_min,R_max,R_max/R_min);
R2_vector=[NaN R2_vec(1:end)]; R3_vector=[NaN R3_vec(1:end)];
n2=length(R2_vector); n3=length(R3_vector);
NB2_storage_R=NaN(n2,n3); NB3_storage_R=NaN(n2,n3); TNB_storage_R=NaN(n2,n3);
N_efficiency_storage=NaN(n2,n3);
share2_storage=NaN(n2,n3); share3_storage=NaN(n2,n3);
behavior_storage=NaN(n2,n3);
VL_lost_storage=NaN(n2,n3); VC_lost_storage=NaN(n2,n3); V_lost_storage=NaN(n2,n3); V_nourish_storage=NaN(n2,n3);
TnE2_storage=NaN(n2,n3); TnE3_storage=NaN(n2,n3); TnE_storage=NaN(n2,n3);

parfor iR2=1:numel(R2_vector)
    NB2_vector_R=NaN(1,n3); NB3_vector_R=NaN(1,n3); TNB_vector_R=NaN(1,n3);
    N_efficiency_vector=NaN(1,n3);
    share2_vector=NaN(1,n3); share3_vector=NaN(1,n3);
    behavior_vector=NaN(1,n3);
    VL_lost_vector=NaN(1,n3); VC_lost_vector=NaN(1,n3); V_lost_vector=NaN(1,n3); V_nourish_vector=NaN(1,n3);
    TnE2_vector=NaN(1,n3); TnE3_vector=NaN(1,n3); TnE_vector=NaN(1,n3);
    
    for iR3=1:numel(R3_vector)
        %% Input Physical Parameters %%
        lot_size=25; 
        w_init=25;
        beta=0.25; %beach width hedonic parameter 
        alpha2=1e6/(w_init^beta);
        alpha3=1e6/(w_init^beta); 
        s=1000; %alongshore compartment length (m)
        rows_along=s/lot_size;
        rows_cross=1; %# of cross-shore proeprty rows
        properties_total=rows_cross*rows_along;
        comm_width=rows_cross*lot_size; %initial Community Width (m)
        psi=0.2;
        D=10; %depth of closure (m)
        gamma=1; %erosion rate (m/yr)
        d=10000; %alongshore flux coeff
        K=2000; %cross-shore flux coeff
        phi=20; %sand cost ($/m^3) 
        c=2e6; %fixed nourishment cost ($)
        rho=0.06; %discount rate
        R2=11.2; %R2_vector(iR2); %Rotation Length
        R3=7.7; %R3_vector(iR3); %Rotation Length
        xN=50; %nourishment magnitude (m)
        nu=0; %beach width decline beyond max threshold
        theta_eq=0.02; %equilibrium shoreface slope
        k2=0;
        k3=0;
        %% Computational Parameters %%
        tmax=100; dt=0.1; t=0:dt:tmax; n=length(t); 
        Smax=3; ds=1; S=1:ds:Smax; m=length(S);
        A2=alpha2*rho; A3=alpha3*rho; 

        theta=zeros(n,m); qL=zeros(n,m); qC=zeros(n,m); fS=zeros(n,m); fT=zeros(n,m); 
        xS=zeros(n,m); xT=zeros(n,m); xH=zeros(n,m); w=zeros(n,m);
        B=zeros(n,m); B_disc=zeros(n,m); TB=zeros(n,m);
        C=zeros(n,m); TC=zeros(n,m);
        nb=NaN(n,m); N_h=zeros(n,m);        
        volume=zeros(n,m); nE=zeros(n,m); Vtotal=zeros(n,m);
        TQL_out=zeros(1,n); TQC_out=zeros(1,n); TQ_out=zeros(1,n);
        share2=NaN(1); share3=NaN(1); behavior=NaN(1); NB2=NaN(1); NB3=NaN(1); TNB=NaN(1);

        %% Initial Conditions %%
        xS(1,:)=comm_width+w_init; 
        xT(1,:)=xS(1,:)+(D/theta_eq);
        xH(1,:)=comm_width;
        N_h(1,:)=rows_cross;
        w(1,:)=xS(1,:)-xH(1,:);

        %% Main Code %%
        for i=1:n-1
            for j=2:m-1                           
                %% Nourishment Initiation + Volume
                ff2=round(t(i)-k2*R2,4);
                if ff2==0
                    k2=k2+1; volume(i,2)=0.5*s*D*xN; nE(i,2)=1;
                end
                Vtotal(i+1,2)=Vtotal(i,2)+volume(i,2);
                ff3=round(t(i)-k3*R3,4);
                if ff3==0
                    k3=k3+1; volume(i,3)=0.5*s*D*xN; nE(i,3)=1;
                end
                Vtotal(i+1,3)=Vtotal(i,3)+volume(i,3);
                
                %% Fluxes (Along/Cross-shore) and Shoreface Dynamics
                qL(i,j)=d*D*((xS(i,j-1)-xS(i,j))/s); qL(i,1)=d*D*((xS(i,m-1)-xS(i,1))/s); qL(i,m)=d*D*((xS(i,m-1)-xS(i,m))/s);
                theta(i,j)=D/(xT(i,j)-xS(i,j)); theta(i,1)=D/(xT(i,1)-xS(i,1)); theta(i,m)=D/(xT(i,m)-xS(i,m));
                qC(i,j)=4*s*K*(theta(i,j)-theta_eq); qC(i,1)=4*s*K*(theta(i,1)-theta_eq); qC(i,m)=4*s*K*(theta(i,m)-theta_eq);
                TQC_out(i)=-sum(qC(:,2))-sum(qC(:,3)); TQL_out(i)=sum(qL(:,2))-sum(qL(:,m)); TQ_out(i)=TQC_out(i)+TQL_out(i);
                
                %% ODE's
                fT(i,j)=(4*K*(theta(i,j)-theta_eq)/D)-gamma; fT(i,1)=(4*K*(theta(i,1)-theta_eq)/D)-gamma; fT(i,m)=(4*K*(theta(i,m)-theta_eq)/D)-gamma;
                fS(i,j)=((2*d/(s^2))*(xS(i,j-1)-2*xS(i,j)+xS(i,j+1)))-(4*K*(theta(i,j)-theta_eq)/D)-gamma; fS(i,1)=((2*d/(s^2))*(xS(i,m-1)-2*xS(i,1)+xS(i,2)))-(4*K*(theta(i,1)-theta_eq)/D)-gamma; fS(i,m)=fS(i,1); 

                %% Numerical Approximations

                xT(i+1,j)=xT(i,j)+dt*fT(i,j); xT(i+1,1)=xT(i,1)+dt*fT(i,1); xT(i+1,m)=xT(i,m)+dt*fT(i,m);
                xS(i+1,j)=xS(i,j)+dt*fS(i,j); xS(i+1,1)=xS(i,1)+dt*fS(i,1); xS(i+1,m)=xS(i,m)+dt*fS(i,m); 
                if volume(i,2)~=0
                    xS(i+1,2)=xS(i,2)+xN;
                end
                if volume(i,3)~=0
                    xS(i+1,3)=xS(i,3)+xN;
                end
                if xS(i,j)<=lot_size
                    N_h(i,j)=0; xH(i,j)=0; xS(i,j)=0; N_h(i+1,j)=0; xH(i+1,j)=0;
                elseif xS(i,j)<=xH(i,j)+0.5
                    xH(i+1,j)=xH(i,j)-lot_size; N_h(i+1,j)=N_h(i,j)-1;
                else
                    xH(i+1,j)=xH(i,j);
                end
                if xS(i,m)<=lot_size
                    N_h(i,m)=0; xH(i,m)=0; xS(i,m)=0; N_h(i+1,m)=0; xH(i+1,m)=0;
                elseif xS(i,j)<=xH(i,m)+0.5
                    xH(i+1,m)=xH(i,m)-lot_size; N_h(i+1,m)=N_h(i,m)-1;
                else
                    xH(i+1,m)=xH(i,m);
                end
                w(i,j)=xS(i,j)-xH(i,j); 

                %% Housing Lines 
                N_h(i,j)=xH(i,j)/lot_size; 
                
                %% Benefit 
                B(i,2)=rows_along*A2*((w(i,2)).^beta)*((N_h(i,2)).^psi)-(nu*(w(i,2).^2)); B(i,3)=rows_along*A3*((w(i,3)).^beta)*((N_h(i,3)).^psi)-(nu*(w(i,3).^2));
                B_disc(i,j)=B(i,j)*exp(-rho*t(i));
                TB(i,j)=(dt/3)*(B_disc(1,j)+4*(sum(B_disc(2:2:end,j)))+2*sum(B_disc(2:1:end,j))+B_disc(end,j));

                %% Cost 
                C(i,j)=(c+volume(i,j)*phi)*exp(-rho*t(i)); 
                if volume(i,j)==0
                    C(i,j)=0;
                end
                TC(i,j)=sum(C(:,j));

                %% marginal net benefit 
                nb(i,j)=TB(i,j)-TC(i,j);

                %% Net Benefit
                NB2=nb(i,2); NB3=nb(i,3); TNB=NB2+NB3;
                share2=NB2/TNB; share3=NB3/TNB;
            end
        end
        
        tol_eq=0.5;
        %% Identify Behavior
        if (max(xS(501:end,2))>xS(2,2)+w_init || max(xS(501:end,3))>xS(2,3)+w_init) && abs(R2-R3)<=tol_eq && isnan(R2)==0 && isnan(R3)==0
            behavior=0; %equal effort and seaward growth
        elseif (max(xS(501:end,2))>xS(2,2)+w_init || max(xS(501:end,3))>xS(2,3)+w_init) && abs(R2-R3)>tol_eq && isnan(R2)==0 && isnan(R3)==0
            behavior=1; %unequal effort and seaward growth
        elseif (max(xS(501:end,2))>xS(2,2)+w_init || max(xS(501:end,3))>xS(2,3)+w_init) && ((isnan(R2)==0 && isnan(R3)==1) || (isnan(R2)==1 && isnan(R3)==0))
            behavior=2; %free riding and seaward growth
        elseif max(xS(501:end,2))<=xS(2,2)+w_init && max(xS(501:end,3))<=xS(2,3)+w_init && xH(end,2)-xH(1,2)==0 && xH(end,3)-xH(1,3)==0 && abs(R2-R3)<=tol_eq && isnan(R2)==0 && isnan(R3)==0
            behavior=3; %equal effort and maintain position
        elseif max(xS(501:end,2))<=xS(2,2)+w_init && max(xS(501:end,3))<=xS(2,3)+w_init && xH(end,2)-xH(1,2)==0 && xH(end,3)-xH(1,3)==0 && abs(R2-R3)>tol_eq && isnan(R2)==0 && isnan(R3)==0
            behavior=4; %unequal effort and maintain position
        elseif max(xS(501:end,2))<=xS(2,2)+w_init && max(xS(501:end,3))<=xS(2,3)+w_init && xH(end,2)-xH(1,2)==0 && xH(end,3)-xH(1,3)==0 && ((isnan(R2)==0 && isnan(R3)==1) || (isnan(R2)==1 && isnan(R3)==0))
            behavior=5; %free riding and maintain position
        elseif  (xH(end,2)-xH(1,2)<0 || xH(end,3)-xH(1,3)<0) && abs(R2-R3)<=tol_eq && isnan(R2)==0 && isnan(R3)==0
            behavior=6; %equal effort and slow retreat
        elseif (xH(end,2)-xH(1,2)<0 || xH(end,3)-xH(1,3)<0) && abs(R2-R3)>tol_eq && isnan(R2)==0 && isnan(R3)==0
            behavior=7; %unequal effort and slow retreat
        elseif (xH(end,2)-xH(1,2)<0 || xH(end,3)-xH(1,3)<0) && ((isnan(R2)==0 && isnan(R3)==1) || (isnan(R2)==1 && isnan(R3)==0))
            behavior=8; %free riding and slow retreat
        elseif xS(end,2)<=xH(1,2) && xS(end,3)<=xH(1,3) && isnan(R2)==1 && isnan(R3)==1
            behavior=9; %max retreat
        end

        %% Sediment Loss
        I90=round(log(1/0.1)/rho/dt);
        QN=Vtotal(end-1,2)+Vtotal(end-1,3)
        QL=TQL_out(end-1)
        QC=TQC_out(end-1)
        Q=TQ_out(end-1)
        Eff=Q/QN

        %% V_gamma vs. V_n
        V_gamma=2*gamma*tmax*s*D;
        V_n=Vtotal(end,2)+Vtotal(end,3);
        ratio=V_n/V_gamma;

        %% # nourishment episodes
        TnE2=sum(nE(:,2));
        TnE3=sum(nE(:,3));

        %% Output Storage
        VL_lost_vector(iR3)=QL;
        VC_lost_vector(iR3)=QC;
        V_lost_vector(iR3)=QL+QC;
        V_nourish_vector(iR3)=QN;
        N_efficiency_vector(iR3)=Eff;
        TnE2_vector(iR3)=TnE2;
        TnE3_vector(iR3)=TnE3;
        TnE_vector(iR3)=TnE2+TnE3;
        NB2_vector_R(iR3)=NB2;
        NB3_vector_R(iR3)=NB3;
        TNB_vector_R(iR3)=TNB;
        share2_vector(iR3)=share2;
        share3_vector(iR3)=share3;
        behavior_vector(iR3)=behavior;
    end
    VL_lost_storage(iR2,:)=VL_lost_vector;
    VC_lost_storage(iR2,:)=VC_lost_vector;
    V_lost_storage(iR2,:)=V_lost_vector;    
    V_nourish_storage(iR2,:)=V_nourish_vector;
    N_efficiency_storage(iR2,:)=N_efficiency_vector;
    TnE2_storage(iR2,:)=TnE2_vector;
    TnE3_storage(iR2,:)=TnE3_vector;
    TnE_storage(iR2,:)=TnE2_vector+TnE3_vector;
    NB2_storage_R(iR2,:)=NB2_vector_R;
    NB3_storage_R(iR2,:)=NB3_vector_R;
    TNB_storage_R(iR2,:)=TNB_vector_R;
    share2_storage(iR2,:)=share2_vector;
    share3_storage(iR2,:)=share3_vector;
    behavior_storage(iR2,:)=behavior_vector;
end
%% Clear Variables
% clearvars -except TNB_storage_R R2_vector R3_vector col_coord row_coord col3 row2 col3_cons row2_cons

% find optimum
% coordinated
M_coord=max(TNB_storage_R(:));
[row_coord,col_coord]=find(TNB_storage_R==M_coord);
R2coord_star=R2_vector(row_coord);
R3coord_star=R3_vector(col_coord);
Vratio_coord=N_efficiency_storage(row_coord,col_coord);
share2_coord=share2_storage(row_coord,col_coord);
share3_coord=share3_storage(row_coord,col_coord);
behavior_coord=behavior_storage(row_coord,col_coord);
Vlost_coord=VL_lost_storage(row_coord(1),col_coord(1));
V_diag=diag(VL_lost_storage);
V_equaleffort_coord=(V_diag(row_coord(1))+V_diag(col_coord(1)))/2;

%uncoordinated
M2=max(NB2_storage_R(:));
[row2,col2]=find(NB2_storage_R==M2);
R2uncoord_star=R2_vector(row2);
R2_uncoord_assump=R3_vector(col2);
M2_cons=max(NB2_storage_R(:,1));
[row2_cons,col2_cons]=find(NB2_storage_R==M2_cons);
R2uncoord_star_cons=R2_vector(row2_cons);
R2_uncoord_assump_cons=R3_vector(col2_cons);

M3=max(NB3_storage_R(:));
[row3,col3]=find(NB3_storage_R==M3);
R3uncoord_star=R3_vector(col3);
R3_uncoord_assump=R2_vector(row3);
M3_cons=max(NB3_storage_R(1,:));
[row3_cons,col3_cons]=find(NB3_storage_R==M3_cons);
R3uncoord_star_cons=R3_vector(col3_cons);
R3_uncoord_assump_cons=R2_vector(row3_cons);

Vratio_uncoord=N_efficiency_storage(row2,col3);
Vratio_uncoord_cons=N_efficiency_storage(row2_cons,col3_cons);
M_uncoord=TNB_storage_R(row2,col3);
M_uncoord_cons=TNB_storage_R(row2_cons,col3_cons);
share2_uncoord=share2_storage(row2,col3);
share2_uncoord_cons=share2_storage(row2_cons,col3_cons);
share3_uncoord=share3_storage(row2,col3);
share3_uncoord_cons=share3_storage(row2_cons,col3_cons);
behavior_uncoord_cons=behavior_storage(row2_cons,col3_cons);
behavior_uncoord=behavior_storage(row2,col3);
VLlost_uncoord=VL_lost_storage(row2,col3);
VLlost_uncoord_cons=VL_lost_storage(row2_cons,col3_cons);


% %% figures
% figure (1)
% hold on
% pcolor(TNB_storage_R/(1e9))
% colormap(flipud(hot))
% shading flat
% xlim([1 length(R3_vector)])
% ylim([1 length(R3_vector)])
% xlabel('R_2 (yr)')
% ylabel('R_1 (yr)')
% xticks([1 6 11 16 21 26 31])
% xticklabels({'0','5','10','15','20','25','30'})
% yticks([1 6 11 16 21 26 31])
% yticklabels({'0','5','10','15','20','25','30'})
% title({'Optimal Rotation Lengths';'PV_1 = $1.9 million; PV_2 = $2.5 million'})
% cbar = colorbar;
% cbar.Label.String = 'TNB ($10^9)';
% pbaspect([1 1 1])
% set(gca,'FontSize',20)
% plot(col_coord,row_coord,'co','linewidth',6)
% plot(col3,row2,'go','linewidth',6)
% plot(col3_cons,row2_cons,'bo','linewidth',6)
% pbaspect([1 1 1])

% figure (2)
% subplot(1,2,1)
% hold on
% contourf(NB2_storage_R)
% colormap(flipud(hot))
% plot(col_coord,row_coord,'co','linewidth',6)
% plot(col2,row2,'bo','linewidth',6)
% xlabel('R_2 (yr)')
% ylabel('R_1 (yr)')
% title('R_1^*')
% cbar = colorbar;
% cbar.Label.String = 'NB_1 ($)';
% pbaspect([1 1 1])
% set(gca,'FontSize',20)
% subplot(1,2,2)
% hold on
% contourf(NB3_storage_R)
% colormap(flipud(hot))
% plot(col_coord,row_coord,'co','linewidth',6)
% plot(col3,row3,'bo','linewidth',6)
% xlabel('R_2 (yr)')
% ylabel('R_1 (yr)')
% title('R_2^*')
% cbar = colorbar;
% cbar.Label.String = 'NB_2 ($)';
% pbaspect([1 1 1])
% set(gca,'FontSize',20)