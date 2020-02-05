gamma_vec=linspace(0.5,2.5,30);
phi_vec=linspace(5,20,30);
PV=linspace(1e6,4e6,30);
nn=length(gamma_vec);
mm=length(phi_vec);
R2_coord=zeros(nn,mm);
R3_coord=zeros(nn,mm);
R2_uncoord_cons=zeros(nn,mm);
R3_uncoord_cons=zeros(nn,mm);
ratio_coord=zeros(nn,mm);
ratio_uncoord_cons=zeros(nn,mm);
TNB_coord=zeros(nn,mm);
TNB_uncoord_cons=zeros(nn,mm);
share2coord=zeros(nn,mm);
share3coord=zeros(nn,mm);
share2uncoord_cons=zeros(nn,mm);
share3uncoord_cons=zeros(nn,mm);
behaviorcoord_storage=zeros(nn,mm);
behavioruncoord_storage=zeros(nn,mm);
TNB_diff=zeros(nn,mm);

%% Main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor ii=1:numel(gamma_vec) 
    for jj=1:numel(phi_vec)
        gamma=gamma_vec(ii);
        phi=phi_vec(jj);
        [Vratio_coord,Vratio_uncoord_cons,Vratio_uncoord,R2coord_star,R3coord_star,R2uncoord_star_cons,R3uncoord_star_cons,R2uncoord_star,R3uncoord_star,M_coord,M_uncoord_cons,M_uncoord,share2_coord,share3_coord,share2_uncoord_cons,share3_uncoord_cons,share2_uncoord,share3_uncoord,behavior_coord,behavior_uncoord_cons,behavior_uncoord,VLlost_coord,VLlost_uncoord_cons,VLlost_uncoord,VL_equaleffort_coord]=maincode(alpha2,alpha3);
        %% Storage %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ratio_coord(ii,jj)=Vratio_coord(1);
        ratio_uncoord_cons(ii,jj)=Vratio_uncoord_cons;
        R2_coord(ii,jj)=R2coord_star(1);
        R3_coord(ii,jj)=R3coord_star(1);
        R2_uncoord_cons(ii,jj)=R2uncoord_star_cons;
        R3_uncoord_cons(ii,jj)=R3uncoord_star_cons;
        TNB_coord(ii,jj)=M_coord(1);
        TNB_uncoord_cons(ii,jj)=M_uncoord_cons;
        share2coord(ii,jj)=share2_coord(1);
        share3coord(ii,jj)=share3_coord(1);
        share2uncoord_cons(ii,jj)=share2_uncoord_cons;
        share3uncoord_cons(ii,jj)=share3_uncoord_cons;
        behaviorcoord_storage(ii,jj)=behavior_coord(1);
        behavioruncoord_storage(ii,jj)=behavior_uncoord_cons;
        TNB_diff(ii,jj)=TNB_coord(ii,jj)-TNB_uncoord_cons(ii,jj);
    end
end

%% figures
% figure (2)
% subplot(1,2,1)
% pcolor(phi_vec,gamma_vec,TNB_diff_cons/(1e6))
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% set(gca,'FontSize',20)
% c1=colorbar;
% c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^6)';
% title({'Benefit of Coordination';'(cautionary)'})
% subplot(1,2,2)
% pcolor(phi_vec,gamma_vec,TNB_diff_lib/(1e7))
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% set(gca,'FontSize',20)
% c1=colorbar;
% c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^7)';
% title({'Benefit of Coordination';'(risky)'})
% 
% figure (4)
% subplot(2,3,1)
% pcolor(phi_vec,gamma_vec,R2_coord)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel({'Community 1';'\gamma (m/yr)'})
% set(gca,'FontSize',15)
% c8=colorbar;
% c8.Label.String='R_1 (yrs)';
% set(gca,'FontSize',15)
% title('Coordination')
% caxis([0 55])
% subplot(2,3,4)
% pcolor(phi_vec,gamma_vec,R3_coord)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel({'Community 2';'\gamma (m/yr)'})
% set(gca,'FontSize',15)
% c9=colorbar;
% c9.Label.String='R_2 (yrs)';
% caxis([0 55])
% subplot(2,3,2)
% pcolor(phi_vec,gamma_vec,R2_uncoord_cons)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% title({'Non-coordination';'neighbor does nothing'})
% set(gca,'FontSize',15)
% c10=colorbar;
% c10.Label.String='R_1 (yrs)';
% caxis([0 55])
% subplot(2,3,5)
% pcolor(phi_vec,gamma_vec,R3_uncoord_cons)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% set(gca,'FontSize',15)
% c11=colorbar;
% c11.Label.String='R_2 (yrs)';
% caxis([0 55])
% subplot(2,3,3)
% pcolor(phi_vec,gamma_vec,R2_uncoord_lib)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% title({'Non-coordination';'neighbor nourishes'})
% set(gca,'FontSize',15)
% c12=colorbar;
% c12.Label.String='R_1 (yrs)';
% caxis([0 55])
% subplot(2,3,6)
% pcolor(phi_vec,gamma_vec,R3_uncoord_lib)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% set(gca,'FontSize',15)
% c13=colorbar;
% c13.Label.String='R_2 (yrs)';
% caxis([0 55])
% 
% figure (5)
% subplot(1,3,1)
% pcolor(phi_vec,gamma_vec,behaviorcoord_storage)
% colormap hot
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% set(gca,'FontSize',20)
% title('Coordinated Behaviors')
% caxis([0 9])
% subplot(1,3,2)
% pcolor(phi_vec,gamma_vec,behavioruncoord_cons_storage)
% colormap hot
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% set(gca,'FontSize',20)
% title({'Uncoordinated Behaviors';'neighbor does nothing'})
% caxis([0 9])
% subplot(1,3,3)
% pcolor(phi_vec,gamma_vec,behavioruncoord_lib_storage)
% colormap hot
% shading flat
% pbaspect([1 1 1])
% xlabel('\phi_N ($/m^3)')
% ylabel('\gamma (m/yr)')
% set(gca,'FontSize',20)
% title({'Uncoordinated Behaviors';'neighbor nourishes'})
% caxis([0 9])

figure (1)
subplot(1,3,1)
pcolor(phi_vec,gamma_vec,behaviorcoord_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
title('Coordination')
caxis([0 10.5])
subplot(1,3,2)
pcolor(phi_vec,gamma_vec,behavioruncoord_cons_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
title({'Future Behaviors (Wealth Symmetry)';'Non-coordination (cautionary)'})
caxis([0 10.5])
subplot(1,3,3)
pcolor(phi_vec,gamma_vec,behavioruncoord_lib_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
title('Non-coordination (risky)')
caxis([0 10.5])

figure (2)
subplot(1,2,1)
pcolor(phi_vec,gamma_vec,TNB_diff_cons/(1e7))
colormap(flipud(hot))
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
c1=colorbar;
c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^7)';
title({'Benefit of Coordination';'(cautionary)'})
subplot(1,2,2)
pcolor(phi_vec,gamma_vec,TNB_diff_lib/(1e8))
colormap(flipud(hot))
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
c1=colorbar;
c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^8)';
title({'Benefit of Coordination';'(risky)'})


