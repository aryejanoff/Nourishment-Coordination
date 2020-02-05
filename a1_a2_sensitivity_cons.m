a2_vec=linspace(1e6/25^0.6,4e6/25^0.6,30);
a3_vec=linspace(1e6/25^0.6,4e6/25^0.6,30);
PV=linspace(0.1e6,1.5e6,20);
nn=length(a2_vec);
mm=length(a3_vec);
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
parfor ii=1:numel(a2_vec) 
    R2_coord_vector=zeros(1,mm);
    R3_coord_vector=zeros(1,mm);
    R2_uncoord_cons_vector=zeros(1,mm);
    R3_uncoord_cons_vector=zeros(1,mm);
    ratio_coord_vector=zeros(1,mm);
    ratio_uncoord_cons_vector=zeros(1,mm);
    TNB_coord_vector=zeros(1,mm);
    TNB_uncoord_cons_vector=zeros(1,mm);
    share2coord_vector=zeros(1,mm);
    share3coord_vector=zeros(1,mm);
    share2uncoord_cons_vector=zeros(1,mm);
    share3uncoord_cons_vector=zeros(1,mm);
    behaviorcoord_storage_vector=zeros(1,mm);
    behavioruncoord_storage_vector=zeros(1,mm);
    TNB_diff_vector=zeros(1,mm);
    for jj=1:numel(a3_vec)
        alpha2=a2_vec(ii);
        alpha3=a3_vec(jj);
        [Vratio_coord,Vratio_uncoord_cons,R2coord_star,R3coord_star,R2uncoord_star_cons,R3uncoord_star_cons,M_coord,M_uncoord_cons,share2_coord,share3_coord,share2_uncoord_cons,share3_uncoord_cons,behavior_coord,behavior_uncoord_cons]=maincode(alpha2,alpha3);
        %% Storage %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ratio_coord_vector(jj)=Vratio_coord(1);
        ratio_uncoord_cons_vector(jj)=Vratio_uncoord_cons;
        R2_coord_vector(jj)=R2coord_star(1);
        R3_coord_vector(jj)=R3coord_star(1);
        R2_uncoord_cons_vector(jj)=R2uncoord_star_cons;
        R3_uncoord_cons_vector(jj)=R3uncoord_star_cons;
        TNB_coord_vector(jj)=M_coord(1);
        TNB_uncoord_cons_vector(jj)=M_uncoord_cons;
        share2coord_vector(jj)=share2_coord(1);
        share3coord_vector(jj)=share3_coord(1);
        share2uncoord_cons_vector(jj)=share2_uncoord_cons;
        share3uncoord_cons_vector(jj)=share3_uncoord_cons;
        behaviorcoord_storage_vector(jj)=behavior_coord(1);
        behavioruncoord_storage_vector(jj)=behavior_uncoord_cons;
        TNB_diff_vector(jj)=TNB_coord_vector(jj)-TNB_uncoord_cons_vector(jj);
    end
    ratio_coord(ii,:)=ratio_coord_vector;
    ratio_uncoord_cons(ii,:)=ratio_uncoord_cons_vector;
    R2_coord(ii,:)=R2_coord_vector;
    R3_coord(ii,:)=R3_coord_vector;
    R2_uncoord_cons(ii,:)=R2_uncoord_cons_vector;
    R3_uncoord_cons(ii,:)=R3_uncoord_cons_vector;
    TNB_coord_(ii,:)=TNB_coord_vector;
    TNB_uncoord_cons(ii,:)=TNB_uncoord_cons_vector;
    share2coord(ii,:)=share2coord_vector;
    share3coord(ii,:)=share3coord_vector;
    share2uncoord_cons(ii,:)=share2uncoord_cons_vector;
    share3uncoord_cons(ii,:)=share3uncoord_cons_vector;
    behaviorcoord_storage(ii,:)=behaviorcoord_storage_vector;
    behavioruncoord_storage(ii,:)=behavioruncoord_storage_vector;
    TNB_diff(ii,:)=TNB_diff_vector;
end

%% figures

figure (1)
subplot(1,2,1)
pcolor(PV/1e6,PV/1e6,TNB_diff_cons/(1e7))
colormap(flipud(hot))
shading flat
pbaspect([1 1 1])
xlabel('Initial PV_2 ($10^6)')
ylabel('Initial PV_1 ($10^6)')
set(gca,'FontSize',20)
c1=colorbar;
c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^7)';
title({'Benefit of Coordination';'(cautionary)'})
subplot(1,2,2)
pcolor(PV/1e6,PV/1e6,TNB_diff_lib/(1e8))
colormap(flipud(hot))
shading flat
pbaspect([1 1 1])
xlabel('Initial PV_2 ($10^6)')
ylabel('Initial PV_1 ($10^6)')
set(gca,'FontSize',20)
c1=colorbar;
c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^8)';
title({'Benefit of Coordination';'(risky)'})

figure (2)
subplot(1,3,1)
pcolor(PV/1e6,PV/1e6,behaviorcoord_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('Initial PV_2 ($10^6)')
ylabel('Initial PV_1 ($10^6)')
set(gca,'FontSize',20)
title('Coordination')
caxis([0 10.5])
subplot(1,3,2)
pcolor(PV/1e6,PV/1e6,behavioruncoord_cons_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('Initial PV_2 ($10^6)')
ylabel('Initial PV_1 ($10^6)')
set(gca,'FontSize',20)
title({'Key Behaviors';'Non-coordination (cautionary)'})
caxis([0 10.5])
subplot(1,3,3)
pcolor(PV/1e6,PV/1e6,behavioruncoord_lib_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('Initial PV_2 ($10^6)')
ylabel('Initial PV_1 ($10^6)')
set(gca,'FontSize',20)
title('Non-coordination (risky)')
caxis([0 10.5])

% figure (3)
% subplot(2,3,1)
% pcolor((a3_vec*25^0.25)/1e6,(a2_vec*25^0.25)/1e6,R2_coord)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('PV_2 ($10^6)')
% ylabel({'Community 1';'PV_1 ($10^6)'})
% set(gca,'FontSize',15)
% c8=colorbar;
% c8.Label.String='R_1 (yrs)';
% set(gca,'FontSize',15)
% title('Coordination')
% caxis([0 30])
% subplot(2,3,4)
% pcolor((a3_vec*25^0.25)/1e6,(a2_vec*25^0.25)/1e6,R3_coord)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('PV_2 ($10^6)')
% ylabel({'Community 2';'PV_1 ($10^6)'})
% set(gca,'FontSize',15)
% c9=colorbar;
% c9.Label.String='R_2 (yrs)';
% caxis([0 30])
% subplot(2,3,2)
% pcolor((a3_vec*25^0.25)/1e6,(a2_vec*25^0.25)/1e6,R2_uncoord_cons)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('PV_2 ($10^6)')
% ylabel('PV_1 ($10^6)')
% title({'Non-coordination';'neighbor does nothing'})
% set(gca,'FontSize',15)
% c10=colorbar;
% c10.Label.String='R_1 (yrs)';
% caxis([0 30])
% subplot(2,3,5)
% pcolor((a3_vec*25^0.25)/1e6,(a2_vec*25^0.25)/1e6,R3_uncoord_cons)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('PV_2 ($10^6)')
% ylabel('PV_1 ($10^6)')
% set(gca,'FontSize',15)
% c11=colorbar;
% c11.Label.String='R_2 (yrs)';
% caxis([0 30])
% subplot(2,3,3)
% pcolor((a3_vec*25^0.25)/1e6,(a2_vec*25^0.25)/1e6,R2_uncoord_lib)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('PV_2 ($10^6)')
% ylabel('PV_1 ($10^6)')
% title({'Non-coordination';'neighbor nourishes'})
% set(gca,'FontSize',15)
% c12=colorbar;
% c12.Label.String='R_1 (yrs)';
% caxis([0 30])
% subplot(2,3,6)
% pcolor((a3_vec*25^0.25)/1e6,(a2_vec*25^0.25)/1e6,R3_uncoord_lib)
% colormap(flipud(hot))
% shading flat
% pbaspect([1 1 1])
% xlabel('PV_2 ($10^6)')
% ylabel('PV_1 ($10^6)')
% set(gca,'FontSize',15)
% c13=colorbar;
% c13.Label.String='R_2 (yrs)';
% caxis([0 30])
% 
