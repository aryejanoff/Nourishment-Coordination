load('future_erosion_sandcost_DATA.mat')

figure (1) %figure 10a-b in paper
subplot(1,2,1)
pcolor(phi_vec,gamma_vec,behaviorcoord_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
title('Coordination')
caxis([0 10.5])
subplot(1,2,2)
pcolor(phi_vec,gamma_vec,behavioruncoord_cons_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
title('Non-coordination')
caxis([0 10.5])

figure (2) %figure A2a-b in paper
subplot(1,2,1)
pcolor(phi_vec,gamma_vec,behaviorcoord_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
title('Coordination')
caxis([0 10.5])
subplot(1,2,2)
pcolor(phi_vec,gamma_vec,behavioruncoord_risk_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
title('Non-coordination')
caxis([0 10.5])

figure (3) %figure A2c in paper
pcolor(phi_vec,gamma_vec,TNB_diff_risk/(1e6))
colormap(flipud(hot))
shading flat
pbaspect([1 1 1])
xlabel('\phi_N ($/m^3)')
ylabel('\gamma (m/yr)')
set(gca,'FontSize',20)
c1=colorbar;
c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^6)';
title('Benefit of Coordination')
