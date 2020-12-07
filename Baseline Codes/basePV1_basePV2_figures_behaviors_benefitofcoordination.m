load('basePV1_basePV2_sensitivity_analyses_DATA.mat')

figure (1) %figure 5a-b in paper
subplot(1,2,1)
pcolor(a3_vec/1e3,a2_vec/1e3,behaviorcoord_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^5)')
ylabel('Baseline PV_1 ($10^5)')
set(gca,'FontSize',20)
title('Coordination')
caxis([0 10.5])
subplot(1,2,2)
pcolor(a3_vec/1e3,a2_vec/1e3,behavioruncoord_cons_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^5)')
ylabel('Baseline PV_1 ($10^5)')
set(gca,'FontSize',20)
title('Non-coordination')
caxis([0 10.5])

figure (3) %figure 5e in paper
pcolor(a3_vec/1e3,a2_vec/1e3,TNB_diff_cons/(1e6))
colormap(flipud(hot))
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^5)')
ylabel('Baseline PV_1 ($10^5)')
set(gca,'FontSize',20)
c1=colorbar;
c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^6)';
title('Benefit of Coordination')

figure (4) %figure A1c in paper
pcolor(a3_vec/1e3,a2_vec/1e3,TNB_diff_risk/(1e6))
colormap(flipud(hot))
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^5)')
ylabel('Baseline PV_1 ($10^5)')
set(gca,'FontSize',20)
c1=colorbar;
c1.Label.String='TNB_{coord} - TNB_{non-coord} ($10^6)';
title('Benefit of Coordination')

figure (5) %figure A1a-b in paper
subplot(1,2,1)
pcolor(a3_vec/1e3,a2_vec/1e3,behaviorcoord_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^5)')
ylabel('Baseline PV_1 ($10^5)')
set(gca,'FontSize',20)
title('Coordination')
caxis([0 10.5])
subplot(1,2,2)
pcolor(a3_vec/1e3,a2_vec/1e3,behavioruncoord_risk_storage)
colormap jet
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^5)')
ylabel('Baseline PV_1 ($10^5)')
set(gca,'FontSize',20)
title('Non-coordination')
caxis([0 10.5])
