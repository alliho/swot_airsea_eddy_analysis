

%% [save] Umag vs theta | quick test just for pcolor
nmin = 5;
Lavg = 6; 
Rmax = 200;

% X = TH_
rotfldstr = '';
run load_vars_into_memory_universal
msk = ones(size(X)); msk(mag(X,Y)>Rmax) = NaN;
TH = cart2pol(U,V); TH = rad2deg(TH);
DTH = TH - TH_wv; %DTH = mod(DTH,360);
DTH = abs(DTH); DTH(DTH>360) = DTH(DTH>360) - 360;
DTH(DTH>180) = 360 - DTH(DTH>180);
x = DTH; y = mag(U,V);
z = SWHa./nanmean(SWH,[1 2]); z = remove_planefit(X,Y,z);
% z = SWHa; z = remove_planefit(X,Y,z);

A = {x,y,z}; 
A = cellfun(@(Z) cell2matn(arrayfun(@(i)    gaussfilt_2d(Z(:,:,i),eddygrid.drkm, Lavg), 1:size(Z,3), 'Un', 0),3), A, 'Un', 0);
A = cellfun(@(Z) Z.*msk, A, 'Un', 0);
A = cellfun(@(Z) Z(:), A, 'Un', 0);
[x y z] = A{:};
nn = find(~isnan(x) & ~isnan(y) & ~isnan(z));
A = cellfun(@(Z) Z(nn), A, 'Un', 0);
[x y z] = A{:};



figure(232); clf;
% setfigsize(gcf, [500 600])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 2, [0.025 0.025], [0.14 0.05], [0.12 0.19]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
format_fig4print(ha);
% setaxes(ha, 'DataAspectRatio', [1 1 1]); 
setaxes(ha, 'XTick', [0:30:180]); 
setaxes(ha, 'YTick', [0:0.5:3]); 
setaxes(ha, 'CLim', [-1 1].*0.25); 
setaxes(ha, 'Colormap', niceredblue); 
setaxes(ha, 'YLim', [0 2.4]); 
setaxes(ha, 'XLim', [0 180]); 
setaxes(ha, 'XLabel', '$|\Delta \theta|$'); 
setaxes(ha, 'YLabel', '$|U_{geo}|$ [m/s]'); 
% linkaxes(ha)
quietaxes(ha(2:end), 'y')
format_mapticks(ha, 'x');
% -------------------------------------------------------------------------
xbins = [0:2.5:180]; ybins = [0:0.05:3];
[M, N, xmid, ymid] = bin2D( xbins, ybins, x, y, z , 'stat', 'median');
[Xmid Ymid] = ndgrid(xmid,ymid); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
Z = M; 
Z(N<=nmin) = NaN;
pcolor(Xmid,Ymid,Z); shading flat;
Z = M; Z(N<=nmin) = NaN;
Z = movmean(Z,3,1, 'omitnan'); Z = movmean(Z,3,2, 'omitnan'); 
Z = movmean(Z,5,1, 'omitnan'); Z = movmean(Z,5,2, 'omitnan'); 
Z(N<=nmin) = NaN;
contour(Xmid,Ymid,Z, [-1:0.025:1], 'k-','EdgeAlpha', 0.3, 'LineWidth',1)

% %%% PHILLIPS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
% 
% contour(Xmid,Ymid,Z, [-1:0.025:1], '-','LineWidth',1)

%%% phillips
omp = 2*pi*0.1; 
kp = findk_u(omp, 4000,0,0);
cp = omp./kp;


% umagbins = [0:0.1:3];
% thbins = [0:180];
% [umag th] = ndgrid(umagbins, thbins);
[Ubins] = pol2cart(deg2rad(Xmid), Ymid);

cp_new = [0.5+0.5.*sqrt(1+4.*Ubins./cp)].*cp;

var = cp./sqrt(cp_new.*(cp_new + 2.*Ubins));
var(isnan(Z)) = NaN;
% contourf(Xmid,Ymid,var-1, [-1:0.05:1].*0.5, 'LineWidth',2)
contour(Xmid,Ymid,var-1, [-1:0.05:1].*0.5, 'LineWidth',2)

posx = thisha.Position(1) + 0.01; posy = sum(thisha.Position([2 4])) - 0.02;
textbypos(posx, posy, '$\frac{a-a_0}{a_0}$', 'Interpreter', 'latex')


% a/a0 = c0/sqrt(c*(c + 2U))
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c = fixedcolorbar(gca);
c.Position(1) = sum(thisha.Position([1 3])) + 0.03;
ylabel(c, 'SWHa/SWH', 'Interpreter', 'latex')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


savejpg(gcf, ['theta_vs_mag_vs_swha_Rmax' num2str(Rmax)], savepath, saveopt)




%% scatter variation over time

allera5 = 0;
errbopt = 1;
fs = 13;
prcnan = 0.25;
legvarargin = {'southwest','FontSize', max([fs-5 10]),'Interpreter', 'latex', 'NumColumns',3, 'Box', 'off'};
prcs = [5 95];
thwid = 35; rwid = [25 50];
thwid = 45; rwid = [25 50];
% -------------------------------------------------------------------------
wspd = squeeze(nanmean(WSP_model, [1 2])); 
swh = squeeze(nanmean(SWH_model, [1 2])); 

% wspd = naninterp1(era5.time, (squeeze(nanmean(mag(era5.u10, era5.v10),[1 2]))), t);
% swh = naninterp1(era5.time, (squeeze(nanmean(era5.swh,[1 2]))), t);



dth = squeeze(dirnanmean(abs(TH_wv_outer - TH_wnd_outer), [1 2])); %dth = abs(dth);
dth = squeeze(dirnanmean((TH_wv_outer - TH_wnd_outer), [1 2])); %dth = abs(dth);
wdw = dirnaninterp1(era5.time, rad2deg(squeeze(nanmean(era5.wdw,[1 2]))), t);

mwp = naninterp1(era5.time, (squeeze(nanmean(era5.mwp,[1 2]))), t);
% dswh = squeeze(nanmean(era5.swh,[1 2])); dswh = movmean(dswh,12); dswh = diff(dswh); 
% dswh = naninterp1(era5.time(2:end), dswh, t);
% dwspd = squeeze(nanmean(mag(era5.u10, era5.v10),[1 2])); dwspd = movmean(dwspd,12); dwspd = diff(dwspd); 
% dwspd = naninterp1(era5.time(2:end), dwspd, t);
% cp = om/k; om^2 = gk; k = om^2/g; cp = om/(om^2/g); cp = g./om; 
cp = g./(2*pi./mwp);
wag = cp./wspd;
vwag = cp./(wspd.*cosd(dth));

% -------------------------------------------------------------------------


rotfldstr = '_rot_wv';
run load_vars_into_memory
VARs = {U_wv, SWH,SWH./nanmean(SWH_model.*cmsk, [1 2])};
rotfldstr = '_rot_wnd';
run load_vars_into_memory       
VARs = [VARs {U_wnd, WSP,WSP./nanmean(WSP_model.*cmsk, [1 2])}];
rotfldstr = '';
run load_vars_into_memory       
VARs = [VARs {WSP,SWH}];



% VARs = {U_wv, SWH./nanmean(SWH_model.*cmsk, [1 2])};
VAR_anoms = cell(size(VARs));
tmp = NaN(2, length(t)); 
VAR_anoms(:) = {tmp}
for vi=1:length(VARs)
    VAR = VARs{vi};
    msk = ones(size(X)); msk(~isin(mag(X,Y), [75 200])) = NaN;
    VAR = VAR - nanmedian(VAR.*msk, [1 2]); 
    % VAR = remove_planefit(X,Y,VAR, 'smooth', 50, 'drkm', eddygrid.drkm);
    
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
    VAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;
    
    
    msk = ones(size(VAR)); msk(isnan(VAR)) = NaN; 
    VAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(VAR(:,:,i),eddygrid.drkm, Lavg), 1:size(VAR,3), 'Un', 0),3);
    VAR = cell2matn(arrayfun(@(i)    filt2(VAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(VAR,3), 'Un', 0),3);
    VAR = VAR.*msk;

    VARs{vi} = VAR;

    

    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = 90 + [-1 1].*thwid - 25;
    % thlims = 90 + [-1 1].*thwid - 15;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';
    
    VAR_anoms{vi}(1,:) = y;

    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = -90 + [-1 1].*thwid + 25;
    % thlims = -90 + [-1 1].*thwid + 15;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';

    VAR_anoms{vi}(2,:) = y;
end
    
figure(5); clf;
% setfigsize(gcf, [500 600])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 3, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
format_fig4print(ha);
% setaxes(ha, 'DataAspectRatio', [1 1 1]); 
% setaxes(ha, 'XLim', lonlims); 
setaxes(ha(2:3), 'YLim', [0 0.6]); 
setaxes(ha(3), 'XLim', [0 0.6]); 
% linkaxes(ha)
% format_mapticks(ha, 'xy');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
x = swh'; y = range(VAR_anoms{2});z = range(VAR_anoms{1});
nn = find(~isnan(x + y + z)); 
x = x(nn); y = y(nn); z = z(nn);
scatter(x,y,100,z,'filled', 'MarkerEdgeColor', 'k')
colorbar; colormap(flipud(cmocean('deep')))
plotBestOptFit(gca, x,y); 
cleanLegend(gca)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
x = range(VAR_anoms{1}); y = range(VAR_anoms{3}); z = cp';
% x = range(VAR_anoms{1}); y = range(VAR_anoms{2})./swh'; z = abs(dth');
% x = range(VAR_anoms{1}); y = range(VAR_anoms{2})./swh'; z = range(VAR_anoms{6});
nn = find(~isnan(x + y + z)); 
x = x(nn); y = y(nn); z = z(nn);
scatter(x,y,100,z,'filled', 'MarkerEdgeColor', 'k')
colorbar; colormap(flipud(cmocean('deep')))
plotBestOptFit(gca, x,y); 
cleanLegend(gca)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(3); axes(thisha); hold on; 



%%% phillips
% omp = 2*pi./[mwp+5]; 
omp = 2*pi./[11]; 
kp = findk_u(omp, 4000,0,0);
cp = omp./kp;


% umagbins = [0:0.1:3];
% thbins = [0:180];
% [umag th] = ndgrid(umagbins, thbins);
[Ubins] = pol2cart(deg2rad(Xmid), Ymid);

urng = range(VAR_anoms{1})'; 
urng = (VAR_anoms{1}(1,:) - VAR_anoms{1}(2,:))'; 


cp_new = [0.5+0.5.*sqrt(1+4.*urng./cp)].*cp;

var = cp./sqrt(cp_new.*(cp_new + 2.*urng));

x = var'-1; y = (VAR_anoms{3}(1,:) - VAR_anoms{3}(2,:)); z = range(VAR_anoms{5});

msk = ones(size(X)); msk(~isin(mag(X,Y), [0 100])) = NaN;
% z = arrayfun(@(ti) ccfxn(msk.*VARs{1}(:,:,ti), msk.*VARs{2}(:,:,ti)), 1:size(VARs{1},3)); caxis([-1 -0.3])
z = arrayfun(@(ti) ccfxn(msk.*VARs{7}(:,:,ti), msk.*VARs{8}(:,:,ti)), 1:size(VARs{1},3)); %caxis([-1 -0.3])
z = arrayfun(@(ti) ccfxn(msk.*VARs{4}(:,:,ti), msk.*VARs{5}(:,:,ti)), 1:size(VARs{1},3)); %caxis([-1 -0.3])
z = abs(z); caxis([0 1])

z = 1./vwag'; caxis([-0.5 1].*1.5)
nn = find(~isnan(x + y + z)); 
x = x(nn); y = y(nn); z = z(nn);
scatter(x,y,100,z,'filled', 'MarkerEdgeColor', 'k')
colorbar; colormap(gca,(nicejet))
plotBestOptFit(gca, x,y); 
cleanLegend(gca)
plotone2one(gca)
% x = range(VAR_anoms{1}); y = range(VAR_anoms{3}); z = cp';
% x = range(VAR_anoms{4}); y = range(VAR_anoms{5}); z = abs(dth)';
% % x = range(VAR_anoms{1}); y = range(VAR_anoms{2})./swh'; z = range(VAR_anoms{6});
% nn = find(~isnan(x + y + z)); 
% x = x(nn); y = y(nn); z = z(nn);
% scatter(x,y,100,z,'filled', 'MarkerEdgeColor', 'k')
% colorbar; colormap(flipud(cmocean('deep')))




% 
% col = scaledcmap(cmap,nanmedian(y), clims);
% y = y + y0; lims = lims + y0;
% shaded_limits(x,y,lims, col, 'LineWidth', 0);
% n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])
% 
% col = [0.2 0.1 0.6]; 
% msk = NaN(size(X)); msk(isin(mag(X,Y), rwid+10)) = 1; 
% thlims = -90 + [-1 1].*thwid + 25;
% msk(~isin(atan2d(Y,X), thlims)) = NaN;









%% [save] combine scatter vs pcolor

nmin = 5;
Lavg = 6; 
Rmax = 200;
daterange = [datenum(2023,4,13) datenum(2023,5,12)];

% X = TH_
rotfldstr = '';
run load_vars_into_memory_universal
msk = ones(size(X)); msk(mag(X,Y)>Rmax) = NaN;
TH = cart2pol(U,V); TH = rad2deg(TH);
DTH = TH - TH_wv; %DTH = mod(DTH,360);
DTH = abs(DTH); DTH(DTH>360) = DTH(DTH>360) - 360;
DTH(DTH>180) = 360 - DTH(DTH>180);
x = DTH; y = mag(U,V);
z = SWHa./nanmean(SWH,[1 2]); z = remove_planefit(X,Y,z);
% z = SWHa; z = remove_planefit(X,Y,z);

A = {x,y,z}; 
A = cellfun(@(Z) cell2matn(arrayfun(@(i)    gaussfilt_2d(Z(:,:,i),eddygrid.drkm, Lavg), 1:size(Z,3), 'Un', 0),3), A, 'Un', 0);
A = cellfun(@(Z) Z.*msk, A, 'Un', 0);
A = cellfun(@(Z) Z(:), A, 'Un', 0);
[x y z] = A{:};
nn = find(~isnan(x) & ~isnan(y) & ~isnan(z));
A = cellfun(@(Z) Z(nn), A, 'Un', 0);
[x y z] = A{:};



figure(234); clf;
% setfigsize(gcf, [500 600])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 3, [0.025 0.015], [0.17 0.05], [0.05 0.19]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
ha(end).Position(1) = ha(end).Position(1) + 0.14;
% -------------------------------------------------------------------------
format_fig4print(ha, 'Fontsize',12);
setaxes(ha(1:2), 'XTick', [0:30:180]); 
setaxes(ha(1:2), 'YTick', [0:0.5:3]); 
setaxes(ha, 'CLim', [-1 1].*0.3); 
setaxes(ha, 'Colormap', niceredblue); 
setaxes(ha, 'YLim', [0 2.4]); 
setaxes(ha, 'XLim', [0 180]); 
setaxes(ha, 'XLabel', '$|\Delta \theta|$'); 
setaxes(ha, 'YLabel', '$|\vec{u}_{geo}|$ [m/s]'); 

% setaxes(ha(end), 'DataAspectRatio', [1 1 1]); 


% linkaxes(ha)
quietaxes(ha(2:end-1), 'y')
format_mapticks(ha(1:2), 'x');
% -------------------------------------------------------------------------
xbins = [0:2.5:180]; ybins = [0:0.05:3];
[M, N, xmid, ymid] = bin2D( xbins, ybins, x, y, z , 'stat', 'median');
[Xmid Ymid] = ndgrid(xmid,ymid); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
Z = M; 
Z(N<=nmin) = NaN;
pcolor(Xmid,Ymid,Z); shading flat;
Z = M; Z(N<=nmin) = NaN;
Z = movmean(Z,3,1, 'omitnan'); Z = movmean(Z,3,2, 'omitnan'); 
Z = movmean(Z,5,1, 'omitnan'); Z = movmean(Z,5,2, 'omitnan'); 
Z(N<=nmin) = NaN;
contour(Xmid,Ymid,Z, [-1:0.025:1], 'k-','EdgeAlpha', 0.3, 'LineWidth',1)

% %%% PHILLIPS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
% 
% contour(Xmid,Ymid,Z, [-1:0.025:1], '-','LineWidth',1)

%%% phillips
omp = 2*pi*(1/11); 
kp = findk_u(omp, 4000,0,0);
cp = omp./kp;


% umagbins = [0:0.1:3];
% thbins = [0:180];
% [umag th] = ndgrid(umagbins, thbins);
[Ubins] = pol2cart(deg2rad(Xmid), Ymid);

cp_new = [0.5+0.5.*sqrt(1+4.*Ubins./cp)].*cp;

var = cp./sqrt(cp_new.*(cp_new + 2.*Ubins));
var(isnan(Z)) = NaN;
% contourf(Xmid,Ymid,var-1, [-1:0.05:1].*0.5, 'LineWidth',2)
contour(Xmid,Ymid,var-1, [-1:0.05:1].*0.5, 'LineWidth',2)

posx = thisha.Position(1) + 0.01; posy = sum(thisha.Position([2 4])) - 0.02;
textbypos(posx, posy, '$\frac{a-a_0}{a_0}$', 'Interpreter', 'latex')


% a/a0 = c0/sqrt(c*(c + 2U))
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c = fixedcolorbar(gca);
c.Position(1) = sum(thisha.Position([1 3])) + 0.013;
ylabel(c, 'SWHa/$\langle$SWH$\rangle$', 'Interpreter', 'latex')




% -------------------------------------------------------------------------
prcs = [5 95];
thwid = 35; rwid = [25 50];
% thwid = 45; rwid = [25 50]; thoff = 25; 
thwid = 45; rwid = [25 50]; thoff = 5; 
% -------------------------------------------------------------------------
wspd = squeeze(nanmean(WSP_model, [1 2])); 
swh = squeeze(nanmean(SWH_model, [1 2])); 
dth = squeeze(dirnanmean((TH_wv_outer - TH_wnd_outer), [1 2])); %dth = abs(dth);
wdw = dirnaninterp1(era5.time, rad2deg(squeeze(nanmean(era5.wdw,[1 2]))), t);
mwp = naninterp1(era5.time, (squeeze(nanmean(era5.mwp,[1 2]))), t);
cp = g./(2*pi./mwp);
wag = cp./wspd;
vwag = cp./(wspd.*cosd(dth));
% -------------------------------------------------------------------------
rotfldstr = '_rot_wv';
run load_vars_into_memory
VARs = {U_wv, SWH,SWH./nanmean(SWH_model.*cmsk, [1 2])};
% rotfldstr = '_rot_wnd';
run load_vars_into_memory       
VARs = [VARs {U_wnd, WSP,WSP./nanmean(WSP_model.*cmsk, [1 2])}];
% rotfldstr = '';
% run load_vars_into_memory       
% VARs = [VARs {WSP,SWH}];


omp = 2*pi*(1./mwp); omp = permute(omp, [2 3 1]);
kp = findk_u(omp, 4000,0,0);
cp = omp./kp;
cp_new = [0.5+0.5.*sqrt(1+4.*U_wv./cp)].*cp;
x = cp./sqrt(cp_new.*cp_new + 2.*U_wv);
x = x - 1;
VARs = [VARs, {x}];

VAR_anoms = cell(size(VARs));
tmp = NaN(2, length(t)); 
VAR_anoms(:) = {tmp}
for vi=1:length(VARs)
    VAR = VARs{vi};
    msk = ones(size(X)); msk(~isin(mag(X,Y), [75 200])) = NaN;
    VAR = VAR - nanmedian(VAR.*msk, [1 2]); 
    % VAR = remove_planefit(X,Y,VAR, 'smooth', 50, 'drkm', eddygrid.drkm);
    
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
    VAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;
    
    
    msk = ones(size(VAR)); msk(isnan(VAR)) = NaN; 
    VAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(VAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
    VAR = cell2matn(arrayfun(@(i)    filt2(VAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(aVAR,3), 'Un', 0),3);
    VAR = VAR.*msk;

    % no out of time bounds
    tt = find(~isin(t, daterange)); 
    VAR(:,:,tt) = NaN;

    VARs{vi} = VAR;

    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = 90 + [-1 1].*thwid - thoff;
    % thlims = 90 + [-1 1].*thwid - 15;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';
    
    VAR_anoms{vi}(1,:) = y;

    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = -90 + [-1 1].*thwid + thoff;
    % thlims = -90 + [-1 1].*thwid + 15;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';

    VAR_anoms{vi}(2,:) = y;
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(end); axes(thisha); cla; hold on; 
xlim([-2 2]); ylim([-1 1].*0.3); 
xlabel('$\frac{a-a_0}{a_0}$'); 
xlabel('$\hat{\theta}_{wv}\cdot\vec{u}_{geo}$ [m/s]'); 
ylabel('SWHa/$\langle$SWH$\rangle$')


x = [VAR_anoms{1}]';
y = [VAR_anoms{3}]'; 
z = [VAR_anoms{5}]'; 
% z = [mwp, mwp] + 2;
x = x(:);y = y(:); z = z(:);
scatter(x,y,20, 'k', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
% scatter(x,y,50, z, 'filled', 'MarkerEdgeColor', 'k')
% caxis([-2 2]); colormap(gca, flipud(cmocean('speed')))
% caxis([6 12]); colormap(gca, cmocean('gray'))
% 
% msk = ones(size(X)); 
% msk(~[isin(mag(X,Y), rwid) & isin(atan2d(Y,X), 90 + [-1 1].*thwid)]) = NaN;
% x = nanmean(VARs{1}.*msk,3); y = nanmean(VARs{3}.*msk,3);
% x = x(:);y = y(:);
% scatter(x,y,5, 'k', 'filled')


msk = ones(size(X)); 
msk(~isin(mag(X,Y), [0 100])) = NaN;
x = nanmean(VARs{1}.*msk,3); y = nanmean(VARs{3}.*msk,3); z = atan2d(Y,X);
x = nanmean(VARs{1}.*msk,3); y = nanmean(VARs{3}.*msk,3); z = X;
% x = VARs{1}.*msk; y = VARs{3}.*msk; z = repmat(X, [1 1 size(y,3)]); 
x = x(:);y = y(:); z = z(:);


% cmap = [buildcmap(niceredblue([end 1],:), 80)]; cmap = [cmap(end/2:end,:); niceredblue; cmap(1:end/2,:); ]; 
% colormap(gca, cmap); caxis([-180 180]);
% z = abs(z); 
% cmap = subsetcmap(nicejet, 12); colormap(gca, cmap); caxis([0 180]);
% cmap = subsetcmap(cmocean('thermo'), 12); colormap(gca, cmap); caxis([0 180]);
% ii = find(isin(z, [90 + [-1 1].*15]));
ii = find(isin(z, [-1 1].*10));
scatter(x(ii),y(ii),5, 'r', 'filled', 'MarkerFaceAlpha',0.3);  c = fixedcolorbar(gca);
scatter(x(ii),y(ii),5, 'r', 'filled', 'MarkerFaceAlpha',0.3, 'DisplayName', '$-10$km$< x_{wv} < 10$km');  c = fixedcolorbar(gca);
% plotBestOptFit(gca, x(ii), y(ii))
% ii = find(isin(z, [15 + [-1 1].*15]));
ii = find(isin(z, [-1 1].*10 + 50));
scatter(x(ii),y(ii),5, 'b', 'filled', 'MarkerFaceAlpha',0.3, 'DisplayName', '$40$km$< x_{wv} < 60$km');  c = fixedcolorbar(gca);
% plotBestOptFit(gca, x(ii), y(ii))
% ii = find(isin(z, [-1 1].*10 - 50));
% scatter(x(ii),y(ii),8, 'g', 'filled');  c = fixedcolorbar(gca);


ubins = [-2:0.1:2];
tp = [5 10 15]';
tp = 10;
omp = 2*pi*(1./tp); 
kp = findk_u(omp, 4000,0,0);
cp = omp./kp;
cp_new = [0.5+0.5.*sqrt(1+4.*ubins./cp)].*cp;
var = cp./sqrt(cp_new.*(cp_new + 2.*ubins)); var = var - 1;

cmap = scaledcmap(cmocean('gray'), tp, [4 16]); 
for tpi=1:length(tp)
    plot(ubins,var(tpi,:),'k-', 'Linewidth',2,'DisplayName','$\frac{a-a_0}{a_0}$', 'Color', cmap(tpi,:))
end
cleanLegend(gca)

savejpg(gcf, 'u_vs_th_vs_swh_plus_scatter_withphillips', savepath, saveopt)


%% [save] SWH 
Rcc = 100;

% composite + phillips

cgopt = 0;
nmin = 5;
Lavg = 6; 
Rmax = 200;
prcnan = 0.25;
cclim = 0.25;
cclim = 0.0;
% cclim = 0.5;
% cclim = 0.625;
thwid = 45; rwid = [25 50]; thoff = 5; 
thwid = 45; rwid = [25 50]; thoff = 15; 
thwid = 45; rwid = [25 75]; thoff = 5; 
thwid = 45; rwid = [30 60]; thoff = 5; 
% thwid = 30; rwid = [25 50]; thoff = 5; 

% thwid = 30; rwid = [25 50]; thoff = 5; 

daterange = [datenum(2023,4,13) datenum(2023,5,12)];
% daterange = [datenum(2023,4,13) datenum(2023,4,17)]; prcnan = 0.05;
rotfldstr = '_rot_wv';
run load_vars_into_memory_universal
% -------------------------------------------------------------------------
wspd = squeeze(nanmean(WSP_model, [1 2])); 
swh = squeeze(nanmean(SWH_model, [1 2])); 
dth = squeeze(dirnanmean((TH_wv_outer - TH_wnd_outer), [1 2])); %dth = abs(dth);
wdw = dirnaninterp1(era5.time, rad2deg(squeeze(nanmean(era5.wdw,[1 2]))), t);
mwp = naninterp1(era5.time, (squeeze(nanmean(era5.mwp,[1 2]))), t);


VAR = SSTa; %VAR = movmean(VAR,3,3,'omitnan'); 
cmsk = ones(size(X)); cmsk(mag(X,Y) > 20) = NaN; emsk = ones(size(X)); emsk(~isin(mag(X,Y), [120 200])) = NaN; 
dsst = squeeze(nanmedian(VAR.*cmsk, [1 2]) - nanmedian(VAR.*emsk, [1 2])); 
n = [squeeze(sum(~isnan(SSTa.*emsk),[1 2])./sum(~isnan(emsk),[1 2])) squeeze(sum(~isnan(SSTa.*cmsk),[1 2])./sum(~isnan(cmsk),[1 2]))];
n = min(n,[],2);
dsst(n<0.05) = NaN;

cp = g./(2*pi./mwp);
wag = cp./wspd;
vwag = cp./(wspd.*cosd(dth));
% -------------------------------------------------------------------------
% VARy = SWH./nanmean(SWH_model.*cmsk, [1 2]);
msk = ones(size(X)); msk(mag(X,Y)>Rmax) = NaN;
VARy = SWH./nanmean(SWH.*msk, [1 2]);
if cgopt; VARxs = U_wv_cg; else; VARxs = U_wv; end


tp = 10;
% tp = mwp + 4; tp = permute(tp, [2 3 1]);
omp = 2*pi*(1./tp); 
kp = findk_u(omp, 4000,0,0);
kp = omp^2/g;
cp = omp./kp;
% cp = omp./(omp^2/g) = g/omp = g/(2*pi*1/tp) = tp,*g./(2*pi);
cp = tp.*g./(2*pi);

cp_new = [0.5+0.5.*sqrt(1+4.*U_wv./cp)].*cp;
var = cp./sqrt(cp_new.*(cp_new + 2.*U_wv)); 
var = var - 1;

% omp = 2*pi*(1./mwp); omp = permute(omp, [2 3 1]);
% omp = 2*pi*(1./(mwp-1)); omp = permute(omp, [2 3 1]);
% omp = 2*pi*(1./(10)); omp = permute(omp, [2 3 1]);
% kp = findk_u(omp, 4000,0,0);
% cp = omp./kp;
% cp_new = [0.5+0.5.*sqrt(1+4.*U_wv./cp)].*cp;
% var = cp./sqrt(cp_new.*cp_new + 2.*U_wv);
% var = var - 1;
VARxs = [VARxs, {var}];
VARs = [VARy VARxs];
tmp = NaN(2, length(t)); 
VAR_anoms(:) = {tmp};
Nt = sum(isin(t, daterange));
for vi=1:length(VARs)
    VAR = VARs{vi};
    msk = ones(size(X)); msk(~isin(mag(X,Y), [75 200])) = NaN;
    VAR = VAR - nanmedian(VAR.*msk, [1 2]); % remove pass mean
    VAR = remove_planefit(X,Y,VAR, 'smooth', 50, 'drkm', eddygrid.drkm);
    
    % make sure even coverage on top and bottom
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
    VAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;
    
    % smoothing
    msk = ones(size(VAR)); msk(isnan(VAR)) = NaN; 
    VAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(VAR(:,:,i),eddygrid.drkm, Lavg), 1:size(VAR,3), 'Un', 0),3);
    VAR = cell2matn(arrayfun(@(i)    filt2(VAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(VAR,3), 'Un', 0),3);
    VAR = VAR.*msk;
    
    % no out of time bounds
    tt = find(~[isin(t, daterange)]); 
    % tt = find(~[isin(t, daterange) & abs(dth)' < 35]); 
    % tt = find(~[isin(t, daterange) & mwp'<8.5 & wdw' < 35]); Nt = sum(isin(t, daterange)) - length(tt);
    % tt = find(~[isin(t, daterange) & mwp'>8.5 & wdw' > 35]); Nt = sum(isin(t, daterange)) - length(tt);
    % tt = find(~[isin(t, daterange) & mwp'> 8.25 ]); Nt = sum(isin(t, daterange)) - length(tt);
    % tt = find(~[isin(t, daterange) & mwp'< 8.25 ]); Nt = sum(isin(t, daterange)) - length(tt);
    tt = find(~[isin(t, daterange) & mwp' > 8 ]); Nt = sum(isin(t, daterange)) - length(tt);
    % 
    VAR(:,:,tt) = NaN;

    VARs{vi} = VAR;


    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = 90 + [-1 1].*thwid - thoff;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    % lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';
    
    VAR_anoms{vi}(1,:) = y;

    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = -90 + [-1 1].*thwid + thoff;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    % lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';

    VAR_anoms{vi}(2,:) = y;


    % no out of time bounds
    tt = find(~isin(t, daterange)); 
    % VAR_anoms{vi}(:,tt) = NaN;

end
VARy = VARs{1};
VARxs = VARs(2:end);
VARy_anom = VAR_anoms{1};
VARxs_anom = VAR_anoms(2:end);


%%% choice
VARx = VARxs{1};
msk = ones(size(VARy)); msk(isnan(VARx + VARy)) = NaN;
VARx = VARx.*msk; VARy = VARy.*msk;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

figure(21355); clf;
setfigsize(gcf, [1148         303])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 4, [0.05 0.015], [0.23 0.05], [0.06 0.15]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
dx = 0.05; 
hi = 3; 
ha(hi).Position(1) =ha(hi).Position(1) + dx;
ha(hi).Position(3) =ha(hi).Position(3) + dx/2;
hi = 4;
ha(hi).Position(1) =ha(hi).Position(1) + dx*2.5;

% -------------------------------------------------------------------------
format_fig4print(ha, 'FontSize',13);
setaxes(ha, 'Colormap', niceredblue); 
setaxes(ha([1 2 4]), 'YLabel', '$y_{wv}$ [km]'); 
setaxes(ha([1 2 4]), 'XLabel', '$x_{wv}$ [km]'); 
setaxes(ha([1 2 4]), 'DataAspectRatio', [1 1 1]); 
setaxes(ha([1 2 4]), 'XLim', [-1 1].*130); 
setaxes(ha([1 2 4]), 'YLim', [-1 1].*130); 
setaxes(ha, 'CLim', [-1 1].*0.22); 
setaxes(ha(2), 'CLim', [-1 1].*1.3); 
setaxes(ha(2), 'Colormap', flipud(niceredblue)); 
linkaxes(ha([1 2 4]), 'xy')

setaxes(ha(3), 'YLim', [-1 1].*0.4); 
setaxes(ha(3), 'YLim', [-0.3 0.4]); 
setaxes(ha([3]), 'XLabel', '$\hat{\theta}_{wv} \cdot\vec{u}_{geo}$ [m/s]'); 
setaxes(ha([3]), 'YLabel', 'SWHa/$\langle$SWH$\rangle$'); 

quietaxes(ha(2), 'y')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
Z = VARy; 
% cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
Z(:,:,abs(cc) < cclim) = NaN;
N = sum(~isnan(Z),3); 

Z = nanmean(Z,3); 
Z = gaussfilt_2d(Z,drkm, Lavg); 
Z(N < Nt.*prcnan./2) = NaN;

pcolor(X,Y,Z); shading flat;

clims = get(gca, 'CLim'); 
msk = ones(size(Z)); msk(isnan(Z)) = NaN;
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg.*1.75); Z = Z.*msk; 
% contourf(X,Y,Z,[-1:0.075:1].*max(clims),'Edgealpha',0.1)
% contour(X,Y,Z,[-1:0.075:1].*max(clims),'Edgealpha',0.2, 'EdgeColor', 'k')
contour(X,Y,Z,linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha',0.3)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
Z = VARx; 
msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
Z(:,:,abs(cc) < cclim) = NaN;
N = sum(~isnan(Z),3); 

Z = nanmean(Z,3); 
Z = gaussfilt_2d(Z,drkm, Lavg); 
Z(N < Nt.*prcnan./2) = NaN;

pcolor(X,Y,Z); shading flat;

clims = get(gca, 'CLim'); 
msk = ones(size(Z)); msk(isnan(Z)) = NaN;
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg.*1.75); Z = Z.*msk; 
% contourf(X,Y,Z,[-1:0.075:1].*max(clims),'Edgealpha',0.1)
contour(X,Y,Z,linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha',0.3)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(4); axes(thisha); hold on; 
% Z = nanmean(VARy,3) - nanmean(VARxs{2},3); 
% Z = nanmean(VARy - VARxs{2},3); 


Z = VARy - VARxs{2}; 
msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
Z(:,:,abs(cc) < cclim) = NaN;
N = sum(~isnan(Z),3); 

Z = nanmean(Z,3); 
Z = gaussfilt_2d(Z,drkm, Lavg); 
Z(N < Nt.*prcnan./2) = NaN;


% N = sum(~isnan(VARy + VARx),3); 
% Z(N < Nt.*prcnan./2) = NaN;
% cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
% Z(abs(cc) < 0.5) = NaN;

pcolor(X,Y,Z); shading flat; 


clims = get(gca, 'CLim'); 
msk = ones(size(Z)); msk(isnan(Z)) = NaN;
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg.*1.75); Z = Z.*msk; 
% contourf(X,Y,Z,[-1:0.075:1].*max(clims),'Edgealpha',0.1)
contour(X,Y,Z,linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha',0.3)



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(3); axes(thisha); cla; hold on; 
x = VARx;
y = VARy; 

msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3)); cc = abs(cc);

% cc = arrayfun(@(i) ccfxn(x(:,:,i), y(:,:,i)), 1:size(x,3)); cc = abs(cc);
msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% msk(~isin(X, 5.*[-1 1]+0)) = NaN;
x(:,:,cc<cclim) = NaN; y(:,:,cc<cclim) = NaN;
nn = find(~isnan(x.*msk + y.*msk));
scatter(x(nn),y(nn),5,[1 1 1].*0.4, 'filled', 'MarkerFaceAlpha',0.1);

binx = [-2:0.05:2]; biny = [-1:0.01:1];
binx = [-2:0.1:2]; biny = [-1:0.02:1];
binx = [-2:0.15:2]; biny = [-1:0.025:1];
% binx = [-2:0.2:2]; biny = [-1:0.04:1];
n = hist2d(x(nn), y(nn), binx, biny); %n(n==0) = NaN;

% n = n./nansum(n(:)); %n = log10(n); 

pcolor(binx, biny, n'); shading flat; 
% caxis(10.^[-4.5 -1.3]); colormap(gca, flipud(cmocean('gray')))
caxis([0 max(n(:)).*1.1]); colormap(gca, flipud(cmocean('gray')))
% caxis([0 max(n(:)).*1.2]); colormap(gca, (cmocean('algae')))
[bX bY] = ndgrid(binx,biny); 
% contourf(bX,bY,n,10.^[-6:0.25:1])

% contourf(bX,bY,n,[0:0.0025:1], 'EdgeColor', 'none')
% c = fixedcolorbar(gca);


% % 
% % nn = find(~isnan(x.*msk + y.*msk) & mwp'<8.5);
% % scatter(x(nn),y(nn),5,[1 1 1].*0.4, 'filled', 'MarkerFaceAlpha',0.1)


% rbins = [-15:30:80]; 
% ubins = [-2:0.25:2];
% cmap = scaledcmap(nicejet, rbins(1:end-1));
% for ri=1:length(rbins)-1
%     x = nanmean(VARx,3); y = nanmean(VARy,3); 
%     msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
%     msk(~isin(X, rbins(ri:ri+1))) = NaN;
%     nn = find(~isnan(x.*msk + y.*msk));
%     % scatter(x(nn),y(nn),5,'k', 'filled', 'MarkerFaceColor', cmap(ri,:));
% 
%     x = VARx; y = VARy; 
%     msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
%     msk(~isin(X, rbins(ri:ri+1))) = NaN;
%     nn = find(~isnan(x.*msk + y.*msk));
%     [sts] = boxplot_bybin(gca, y(nn), x(nn), ubins, 'Color', cmap(ri,:), 'width', mode(diff(ubins)).*0.4, 'outliars', 0, 'plot', 0);
%     plot(sts.x0, sts.q2, 'k.-', 'Color', cmap(ri,:), 'MarkerSize',10)
% end
% 
% x = nanmean(VARx,3); y = nanmean(VARy,3); 
% msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% nn = find(~isnan(x.*msk + y.*msk));
% scatter(x(nn),y(nn),5,'k', 'filled')
% scatter(x(nn),y(nn),20, 'k', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'DisplayName', 'SWHa/$\langle$ SWH$\rangle$')

% 
varx = VARx; vary = VARy;
msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
cc = repmat(cc, [2 1]);
y = VARy_anom; 
x = VARxs_anom{1}; 
% plot(x,y,'k-')
x = x(:); y = y(:); cc = cc(:); 
scatter(x,y,14, 'x', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'linewidth',1)
nn = find(abs(cc)>cclim); 
scatter(x(nn),y(nn),20, 'k', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'linewidth',1, 'DisplayName', 'SWHa(t)/$\langle$ SWH$\rangle$')





% plot(nanmean(VARx,3), nanmean(VARxs{2},3), 'b.')


ubins = [-2:0.1:2];
tp = [7.5 10 12.5]';
tp = [5 10 15]';
omp = 2*pi*(1./tp); 
kp = findk_u(omp, 4000,0,0);
cp = omp./kp;
cp_new = [0.5+0.5.*sqrt(1+4.*ubins./cp)].*cp;
var = cp./sqrt(cp_new.*(cp_new + 2.*ubins)); 
var = var - 1; var = real(var);
% plot(ubins, var, 'k:', 'LineWidth',2)
shaded_limits(ubins, var(2,:), var([1 3],:), 'k', 'Linestyle', ':', 'EdgeColor', 'none', 'DisplayName', '$\Delta a/a_0$', 'facealpha',0.15);


% 
% ubins = [-2:0.25:2];
% % col = [0.7 0.8 0.3];
% col = [0.26 0.5 0.05];
% x = nanmean(VARx,3); y = nanmean(VARy,3); 
% % x = VARx; y = VARy; 
% msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% msk(~isin(X, 50 + [-1 1].*5)) = NaN;
% % msk(~isin(mag(X,Y), [45 55])) = NaN;
% nn = find(~isnan(x.*msk + y.*msk));
% % scatter(x(nn),y(nn),8,col, '^', 'filled', 'Displayname', 'composite @ $x_{wv}=50$km');
% [sts] = boxplot_bybin(gca, y(nn), x(nn), ubins, 'Color', col, 'width', mode(diff(ubins)).*0.4, 'outliars', 0, 'plot', 0, 'FaceAlpha', 0.2, 'nmin',3);
% % errorbar(sts.x0, sts.q2,sts.q2 - sts.q1, sts.q3 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% errorbar(sts.x0, sts.q2,sts.q2 - sts.p10, sts.p90 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% 
% 
% col = [0.26 0.5 0.05]./1.5;
% x = nanmean(VARx,3); y = nanmean(VARy,3); 
% % x = VARx; y = VARy; 
% msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% msk(~isin(X, [-1 1].*5)) = NaN;
% % msk(~isin(mag(X,Y), [-5 5])) = NaN;
% nn = find(~isnan(x.*msk + y.*msk));
% % scatter(x(nn),y(nn),5,'k', 'filled');
% % scatter(x(nn),y(nn),8 ,col,'^', 'filled', 'Displayname', 'composite @ $x_{wv}=0$km');
% [sts] = boxplot_bybin(gca, y(nn), x(nn), ubins, 'Color', col, 'width', mode(diff(ubins)).*0.4, 'outliars', 0, 'plot', 0, 'FaceAlpha', 0.2, 'nmin',3);
% % errorbar(sts.x0, sts.q2,sts.q2 - sts.q1, sts.q3 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% errorbar(sts.x0, sts.q2,sts.q2 - sts.p10, sts.p90 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% 



% cleanLegend(gca, 'northeast', 'FontSize',8, 'Color', 'none', 'EdgeColor', 'none', 'NumColumns',2)
leg = cleanLegend(gca, 'southwest', 'FontSize',8, 'Color', 'none', 'EdgeColor', 'none', 'NumColumns',2);
leg.ItemTokenSize = leg.ItemTokenSize./2;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x0 = -1.5; y0 = 0.2; 
x0 = thisha.Position(1) + thisha.Position(3).*0.05; y0 = sum(thisha.Position([2 4])) - thisha.Position(4).*0.05; 
varx = VARx; vary = VARy;  
varx(isnan(varx + vary)) = NaN; vary(isnan(varx + vary)) = NaN;
msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
varx = varx.*msk; vary = vary.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3)); cc = abs(cc);
% varx(:,:,abs(cc) < cclim) = NaN; vary(:,:,abs(cc) < cclim) = NaN;
% cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
cc = cc.^2;
cclims = prctile(cc, [25 75]); 
N = sum(~isnan(varx + vary),3); 
msk = ones(size(X)); msk(N < Nt.*prcnan./2) = NaN;
varx = nanmedian(varx,3).*msk; vary = nanmedian(vary,3).*msk; 
cc = ccfxn(varx,vary);
cc = cc.^2;
textbypos(x0,y0,['$r^2$=' num2str(round(cc,2)) ' (' num2str(round(cclims(1),2)) '-'  num2str(round(cclims(2),2)) ')'], 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize',8)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% for hi=1:2; 
%     axes(ha(hi)); hold on; 
%     xline([0 50], 'k--')
% end


for hi=[1 2 4];
    thisha = ha(hi); axes(thisha); hold on; 
    % xline([0 50], 'k--')
    thbins = [0:5:360];
    for r=[50 100 150]
        [xb yb] = pol2cart(deg2rad(thbins), ones(size(thbins)).*r);
        plot(xb,yb,'k:', 'Color', [0 0 0 0.6]);
    end
    % [xb yb] = pol2cart(deg2rad(thbins), ones(size(thbins)).*25);
    % plot(xb,yb,'k:', 'Color', [0 0 0 0.6]);
    % [xb yb] = pol2cart(deg2rad(thbins), ones(size(thbins)).*50);
    % plot(xb,yb,'k:', 'Color', [0 0 0 0.6]);

    c = fixedcolorbar(gca, 'Location','northoutside', 'AxisLocation', 'in', 'FontSize',8);
    c.Position([3 4]) = c.Position([3 4])./[5 2];
    c.Position([1 2]) = thisha.Position([1 2]) + thisha.Position([3 4]).*[0 1];
    c.Position([1 2]) = c.Position([1 2]) - c.Position([3 4]).*[0 1] + c.Position([3 4]).*0.2.*[1 -5];
    c.TickLength = c.TickLength.*6;
    c.Ruler.TickLabelRotation = 0;
    % c.AxisLocation = 'in';
    if hi==1; lab = get(ha(3),'YLabel'); lab = lab.String;
    elseif hi==2; lab = get(ha(3),'XLabel'); lab = lab.String;
    elseif hi==4; lab = 'residual';
    end
    textbypos(sum(c.Position([1 3])) + 0.005, c.Position(2) + c.Position(4)/2, lab, 'Interpreter', 'latex', 'FontSize',get(c, 'FontSize')+2, 'VerticalAlignment', 'middle')
end

for hi=[3]
    thisha = ha(hi); axes(thisha); hold on; 

    c = fixedcolorbar(gca, 'Location','northoutside', 'AxisLocation', 'in', 'FontSize',8);
    c.Position([3 4]) = c.Position([3 4])./[8 2];
    c.Position([1 2]) = thisha.Position([1 2]) + thisha.Position([3 4]).*[1 1];
    
    c.Position([1 2]) = c.Position([1 2]) - c.Position([3 4]).*[1 1] + c.Position([3 4]).*0.2.*[-2 -5];
    c.TickLength = c.TickLength.*6;
    c.Ticks = [0 500]; c.Ruler.TickLabelRotation = 0;
    % c.AxisLocation = 'in';
    if hi==3; lab = 'n';
    end
    textbypos(sum(c.Position([1])) - 0.005, c.Position(2) + c.Position(4)/2, lab,...
        'Interpreter', 'latex', 'FontSize',get(c, 'FontSize')+2, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
end


optstr = '';
if cgopt; optstr = [optstr '_cyclogeostrophy']; end
savejpg(gcf, ['fig_swh_explainer_phillips' optstr], savepath, saveopt)






%% [save] WIND composite 
Rcc = 100;

cgopt = 0;
nmin = 5;
Lavg = 8; 
Rmax = 200;
prcnan = 0.25;
% cclim = 0.75; prcnan = 0.2;
cclim = 0.25;
cclim = 0.25;
% cclim = 0.5;
% cclim = 0.625;
thwid = 45; rwid = [25 50]; thoff = 5; 
thwid = 45; rwid = [25 50]; thoff = 15; 
thwid = 45; rwid = [25 75]; thoff = 5; 
% thwid = 45; rwid = [30 60]; thoff = 5; 
thwid = 45; rwid = [30 60]; thoff = 0; 
% thwid = 30; rwid = [25 50]; thoff = 5; 

% thwid = 30; rwid = [25 50]; thoff = 5; 

daterange = [datenum(2023,4,13) datenum(2023,5,12)];
daterange = [datenum(2023,4,15) datenum(2023,5,12)];
% daterange = [datenum(2023,4,14) datenum(2023,4,26)];
% daterange = [datenum(2023,4,25) datenum(2023,5,3)];
% daterange = [datenum(2023,4,12) datenum(2023,4,18)];
% % daterange = [datenum(2023,4,17) datenum(2023,4,21)];
% daterange = [datenum(2023,4,14) datenum(2023,4,17)];
% daterange = [datenum(2023,4,13) datenum(2023,4,17)]; prcnan = 0.05;
rotfldstr = '_rot_wnd';
run load_vars_into_memory
% -------------------------------------------------------------------------
wspd = squeeze(nanmean(WSP_model, [1 2])); 
swh = squeeze(nanmean(SWH_model, [1 2])); 
dth = squeeze(dirnanmean((TH_wv_outer - TH_wnd_outer), [1 2])); %dth = abs(dth);
wdw = dirnaninterp1(era5.time, rad2deg(squeeze(nanmean(era5.wdw,[1 2]))), t);
mwp = naninterp1(era5.time, (squeeze(nanmean(era5.mwp,[1 2]))), t);
cp = g./(2*pi./mwp);
wag = cp./wspd;
vwag = cp./(wspd.*cosd(dth));
% -------------------------------------------------------------------------
% VARy = SWH./nanmean(SWH_model.*cmsk, [1 2]);
% VARy = SWH./nanmean(SWH.*cmsk, [1 2]);
% VARxs = U_wv;

VARy = WSP;
if cgopt; VARxs = U_wnd_cg; else VARxs = U_wnd; end

clear i; 
var = [nanmean(WSP_U_model.*cntrmsk, [1 2]) + i.*nanmean(WSP_V_model.*cntrmsk, [1 2])] - [U + i.*V]; var = abs(var);

if cgopt; 
    var = [WSP_U_model + i.*WSP_V_model] - [U_cg + i.*V_cg]; var = abs(var);
else 
    var = [WSP_U_model + i.*WSP_V_model] - [U + i.*V]; var = abs(var);
end


VARxs = [VARxs, {var}, SSTa];
VARs = [VARy VARxs];
tmp = NaN(2, length(t)); 
VAR_anoms(:) = {tmp};
Nt = sum(isin(t, daterange));
for vi=1:length(VARs)
    VAR = VARs{vi};
    msk = ones(size(X)); msk(~isin(mag(X,Y), [75 200])) = NaN;
    VAR = VAR - nanmedian(VAR.*msk, [1 2]); % remove pass mean
    VAR = remove_planefit(X,Y,VAR, 'smooth', 50, 'drkm', eddygrid.drkm);
    
    % make sure even coverage on top and bottom
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
    msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
    n = sum(~isnan(msk.*[VAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
    VAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;
    
    % smoothing
    msk = ones(size(VAR)); msk(isnan(VAR)) = NaN; 
    VAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(VAR(:,:,i),eddygrid.drkm, Lavg), 1:size(VAR,3), 'Un', 0),3);
    VAR = cell2matn(arrayfun(@(i)    filt2(VAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(VAR,3), 'Un', 0),3);
    VAR = VAR.*msk;
    
    % no out of time bounds
    tt = find(~[isin(t, daterange)]); 
    % tt = find(~[isin(t, daterange) & abs(dth)' < 35]); 
    % tt = find(~[isin(t, daterange) & mwp'<8.5 & wdw' < 35]); Nt = sum(isin(t, daterange)) - length(tt);
    % tt = find(~[isin(t, daterange) & mwp'>8.5 & wdw' > 35]); Nt = sum(isin(t, daterange)) - length(tt);
    % tt = find(~[isin(t, daterange) & mwp'> 8.25 ]); Nt = sum(isin(t, daterange)) - length(tt);
    % tt = find(~[isin(t, daterange) & mwp'< 8.25 ]); Nt = sum(isin(t, daterange)) - length(tt);
    % tt = find(~[isin(t, daterange) & mwp'> 8.5 ]); Nt = sum(isin(t, daterange)) - length(tt);
    % 
    VAR(:,:,tt) = NaN;

    VARs{vi} = VAR;


    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = 90 + [-1 1].*thwid - thoff;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    % lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';
    
    VAR_anoms{vi}(1,:) = y;

    msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
    thlims = -90 + [-1 1].*thwid + thoff;
    msk(~isin(atan2d(Y,X), thlims)) = NaN;
    y = VAR.*msk;
    y = squeeze(nanmedian(y,[1 2])); 
    % lims = prctile(y, prcs, [1 2]); lims = squeeze(lims)';

    VAR_anoms{vi}(2,:) = y;


    % no out of time bounds
    tt = find(~isin(t, daterange)); 
    % VAR_anoms{vi}(:,tt) = NaN;

end
VARy = VARs{1};
VARxs = VARs(2:end);
VARy_anom = VAR_anoms{1};
VARxs_anom = VAR_anoms(2:end);


%%% choice
VARx = VARxs{2}; 
msk = ones(size(VARy)); msk(isnan(VARx + VARy)) = NaN;
VARx = VARx.*msk; VARy = VARy.*msk;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

figure(21356); clf;
setfigsize(gcf, [1148         293])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 4, [0.05 0.015], [0.23 0.05], [0.06 0.15]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
dx = 0.05; 
hi = 3; 
ha(hi).Position(1) =ha(hi).Position(1) + dx;
ha(hi).Position(3) =ha(hi).Position(3) + dx/2;
hi = 4;
ha(hi).Position(1) =ha(hi).Position(1) + dx*2.5;

% -------------------------------------------------------------------------
format_fig4print(ha, 'FontSize',13);
setaxes(ha, 'Colormap', niceredblue); 
setaxes(ha([1 2 4]), 'YLabel', '$y_{wnd}$ [km]'); 
setaxes(ha([1 2 4]), 'XLabel', '$x_{wnd}$ [km]'); 
setaxes(ha([1 2 4]), 'DataAspectRatio', [1 1 1]); 
setaxes(ha([1 2 4]), 'XLim', [-1 1].*130); 
setaxes(ha([1 2 4]), 'YLim', [-1 1].*130); 
setaxes(ha, 'CLim', [-1 1].*1.5); 
% setaxes(ha(2), 'CLim', [-1 1].*1.5); 
linkaxes(ha([1 2 4]), 'xy')

setaxes(ha(3), 'XLim', [-1 1].*2.8); 
setaxes(ha(3), 'XLim', [-1 1].*2); 
setaxes(ha(3), 'YLim', [-1 1].*3); 
% setaxes(ha([3]), 'XLabel', '$\hat{\theta}_{wv} \cdot\vec{u}_{geo}$ [m/s]'); 
% setaxes(ha([3]), 'YLabel', 'SWHa/$\langle$SWH$\rangle$'); 
setaxes(ha([3]), 'YLabel', 'WSPa [m/s]'); 
setaxes(ha([3]), 'XLabel', '$|\vec{U}_{10} - \vec{u}_{geo}|_a$ [m/s]'); 
% setaxes(ha([3]), 'XLabel', '$\hat{\theta}_{wnd} \cdot\vec{u}_{geo}$ [m/s]'); 
% setaxes(ha(2), 'Colormap', flipud(niceredblue)); 

quietaxes(ha(2), 'y')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
Z = VARy; 
% cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
Z(:,:,abs(cc) < cclim) = NaN;
N = sum(~isnan(Z),3); 

Z = nanmedian(Z,3); 
Z = gaussfilt_2d(Z,drkm, Lavg); 
Z(N < Nt.*prcnan./2) = NaN;

pcolor(X,Y,Z); shading flat;

clims = get(gca, 'CLim'); 
msk = ones(size(Z)); msk(isnan(Z)) = NaN;
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg.*1.75); Z = Z.*msk; 
% contourf(X,Y,Z,[-1:0.075:1].*max(clims),'Edgealpha',0.1)
contour(X,Y,Z,linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha',0.3)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
Z = VARx; 
msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
Z(:,:,abs(cc) < cclim) = NaN;
N = sum(~isnan(Z),3); 

Z = nanmedian(Z,3); 
Z = gaussfilt_2d(Z,drkm, Lavg); 
Z(N < Nt.*prcnan./2) = NaN;

pcolor(X,Y,Z); shading flat;

clims = get(gca, 'CLim'); 
msk = ones(size(Z)); msk(isnan(Z)) = NaN;
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg.*1.75); Z = Z.*msk; 
% contourf(X,Y,Z,[-1:0.075:1].*max(clims),'Edgealpha',0.1)
contour(X,Y,Z,linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha',0.3)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(4); axes(thisha); hold on; 


Z = VARy - VARxs{2}; 
msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
Z(:,:,abs(cc) < cclim) = NaN;
N = sum(~isnan(Z),3); 

Z = nanmedian(Z,3); 
Z = gaussfilt_2d(Z,drkm, Lavg); 
Z(N < Nt.*prcnan./2) = NaN;


% N = sum(~isnan(VARy + VARx),3); 
% Z(N < Nt.*prcnan./2) = NaN;
% cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
% Z(abs(cc) < 0.5) = NaN;

pcolor(X,Y,Z); shading flat; 


clims = get(gca, 'CLim'); 
msk = ones(size(Z)); msk(isnan(Z)) = NaN;
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg.*1.75); Z = Z.*msk; 
% contourf(X,Y,Z,[-1:0.075:1].*max(clims),'Edgealpha',0.1)
contour(X,Y,Z,linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha',0.3)



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(3); axes(thisha); cla; hold on; 
x = VARx;
y = VARy; 

msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3)); cc = abs(cc);

% cc = arrayfun(@(i) ccfxn(x(:,:,i), y(:,:,i)), 1:size(x,3)); cc = abs(cc);
msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% msk(~isin(X, 5.*[-1 1]+0)) = NaN;
x(:,:,cc<cclim) = NaN; y(:,:,cc<cclim) = NaN;
nn = find(~isnan(x.*msk + y.*msk));
scatter(x(nn),y(nn),5,[1 1 1].*0.4, 'filled', 'MarkerFaceAlpha',0.1);


binx = [-3:0.18:3]; biny = binx;
% binx = [-2:0.2:2]; biny = [-1:0.04:1];
n = hist2d(x(nn), y(nn), binx, biny); %n(n==0) = NaN;
% n = n./nansum(n(:)); %n = log10(n); 
pcolor(binx, biny, n'); shading flat; 
% caxis(10.^[-4.5 -1.3]); colormap(gca, flipud(cmocean('gray')))
caxis([0 max(n(:)).*1.1]); colormap(gca, flipud(cmocean('gray')))
% caxis([0 max(n(:)).*1.2]); colormap(gca, (cmocean('algae')))
[bX bY] = ndgrid(binx,biny); 
% contourf(bX,bY,n,10.^[-6:0.25:1])

% contourf(bX,bY,n,[0:0.0025:1], 'EdgeColor', 'none')
% c = fixedcolorbar(gca);


% % 
% % nn = find(~isnan(x.*msk + y.*msk) & mwp'<8.5);
% % scatter(x(nn),y(nn),5,[1 1 1].*0.4, 'filled', 'MarkerFaceAlpha',0.1)


% rbins = [-15:30:80]; 
% ubins = [-2:0.25:2];
% cmap = scaledcmap(nicejet, rbins(1:end-1));
% for ri=1:length(rbins)-1
%     x = nanmean(VARx,3); y = nanmean(VARy,3); 
%     msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
%     msk(~isin(X, rbins(ri:ri+1))) = NaN;
%     nn = find(~isnan(x.*msk + y.*msk));
%     % scatter(x(nn),y(nn),5,'k', 'filled', 'MarkerFaceColor', cmap(ri,:));
% 
%     x = VARx; y = VARy; 
%     msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
%     msk(~isin(X, rbins(ri:ri+1))) = NaN;
%     nn = find(~isnan(x.*msk + y.*msk));
%     [sts] = boxplot_bybin(gca, y(nn), x(nn), ubins, 'Color', cmap(ri,:), 'width', mode(diff(ubins)).*0.4, 'outliars', 0, 'plot', 0);
%     plot(sts.x0, sts.q2, 'k.-', 'Color', cmap(ri,:), 'MarkerSize',10)
% end
% 
% x = nanmean(VARx,3); y = nanmean(VARy,3); 
% msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% nn = find(~isnan(x.*msk + y.*msk));
% scatter(x(nn),y(nn),5,'k', 'filled')
% scatter(x(nn),y(nn),20, 'k', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'DisplayName', 'SWHa/$\langle$ SWH$\rangle$')

% 
varx = VARx; vary = VARy;
msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
varx = VARx.*msk; vary = VARy.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
cc = repmat(cc, [2 1]);
y = VARy_anom; 
x = VARxs_anom{2}; 
% plot(x,y,'k-')
x = x(:); y = y(:); cc = cc(:); 
scatter(x,y,14, 'x', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'linewidth',1)
nn = find(abs(cc)>cclim); 
scatter(x(nn),y(nn),20, 'k', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'linewidth',1, 'DisplayName', 'WSPa(t)')





% plot(x(nn), y(nn), 'r.')
% [pl, sts] = plotBestOptFit(gca, x(nn), y(nn), 'color', [0.4 0.2 0.2]);
% pl.DisplayName = ['$r^2=$' num2str(round(sts.rsq,2))];


% plot(nanmean(VARx,3), nanmean(VARxs{2},3), 'b.')


ubins = [-3:0.1:3];
% tp = [7.5 10 12.5]';
% tp = [5 10 15]';
% omp = 2*pi*(1./tp); 
% kp = findk_u(omp, 4000,0,0);
% cp = omp./kp;
% cp_new = [0.5+0.5.*sqrt(1+4.*ubins./cp)].*cp;
% var = cp./sqrt(cp_new.*(cp_new + 2.*ubins)); 
% var = var - 1; var = real(var);
plot(ubins, ubins, 'k:', 'LineWidth',2)
% shaded_limits(ubins, var(2,:), var([1 3],:), 'k', 'Linestyle', ':', 'EdgeColor', 'none', 'DisplayName', '$\Delta a/a_0$', 'facealpha',0.15);


% 
% ubins = [-2:0.25:2];
% % col = [0.7 0.8 0.3];
% col = [0.26 0.5 0.05];
% x = nanmean(VARx,3); y = nanmean(VARy,3); 
% % x = VARx; y = VARy; 
% msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% msk(~isin(X, 50 + [-1 1].*5)) = NaN;
% % msk(~isin(mag(X,Y), [45 55])) = NaN;
% nn = find(~isnan(x.*msk + y.*msk));
% % scatter(x(nn),y(nn),8,col, '^', 'filled', 'Displayname', 'composite @ $x_{wv}=50$km');
% [sts] = boxplot_bybin(gca, y(nn), x(nn), ubins, 'Color', col, 'width', mode(diff(ubins)).*0.4, 'outliars', 0, 'plot', 0, 'FaceAlpha', 0.2, 'nmin',3);
% % errorbar(sts.x0, sts.q2,sts.q2 - sts.q1, sts.q3 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% errorbar(sts.x0, sts.q2,sts.q2 - sts.p10, sts.p90 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% 
% 
% col = [0.26 0.5 0.05]./1.5;
% x = nanmean(VARx,3); y = nanmean(VARy,3); 
% % x = VARx; y = VARy; 
% msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
% msk(~isin(X, [-1 1].*5)) = NaN;
% % msk(~isin(mag(X,Y), [-5 5])) = NaN;
% nn = find(~isnan(x.*msk + y.*msk));
% % scatter(x(nn),y(nn),5,'k', 'filled');
% % scatter(x(nn),y(nn),8 ,col,'^', 'filled', 'Displayname', 'composite @ $x_{wv}=0$km');
% [sts] = boxplot_bybin(gca, y(nn), x(nn), ubins, 'Color', col, 'width', mode(diff(ubins)).*0.4, 'outliars', 0, 'plot', 0, 'FaceAlpha', 0.2, 'nmin',3);
% % errorbar(sts.x0, sts.q2,sts.q2 - sts.q1, sts.q3 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% errorbar(sts.x0, sts.q2,sts.q2 - sts.p10, sts.p90 - sts.q2, 'k.-', 'Color', col, 'MarkerSize',10, 'LineWidth',1)
% 



% cleanLegend(gca, 'northeast', 'FontSize',8, 'Color', 'none', 'EdgeColor', 'none', 'NumColumns',2)
leg = cleanLegend(gca, 'southeast', 'FontSize',8, 'Color', 'none', 'EdgeColor', 'none', 'NumColumns',2);
leg.ItemTokenSize = leg.ItemTokenSize./2;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x0 = -2.5; y0 = 2.5; 
x0 = thisha.Position(1) + thisha.Position(3).*0.05; y0 = sum(thisha.Position([2 4])) - thisha.Position(4).*0.05; 
varx = VARx; vary = VARy;  
varx(isnan(varx + vary)) = NaN; vary(isnan(varx + vary)) = NaN;
msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
varx = varx.*msk; vary = vary.*msk;
cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3)); cc = abs(cc);
cc = cc.^2;
mm = arrayfun(@(i) slpfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3)); mm(isnan(cc)) =NaN;
% varx(:,:,abs(cc) < cclim) = NaN; vary(:,:,abs(cc) < cclim) = NaN;
% cc = arrayfun(@(i) ccfxn(varx(:,:,i), vary(:,:,i)), 1:size(vary,3));
cclims = prctile(cc, [25 75]); 
mmlims = prctile(mm, [25 75]); 
N = sum(~isnan(varx + vary),3); 
msk = ones(size(X)); msk(N < Nt.*prcnan./2) = NaN;
varx = nanmedian(varx,3).*msk; vary = nanmedian(vary,3).*msk; 
cc = ccfxn(varx,vary);
cc = cc.^2;
mm = slpfxn(varx,vary);
textbypos(x0,y0,['$r^2$=' num2str(round(cc,2)) ' (' num2str(round(cclims(1),2)) '-'  num2str(round(cclims(2),2)) ')'], 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize',8)
textbypos(x0,y0-thisha.Position(4).*0.055,['$m$=' num2str(round(mm,2)) ' (' num2str(round(mmlims(1),2)) '-'  num2str(round(mmlims(2),2)) ')'], 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize',8)



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% for hi=1:2; 
%     axes(ha(hi)); hold on; 
%     xline([0 50], 'k--')
% end


for hi=[1 2 4];
    thisha = ha(hi); axes(thisha); hold on; 
    % xline([0 50], 'k--')
    thbins = [0:5:360];
    for r=[50 100 150]
        [xb yb] = pol2cart(deg2rad(thbins), ones(size(thbins)).*r);
        plot(xb,yb,'k:', 'Color', [0 0 0 0.6]);
    end
    % [xb yb] = pol2cart(deg2rad(thbins), ones(size(thbins)).*25);
    % plot(xb,yb,'k:', 'Color', [0 0 0 0.6]);
    % [xb yb] = pol2cart(deg2rad(thbins), ones(size(thbins)).*50);
    % plot(xb,yb,'k:', 'Color', [0 0 0 0.6]);

    c = fixedcolorbar(gca, 'Location','northoutside', 'AxisLocation', 'in', 'FontSize',8);
    c.Position([3 4]) = c.Position([3 4])./[5 2];
    c.Position([1 2]) = thisha.Position([1 2]) + thisha.Position([3 4]).*[0 1];
    c.Position([1 2]) = c.Position([1 2]) - c.Position([3 4]).*[0 1] + c.Position([3 4]).*0.2.*[1 -5];
    c.TickLength = c.TickLength.*6;
    c.Ruler.TickLabelRotation = 0;
    % c.AxisLocation = 'in';
    if hi==1; lab = get(ha(3),'YLabel'); lab = lab.String;
    elseif hi==2; lab = get(ha(3),'XLabel'); lab = lab.String;
    elseif hi==4; lab = 'residual';
    end
    textbypos(sum(c.Position([1 3])) + 0.005, c.Position(2) + c.Position(4)/2, lab, 'Interpreter', 'latex', 'FontSize',get(c, 'FontSize')+2, 'VerticalAlignment', 'middle')
end

for hi=[3]
    thisha = ha(hi); axes(thisha); hold on; 

    c = fixedcolorbar(gca, 'Location','northoutside', 'AxisLocation', 'in', 'FontSize',8);
    c.Position([3 4]) = c.Position([3 4])./[8 2];
    c.Position([1 2]) = thisha.Position([1 2]) + thisha.Position([3 4]).*[1 1];
    
    c.Position([1 2]) = c.Position([1 2]) - c.Position([3 4]).*[1 1] + c.Position([3 4]).*0.2.*[-2 -5];
    c.TickLength = c.TickLength.*6;
    c.Ticks = [0 500]; c.Ruler.TickLabelRotation = 0;

    % c.AxisLocation = 'in';
    if hi==3; lab = 'n';
    end
    textbypos(sum(c.Position([1])) - 0.005, c.Position(2) + c.Position(4)/2, lab,...
        'Interpreter', 'latex', 'FontSize',get(c, 'FontSize')+2, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
end


optstr = '';
if cgopt; optstr = [optstr '_cyclogeostrophy']; end
savejpg(gcf, ['fig_wsp_explainer_phillips' optstr], savepath, saveopt)
