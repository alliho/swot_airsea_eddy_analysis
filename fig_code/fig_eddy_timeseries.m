%% setup
[basepath repopath swotpath datapath] = checksystem();

run([repopath 'code/startup_append.m'])
savepath = [repopath 'figs/']
saveopt = 'time'

addpath([repopath 'code/swot_airsea_eddy_analysis/preprocess/'])


%% [save] [4x1]       plot EDDY and TIME SERIES no vorticity; add quiver
vortopt = 0; 
daterange = [datenum(2023,4,12) datenum(2023,5,12)];
% daterange = [datenum(2023,4,26) datenum(2023,4,30)];
% daterange = [datenum(2023,4,10) datenum(2023,5,12)];
% daterange = [datenum(2023,5,4) datenum(2023,5,7)];
% daterange = [datenum(2023,4,12) datenum(2023,4,20)];
tt = find(isin(eddygrid.time, daterange));
Lavg = 6;
allera5 = 0;
errbopt = 1;
% errbopt = 0;
fs = 10;
prcnan = 0.25;
legvarargin = {'southwest','FontSize', max([fs-3]),'Interpreter', 'latex', 'NumColumns',5, 'Box', 'off'};
prcs = [5 95];
% prcs = [10 90];
% prcs = [25 75];
% prcs = [33 66];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

figure(2315); clf; Nh = 4;
setfigsize(gcf, [947         917])
% setfigsize(gcf, [1         0.96])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
oldha = tight_subplot(Nh, 2, [0.01 0.065], [0.05 0.025], [0.08 0.02]);
[ha] = reorg_tightsubplot(oldha, {1,3,5,7,2,4,6,8});
% -------------------------------------------------------------------------
dx = 0.2; 
for hi=1:4
    ha(hi).Position(3) = ha(hi).Position(3) + dx; 
    ha(hi+Nh).Position(1) = ha(hi+Nh).Position(1) + dx; 
    ha(hi+Nh).Position(3) = ha(hi+Nh).Position(3) - dx; 
end


% -------------------------------------------------------------------------
format_fig4print(ha, 'FontSize',fs);
setaxes(ha(1:Nh), 'XLim', [datenum(2023,4,13) datenum(2023,5,15)]); 
setaxes(ha(1:Nh), 'XLim', [datenum(2023,4,12.5) datenum(2023,5,9.35)]); 

linkaxes(ha(1:Nh), 'x')
for hi=1:Nh; autodatetick(ha(hi), 'x'); end


setaxes(ha((Nh+1):end), 'XLim', [-1 1]*110); 
setaxes(ha((Nh+1):end), 'YLim', [-1 1]*110); 
setaxes(ha((Nh+1):end), 'DataAspectRatio', [1 1 1]); 
setaxes(ha((Nh+1):end), 'YLabel', '$y$ [km]')
% 
% setaxes(ha(5), 'YLabel', '$y_{wnd}$ [km]')
% setaxes(ha(6), 'YLabel', '$y_{wv}$ [km]')

setaxes(ha((Nh+1):end), 'XLabel', '$x$ [km]')
% setaxes(ha(6), 'YLabel', '$y_{wnd}$ [km]')
% setaxes(ha(7), 'YLabel', '$y_{wv}$ [km]')


quietaxes(ha(setxor(1:8,[4 8])), 'x')


setaxes(ha(5), 'Colormap', flipud(cmocean('speed')))
setaxes(ha(5), 'CLim', [-1 1]*1.5)
setaxes(ha(6), 'Colormap', nicejet)
setaxes(ha(6), 'CLim', [-1 1]*0.8)
setaxes(ha(7), 'Colormap', cmocean('thermo'))
setaxes(ha(7), 'CLim', [16.85 19.5])
% setaxes(ha(10), 'Colormap', niceredblue)
% setaxes(ha(10), 'CLim', [-1 1]*1.5)
% setaxes(ha(10), 'CLim', [-1 1]*0.7)

% cmap = cmocean('thermo') + [0 0 0.2]; 
% cmap = cmap.*0.5; 
% cmap = cmap + [0.1 0 0.2];
% cmap(cmap>1) = 1; 
% cmap =  cmap* (1./max(cmap(:))); cmap = cmap*0.8;
% setaxes(ha(10), 'Colormap', cmap)
cmap = flipud(cmocean('deep')); cmap = cmap(20:end-5,:);
setaxes(ha(8), 'Colormap', cmap)
setaxes(ha(8), 'CLim', [0 1.6])
if vortopt; 
    setaxes(ha(8), 'Colormap', niceredblue)
    setaxes(ha(8), 'CLim', [-1 1])
end




% -------------------------------------------------------------------------

msk = NaN(size(X)); msk(mag(X,Y) < 60) = 1;
% -------------------------------------------------------------------------
rotfldstr = '_rot_wnd';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 1;

bVAR = WSP; aVAR = WSPa;
aVAR = remove_planefit(X,Y,aVAR, 'smooth', 50, 'drkm', eddygrid.drkm);


msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
bVAR(:,:,[n_do<prcnan | n_up<prcnan/2]) = NaN;
aVAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;


msk = ones(size(aVAR)); msk(isnan(aVAR)) = NaN; 
aVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(aVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
aVAR = cell2matn(arrayfun(@(i)    filt2(aVAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(aVAR,3), 'Un', 0),3);
aVAR = aVAR.*msk;

msk = ones(size(bVAR)); msk(isnan(bVAR)) = NaN; 
bVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(bVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
bVAR = bVAR.*msk;
% clear i; VAR =  (WSP_U + i.*WSP_V)-[ - (U + i.*V)]; VAR = abs(VAR); VAR = VAR - nanmean(VAR.*edgemsk, [1 2]); hVAR = VAR;




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim');clims = nanmean(clims) + [-1 1].*range(clims)*0.4;
ylabel('wind [m/s]'); ylim([0 16.5]); col = [0 0 0] + 0.2; 

x = era5.time; y = mag(era5.u10, era5.v10); 
y = squeeze(nanmean(y, [1 2])); 
if ~allera5
    x = x(1:4:end); y = y(1:4:end);
end
plot(x,y, 'k.', 'Color', col*0.8, 'DisplayName', 'ERA5')

msk = NaN(size(X)); msk(isin(mag(X,Y), [75 125])) = 1; 
VAR = bVAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
if errbopt; errorbar(x,y,abs(y-lims(:,1)),abs(y-lims(:,2)), 'Color', col, 'LineStyle',':', 'LineWidth',1); 
% else; plot(x,y,'^:','Color', col, 'LineStyle',':', 'LineWidth',1.5, 'MarkerSize',6, 'MarkerFaceColor', 'w'); 
end
    

y0 = y;

col = [0.6 0.1 0.2]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1; 
thlims = 90 + [-1 1].*45 - 20;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth',0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% lab = 
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wnd}>' num2str(thlims(2))  '^o$'])

col = [0.2 0.1 0.6]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1; 
thlims = -90 + [-1 1].*45 + 10;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth',0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wnd}>' num2str(thlims(2))  '^o$'])

leg = cleanLegend(gca,legvarargin{:});




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pos = thisha.Position; xlims = get(gca, 'XLim'); 
pos = thisha.Position; xlims = get(gca, 'XLim'); 
thisha = axes;  
thisha.Position = pos; thisha.Position(4) = thisha.Position(4)/8; thisha.Position(2) = sum(pos([2 4])) - thisha.Position(4);
set(gca, 'XColor', 'none'); set(gca, 'YColor', 'none');  set(gca, 'Color', 'none'); set(gca, 'XTick', ''); set(gca, 'YTick', ''); xlim(xlims); ylim([0 2]);
% box on; 
hold on; 

th = TH_wnd; 
th = squeeze(dirnanmean(th, [1 2]))';  spd = squeeze(nanmean(WSP_model, [1 2]))' ; 
[u v] = pol2cart(deg2rad(th), spd.*0.05);

dar = daspect;
pdar = get(gca, 'PlotBoxAspectRatio');
set(thisha, 'DataAspectRatioMode', 'manual')
darrat = dar(1)./dar(2); pdarrat = pdar(1)/pdar(2);
resetvar = pdarrat./darrat;
ndar = [1 resetvar 1]; rat = ndar(1)./ndar(2);
set(thisha, 'DataAspectRatio', ndar)

quiver(t, ones(size(t)),u.*rat,v, 'ShowArrowHead', 'on', 'Color', 'k', 'LineWidth',1, ...
    'AutoScale', 'off', 'MaxHeadSize', 1, 'DisplayName', 'buoy drift');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha); cla; hold on; 
VAR = aVAR; VAR = VAR(:,:,tt);
N = sum(~isnan(VAR),[3])./size(VAR,3); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(N < prcnan/2) = NaN;
pcolor(X,Y,Z); shading flat;

Z = gaussfilt_2d(Z,mode(mode(diff(X))),10);
% contour(X,Y,Z, [-2:0.4:2], 'k-', 'EdgeAlpha', 0.25)
clims = get(gca, 'CLim');
contour(X,Y,Z, linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha', 0.25)



thlims = 90 + [-1 1].*45 - 10;
rlims = [25 50];
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')

thlims = -90 + [-1 1].*45 + 10;
rlims = [25 50];
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')



x0 = -100+25; y0 = -100; 
sc = 20;
quiver(x0,y0,sc.*cosd(0), sc.*sind(0), 'k-', 'AutoScale','off', 'MaxHeadSize',1, 'LineWidth',1);
text(x0-25,y0,'$\theta_{wnd}$', 'Interpreter', 'latex')

c = cbarfun(thisha, 'WSPa [m/s]');


% -------------------------------------------------------------------------
rotfldstr = '_rot_wv';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 2;
bVAR = SWH; aVAR = SWHa;
% aVAR = remove_planefit(X,Y,aVAR);
aVAR = remove_planefit(X,Y,aVAR, 'smooth', 50, 'drkm', eddygrid.drkm);


msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
bVAR(:,:,[n_do<prcnan | n_up<prcnan/2]) = NaN;
aVAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;


msk = ones(size(aVAR)); msk(isnan(aVAR)) = NaN; 
aVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(aVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
aVAR = cell2matn(arrayfun(@(i)    filt2(aVAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(aVAR,3), 'Un', 0),3);
aVAR = aVAR.*msk;

msk = ones(size(bVAR)); msk(isnan(bVAR)) = NaN; 
bVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(bVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
bVAR = bVAR.*msk;

thwid = 35; rwid = [25 50];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim'); clims = nanmean(clims) + [-1 1].*range(clims)*0.3;
ylabel('SWH [m]'); ylim([0 7.9])

col = [0 0 0];
x = era5.time; y = era5.swh; 
y = squeeze(nanmean(y, [1 2])); 
% x = x(1:6:end); y = y(1:6:end);
if ~allera5
    x = x(1:4:end); y = y(1:4:end);
else
    x = x(1:2:end); y = y(1:2:end);
end

plot(x,y, 'k.', 'Color', col*0.8, 'DisplayName', 'ERA5')

msk = NaN(size(X)); msk(isin(mag(X,Y), [75 125])) = 1; 
VAR = bVAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
if errbopt; errorbar(x,y,abs(y-lims(:,1)),abs(y-lims(:,2)), 'Color', col, 'LineStyle',':', 'LineWidth',1); 
% else; plot(x,y,'^:','Color', col, 'LineStyle',':', 'LineWidth',1.5, 'MarkerSize',6, 'MarkerFaceColor', 'w'); 
end
y0 = y; 


col = [0.6 0.1 0.2]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
thlims = 90 + [-1 1].*thwid - 25;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])

col = [0.2 0.1 0.6]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), rwid+10)) = 1; 
thlims = -90 + [-1 1].*thwid + 25;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
y = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])

% 
% col = [0.2 0.1 0.6]; 
% msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50]+50)) = 1; 
% thlims = -90 + [-1 1].*20 + 50;
% msk(~isin(atan2d(Y,X), thlims)) = NaN;
% VAR = SWHa; VAR = VAR.*msk;
% x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
% col = scaledcmap(cmap,nanmedian(y), clims);
% y = y + y0; lims = lims + y0;
% shaded_limits(x,y,lims, col, 'LineWidth', 0);
% n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% plot(x,y,'k.-', 'color', col, 'LineWidth',2, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])
% 
% col = [0.2 0.1 0.6]; 
% msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50]+60)) = 1; 
% thlims = 90 + [-1 1].*20 - 50;
% % thlims = 180 + [-1 1].*20 - 50;
% msk(~isin(atan2d(Y,X), thlims)) = NaN;
% VAR = SWHa; VAR = VAR.*msk;
% x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
% col = scaledcmap(cmap,nanmedian(y), clims);
% y = y + y0; lims = lims + y0;
% shaded_limits(x,y,lims, col, 'LineWidth', 0);
% n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% plot(x,y,'k.-', 'color', col, 'LineWidth',2, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])

leg = cleanLegend(gca,legvarargin{:});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pos = thisha.Position; xlims = get(gca, 'XLim'); 
thisha = axes;  
thisha.Position = pos; thisha.Position(4) = thisha.Position(4)/8; thisha.Position(2) = sum(pos([2 4])) - thisha.Position(4);
set(gca, 'XColor', 'none'); set(gca, 'YColor', 'none');  set(gca, 'Color', 'none'); set(gca, 'XTick', ''); set(gca, 'YTick', ''); xlim(xlims); ylim([0 2]);
% box on; 
hold on; 

th = TH_wv; 
th = squeeze(dirnanmean(th, [1 2]))';  
spd = squeeze(nanmean(TP, [1 2]))' ; spd = spd.^4; spd = spd.*0.00028;
[u v] = pol2cart(deg2rad(th), spd);

dar = daspect;
pdar = get(gca, 'PlotBoxAspectRatio');
set(thisha, 'DataAspectRatioMode', 'manual')
darrat = dar(1)./dar(2); pdarrat = pdar(1)/pdar(2);
resetvar = pdarrat./darrat;
ndar = [1 resetvar 1]; rat = ndar(1)./ndar(2);
set(thisha, 'DataAspectRatio', ndar)

quiver(t, ones(size(t)),u.*rat,v, 'ShowArrowHead', 'on', 'Color', 'k', 'LineWidth',1, ...
    'AutoScale', 'off', 'MaxHeadSize', 1, 'DisplayName', 'buoy drift');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha); cla; hold on; 


VAR = aVAR; VAR = VAR(:,:,tt);
N = sum(~isnan(VAR),[3])./size(VAR,3); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(N < prcnan/2) = NaN;
pcolor(X,Y,Z); shading flat;


Z = gaussfilt_2d(Z,mode(mode(diff(X))),10);
% contour(X,Y,Z, [-1:0.2:1], 'k-', 'EdgeAlpha', 0.25)
clims = get(gca, 'CLim');
contour(X,Y,Z, linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha', 0.25)

thlims = 90 + [-1 1].*thwid - 25;
rlims = rwid;
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')

thlims = -90 + [-1 1].*thwid + 25;
rlims = rwid+10;
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')

x0 = -100+25; y0 = -100; 
sc = 20;
quiver(x0,y0,sc.*cosd(0), sc.*sind(0), 'k-', 'AutoScale','off', 'MaxHeadSize',1, 'LineWidth',1);
text(x0-25,y0,'$\theta_{wv}$', 'Interpreter', 'latex')



c = cbarfun(thisha, 'SWHa [m]');
c.Ticks = [-0.5 0 0.5];
c.Ruler.TickLabelRotation = 0;
c.FontSize = 6;


% -------------------------------------------------------------------------
rotfldstr = '';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 3;
bVAR = SST; aVAR = SSTa; oVAR = SST_L4; Rc = 25;
aVAR = remove_planefit(X,Y,aVAR);

% bVAR = movmean(bVAR, 3, 3, 'omitnan');
% aVAR = movmean(aVAR, 3, 3, 'omitnan');
% 
msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n);
bVAR(:,:,[n<prcnan./2/2]) = NaN;
aVAR(:,:,[n<prcnan./2]) = NaN;
% oVAR(:,:,[n<prcnan./2]) = NaN;

msk = ones(size(aVAR)); msk(isnan(aVAR)) = NaN; 
aVAR = gaussfilt_2d(aVAR, eddygrid.drkm, Lavg);
aVAR = aVAR.*msk;

msk = ones(size(bVAR)); msk(isnan(bVAR)) = NaN; 
bVAR = gaussfilt_2d(bVAR, eddygrid.drkm, Lavg);
bVAR = bVAR.*msk;

msk = ones(size(oVAR)); msk(isnan(oVAR)) = NaN; 
oVAR = gaussfilt_2d(oVAR, eddygrid.drkm, Lavg);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); cla; hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim'); 
ylabel('SST [$^o$C]'); ylim([15 21.2])

col = [0 0 0]+0.2;
x = era5.time; 
y = era5.sst - 273.15;
y = squeeze(nanmean(y, [1 2])); 
x = x(1:12:end); y = y(1:12:end);
% y = movmean(y,9); y = movmean(y,15); 
plot(x,y, 'k.', 'Color', col*0.8, 'DisplayName', 'ERA5')



%%% PLOT L4
msk = ones(size(X)); msk(~isin(mag(X,Y), [0 Rc])) = NaN; 
VAR = oVAR.*msk; 
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims, col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',2, 'MarkerSize',10, 'DisplayName', ['GHRSST L4 $R<' num2str(Rc) '$km'])

%%% PLOT L4
msk = ones(size(X)); msk(~isin(mag(X,Y), [75 150])) = NaN; 
VAR = oVAR.*msk; 
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims, col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',2, 'MarkerSize',10, 'DisplayName', ['GHRSST L4 $R>75$km'])




%%% DIRECTLY PLOT GHRSST
gtt = find(isin(ghrsst_l3.t0, [datenum(2023,4,12) datenum(2023,5,12)] + [-1 1].*5));

x = ghrsst_l3.t0(gtt); y = NaN(size(x)); ylims = NaN(2,length(x));
for i = 1:length(gtt)
    gti = gtt(i); var = ghrsst_l3.sst_proc_patch(:,:,gti);
    [~, eti] = min(abs(ghrsst_l3.t0(gti) - eddygrid.time));
    [glo gla] = ndgrid(ghrsst_l3.lon, ghrsst_l3.lat);
    dist = distance(eddygrid.center_lat(eti), eddygrid.center_lon(eti), gla, glo, referenceEllipsoid('wgs84'));
    dist = dist./1000; 
    msk = NaN(size(dist)); msk(isin(dist, [0 30])) = 1;
    % ii = find(isin(dist, [75 150]));
    if sum(~isnan(var.*msk),[1 2])./sum(~isnan(msk),[1 2]) > 0.05 | sum(~isnan(var.*msk),[1 2]) > 20
        y(i) = nanmedian(var.*msk, [1 2]);
        ylims(:,i) = squeeze(prctile(var.*msk, prcs, [1 2]));
    end
end
ylims = movmean(ylims,12,2, 'omitnan');
y = movmean(y,3, 'omitnan');
col = scaledcmap(cmap,nanmedian(y), clims);
nn = find(abs(diff(isnan(ylims(1,:))))); 
for i=1:length(nn)-1
    ii = nn(i):nn(i+1)-1;
    shaded_limits(x(ii),y(ii),ylims(:,ii)', col, 'LineWidth',0);
end
plot(x,y,'k.', 'color', col, 'LineWidth',3, 'MarkerSize',5, 'DisplayName', ['GOES16 $R<'  num2str(Rc) '$km'])

x = ghrsst_l3.t0(gtt); y = NaN(size(x)); ylims = NaN(2,length(x));
for i = 1:length(gtt)
    gti = gtt(i); var = ghrsst_l3.sst_proc_patch(:,:,gti);
    [~, eti] = min(abs(ghrsst_l3.t0(gti) - eddygrid.time));
    [glo gla] = ndgrid(ghrsst_l3.lon, ghrsst_l3.lat);
    dist = distance(eddygrid.center_lat(eti), eddygrid.center_lon(eti), gla, glo, referenceEllipsoid('wgs84'));
    dist = dist./1000; 
    msk = NaN(size(dist)); msk(isin(dist, [75 150])) = 1;
    % ii = find(isin(dist, [75 150]));
    if sum(~isnan(var.*msk),[1 2])./sum(~isnan(msk),[1 2]) > 0.05 | sum(~isnan(var.*msk),[1 2]) > 50
        y(i) = nanmedian(var.*msk, [1 2]);
        ylims(:,i) = squeeze(prctile(var.*msk, prcs, [1 2]));
    end
end
ylims = movmean(ylims,12,2, 'omitnan');
y = movmean(y,3, 'omitnan');
col = scaledcmap(cmap,nanmedian(y), clims);
nn = find(abs(diff(isnan(ylims(1,:))))); 
for i=1:length(nn)-1
    ii = nn(i):nn(i+1)-1;
    shaded_limits(x(ii),y(ii),ylims(:,ii)', col, 'LineWidth',0);
end
plot(x,y,'k.', 'color', col, 'LineWidth',3, 'MarkerSize',5, 'DisplayName', 'GOES16 $R>75$km')


% 
% msk = NaN(size(X)); msk(mag(X,Y) < 30) = 1;
% VAR = bVAR; VAR = VAR.*msk; 
% x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
% col = scaledcmap(cmap,nanmedian(y), clims);
% % shaded_limits(x,y,lims', col, 'LineWidth',0);
% n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% plot(x,y,'k.', 'color', col, 'LineWidth',3, 'MarkerSize',10)
% % errorbar(x,y,abs(y'-lims(1,:)),abs(y'-lims(2,:)), 'color',col.*0.9, 'LineWidth',0.2, 'LineStyle', 'none');
% 
% msk = NaN(size(X)); msk(isin(mag(X,Y), [75 150])) = 1; 
% % msk = NaN(size(X)); msk(isin(mag(X,Y), [75 90])) = 1; 
% VAR = bVAR; VAR = VAR.*msk;
% x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
% col = scaledcmap(cmap,nanmedian(y), clims);
% % shaded_limits(x,y,lims', col, 'LineWidth',0);
% n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% plot(x,y,'k.', 'color', col, 'LineWidth',3, 'MarkerSize',10)
% % errorbar(x,y,abs(y'-lims(1,:)),abs(y'-lims(2,:)), 'color',col.*0.9, 'LineWidth',0.4, 'LineStyle', 'none');
% 

leg = cleanLegend(gca,legvarargin{:});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha); hold on; 
% VAR = SST - nanmean(SST,[1 2]);
VAR = bVAR; VAR = VAR(:,:,tt);
nmsk = sum(~isnan(VAR),[3]); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(nmsk<length(tt).*0.1) = NaN;
pcolor(X,Y,Z); shading flat;
% c = fixedcolorbar(gca, 'Location', 'northoutside');

clims = get(gca, 'CLim');
contour(X,Y,Z, linspace(clims(1), clims(2), 20), 'k-', 'EdgeAlpha', 0.25)

rlims = [Rc 75];
thtx = 0:5:360;
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')

c = cbarfun(thisha, 'SSTa [$^o$C]');
c.Ticks = [17 19];
c.Ruler.TickLabelRotation = 0;
c.FontSize = 7;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 4;
bVAR = mag(U,V); aVAR = mag(U,V); 
if vortopt
    bVAR = relvort; 
    aVAR = relvort; 
end

msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n);
bVAR(:,:,[n<prcnan/2]) = NaN;
aVAR(:,:,[n<prcnan]) = NaN;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); cla; hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim');
ylabel('$|U_{geo}|$ [m/s]'); ylim([0 1.8])
if vortopt; ylabel('$\zeta/f$'); ylim([-1 1.5].*1.2); end
msk = NaN(size(X)); msk(isin(mag(X,Y), [75 125])) = 1; 
if vortopt; 
    msk = NaN(size(X)); msk(isin(mag(X,Y), [50 75])) = 1; 
end
VAR = bVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims', col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', '$R>75$km')

msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1;
if vortopt; 
    msk = NaN(size(X)); msk(isin(mag(X,Y), [0 45])) = 1; 
end
VAR = bVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims', col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', '$25$km $<R<50$km')

leg = cleanLegend(gca,legvarargin{:});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha);cla; hold on; 
VAR = aVAR; VAR = VAR(:,:,tt);
N = sum(~isnan(VAR),[3])./size(VAR,3); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(N < prcnan/2) = NaN;
pcolor(X,Y,Z); shading flat;
% c = fixedcolorbar(gca, 'Location', 'northoutside');

rlims = [25 50]; if vortopt; rlims = [0 45 50 75]; end
thtx = 0:5:360;
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')


qU = nanmedian(U,3); qU(isnan(Z)) = NaN;
qV = nanmedian(V,3); qV(isnan(Z)) = NaN;
di = 4; xx = 2:di:size(X,1); yy = 2:di:size(X,2); 
sc = 7;
quiver(X(xx,yy),Y(xx,yy),sc.*qU(xx,yy),sc.*qV(xx,yy), 'k-', 'AutoScale', 'off', 'Color', [1 1 1].*0.2)

c = cbarfun(thisha, '|$u_{geo}|$');


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
for hi=(Nh+1):length(ha)
    axes(ha(hi)); 
    rlims = [0:25:105];
    thtx = 0:5:360;
    [xtx ytx] = pol2cart(deg2rad(thtx), rlims');
    plot(xtx',ytx','k-', 'Color', [0 0 0 0.08])


end

% for hi=(1):Nh
%     axes(ha(hi)); 
%     xline(daterange)
% 
% 
% end


xlims = get(ha(1), 'XLim');
alignyaxes(ha(1:Nh), 'offset',-1)

optstr = '';

if ~vortopt; optstr = ['novort' optstr]; end
if errbopt; optstr = [optstr '_werrorbars']; end

savejpg(gcf, ['timeseries_andxy' optstr '_w_quiv_' savestr], savepath, saveopt)



%% [save] basic [4x1]        plot EDDY and TIME SERIES no vorticity; add quiver
vortopt = 0; 
daterange = [datenum(2023,4,12) datenum(2023,5,12)];
% daterange = [datenum(2023,4,26) datenum(2023,4,30)];
% daterange = [datenum(2023,4,10) datenum(2023,5,12)];
% daterange = [datenum(2023,5,4) datenum(2023,5,7)];
% daterange = [datenum(2023,4,12) datenum(2023,4,20)];
tt = find(isin(eddygrid.time, daterange));
Lavg = 6;
allera5 = 0;
errbopt = 1;
fs = 10;
prcnan = 0.25;
legvarargin = {'southwest','FontSize', max([fs-3]),'Interpreter', 'latex', 'NumColumns',3, 'Box', 'off'};
prcs = [5 95];
% prcs = [10 90];
% prcs = [25 75];
% prcs = [33 66];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

figure(2315); clf; Nh = 4;
setfigsize(gcf, [947         917]*0.8)
% setfigsize(gcf, [1         0.96])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
oldha = tight_subplot(Nh, 2, [0.01 0.065], [0.05 0.025], [0.08 0.02]);
[ha] = reorg_tightsubplot(oldha, {1,3,5,7,2,4,6,8});
% -------------------------------------------------------------------------
dx = 0.2; 
for hi=1:4
    ha(hi).Position(3) = ha(hi).Position(3) + dx; 
    ha(hi+Nh).Position(1) = ha(hi+Nh).Position(1) + dx; 
    ha(hi+Nh).Position(3) = ha(hi+Nh).Position(3) - dx; 
end


% -------------------------------------------------------------------------
format_fig4print(ha, 'FontSize',fs);
setaxes(ha(1:Nh), 'XLim', [datenum(2023,4,13) datenum(2023,5,15)]); 
linkaxes(ha(1:Nh), 'x')
for hi=1:Nh; autodatetick(ha(hi), 'x'); end


setaxes(ha((Nh+1):end), 'XLim', [-1 1]*110); 
setaxes(ha((Nh+1):end), 'YLim', [-1 1]*110); 
setaxes(ha((Nh+1):end), 'DataAspectRatio', [1 1 1]); 
setaxes(ha((Nh+1):end), 'YLabel', '$y$ [km]')
% 
% setaxes(ha(5), 'YLabel', '$y_{wnd}$ [km]')
% setaxes(ha(6), 'YLabel', '$y_{wv}$ [km]')

setaxes(ha((Nh+1):end), 'XLabel', '$x$ [km]')
% setaxes(ha(6), 'YLabel', '$y_{wnd}$ [km]')
% setaxes(ha(7), 'YLabel', '$y_{wv}$ [km]')


quietaxes(ha(setxor(1:8,[4 8])), 'x')


setaxes(ha(5), 'Colormap', flipud(cmocean('speed')))
setaxes(ha(5), 'CLim', [-1 1]*1.5)
setaxes(ha(6), 'Colormap', nicejet)
setaxes(ha(6), 'CLim', [-1 1]*0.8)
setaxes(ha(7), 'Colormap', niceredblue)
setaxes(ha(7), 'CLim', [-1 1]*0.6)
cmap = flipud(cmocean('deep')); cmap = cmap(20:end-5,:);
setaxes(ha(8), 'Colormap', cmap)
setaxes(ha(8), 'CLim', [0 1.6])
if vortopt; 
    setaxes(ha(8), 'Colormap', niceredblue)
    setaxes(ha(8), 'CLim', [-1 1])
end




% -------------------------------------------------------------------------

msk = NaN(size(X)); msk(mag(X,Y) < 60) = 1;
% -------------------------------------------------------------------------
rotfldstr = '_rot_wnd';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 1;

bVAR = WSP; aVAR = WSPa;
aVAR = remove_planefit(X,Y,aVAR, 'smooth', 50, 'drkm', eddygrid.drkm);


msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
bVAR(:,:,[n_do<prcnan | n_up<prcnan/2]) = NaN;
aVAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;


msk = ones(size(aVAR)); msk(isnan(aVAR)) = NaN; 
aVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(aVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
aVAR = cell2matn(arrayfun(@(i)    filt2(aVAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(aVAR,3), 'Un', 0),3);
aVAR = aVAR.*msk;

msk = ones(size(bVAR)); msk(isnan(bVAR)) = NaN; 
bVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(bVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
bVAR = bVAR.*msk;
% clear i; VAR =  (WSP_U + i.*WSP_V)-[ - (U + i.*V)]; VAR = abs(VAR); VAR = VAR - nanmean(VAR.*edgemsk, [1 2]); hVAR = VAR;




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim');clims = nanmean(clims) + [-1 1].*range(clims)*0.4;
ylabel('wind [m/s]'); ylim([0 16.5]); col = [0 0 0] + 0.2; 

x = era5.time; y = mag(era5.u10, era5.v10); 
y = squeeze(nanmean(y, [1 2])); 
if ~allera5
    x = x(1:4:end); y = y(1:4:end);
end
plot(x,y, 'k.', 'Color', col*0.8, 'DisplayName', 'ERA5')

msk = NaN(size(X)); msk(isin(mag(X,Y), [75 125])) = 1; 
VAR = bVAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
if errbopt; errorbar(x,y,abs(y-lims(:,1)),abs(y-lims(:,2)), 'Color', col, 'LineStyle',':', 'LineWidth',1); 
else; plot(x,y,'^:','Color', col, 'LineStyle',':', 'LineWidth',1.5, 'MarkerSize',6, 'MarkerFaceColor', 'w'); end
    

y0 = y;

col = [0.6 0.1 0.2]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1; 
thlims = 90 + [-1 1].*45 - 20;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth',0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% lab = 
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wnd}>' num2str(thlims(2))  '^o$'])

col = [0.2 0.1 0.6]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1; 
thlims = -90 + [-1 1].*45 + 10;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth',0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wnd}>' num2str(thlims(2))  '^o$'])

leg = cleanLegend(gca,legvarargin{:});




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pos = thisha.Position; xlims = get(gca, 'XLim'); 
pos = thisha.Position; xlims = get(gca, 'XLim'); 
thisha = axes;  
thisha.Position = pos; thisha.Position(4) = thisha.Position(4)/8; thisha.Position(2) = sum(pos([2 4])) - thisha.Position(4);
set(gca, 'XColor', 'none'); set(gca, 'YColor', 'none');  set(gca, 'Color', 'none'); set(gca, 'XTick', ''); set(gca, 'YTick', ''); xlim(xlims); ylim([0 2]);
% box on; 
hold on; 

th = TH_wnd; 
th = squeeze(dirnanmean(th, [1 2]))';  spd = squeeze(nanmean(WSP_model, [1 2]))' ; 
[u v] = pol2cart(deg2rad(th), spd.*0.05);

dar = daspect;
pdar = get(gca, 'PlotBoxAspectRatio');
set(thisha, 'DataAspectRatioMode', 'manual')
darrat = dar(1)./dar(2); pdarrat = pdar(1)/pdar(2);
resetvar = pdarrat./darrat;
ndar = [1 resetvar 1]; rat = ndar(1)./ndar(2);
set(thisha, 'DataAspectRatio', ndar)

quiver(t, ones(size(t)),u.*rat,v, 'ShowArrowHead', 'on', 'Color', 'k', 'LineWidth',1, ...
    'AutoScale', 'off', 'MaxHeadSize', 1, 'DisplayName', 'buoy drift');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha); cla; hold on; 
VAR = aVAR; VAR = VAR(:,:,tt);
N = sum(~isnan(VAR),[3])./size(VAR,3); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(N < prcnan/2) = NaN;
pcolor(X,Y,Z); shading flat;

Z = gaussfilt_2d(Z,mode(mode(diff(X))),10);
contour(X,Y,Z, [-2:0.4:2], 'k-', 'EdgeAlpha', 0.25)



thlims = 90 + [-1 1].*45 - 10;
rlims = [25 50];
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')

thlims = -90 + [-1 1].*45 + 10;
rlims = [25 50];
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')



x0 = -100+25; y0 = -100; 
sc = 20;
quiver(x0,y0,sc.*cosd(0), sc.*sind(0), 'k-', 'AutoScale','off', 'MaxHeadSize',1, 'LineWidth',1);
text(x0-25,y0,'$\theta_{wnd}$', 'Interpreter', 'latex')

c = cbarfun(thisha, 'WSPa [m/s]');


% -------------------------------------------------------------------------
rotfldstr = '_rot_wv';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 2;
bVAR = SWH; aVAR = SWHa;
% aVAR = remove_planefit(X,Y,aVAR);
aVAR = remove_planefit(X,Y,aVAR, 'smooth', 50, 'drkm', eddygrid.drkm);


msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [0 180])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_up = n;
msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100]) & isin(atan2d(Y,X), [-180 0])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n); n_do = n;
bVAR(:,:,[n_do<prcnan | n_up<prcnan/2]) = NaN;
aVAR(:,:,[n_do<prcnan | n_up<prcnan]) = NaN;


msk = ones(size(aVAR)); msk(isnan(aVAR)) = NaN; 
aVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(aVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
aVAR = cell2matn(arrayfun(@(i)    filt2(aVAR(:,:,i),eddygrid.drkm, 500, 'hp'), 1:size(aVAR,3), 'Un', 0),3);
aVAR = aVAR.*msk;

msk = ones(size(bVAR)); msk(isnan(bVAR)) = NaN; 
bVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(bVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
bVAR = bVAR.*msk;

thwid = 35; rwid = [25 50];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim'); clims = nanmean(clims) + [-1 1].*range(clims)*0.3;
ylabel('SWH [m]'); ylim([0 7.9])

col = [0 0 0];
x = era5.time; y = era5.swh; 
y = squeeze(nanmean(y, [1 2])); 
% x = x(1:6:end); y = y(1:6:end);
if ~allera5
    x = y(1:4:end); y = y(1:4:end);
else
    x = y(1:2:end); y = y(1:2:end);
end

plot(x,y, 'k.', 'Color', col*0.8, 'DisplayName', 'ERA5')

msk = NaN(size(X)); msk(isin(mag(X,Y), [75 125])) = 1; 
VAR = bVAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
if errbopt; errorbar(x,y,abs(y-lims(:,1)),abs(y-lims(:,2)), 'Color', col, 'LineStyle',':', 'LineWidth',1); 
else; plot(x,y,'^:','Color', col, 'LineStyle',':', 'LineWidth',1.5, 'MarkerSize',6, 'MarkerFaceColor', 'w'); end
y0 = y; 


col = [0.6 0.1 0.2]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), rwid)) = 1; 
thlims = 90 + [-1 1].*thwid - 25;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])

col = [0.2 0.1 0.6]; 
msk = NaN(size(X)); msk(isin(mag(X,Y), rwid+10)) = 1; 
thlims = -90 + [-1 1].*thwid + 25;
msk(~isin(atan2d(Y,X), thlims)) = NaN;
VAR = aVAR; VAR = VAR.*msk;
y = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
col = scaledcmap(cmap,nanmedian(y), clims);
y = y + y0; lims = lims + y0;
shaded_limits(x,y,lims, col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])

% 
% col = [0.2 0.1 0.6]; 
% msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50]+50)) = 1; 
% thlims = -90 + [-1 1].*20 + 50;
% msk(~isin(atan2d(Y,X), thlims)) = NaN;
% VAR = SWHa; VAR = VAR.*msk;
% x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
% col = scaledcmap(cmap,nanmedian(y), clims);
% y = y + y0; lims = lims + y0;
% shaded_limits(x,y,lims, col, 'LineWidth', 0);
% n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% plot(x,y,'k.-', 'color', col, 'LineWidth',2, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])
% 
% col = [0.2 0.1 0.6]; 
% msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50]+60)) = 1; 
% thlims = 90 + [-1 1].*20 - 50;
% % thlims = 180 + [-1 1].*20 - 50;
% msk(~isin(atan2d(Y,X), thlims)) = NaN;
% VAR = SWHa; VAR = VAR.*msk;
% x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims)';
% col = scaledcmap(cmap,nanmedian(y), clims);
% y = y + y0; lims = lims + y0;
% shaded_limits(x,y,lims, col, 'LineWidth', 0);
% n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
% plot(x,y,'k.-', 'color', col, 'LineWidth',2, 'MarkerSize',10, 'DisplayName', ['$' num2str(thlims(1)) '^o>\theta_{wv}>' num2str(thlims(2))  '^o$'])

leg = cleanLegend(gca,legvarargin{:});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pos = thisha.Position; xlims = get(gca, 'XLim'); 
thisha = axes;  
thisha.Position = pos; thisha.Position(4) = thisha.Position(4)/8; thisha.Position(2) = sum(pos([2 4])) - thisha.Position(4);
set(gca, 'XColor', 'none'); set(gca, 'YColor', 'none');  set(gca, 'Color', 'none'); set(gca, 'XTick', ''); set(gca, 'YTick', ''); xlim(xlims); ylim([0 2]);
% box on; 
hold on; 

th = TH_wv; 
th = squeeze(dirnanmean(th, [1 2]))';  
spd = squeeze(nanmean(TP, [1 2]))' ; spd = spd.^4; spd = spd.*0.00028;
[u v] = pol2cart(deg2rad(th), spd);

dar = daspect;
pdar = get(gca, 'PlotBoxAspectRatio');
set(thisha, 'DataAspectRatioMode', 'manual')
darrat = dar(1)./dar(2); pdarrat = pdar(1)/pdar(2);
resetvar = pdarrat./darrat;
ndar = [1 resetvar 1]; rat = ndar(1)./ndar(2);
set(thisha, 'DataAspectRatio', ndar)

quiver(t, ones(size(t)),u.*rat,v, 'ShowArrowHead', 'on', 'Color', 'k', 'LineWidth',1, ...
    'AutoScale', 'off', 'MaxHeadSize', 1, 'DisplayName', 'buoy drift');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha); cla; hold on; 


VAR = aVAR; VAR = VAR(:,:,tt);
N = sum(~isnan(VAR),[3])./size(VAR,3); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(N < prcnan/2) = NaN;
pcolor(X,Y,Z); shading flat;


Z = gaussfilt_2d(Z,mode(mode(diff(X))),10);
contour(X,Y,Z, [-1:0.2:1], 'k-', 'EdgeAlpha', 0.25)

thlims = 90 + [-1 1].*thwid - 25;
rlims = rwid;
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')

thlims = -90 + [-1 1].*thwid + 25;
rlims = rwid+10;
thtx = thlims(1):5:thlims(2);
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')
[xtx ytx] = pol2cart(deg2rad(thlims), rlims');
plot(xtx,ytx,'k:')

x0 = -100+25; y0 = -100; 
sc = 20;
quiver(x0,y0,sc.*cosd(0), sc.*sind(0), 'k-', 'AutoScale','off', 'MaxHeadSize',1, 'LineWidth',1);
text(x0-25,y0,'$\theta_{wv}$', 'Interpreter', 'latex')



c = cbarfun(thisha, 'SWHa [m]');

% -------------------------------------------------------------------------
rotfldstr = '';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 3;
bVAR = SSHa; aVAR = SSHa;
aVAR = remove_planefit(X,Y,aVAR);

% bVAR = movmean(bVAR, 3, 3, 'omitnan');
% aVAR = movmean(aVAR, 3, 3, 'omitnan');
% 
msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n);
bVAR(:,:,[n<prcnan./2/2]) = NaN;
aVAR(:,:,[n<prcnan./2]) = NaN;

msk = ones(size(aVAR)); msk(isnan(aVAR)) = NaN; 
aVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(aVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
aVAR = aVAR.*msk;

msk = ones(size(bVAR)); msk(isnan(bVAR)) = NaN; 
bVAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(bVAR(:,:,i),eddygrid.drkm, Lavg), 1:size(aVAR,3), 'Un', 0),3);
bVAR = bVAR.*msk;



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); cla; hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim'); 
ylabel('SSHa [m]'); ylim([-1 1])

msk = NaN(size(X)); msk(mag(X,Y) < 25) = 1;
VAR = bVAR; VAR = VAR.*msk; 
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims', col, 'LineWidth',0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', '$R<25$km')
% errorbar(x,y,abs(y'-lims(1,:)),abs(y'-lims(2,:)), 'color',col.*0.9, 'LineWidth',0.2, 'LineStyle', 'none');

msk = NaN(size(X)); msk(isin(mag(X,Y), [75 150])) = 1; 
% msk = NaN(size(X)); msk(isin(mag(X,Y), [75 90])) = 1; 
VAR = bVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims', col, 'LineWidth',0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', '$75$km $<R<100$km')
% errorbar(x,y,abs(y'-lims(1,:)),abs(y'-lims(2,:)), 'color',col.*0.9, 'LineWidth',0.4, 'LineStyle', 'none');


leg = cleanLegend(gca,legvarargin{:});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha); hold on; 
% VAR = SST - nanmean(SST,[1 2]);
VAR = bVAR; VAR = VAR(:,:,tt);
nmsk = sum(~isnan(VAR),[3]); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(nmsk<length(tt).*0.1) = NaN;
pcolor(X,Y,Z); shading flat;
% c = fixedcolorbar(gca, 'Location', 'northoutside');

rlims = [25 75];
thtx = 0:5:360;
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')

c = cbarfun(thisha, '');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
hi = 4;
bVAR = mag(U,V); aVAR = mag(U,V); 
if vortopt
    bVAR = relvort; 
    aVAR = relvort; 
end

msk = ones(size(X)); msk(~[isin(mag(X,Y), [0 100])] ) = NaN;
n = sum(~isnan(msk.*[aVAR + bVAR]), [1 2])./sum(~isnan(msk), [1 2]); n = squeeze(n);
bVAR(:,:,[n<prcnan/2]) = NaN;
aVAR(:,:,[n<prcnan]) = NaN;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi); axes(thisha); cla; hold on; 
cmap = get(ha(hi+Nh), 'Colormap'); clims = get(ha(hi+Nh), 'Clim');
ylabel('$|U_{geo}|$ [m/s]'); ylim([0 1.8])
if vortopt; ylabel('$\zeta/f$'); ylim([-1 1.5].*1.2); end
msk = NaN(size(X)); msk(isin(mag(X,Y), [75 125])) = 1; 
if vortopt; 
    msk = NaN(size(X)); msk(isin(mag(X,Y), [50 75])) = 1; 
end
VAR = bVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims', col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', '$R>75$km')

msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1;
if vortopt; 
    msk = NaN(size(X)); msk(isin(mag(X,Y), [0 45])) = 1; 
end
VAR = bVAR; VAR = VAR.*msk;
x = t; y = squeeze(nanmedian(VAR,[1 2])); lims = prctile(VAR, prcs, [1 2]); lims = squeeze(lims);
col = scaledcmap(cmap,nanmedian(y), clims);
shaded_limits(x,y,lims', col, 'LineWidth', 0);
n = squeeze(sum(~isnan(VAR), [1 2])); y(n<10) = NaN;
plot(x,y,'k.-', 'color', col, 'LineWidth',3, 'MarkerSize',10, 'DisplayName', '$25$km $<R<50$km')

leg = cleanLegend(gca,legvarargin{:});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(hi + Nh); axes(thisha);cla; hold on; 
VAR = aVAR; VAR = VAR(:,:,tt);
N = sum(~isnan(VAR),[3])./size(VAR,3); 
Z = nanmedian(VAR,3);
Z = gaussfilt_2d(Z,eddygrid.drkm, Lavg); 
Z(N < prcnan/2) = NaN;
pcolor(X,Y,Z); shading flat;
% c = fixedcolorbar(gca, 'Location', 'northoutside');

rlims = [25 50]; if vortopt; rlims = [0 45 50 75]; end
thtx = 0:5:360;
[xtx ytx] = pol2cart(deg2rad(thtx), rlims');
plot(xtx',ytx','k:')


qU = nanmedian(U,3); qU(isnan(Z)) = NaN;
qV = nanmedian(V,3); qV(isnan(Z)) = NaN;
di = 4; xx = 2:di:size(X,1); yy = 2:di:size(X,2); 
sc = 7;
quiver(X(xx,yy),Y(xx,yy),sc.*qU(xx,yy),sc.*qV(xx,yy), 'k-', 'AutoScale', 'off', 'Color', [1 1 1].*0.2)

c = cbarfun(thisha, '');



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
for hi=(Nh+1):length(ha)
    axes(ha(hi)); 
    rlims = [0:25:105];
    thtx = 0:5:360;
    [xtx ytx] = pol2cart(deg2rad(thtx), rlims');
    plot(xtx',ytx','k-', 'Color', [0 0 0 0.08])
end

% for hi=(1):Nh
%     axes(ha(hi)); 
%     xline(daterange)
% 
% 
% end


xlims = get(ha(1), 'XLim');
alignyaxes(ha(1:Nh), xlims(1)-3)


vortstr = 'vort'; if ~vortopt; vortstr = ['no' vortstr]; end
savejpg(gcf, ['timeseries_andxy_' vortstr '_w_quiv_w_SSHa_' savestr], savepath, saveopt)






%% TIME SERIES: CLEAN IT UP CLEAN THE OLD CODE!!!
A = {}; Alabs = {}; Amsks = {};
rotfldstr = '_rot_wnd';
run load_vars_into_memory_universal
A{1} = {WSP, mag(WSP_U_ASCAT,WSP_V_ASCAT) mag(WSP_U_model,WSP_V_model)};
Alabs{1} = 'WSP [m/s]';
Acmaps{1} = flipud(cmocean('speed'));  Aclims{1} = [-1 1]*1.5;

msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1; 
thlims = 90 + [-1 1].*45 ;
mskhi = msk; 
mskhi(~isin(atan2d(Y,X), thlims)) = NaN;
thlims = -90 + [-1 1].*45 ;
msklo = msk; 
msklo(~isin(atan2d(Y,X), thlims)) = NaN;
Amsks{1} = {mskhi, msklo};

rotfldstr = '_rot_wv';
run load_vars_into_memory_universal
A{2} = {SWH, NaN(size(SWH)),SWH_model};
Alabs{2} = 'SWH [m]'
Acmaps{2} = nicejet; Aclims{2} = [-1 1]*0.56;
msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1; 
thlims = 90 + [-1 1].*35 - 25;
mskhi = msk; 
mskhi(~isin(atan2d(Y,X), thlims)) = NaN;

thlims = -90 + [-1 1].*35 + 25;
msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50]+10)) = 1; 
msklo = msk; 
msklo(~isin(atan2d(Y,X), thlims)) = NaN;
Amsks{2} = {mskhi, msklo};

rotfldstr = '';
run load_vars_into_memory_universal
A{3} = {mag(U,V), mag(U_NEUROST,V_NEUROST),mag(U_DUACS,V_DUACS)};
Alabs{3} = '$|\vec{u}_{geo}|$'
Acmaps{3} = flipud(cmocean('deep')); Aclims{3} = [0 1.5];

msk = NaN(size(X)); msk(isin(mag(X,Y), [25 50])) = 1; 
mskhi = msk; 
msk = NaN(size(X)); msk(isin(mag(X,Y), [75 200])) = 1; 
msklo = msk; 
Amsks{3} = {mskhi, msklo};


% A
% -------------------------------------------------------------------------  % -------------------------------------------------------------------------
% -------------------------------------------------------------------------  % -------------------------------------------------------------------------
% -------------------------------------------------------------------------  % -------------------------------------------------------------------------

NA = length(A);

% -------------------------------------------------------------------------  % -------------------------------------------------------------------------
figure(1555); clf; 
setfigsize(gcf, [916   452])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(NA, 4, [0.03 0.01], [0.09 0.05], [0.05 0.4]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
dx = 0.3; dw = 0.05;
for hi=1:4:length(ha);
    ha(hi).Position(3) = ha(hi).Position(3) + dx;
end
for hi=setxor(1:length(ha), [1:4:NA*4])
    ha(hi).Position(1) = ha(hi).Position(1) + dx + dw;
end
% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',10);
setaxes(ha(setxor(1:length(ha), [1:4:NA*4])), 'DataAspectRatio', [1 1 1]); 
% setaxes(ha(2:2:end), 'XLim', [-1 1]*100); 
% setaxes(ha(2:2:end), 'YLim', [-1 1]*100); 
setaxes(ha(setxor(1:length(ha), [1:4:NA*4])), 'XLim', [-1 1]*100); 
setaxes(ha(setxor(1:length(ha), [1:4:NA*4])), 'YLim', [-1 1]*100); 
quietaxes(ha(setxor(1:length(ha), [1:4:NA*4 2:4:NA*4])), 'y')
quietaxes(ha(setxor(1:length(ha), [NA*4:-1:NA*3])), 'x')


setaxes(ha(1:4:end), 'XLim', [datenum(2023,4,15) datenum(2023,5,12)]); 
autodatetick(ha(1:4:end), 'x')
quietaxes(ha(1:4:end-4), 'x')


% linkaxes(ha)
% format_mapticks(ha, 'xy');

for hi=1:NA
    B = A{hi}; lab = Alabs{hi}; mskhi = Amsks{hi}{1}; msklo = Amsks{hi}{2};
    Blnstyl = {'-', '-', '-'};
    Bms = {1,5,5};
    Blw = {2,1,0.5};
    cmap = Acmaps{hi}; clims = Aclims{hi};
    B0 = B;
    if contains(lab, {'SWH', 'WSP'})
        % B = cellfun(@(VAR) VAR - nanmedian(VAR, [1 2]), B, 'Un', 0);
        B = cellfun(@(VAR) remove_planefit(X,Y,VAR), B, 'Un', 0);
        B = cellfun(@(VAR) movmean(VAR,3,3,'omitnan'), B, 'Un', 0);
        
        B0 = cellfun(@(VAR) movmean(remove_planefit(X,Y,VAR),3,3, 'omitnan') + nanmedian(VAR, [1 2]), B0, 'Un', 0);
        % B = cellfun(@(VAR) movmedian(VAR,3,3,'omitnan'), B, 'Un', 0);

    % else
        % B0 = cellfun(@(VAR) remove_planefit(X,Y,VAR) + nanmedian(VAR, [1 2]), B0, 'Un', 0);
    end
    % if contains(lab, {'WSP'})
        % B{2} = movmedian(B{2},5,3,'omitnan');
        % B{2} = remove_planefit(X,Y,B{2});
        % B0{2} = B{2} + nanmedian(B0{2}, [1 2]);

    % end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi*4-3); axes(thisha); hold on;
    ylabel(lab)
    % x = t; y = ygm; xx = find(isin(xgm, [-1 1]*5));
    % Z = B{1}; Z = nanmedian(Z(xx,:,:),1); Z = squeeze(Z);
    % pcolor(x,y,Z); shading flat;

    for bi=length(B):-1:1
        x = t; Z = B0{bi}; Za = B{bi};
        if allnan(Z); continue; end
        col = scaledcmap(cmap,nanmedian(nanmedian(Za.*mskhi,3), [1 2]), clims);
        if bi==3; col = col*0.2; end
        y = Z.*mskhi; 
        lims = prctile(y, [10 90], [1 2]); 
        y = nanmedian(y,[1 2]); y = squeeze(y); lims = squeeze(lims);
        if bi==1
            shaded_limits(x,y,lims, col)
        else
            plot(x,y,'k.', 'Color', col, 'LineWidth',Blw{bi}, 'LineStyle', Blnstyl{bi}, 'MarkerSize',Bms{bi})
        end


        col = scaledcmap(cmap,nanmedian(nanmedian(Za.*msklo,3), [1 2]), clims);
        if bi==3; col = col*0.2; end
        y = Z.*msklo; 
        lims = prctile(y, [10 90], [1 2]); 
        y = nanmedian(y,[1 2]); y = squeeze(y);
        plot(x,y,'k.', 'Color', col, 'LineWidth',Blw{bi}, 'LineStyle', Blnstyl{bi}, 'MarkerSize',Bms{bi})
        if bi==1
            shaded_limits(x,y,lims, col)
        else
            plot(x,y,'k.', 'Color', col, 'LineWidth',Blw{bi}, 'LineStyle', Blnstyl{bi}, 'MarkerSize',Bms{bi})
        end
    end

    for  bi=1:length(B)
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi*4-(3-bi)); axes(thisha); hold on; 
        VAR = B{bi};
        Z = nanmedian(VAR,3);
        pcolor(X,Y,Z); shading flat;
        colormap(gca,cmap); caxis(clims);

    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi*4); axes(thisha); hold on; 
    c = fixedcolorbar(gca); c.Position(1) = sum(thisha.Position([1 3])) + 0.01;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % thisha = ha(hi*4-1); axes(thisha); hold on; 
    % VAR = B{2};
    % Z = nanmedian(VAR,3);
    % pcolor(X,Y,Z); shading flat;
    % colormap(gca,cmap); caxis(clims);

end

% -------------------------------------------------------------------------
for hi=setxor(1:length(ha), [1:4:NA*4])
    axes(ha(hi)); 
    rlims = [0:25:105];
    thtx = 0:5:360;
    [xtx ytx] = pol2cart(deg2rad(thtx), rlims');
    plot(xtx',ytx','k-', 'Color', [0 0 0 0.08])
end


%%
t0s = [datenum(2023,4,15.5), datenum(2023,4,20) datenum(2023,4,25.5) datenum(2023,5,11)];
A = {WSP, SWH,mag(U,V)};
NA = length(A); Nt = length(t0s);
Acmaps = {flipud(cmocean('speed')) nicejet flipud(cmocean('deep'))}
Aclims = {[-1 1]*2 [-1 1]*0.8 [0 1.5]};
% -------------------------------------------------------------------------figure(5); clf;
figure(142); clf;

setfigsize(gcf, [480   452])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(NA, Nt, [0.03 0.02], [0.09 0.05], [0.05 0.05]);

% setfigsize(gcf, [500 600])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
% ha = tight_subplot(NA, Nt, [0.015 0.015], [0.05 0.05], [0.05 0.05]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',10);
setaxes(ha, 'DataAspectRatio', [1 1 1]); 
setaxes(ha, 'XLim', [-1 1]*110); 
setaxes(ha, 'YLim', [-1 1]*110); 
setaxes(ha, 'XLim', [-1 1]*125); 
setaxes(ha, 'YLim', [-1 1]*150); 
setaxes(ha, 'YTick', [-100 0 100]); 
setaxes(ha, 'YTick', [-100 0 100]); 
% linkaxes(ha)
% format_mapticks(ha, 'xy');
quietaxes(ha(1:end-Nt), 'x')
quietaxes(ha(setxor(1:Nt:length(ha), 1:length(ha))), 'y')
% -------------------------------------------------------------------------
for ai=1:length(A)
% -------------------------------------------------------------------------
    VAR = A{ai};
    cmap = Acmaps{ai}; clims = Aclims{ai};
    VARa = VAR;
    if ai<3
        VARa = VARa - nanmedian(VARa,  [1 2]);
    end
    % VARa = remove_planefit(X,Y,VARa);
    % VARa = VARa - nanmedian()
    for hi = 1:Nt
        t0 = t0s(hi);
        [~, ti] = min(abs(t-t0));
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(Nt*ai - (Nt-hi)); axes(thisha); hold on; 
        Z = VARa(:,:,ti);
        pcolor(X,Y,Z); shading flat;
        colormap(gca, cmap);
        caxis(clims)
    
    end

end


