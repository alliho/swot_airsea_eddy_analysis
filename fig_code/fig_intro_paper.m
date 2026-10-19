%% setup
[basepath repopath swotpath datapath] = checksystem();


run([repopath 'code/startup_append.m'])
savepath = [repopath 'figs/']
saveopt = 'time'

%% load
if ~exist('SWOT')
    load([datapath 'SWOT_proc.mat'])
end

swotcube = cube_swot_by_pass(SWOT);
if ~exist('ERA5')
    load([datapath 'ERA5_grid.mat'])
end

if ~exist('DUACS')
    load([datapath 'DUACS.mat'])
end
if ~exist('NEUROST')
    load([datapath 'NEUROST.mat'])
end
if ~exist('ASCAT')
    load([datapath 'ASCAT.mat'])
end
if ~exist('ghrsst_l3')
    load([datapath 'ghrsst_l3.mat'])
end

if ~exist('eddy')
    load([datapath 'eddy.mat'])
end

try; 
    t0 = datenum(2023,6,28);
    longSWOT = load_swot([erase(datapath, 'processed/eddy/') 'SWOT/onedayrepeat_GSextension/'], lonlims + [-1 1].*90, [-1 1].*90, t0 + [-1 1].*0.5, 'search', '*Expert*.nc');

catch
    longSWOT = SWOT(1:2);
end

load([datapath 'eddy_zeta.mat']); rscale = 0.85; rscale_plot = 1.5; 

%%
t0 = datenum(2023,4,16.5);
% t0 = datenum(2023,4,14.5);

[~, si] = min(abs(t0 - [SWOT.t0])); t0 = SWOT(si).t0;
datestr(t0);
bkgd_source = 'NEUROST';
swocol = [0.25 0.29 0.6]+0.1;

% -------------------------------------------------------------------------

figure(124123); clf; 
setfigsize(gcf, round([930   589].*1.2)); 
frat = 589/930;
gap = 0.015;
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(2, 3, [gap*2 gap], [0.05 0.08], [0.05 0.05]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
dy = 0.12; 
for hi=1:3
    ha(hi).Position([2 4]) =     ha(hi).Position([2 4]) + [1 -1].*dy;
    ha(hi).Position([2 4]) =     ha(hi).Position([2 4]) + [0.03 0];
    ha(hi+3).Position([2 4]) =     ha(hi+3).Position([2 4]) + [0 1].*dy;
end
xrng = sum(cellfun(@(x) x(3), {ha(1:3).Position}))+gap*2;
ha(1).Position(3) = ha(1).Position(4).*frat;
ha(2).Position(1) = sum(ha(1).Position([1 3])) + gap;
ha(2).Position(3) = ha(1).Position(3)/1.4;
ha(3).Position(1) = sum(ha(2).Position([1 3])) + gap;
ha(3).Position(3) = xrng - ha(3).Position(1) + ha(1).Position(1);
ha(2).Position(1) = ha(2).Position(1) + gap.*1;


addha = axes(); 
addha.Position = ha(3).Position;
ha = [ha(1:3); addha; ha(4:end)];
ha(3).Position(4) = (ha(3).Position(4)-gap)/2;
ha(4).Position(4) = ha(3).Position(4);
ha(4).Position(2) = sum(ha(3).Position([2 4])) + gap;
dx = gap*3.5;
for hi=3:4
    ha(hi).Position([1 3]) = ha(hi).Position([1 3]) + [1 -1]*dx;
end

% -------------------------------------------------------------------------
format_fig4print(ha, 'FontSize',11);
setaxes(ha(1), 'Colormap', niceredblue); 

setaxes(ha(2), 'FontSize', get(ha(1), 'FontSize')-1); 
setaxes(ha(2), 'DataAspectRatio', [1 1 1]); 
setaxes(ha(2), 'CLim', [4 10]); 
setaxes(ha(2), 'XLim', -0.4 + -48.5 + [-1 1]*2.88.*0.5);  
setaxes(ha(2), 'YLim', -0.7 + 37.25 + [-1 1]*4.*0.5); 
lonlims = get(ha(2), 'XLim'); latlims = get(ha(2), 'YLim');
setaxes(ha(2), 'Colormap', cmocean('ice')); 



setaxes(ha(3:end), 'XLim', [datenum(2023,4,2) datenum(2023,5,30)]); 
setaxes(ha(3:end), 'FontSize', get(ha(1), 'FontSize')-1); 
autodatetick(ha(3:4), 'x', 'dt',7); 
quietaxes(ha(4), 'x')

setaxes(ha(5:end), 'DataAspectRatio', [1 1 1]); ; 
setaxes(ha(5:end), 'XLim', [-55 -42] + [-1 1].*-3.2); 
setaxes(ha(5:end), 'YLim', [30.5 45] + [-1 1].*-3.2); 

setaxes(ha(5:end), 'XLim', -48.5 + [-1 1]*3.4.*0.6); 
setaxes(ha(5:end), 'YLim', 37.25 + [-1 1]*4.*0.6); 

setaxes(ha([2 5:end]), 'YTick', [30:50]); 
setaxes(ha([2 5:end]), 'XTick', [-60:-40]); 



% setaxes(ha(1), 'Colormap', cmocean('thermo')); 
setaxes(ha(5), 'Colormap', flipud(cmocean('deep'))); 
setaxes(ha(5), 'CLim', [0 2]); 
setaxes(ha(6), 'Colormap', flipud(cmocean('speed'))); 
setaxes(ha(6), 'CLim', [2 14]); 
setaxes(ha(7), 'Colormap', nicejet); 
setaxes(ha(7), 'CLim', [0.66 5.75]); 
% setaxes(ha(5), 'Colormap', cmocean('thermo')); 
% setaxes(ha(5), 'CLim', [10 20]); 

format_mapticks(ha(5:end), 'xy')
format_mapticks(ha(2), 'xy')

quietaxes(ha(6:end), 'y')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
m_proj('ortho','lat',25,'long',nanmean(lonlims)'-15);

if strcmp(bkgd_source, 'DUACS');
    [~, ti] = min(abs(DUACS.t0 - t0))
    m_pcolor(DUACS.lon, DUACS.lat, DUACS.sla(:,:,ti)'); shading flat;
    m_pcolor(swot.longitude-360, swot.latitude, swot.ssha_karin_2_proc_patch+ swot.height_cor_xover); shading flat;
elseif strcmp(bkgd_source, 'NEUROST')
    [~, ti] = min(abs(NEUROST.t0 - t0))
    m_pcolor(NEUROST.lon, NEUROST.lat, NEUROST.adt(:,:,ti)'); shading flat;
    m_pcolor(NEUROST.lon, NEUROST.lat, NEUROST.sla(:,:,ti)'); shading flat; caxis([-1 1]*0.6)

end



for si=1:length(longSWOT)
    swotlong = longSWOT(si); 
    xi = 31;
    X = swotlong.longitude; Y = swotlong.latitude; X(X>0) = X(X>0)-360;
    X = X([5 xi],:); Y = Y([1 xi],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
    m_patch(X,Y,swocol+0.1, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2)
    xi = 69-xi + 2;
    X = swotlong.longitude; Y = swotlong.latitude;  X(X>0) = X(X>0)-360;
    X = X([xi end-3],:); Y = Y([xi end],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
    m_patch(X,Y,swocol+0.1, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2)
end




[xpts,ypts] = drawbox(get(ha(2), 'XLim'),get(ha(2), 'YLim'));
col = [0.7 0 0].*0.2;
m_patch(xpts, ypts, 'w', 'EdgeColor', col, 'LineWidth',0.8, 'FaceColor', 'none');
m_grid('linest','-','xticklabels',[],'yticklabels',[]);

m_coast('patch',[1 1 1].*0.7, 'edgecolor', 'none');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(2); axes(thisha); hold on; 
styleargin = {'FaceAlpha', 0.1, 'EdgeAlpha', 0.6, 'LineWidth',0.5};

for si=1:2
    swot = SWOT(si);
    xi = 31;
    X = swot.longitude - 360; Y = swot.latitude; %X(X>360) = X(X>360)-360;
    X = X([5 xi],:); Y = Y([1 xi],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
    % patch(X,Y,swocol+0.1, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)
    patch(X,Y,swocol+0.1, styleargin{:})
    xi = 69-xi + 2;
    X = swot.longitude - 360; Y = swot.latitude;  %X(X>360) = X(X>360)-360;
    X = X([xi end-3],:); Y = Y([xi end],:); X = [X(1,:) fliplr(X(end,:))]; Y = [Y(1,:) fliplr(Y(end,:))]; 
    % patch(X,Y,swocol+0.1, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)
    patch(X,Y,swocol+0.1, styleargin{:})
end


x = eddy.lon0; y = eddy.lat0; u = eddy.dx./eddy.dt; v = eddy.dy./eddy.dt;
z =  mag(u,v); z = z.*(3600*24)./1000;
sc =3;
y(eddy.sigma>120) = NaN;
plot(x,y,'k:', 'LineWidth',2, 'Color', [1 1 1].*0.6)



c = fixedcolorbar(gca, 'Location', 'northoutside', 'FontSize', get(gca, 'FontSize')-4);
c.Position(4) = c.Position(4)/2.4;
c.Position(3) = c.Position(3)/6;
c.Position(2) = sum(thisha.Position([2]));
c.Position(1) = thisha.Position(1);
c.Position([1 2]) = c.Position([1 2]) + thisha.Position([3 4]).*0.03;

c.AxisLocation = 'in';
textbypos(c.Position(1) + c.Position(3) + 0.005, sum(c.Position([2 4]).*[1 0.5])-0.002,  '[km/day]', 'Interpreter', 'latex', 'FontSize',get(c, 'FontSize')+1, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')

clims = get(gca, 'CLim');
cmap = scaledcmap(get(gca, 'Colormap'), movmean(z,3, 'omitnan'), clims);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -]
[Xc Yc] = ndgrid(NEUROST.lon, NEUROST.lat);
cmap = scaledcmap(get(gca, 'Colormap'), mag(eddy.dx./eddy.dt, eddy.dy./eddy.dt).*(24*3600./1000), clims);

di = 2;
ee = find(isin(eddy.t, get(ha(3), 'XLim'), 'inclusive')); ee = ee(1):di:ee(end);

err = nanstd(cell2matn({eddy.centers_gaussfit, eddy.centers_manual},3), [],3); 
err = mag(err(1,:),err(2,:));
err(movmean(eddy.ratio,3)<0.5) = NaN;

err = abs(eddy.pfit_err - movmedian(eddy.pfit_err,5));


for ei=ee
    x0 = eddy.centers_gaussfit(1,ei); y0 = eddy.centers_gaussfit(2,ei);
    pfit = eddy.pfit(:,ei);
    Rmax = pfit(4);
    pfit(4) = pfit(4)*rscale;
    [arc, az] = distance(y0,x0,Yc,Xc, referenceEllipsoid('wgs84'));
    arc = arc./1000; 
    az = oc2pol((az));
    [Xg Yg] = pol2cart(deg2rad(az),arc); 

    [~, ni] = min(abs(eddy.t(ei) - NEUROST.t0));
    Zc = NEUROST.sla(:,:,ni); 
    % [~, fcor] = calcinertialfreq(Yc);
    % Zc = NEUROST.zeta(:,:,ni)./fcor; 
    Rc = deg2km(mag(Xc - eddy.lon0(ei), Yc - eddy.lat0(ei)));   

    xygrid = [Xg(:) Yg(:)];
    Zfit = nan(size(Xg)); % keep NaNs where original Z was NaN
    Zfit(:) = eddy.gaussfxn(pfit, xygrid);

    % Zc(Rc>pfit(4)*10) = NaN;
    Zc(Rc>Rmax*5) = NaN;
    % Zc(Zfit>0.1) = NaN;
    % Zc()
    % [~, ti] = min(abs(eddygrid.time - eddy.t(ei)));
    if isnan(cmap(ei,:)); continue; end
    if err(ei) < 0.01 & ~isnan(err(ei))
        contour(Xc,Yc, Zc, [1 1].*(-0.1), 'k-', 'EdgeColor', [1 1 1]*0.6, 'Color', cmap(ei,:), 'LineWidth',1.4);
        if abs(eddy.time(ei)-(t0+1))<1; 
            tmp = contour(Xc,Yc, Zc, [1 1].*(-0.1), 'k:', 'EdgeColor', [1 1 1]*0.6, 'Color', 'r', 'LineWidth',2);
        end
    else
        contour(Xc,Yc, Zc, [1 1].*(0.02), 'k:', 'EdgeColor', [1 1 1]*0.6, 'LineWidth',1);
    end


end


u = eddy.dx./eddy.dt; v = eddy.dy./eddy.dt;
x = eddy.lon0  ; y = eddy.lat0;
y(eddy.sigma>120) = NaN;
y(err>0.2 & eddy.ratio>0.5) = NaN;
sc = 3.5;
ee = 1:5:length(y)
quiver(x(ee), y(ee), sc.*u(ee), sc.*v(ee), 'k-', 'Color', 'k', 'AutoScale', 'off', 'linewidth',1, 'MaxHeadSize',1)



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(3); axes(thisha); cla; hold on; 
ylabel(['[km/day]']); ylim([0 17])
u = eddy.dx./eddy.dt; v = eddy.dy./eddy.dt;
x = eddy.time; y = mag(u,v).*(3600*24)./1000; 
err = abs(eddy.pfit_err - movmedian(eddy.pfit_err,5));
y(err>0.01) = NaN; %y(isnan(err)) = NaN;
plot(x,y,'k.-', 'MarkerSize',10, 'LineWidth',1)
xline(t0, 'r:', 'LineWidth',2)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(4); axes(thisha); cla; hold on; 
ylabel(['$R_\sigma$ [km]'], 'interpreter', 'latex'); ylim([25 76])


x = eddy.time; y = eddy.sigma; y(isnan(eddy.lon0))  = NaN;
y = y.*rscale_plot;
err = abs(eddy.pfit_err - movmedian(eddy.pfit_err,5));
y(err>0.01) = NaN; %y(isnan(err)) = NaN;
plot(x,y, 'k.-', 'MarkerSize',10, 'LineWidth',1)
xline(t0, 'r:', 'LineWidth',2)



    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(5); axes(thisha); hold on; 
       
    % textbypos(sum(thisha.Position([1 3]).*[1 0.75]), sum(thisha.Position([2 4]).*[1 1.05]),  {datestr(t0, 'mm/dd/yyyy HH:MM:SS'), []}, 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize',13)
    
    if strcmp(bkgd_source, 'DUACS');
        [~, di] = min(abs([DUACS.t0] - t0));
        x = DUACS.lon; y = DUACS.lat; t = DUACS.t0;
        U = DUACS.ugos; V = DUACS.vgos; 
    elseif strcmp(bkgd_source, 'NEUROST')
        [~, di] = min(abs([NEUROST.t0] - t0));
        x = NEUROST.lon; y = NEUROST.lat; t = NEUROST.t0;
        U = NEUROST.ugos; V = NEUROST.vgos; 
    end
    
    [X Y] = ndgrid(x,y); 
    U = U(:,:,di); V = V(:,:,di);
    Z = mag(U,V); 
    pcolor(X,Y,Z); shading interp;
    sc = 0.1;
    di = 0.25./mode(diff(x)); di = round(di);
    iix = 1:di:size(X,1); iiy = 1:di:size(X,2);
    col = [1 1 1 0.1]; col = [0 0 0];
    quiver(X(iix,iiy),Y(iix,iiy),sc.*U(iix,iiy),sc.*V(iix,iiy),'k-','AutoScale','off', 'MaxHeadSize',1, 'Color', col)
    

    [~, si] = min(abs([SWOT.t0]-t0));
    X = SWOT(si).longitude-360; Y = SWOT(si).latitude;
    U = SWOT(si).U; V = SWOT(si).V; 
    Z = mag(U,V); 
    xx = find(~isin(abs(squeeze(nanmean(SWOT(si).cross_track_distance, [2 3]))./1000), [8 60])); Z(xx,:) = NaN;
    pcolor(X,Y,Z); shading interp;
    % di = deg2km(0.25)./2; di = floor(di); 
    di = di * 2; 
    iix = 1:di:size(X,1); iiy = 1:di:size(X,2);
    quiver(X(iix,iiy),Y(iix,iiy),sc.*U(iix,iiy),sc.*V(iix,iiy),'k-','AutoScale','off', 'MaxHeadSize',1, 'Color', col)

    xx = find(sum(isnan(Z),2)<size(Z,2)); X = X(xx,:); Y = Y(xx,:);
    plot(X(1,:), Y(1,:), 'k-'); plot(X(:,1), Y(:,1), 'k-')
    plot(X(end,:), Y(end,:), 'k-'); plot(X(:,end), Y(:,end), 'k-')
    mi = find(diff(xx)>1);
    plot(X(mi,:), Y(mi,:), 'k-'); plot(X(mi+1,:), Y(mi+1,:), 'k-');


    % c = cbarfun_SW(gca, '$|U_{geo}|$ [m/s]', 'Color', 'w');
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(6); axes(thisha); hold on; 
    
    % textbypos(sum(thisha.Position([1 3]).*[1 0.75]), sum(thisha.Position([2 4]).*[1 1.05]),  {datestr(t0, 'mm/dd/yyyy HH:MM:SS'), []}, 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize',13)

    [~, di] = min(abs([ERA5.time] - t0));
    x = ERA5.longitude; y = ERA5.latitude; t = ERA5.time;
    [X Y] = ndgrid(x,y); [~, fcor] = calcinertialfreq(Y);

    % Z = mag(ERA5.u10 , ERA5.v10); Z = Z(:,:,di);
    U = ERA5.u10; V = ERA5.v10; U = U(:,:,di); V = V(:,:,di);
    Z = mag(U,V); 
    pcolor(X,Y,Z); shading interp;
    sc = 0.02;
    di = 0.25./mode(diff(x)); di = round(di);
    iix = 1:di:size(X,1); iiy = 1:di:size(X,2);
    col = [1 1 1 0.1]; col = [0 0 0];
    quiver(X(iix,iiy),Y(iix,iiy),sc.*U(iix,iiy),sc.*V(iix,iiy),'k-','AutoScale','off', 'MaxHeadSize',1, 'Color', col)



    [~, ei] = min(abs(eddy.t - t0));
    elon0 = eddy.lon0(ei); elat0 = eddy.lat0(ei);
    lims = prctile(Z(isin(mag(X-elon0, Y-elat0), [0 0.5])), [1 99], [1 2])';
    caxis(0.5 + nanmedian(lims) + [-1 1].*2.8 );


    [~, si] = min(abs([SWOT.t0]-t0));
    X = SWOT(si).longitude-360; Y = SWOT(si).latitude;
    Z = SWOT(si).wind_speed_karin_2_proc_patch; 
    xx = find(~isin(abs(squeeze(nanmean(SWOT(si).cross_track_distance, [2 3]))./1000), [8 60])); Z(xx,:) = NaN;
    pcolor(X,Y,Z); shading interp;
    

    xx = find(sum(isnan(Z),2)<size(Z,2)); X = X(xx,:); Y = Y(xx,:);
    plot(X(1,:), Y(1,:), 'k-'); plot(X(:,1), Y(:,1), 'k-')
    plot(X(end,:), Y(end,:), 'k-'); plot(X(:,end), Y(:,end), 'k-')
    mi = find(diff(xx)>1);
    plot(X(mi,:), Y(mi,:), 'k-'); plot(X(mi+1,:), Y(mi+1,:), 'k-');


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(7); axes(thisha); hold on; 

    % textbypos(sum(thisha.Position([1 3]).*[1 0.75]), sum(thisha.Position([2 4]).*[1 1.05]),  {datestr(t0, 'mm/dd/yyyy HH:MM:SS'), []}, 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize',13)

    [~, di] = min(abs([ERA5.time] - t0));
    x = ERA5.wave.longitude; y = ERA5.wave.latitude; t = ERA5.time;
    [X Y] = ndgrid(x,y); [~, fcor] = calcinertialfreq(Y);

    Z = ERA5.wave.swh; 
    R = ERA5.wave.mwp; TH = ERA5.wave.mwd; TH = oc2pol(met2oc(TH));
    [U V] = pol2cart(deg2rad(TH), R);
    Z = Z(:,:,di); U = U(:,:,di); V = V(:,:,di);

    pcolor(X,Y,Z); shading interp;
    sc = 0.02;
    di = 0.25./mode(diff(x)); di = round(di);
    iix = 1:di:size(X,1); iiy = 1:di:size(X,2);
    col = [1 1 1 0.1]; col = [0 0 0];
    quiver(X(iix,iiy),Y(iix,iiy),sc.*U(iix,iiy),sc.*V(iix,iiy),'k-','AutoScale','off', 'MaxHeadSize',1, 'Color', col)


    [~, ei] = min(abs(eddy.t - t0));
    elon0 = eddy.lon0(ei); elat0 = eddy.lat0(ei);
    lims = prctile(Z(isin(mag(X-elon0, Y-elat0), [0 0.5])), [1 99], [1 2])';
    caxis(0.5 + nanmedian(lims) + [-1 1].*1.1 );

    % lims = prctile(Z(isin(X, get(gca, 'XLim')) & isin(Y, get(gca, 'YLim'))), [1 99], [1 2])';
    % caxis(lims + [-1 1].*range(lims)*0.25);


    [~, si] = min(abs([SWOT.t0]-t0));
    X = SWOT(si).longitude-360; Y = SWOT(si).latitude;
    Z = SWOT(si).swh_karin_proc_patch  ; 
    xx = find(~isin(abs(squeeze(nanmean(SWOT(si).cross_track_distance, [2 3]))./1000), [8 60])); Z(xx,:) = NaN;
    pcolor(X,Y,Z); shading interp;

    xx = find(sum(isnan(Z),2)<size(Z,2)); X = X(xx,:); Y = Y(xx,:);
    plot(X(1,:), Y(1,:), 'k-'); plot(X(:,1), Y(:,1), 'k-')
    plot(X(end,:), Y(end,:), 'k-'); plot(X(:,end), Y(:,end), 'k-')
    mi = find(diff(xx)>1);
    plot(X(mi,:), Y(mi,:), 'k-'); plot(X(mi+1,:), Y(mi+1,:), 'k-');



    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for hi=5:length(ha)

    thisha = (ha(hi)); axes(thisha);
    c = fixedcolorbar(gca, 'location', 'northoutside', 'FontSize', get(gca, 'FontSize').*0.6);
    c.Position([3 4]) = c.Position([3 4]).*[1/5.6 0.65];
    c.Position([1 2]) = thisha.Position([1 2]) + [1 1].*thisha.Position([3 4]).*0.05;
    c.AxisLocation = 'in';
    if hi==7 | hi==5 | hi==6; set(c, 'Color' , 'w'); end
    
    if hi==7; labstr = 'SWH [m]'; end
    if hi==6; labstr = 'WSP [m/s]'; end
    if hi==5; labstr = '$\vec{u}_{geo}$ [m/s]'; end
    posx = c.Position(1) + c.Position(3)/2; posy = c.Position(2) + c.Position(4)*2.8; 
    textbypos(posx, posy, labstr,...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'times', 'FontSize', get(c, 'FontSize').*1.1, 'Color', get(c, 'Color'), ...
        'Interpreter', 'latex')

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for hi=3:4
    axes(ha(hi));
    xline(datenum(2023,4,10.5), 'LineWidth',40,'Color', [0.5 0 0 0.1], 'Alpha', 0.1)
    xline([datenum(2023,4,13.45) datenum(2023,5,12.8)], 'r-', 'LineWidth',1)
end


savejpg(gcf, ['fig_intro_t0-' datestr(t0, 'yyyymmddTHHMMSS')], savepath, saveopt) 