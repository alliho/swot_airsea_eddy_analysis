% IDENTIFYING EDDY
% USING NEUROST

[basepath repopath swotpath datapath] = checksystem();


run([repopath 'code/startup_append.m'])
savepath = [repopath 'figs/']
saveopt = 'time'

%% load NEUROST

if ~exist('NEUROST'); load([datapath 'NEUROST.mat']); end
if ~exist('SWOT'); load([datapath 'SWOT.mat']); end


%% initial guess (manual)
[X Y] = ndgrid(NEUROST.lon, NEUROST.lat);
t = NEUROST.t0; 
VAR = NEUROST.adt;
[~, fcor] = calcinertialfreq(Y);
VAR = - NEUROST.zeta./fcor;

clear eddy; 
eddy.t = t; 

figure(80123); clf; hold on; 

eddy.centers_manual = NaN(2,length(t));

for ni=1:length(t)
    t0 = t(ni);
    Z = VAR(:,:,ni); 
    pcolor(X,Y,Z); shading flat; 
    [dt, si] = min(abs([SWOT.t0] - t0));
    swot = SWOT(si); 
    xlim(-50+ [-1 1]*5); ylim(36 + [-1 1]*5); set(gca, 'DataAspectRatio', [1 1 1]); 
    caxis([-1 1]); colormap(niceredblue)


    plot(eddy.centers_manual(1,:), eddy.centers_manual(2,:), 'k.-')
    [xg yg] = ginput(1);
    if ~isempty(xg)
        eddy.centers_manual(:,ni) = [xg yg];
        plot(xg,yg, 'ko');
    end
    pause(0.2);
    cla;
end


%% refine
Rmax = 250; % maximum radius
niter = 2; % number iterations for fitting
sigmamax = 150 ; % maximum sigma for gauss fit

% -------------------------------------------------------------------------
x = NEUROST.lon; y = NEUROST.lat;
[X Y] = ndgrid(x,y);
% field used to identify eddy; either SSHA or vorticty
VAR = NEUROST.adt; Umag = mag(NEUROST.ugosa, NEUROST.vgosa);
[~, fcor] = calcinertialfreq(Y);
VAR = - NEUROST.zeta./fcor; 


%%% MAKE COURSER GRID!
dxy = 0.05; x = x(1):dxy:x(end); y = y(1):dxy:y(end);
nSSHa = NaN(length(x), length(y), size(VAR,3)); 
nUmag = NaN(length(x), length(y), size(Umag,3)); 
[nX nY] = ndgrid(x,y);
for ti=1:size(VAR,3)
    nSSHa(:,:,ti) = interp2(X',Y',VAR(:,:,ti)', nX',nY')';
    nUmag(:,:,ti) = interp2(X',Y',Umag(:,:,ti)', nX',nY')';
end
[X Y] = ndgrid(x,y);
VAR = nSSHa; 
Umag = nUmag;
% -------------------------------------------------------------------------
% Define 2D isotropic Gaussian function
% gauss2d = @(p,xy) p(1) * exp( -((xy(:,1)-p(2)).^2 + (xy(:,2)-p(3)).^2) / (2*p(4)^2) ) + p(5);
gauss2d = @(p,xy) p(1) * exp( ...
    -( ((xy(:,1)-p(2))*cos(p(6)) + (xy(:,2)-p(3))*sin(p(6))).^2 ...
    /(2*p(4)^2) ...
    + (-(xy(:,1)-p(2))*sin(p(6)) + (xy(:,2)-p(3))*cos(p(6))).^2 ...
    /(2*(p(4)*p(5))^2) ) ) + p(7);
% p(1) = amplitude
% p(2) = x center
% p(3) = y center
% p(4) = sigma_x
% p(5) = sigma_y / sigma_x
% p(6) = rotation angle (radians)
% p(7) = offset

% -------------------------------------------------------------------------
centers_firstguess = eddy.centers_manual;
di = 4; 
centers_firstguess = movmean(movmedian(centers_firstguess(:,:),di,2, 'omitnan'),2,di-2);
% -------------------------------------------------------------------------

figure(203); clf;
setfigsize(gcf, [886   420])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(1, 3, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
% [ha] = reorg_tightsubplot(oldha, {[1],[2:3]});
% -------------------------------------------------------------------------
ha(2).Position(1) = ha(2).Position(1) - 0.03;
% -------------------------------------------------------------------------
format_fig4print(ha);
setaxes(ha, 'DataAspectRatio', [1 1 1]); 
setaxes(ha(1:2), 'XLim', nanmean(eddy.centers_manual(1,:))+ [-1 1]*2); 
setaxes(ha(1:2), 'YLim', nanmean(eddy.centers_manual(2,:)) + [-1 1]*2.5); 
linkaxes(ha(1:2))
format_mapticks(ha(1:2), 'xy');
setaxes(ha(end), 'XLim', 0.5.*Rmax*[-1 1]*2); 
setaxes(ha(end), 'YLim', 0.5.*Rmax*[-1 1]*2.5); 
setaxes(ha(:), 'CLim', [-1 1]); 
setaxes(ha(:), 'Colormap', niceredblue); 
quietaxes(ha(2), 'y')

% -------------------------------------------------------------------------
gcol = [0.1 0.5 0.17];
gccol = [0.1 0.3 0.6];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(1); axes(thisha); hold on; 
plot(eddy.centers_manual(1,:), eddy.centers_manual(2,:), 'k.', 'Color', [1 1 1].*0.6);
plot(centers_firstguess(1,:), centers_firstguess(2,:), 'k-', 'Color', [1 1 1].*0.5, 'LineWIdth',1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eddy.centers_gaussfit = NaN(size(centers_firstguess));
eddy.centers_minima = NaN(size(centers_firstguess));
eddy.pfit = NaN(7,length(centers_firstguess));
eddy.pfit_description = 'amplitude, x0 center, y0 center, sigma, sigma_y/sigma_x, offset, rotation'
eddy.pfit_err = NaN(1, length(centers_firstguess));
eddy.pfit_cc = NaN(1, length(centers_firstguess));
xgm = [-Rmax:2:Rmax]; ygm = [-Rmax:2:Rmax];
[Xgm Ygm] = ndgrid(xgm,ygm);
Zgm = NaN([size(Xgm) length(SWOT)]);
[~, fcor] = calcinertialfreq(Y);

% pfit_init = []

for ni=1:length(eddy.t)
    if allnan(centers_firstguess(:,ni)); continue; end
    %%
    Z = VAR(:,:,ni);
    dZ = Umag(:,:,ni);

    [arc, az] = distance(centers_firstguess(2,ni),centers_firstguess(1,ni),Y,X, referenceEllipsoid('wgs84'));
    arc = arc./1000; 
    az = oc2pol((az));
    [Xg Yg] = pol2cart(deg2rad(az),arc); 


    Z0 = gaussfilt_2d(Z, mode(diff(Yg(1,:))), 300);
    % Z0 = nanmedian(Z0,[1 2]);
    % Z0(mag(Xg,Yg)>Rmax) = Z(mag(Xg,Yg)>Rmax); 
    Z = Z - Z0;

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(2); axes(thisha); cla; hold on; 
    pcolor(X,Y,Z); shading flat;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Z = Z - nanmean(Z(isin(mag(Xg,Yg), [100 200])), [1 2]);
    Z(mag(Xg,Yg) > Rmax) = NaN;
    Z = remove_planefit(X,Y,Z); 

    [ix iy] = find(min(Z, [], [1 2])==Z);

    eddy.centers_minima(:,ni) = [X(ix,iy), Y(ix,iy)];
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(3); axes(thisha); cla; hold on; xline(0); yline(0); 
    pcolor(Xg,Yg,Z); shading flat; 
    
    %%% FIT THE GAUSSIAN
    % Flatten the grids into vectors
    % nn = find(~isnan(Z) & isin(abs(Xg), [0 Rmax]) & isin(abs(Yg), [0 Rmax]) );
    nn = find(~isnan(Z));
    x = Xg(nn);
    y = Yg(nn);
    z = Z(nn); z = z - nanmean(z(:));
    
    

    % Initial guesses (based on amplitude at center)
    ii = find(isin(abs(x),[-10 10]) & isin(abs(y),[-10 10])); 
    amp0 = max(z(:)) - median(z(:)); amp0 = amp0.*sign(nanmean(z(ii)));
    x0 = 0; y0   = 0;
    sigma0 = 50;    % in km, adjust based on your data
    if ni>1; sigma = eddy.pfit(4,ni-1); end
    offset0 = min(z(:));
    
    % Fit using nonlinear least squares (ratio with rotation version)
    % Parameters:
    % p = amplitude, x0 center, y0 center, sigma, sigma_y/sigma_x, offset, rotation
    %             1,         2,         3,       4,       5,      6, 7
    xy = [x y];
    opts = optimoptions('lsqcurvefit','Display','iter','MaxFunctionEvaluations',1e5);
    p0 = [amp0,   x0,   y0, sigma0,    1, offset0, 0];
    lb = [ -5,   -25,  -25,      5, 0.75,      -5, -pi];   % lower bounds
    ub = [  5,    25,   25,     sigmamax,    4,       5, pi]; % upper bounds
    pfit = lsqcurvefit(gauss2d, p0, xy, z, lb, ub, opts);
    
    % OPTIONAL: ITERATE
    % for nii=1:niter
    %     pfit_init = pfit; 
    %     p0 = pfit_init;
    %     pfit = lsqcurvefit(gauss2d, p0, xy, z, lb, ub, opts);
    % end



    eddy.pfit(:,ni) = pfit; 
    title(['($\sigma=$' num2str(round(pfit(4))) ', ratio=' num2str(round(pfit(5),2)) ', rotation=' num2str(round(rad2deg(pfit(6)),1))  ')'],'interpreter', 'latex')

    % Reconstruct fitted surface on the original grid using the best pfit
    xygrid = [Xg(:) Yg(:)];
    Zfit = nan(size(Xg)); % keep NaNs where original Z was NaN
    Zfit(:) = gauss2d(pfit, xygrid);
    Zfit0 = nan(size(Xg)); % keep NaNs where original Z was NaN
    Zfit0(:) = gauss2d(p0, xygrid); % for reference
    
    z_ob = Z; z_fit = Zfit; z_fit(isin(z_fit, pfit(7) + [-1 1]*0.05)) = NaN;
    eddy.pfit_err(ni) = nanmean(abs(z_ob-z_fit), [1 2]); 
    eddy.pfit_cc(ni) = ccfxn(z_ob,z_fit);

    %%% PLOT %%%
    contour(Xg,Yg, Zfit,[-1:0.05:1], 'LineWidth',2)
    contour(Xg,Yg,Zfit0, [-1:0.5:1], 'k:')

    
    contour(Xg,Yg,Zfit, [-1 1].*0.01, 'g-', 'EdgeAlpha',0.4, 'Color', gcol, 'LineWidth',2)
    
    plot(pfit(2), pfit(3), 'rx', 'LineWidth',2, 'Color', gccol);
    % Compute distance from every grid point
    dist = hypot(Xg - pfit(2), Yg - pfit(3));
    [~, idx] = min(dist(:));
    % Convert linear index to row/col subscripts
    [irow, icol] = ind2sub(size(Xg), idx);
    
    
    

    eddy.centers_gaussfit(:,ni) = [X(irow,icol), Y(irow,icol)];



    for hi=1:2
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi); axes(thisha); hold on; 
        contour(X,Y,Zfit, [-1 1].*0.01, 'b-', 'EdgeAlpha',0.4, 'EdgeColor', gcol, 'LineWidth',2)
        plot(centers_firstguess(1,ni), centers_firstguess(2,ni) , 'ro', 'LineWidth',1)
        plot(eddy.centers_gaussfit(1,ni), eddy.centers_gaussfit(2,ni), 'bo', 'Color', gccol, 'LineWidth',1);
        plot(eddy.centers_minima(1,ni), eddy.centers_minima(2,ni), 'b.', 'Color', [1 0.6 0.3], 'LineWidth',1);

        if hi==1
            msk = ones(size(X)); 
            msk(mag(Xg,Yg)>100) = NaN;
            contour(X,Y,Z.*msk, [-1 1].*0.01, 'b-', 'EdgeAlpha',0.4, 'EdgeColor', 'k', 'LineWidth',2)
        end
    end

    
    pause(0.01); 
    % keyboard

    
end

eddy.gaussfxn = gauss2d;

%% plot to check 
figure(45); clf; 

subplot(1,3,1); cla; hold on; set(gca, 'DataAspectRatio', [1 1 1])
plot(eddy.centers_manual(1,:), eddy.centers_manual(2,:), 'k.-')
plot(eddy.centers_minima(1,:), eddy.centers_minima(2,:), 'r.-')
plot(eddy.centers_gaussfit(1,:), eddy.centers_gaussfit(2,:), 'b.-', 'Color', gccol, 'LineWidth',2)
% plot(eddy.centers_gaussfit(1,:), eddy.centers_gaussfit(2,:), 'bo-', 'Color', gccol, 'LineWidth',2)

t = eddy.t; sigma = nanmean([eddy.pfit(4,:); eddy.pfit(4,:).*eddy.pfit(5,:)],1); amp = eddy.pfit(1,:);
rat = eddy.pfit(5,:); rat(rat>1) = 1./rat(rat>1);
clear pi;
eddy.area = 2*pi*sigma.^2;
eddy.dt = gradient(t)*24*3600;


eddy.ratio = rat; 
eddy.sigma = sigma; 
eddy.sigma_x = eddy.pfit(4,:);
eddy.sigma_x = eddy.pfit(4,:).*eddy.pfit(5,:);
eddy.amp = amp; 
eddy.offset = eddy.pfit(7,:);
eddy.rotation = eddy.pfit(6,:);
eddy.time = eddy.t; 

loc_err_x = nanstd([eddy.centers_manual(1,:); eddy.centers_gaussfit(1,:);]);
loc_err_y = nanstd([eddy.centers_manual(2,:); eddy.centers_gaussfit(2,:); ]);
eddy.loc_err = deg2km(nanmean([loc_err_x; loc_err_y]));


nn = find(eddy.loc_err > 20 & eddy.sigma > 100);

lats = eddy.centers_gaussfit(2,:); lons = eddy.centers_gaussfit(1,:); 
lons(nn) = NaN; lats(nn) = NaN;
lons = movmean(lons,3, 'omitnan');  lats = movmean(lats,3, 'omitnan');
eddy.lat0 = lats; eddy.lon0 = lons; 


[arc, az] = distance(lats(1:end-1), lons(1:end-1), lats(2:end), lons(2:end),referenceEllipsoid('wgs84'));
lons(arc>1000*20) = NaN; lats(arc>1000*20) = NaN; arc(arc>20*1000) = NaN;


az = oc2pol((az));
[dx dy] = pol2cart(deg2rad(az),arc); 
dx = movmean(dx,3, 'omitnan');  dy = movmean(dy,3, 'omitnan');

eddy.dx = [NaN dx];
eddy.dy = [NaN dy];



plot(lons, lats, 'c-', 'linewidth',2)


subplot(4,3,[2:3]); cla; hold on; 
plot(t, rat, 'k.-');
autodatetick(gca, 'x')
ylabel('ratio')

subplot(4,3,[2:3]+3); cla; hold on; 
plot(t, sigma, 'k.-'); ylim([0 100])
autodatetick(gca, 'x')
ylabel('sigma')

subplot(4,3,[2:3]+6); cla; hold on; 
plot(t, amp, 'k.-');
autodatetick(gca, 'x')
ylabel('amp')

subplot(4,3,[2:3]+9); cla; hold on; 
% plot(t, abs(eddy.pfit_err), 'k.-');
plot(t, mag(eddy.dx./eddy.dt, eddy.dy./eddy.dt)*(3600*24./1000), 'k.-');
ylim([0 10])



% plot(t, eddy.loc_err  , 'k.-');

autodatetick(gca, 'x')

xline([datenum(2023,4,15) datenum(2023,5,12)])


%%

save([datapath 'eddy_zeta.mat'], 'eddy')