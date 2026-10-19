


%% CORRELATION PAIRS (WIND_CURRENT WAVE_CURRENT WIND_WAVE)
Rcc = 100;
cmap = [0.15 0.55 0.85; 0.1 0.43 0.15; 0.8 0.3 0.3;  1 0.6 0.34; 0.92 0.75 0.12 ; 0.2 0.1 0.0 ];
% cmap = [0.15 0.55 0.85; 0.1 0.43 0.15;  1 0.6 0.34; 0.8 0.3 0.3; 0.92 0.75 0.12 ; 0.2 0.1 0.0 ];
% cmap = cmap.*0.8;
A = {}
Alabs = {};

rotfldstr = '';
run load_vars_into_memory_universal

daterange = [datenum(2023,4,12) datenum(2023,5,12) ]
tt = find(isin(t, daterange));
%
% % -------------------------------------------------------------------------
rotfldstr = '_rot_wv';
% rotfldstr = '';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
A = [A {SWHa; -U_wv}];
Alabs = [Alabs, {'SWHa vs $\hat{\theta}_{wv}\cdot\vec{u}_{geo}$'}];
% % -------------------------------------------------------------------------
rotfldstr = '_rot_wnd';
% rotfldstr = '';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
A = [A {WSPa; -U_wnd}];
Alabs = [Alabs, {'WSPa vs $\theta_{wnd}\cdot\vec{u}_{geo}$'}];

% % % -------------------------------------------------------------------------
% rotfldstr = '';
% run load_vars_into_memory_universal
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% A = [A {SWHa; WSPa}];
% Alabs = [Alabs, {'SWHa vs WSPa'}];
% 
% % 
% % -------------------------------------------------------------------------
rotfldstr = '_rot_wnd';
run load_vars_into_memory_universal
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
[dSST_dy dSST_dx] = gradient(SSTa_L4, drkm, drkm, 1);
[dWSP_U_dy dWSP_U_dx] = gradient(WSP_U, drkm, drkm, 1);
[dWSP_V_dy dWSP_V_dx] = gradient(WSP_V, drkm, drkm, 1);

% A = [A {dWSP_V_dx + dWSP_U_dy; dSST_dy}];
A = [A {WSPa; SSTa_L4}];
Alabs = [Alabs, {'WSPa vs SSTa'}];

% % % -------------------------------------------------------------------------
% rotfldstr = '_rot_wnd';
% run load_vars_into_memory_universal
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% msk = ones(size(WSP_V)); msk(isnan(WSP_V + WSP_U)) = NaN;
% WSP_V = gaussfilt_2d(WSP_V,drkm, 6);
% WSP_U = gaussfilt_2d(WSP_U,drkm, 6);
% % WSP_U = movmean(WSP_U,3,3, 'omitnan').*msk;
% % WSP_V = movmean(WSP_V,3,3, 'omitnan').*msk;
% [dWSP_Vdy dWSP_Vdx] = gradient(WSP_V,drkm.*1000,drkm.*1000,1);
% [dWSP_Udy dWSP_Udx] = gradient(WSP_U,drkm.*1000,drkm.*1000,1);
% dWSP_Vdx = dWSP_Vdx.*msk;
% dWSP_Udy = dWSP_Udy.*msk;
% dWSP_Udx = dWSP_Udx.*msk;
% dWSP_Vdy = dWSP_Vdy.*msk;
% % ii = find(isin(t, daterange));

% 
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % A = [A {(dWSP_Udx + dWSP_Vdy); relvort}];
% A = [A {(dWSP_Vdx - dWSP_Udy).*10.*1000; relvort}];
% Alabs = [Alabs, {'$\nabla\times WSPa$ vs $\zeta/f$'}];
% 
% 
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % A = [A {(dWSP_Udx + dWSP_Vdy); dSSTdx_L4_wnd}];
% A = [A {(dWSP_Vdx - dWSP_Udy).*10.*1000; dSSTdy_L4_wnd.*10.*1000}];
% Alabs = [Alabs, {'$\nabla\times \hat{\theta}_{wnd}\cdot WSPa$ vs $\partial SSTa/\partial y$'}];



% % -------------------------------------------------------------------------
msk = ones(size(X)); msk(mag(X,Y)<150) = NaN;
Amsk = cellfun(@(VAR) double(~isnan(VAR)), A, 'Un', 0);
for ai=1:length(Amsk); Amsk{ai}(Amsk{ai}==0) = NaN; end
A = cellfun(@(VAR) VAR - nanmedian(VAR.*msk, [1 2]), A, 'Un', 0);
% A = cellfun(@(VAR) movmedian(VAR,3, 'omitnan'), A, 'Un', 0);
% A = cellfun(@(VAR) remove_planefit(X,Y,VAR), A, 'Un', 0);
A = cellfun(@(VAR) VAR - gaussfilt_2d(VAR,drkm, 300), A, 'Un', 0);
% A = cellfun(@(VAR) movmedian(VAR,3, 'omitnan'), A, 'Un', 0);
A = cellfun(@(VAR) gaussfilt_2d(VAR,drkm, 8), A, 'Un', 0);
% A = cellfun(@(VAR) movmean(VAR,3, 'omitnan'), A, 'Un', 0);
% A = cellfun(@(VAR) VAR - gaussfilt_2d(VAR,drkm, 300), A, 'Un', 0);
% A = cellfun(@(VAR, msk) VAR.*msk, A,Amsk, 'Un', 0);


% A = cellfun(@(VAR) movmean(VAR,5, 'omitnan'), A, 'Un', 0);
% -------------------------------------------------------------------------
A_cc = cell(1, size(A,2));
A_n = cell(1, size(A,2));

%%

ii = find(isin(t, datenum(2023,4,24) + [-1 1]*1));
ii = find(isin(t, datenum(2023,5,7) + [-1 1]*1));
ii = find(isin(t, datenum(2023,4,19) + [-1 1]*2));
figure(1232); clf; 


subplot(211); VAR = A{1,3}; pcolor(X,Y,nanmedian(VAR(:,:,ii),3)); shading flat;
colorbar; caxis([-1 1]*2)
xlim([-1 1]*120)
ylim([-1 1]*120)

subplot(212); VAR = A{2,3}; pcolor(X,Y,nanmedian(VAR(:,:,ii),3)); shading flat;
colorbar; caxis([-1 1]*2)
xlim([-1 1]*120)
ylim([-1 1]*120)

%% basic one at a time.. 
Rcc = 100;
prcnanlim = 0.2;

% -------------------------------------------------------------------------

Na = length(A); 
% if Na>2; Na = Na - 1; end 
% Lavg = 15; 
% -------------------------------------------------------------------------

u10 = squeeze(nanmean(era5.u10,[1 2]));
v10 = squeeze(nanmean(era5.v10,[1 2]));
du10dt = gradient(u10, mode(diff(era5.time))*24);
dv10dt = gradient(v10, mode(diff(era5.time))*24);
du10dt = movmean(du10dt, (6/24)./mode(diff(era5.time)), 'omitnan');
dv10dt = movmean(dv10dt, (6/24)./mode(diff(era5.time)), 'omitnan');

% wsp = squeeze(mag( nanmean(era5.u10,[1 2]), nanmean(era5.v10,[1 2]))); wsp = movmean(wsp, (4/24)./mode(diff(era5.time)));
% dwspdt = gradient(wsp, mode(diff(era5.time))*24);
dwspdt = mag(du10dt, dv10dt);
dwspdt = dwspdt.*sign(gradient(v10, mode(diff(era5.time))*24));
% dwspdt = movmean(dwspdt, (6/24)./mode(diff(era5.time)), 'omitnan');
[~, mi] = min(abs(era5.time -t), [],1);
dwspdt = dwspdt(mi);
% dwspdt = naninterp1(era5.time, dwspdt, t, 2/24);


msk = ones(size(X)); msk(mag(X,Y)>100) = NaN;

th_wv = squeeze(dirnanmean(TH_wv_outer.*msk, [1 2]));
th_wnd = squeeze(dirnanmean(TH_wnd_outer.*msk, [1 2]));
dth = abs(th_wnd - th_wv);

% dth = abs(squeeze(dirnanmean(TH_wnd_outer - TH_wv_outer, [1 2])));
wsp = squeeze(nanmean(WSP, [1 2]));
tp = squeeze(nanmean(TP, [1 2]));
swh = squeeze(nanmean(SWH, [1 2]));
clear pi; cp = [(9.81).*tp]./(2*pi);
wage = cp./wsp;
vwage = cp./wsp.*cosd(dth); % if phase speed of waves is larger than wind speed 
                            % they are old waves (long is fast and swelly)
                            % this means that young waves is equal to small
                            % wave age
msk = ones(size(X)); msk(mag(X,Y)>75) = NaN; 
vortmag = relvort.*msk; vortmag(vortmag<0) = NaN; vortmag = squeeze(nanmedian(vortmag, [1 2])); 

msk = ones(size(X)); msk(mag(X,Y)>75) = NaN; 
currmag = mag(U,V).*msk; currmag(currmag<0) = NaN; currmag = squeeze(nanmedian(currmag, [1 2])); 


% -------------------------------------------------------------------------
figure(512344); clf; 
% setfigsize(gcf, [947         400]*1.2)
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(Na,3, [0.057 0.06], [0.13 0.025], [0.06 0.04]);
% [ha] = reorg_tightsubplot(oldha, [num2cell(2:2:Na*2) {[1:2:Na*2]} ]);

dx = 0.12; 
for hi=3:3:Na*3
    ha(hi).Position([1 3]) = ha(hi).Position([1 3]) + [1 -1].*dx;
    ha(hi-1).Position([1 3]) = ha(hi).Position([1 3]) + [-1 0].*(ha(hi).Position(3)+0.03);

    ha(hi-2).Position([1 3]) = ha(hi-2).Position([1 3]) + [0 1].*(dx.*2+0.06);
end

% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',9);
setaxes(ha(1:3:Na*2), 'YLabel', '$r^2$'); 
% setaxes(ha(1:Na+1), 'YLim', [0 1]); 
% setaxes(ha(1:Na), 'XLabel', '$\theta_{wnd} - \theta_{wv}$ [$^o$]'); 
% setaxes(ha(1:Na), 'XLim', [0 180]); 
% setaxes(ha(1:Na), 'YAxisLocation', 'right'); 
% setaxes(ha(1:Na), 'XTick', [0 45 90 135 ]); 
% quietaxes(ha(1:Na-1), 'x')
setaxes(ha(1:3:end), 'XLim', [datenum(2023,4,12.5) datenum(2023,5,9.35)]); 
linkaxes(ha(1:3:end), 'x')
autodatetick(ha(1:3:end), 'x')
quietaxes(ha(1:3:end-3), 'x')
quietaxes(ha(2:3), 'y')
quietaxes(ha(5:6), 'y')
quietaxes(ha(2:3), 'x')


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for ai=1:length(A)
    hi = ai; 
    % if ai>Na; hi = ai-1; end
    
    VARx = A{2, ai}; VARy = A{1, ai}; lab = Alabs{ai}; col = cmap(ai,:);
    tt = find(~isin(t, get(ha(1), 'XLim'))); 
    VARx(:,:,tt) = NaN; VARy(:,:,tt) = NaN;


    msk = ones(size(VARx)); msk(isnan(VARx + VARy)) = NaN;
    VARx = VARx.*msk; VARy = VARy.*msk;
    

    msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
    cc = arrayfun(@(i) ccfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    p = arrayfun(@(i) pfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    ccmsk = ones(size(X)); ccmsk(~isin(mag(X,Y), [30 50])) = NaN;

    n = sum(~isnan([VARx + VARy].*ccmsk), [1 2])./sum(~isnan(ccmsk), [1 2]); n = squeeze(n);
    cc(n<prcnanlim) = NaN;

    % tt = find(~isin(t, get(ha(end), 'XLim') + [-1 1]*0.4, 'inclusive')); 
    % cc(tt) = NaN;

    
    A_cc{ai} = cc'; A_n{ai} = n';
    % cc = abs(cc);
    cc_raw = cc;
    

    % if hi<=Na
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi*3-1); axes(thisha); hold on; 
    
        
        cc = A_cc{ai}; 
        
        x = dth; y = cc.^2; xlim([0 180]); xticks([0 90 180]);
        if ai>1; xlabel('$|\theta_{wnd} - \theta_{wave}|$ [$^o$]'); end

        % x = dirmovmean(x, 3, 'omitnan'); y = movmean(y, 3, 'omitnan');


        ii = find(~isnan(wsp + cc + dth));
        plot(x(ii), y(ii), 'kx', ...
            'MarkerSize',4, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);

        % ii = intersect(ii, find(wsp>10 & vortmag<0.5));
        % ii = intersect(ii, find(movmean(wsp,3, 'omitnan')>10 & movmean(vortmag,3, 'omitnan')<0.55));
        % ii = intersect(ii, find(movmean(wsp,3, 'omitnan')>10));
        % ii = intersect(ii, find(wsp>10));
        ii = intersect(ii, find(movmean(vortmag,3, 'omitnan')>0.55));
        % ii = intersect(ii, find(isin(wsp, [7 9])));
        plot(x(ii), y(ii), 'k^', ...
            'MarkerSize',4, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);
        % plotBestOptFit(gca, x(ii), y(ii), 'Color', col.*0.7)
        % if ai>2; ylim([0 0.5]); end
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi*3); axes(thisha); hold on; 

        msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
        
        x = VARx.*msk; y = VARy.*msk; 
        ii = find(cc.^2>0.25); 
        if isempty(ii); ii = find(cc.^2>0.125); end
        % x = nanmean(x(:,:,ii), 3); y = nanmean(y(:,:,ii),3);
        x =x(:,:,ii); y = y(:,:,ii);
        ii = find(~isnan(x + y));
        plot(x(ii), y(ii), 'k.', ...
            'MarkerSize',4, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);
        bins = [-1:0.1:1].*2;
        % bins = [-1:0.1:1].*prctile(abs(x(:)), 99);
        boxplot_bybin(gca, y(ii), x(ii), bins, 'Color', col, 'outliars', 0, 'nmin',5);

        % x = VARx.*msk; y = VARy.*msk; 
        % ii = find(cc.^2<0.25); 
        % % x = nanmean(x(:,:,ii), 3); y = nanmean(y(:,:,ii),3);
        % x =x(:,:,ii); y = y(:,:,ii);
        % ii = find(~isnan(x + y));
        % plot(x(ii), y(ii), 'k.', ...
        %     'MarkerSize',4, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));

        % cc = A_cc{ai}; 
        % 
        % % x = wsp; y = cc.^2; xlim([3 16]); xticks([0 5 10 15])
        % % x = dwspdt; y = cc.^2; %xlim([3 16]); xticks([0 5 10 15]);
        % % y = gradient(y);
        % x = vortmag; x = movmean(x,3,'omitnan');
        % % x = currmag; x = movmean(x,5,'omitnan');
        % 
        % x = movmean(x, 3, 'omitnan'); y = movmean(y, 3, 'omitnan');
        % 
        % 
        % ii = find(~isnan(wsp + cc + dth));
        % plot(x(ii), y(ii), 'kx', ...
        %     'MarkerSize',4, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);
        % plotBestOptFit(gca, x(ii), y(ii), 'Color', col.*0.7)
        % 
        % % ii = intersect(ii, find(dth>45 & movmean(vortmag,3, 'omitnan')<0.55));
        % ii = intersect(ii, find(movmean(vortmag,3, 'omitnan')>0.55));
        % plot(x(ii), y(ii), 'k^', ...
        %     'MarkerSize',4, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);

        % ii = find(~isnan(wsp + cc + dth));
        % ii = intersect(ii, find(isin(cc.^2, [0 0.25])));
        % plot(x(ii), y(ii), 'kx', ...
        %     'MarkerSize',4, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));
        % 
        % 
        % ii = intersect(ii, find(isin(cc.^2, [0.25 0.5])));
        % plot(x(ii), y(ii), 'ko', ...
        %     'MarkerSize',4, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));
        % 
        % ii = find(~isnan(wsp + cc + dth));
        % ii = intersect(ii, find(isin(cc.^2, [0.5 1])));
        % plot(x(ii), y(ii), 'k^', ...
        %     'MarkerSize',5, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);

        % 
        % ii = find(~isnan(wsp + cc + dth));
        % % xlim([0 180]); 
        % scatter(dth(ii), dwspdt(ii), 30, cc(ii).^2, 'filled', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2); 
        % caxis([0 1]);
        % colormap(flipud(cmocean('gray')))
        % colormap(nicejet)
        % % ylim([0 2])
    % end

    % if ai>2; ylim([0 0.5]); end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi*3-2); axes(thisha); hold on; 

    x = t; y = cc; ylim([-0.2 1])
    x = t; y = cc.^2; ylim([0 1]); if ai>2; ylim([0 0.5]); end
    plot(x,y, 'kx', 'Color', col, 'LineWidth',1); 


    % y95
    % arrayfun(@())
    % y = y.^2; 
    ystd = movstd(y, 5, 'omitnan'); 
    y0 = movmean(y,5,'omitnan'); 
    y0(isnan(y)) = NaN;
    % lims =
    shaded_limits(x,y0,ystd.*0.5.*[-1 1]+y0, col, 'DisplayName', lab)

    % yyaxis right; plot(t, vortmag, 'r-', 'DisplayName', 'vort')
    
    if hi<3 | ai == 4
        loc = 'northeast'; if contains(lab, 'SWH'); loc = 'southeast'; end
        leg = cleanLegend(gca, loc, 'NumColumns',4', 'Box', 'off', 'FontSize',8);
    end
end

% leg.Position(2) = leg.Position(2) + 0.02;
% leg.Position(1) = leg.Position(1) + 0.01;

savejpg(gcf, 'simple_timeseries_corr', savepath, saveopt)

%% basic one at a time.. 
Rcc = 100;
prcnanlim = 0.25;

% -------------------------------------------------------------------------

% Na = length(A)-1; 
% Lavg = 15; 
% -------------------------------------------------------------------------



wsp = squeeze(mag( nanmean(era5.u10,[1 2]), nanmean(era5.v10,[1 2]))); wsp = movmean(wsp, (4/24)./mode(diff(era5.time)));
dwspdt = gradient(wsp, mode(diff(era5.time))*24);
dwspdt = naninterp1(era5.time, dwspdt, t, 2/24);

dth = abs(squeeze(dirnanmean(TH_wnd_outer - TH_wv_outer, [1 2])));
wsp = squeeze(nanmean(WSP, [1 2]));
tp = squeeze(nanmean(TP, [1 2]));
swh = squeeze(nanmean(SWH, [1 2]));
clear pi; cp = [(9.81).*tp]./(2*pi);
wage = cp./wsp;
vwage = cp./wsp.*cosd(dth); % if phase speed of waves is larger than wind speed 
                            % they are old waves (long is fast and swelly)
                            % this means that young waves is equal to small
                            % wave age
                            



% -------------------------------------------------------------------------
figure(512344); clf; 
% setfigsize(gcf, [947         400]*1.2)
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
ha = tight_subplot(Na,2, [0.05 0.06], [0.09 0.025], [0.06 0.04]);
% [ha] = reorg_tightsubplot(oldha, [num2cell(2:2:Na*2) {[1:2:Na*2]} ]);

dx = 0.25; 
for hi=2:2:Na*2
    ha(hi).Position([1 3]) = ha(hi).Position([1 3]) + [1 -1].*dx;
    ha(hi-1).Position([1 3]) = ha(hi-1).Position([1 3]) + [0 1].*dx;
end

% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',9);
% setaxes(ha(1:Na+1), 'YLabel', '$r^2$'); 
% setaxes(ha(1:Na+1), 'YLim', [0 1]); 
% setaxes(ha(1:Na), 'XLabel', '$\theta_{wnd} - \theta_{wv}$ [$^o$]'); 
% setaxes(ha(1:Na), 'XLim', [0 180]); 
% setaxes(ha(1:Na), 'YAxisLocation', 'right'); 
% setaxes(ha(1:Na), 'XTick', [0 45 90 135 ]); 
% quietaxes(ha(1:Na-1), 'x')
setaxes(ha(1:2:end), 'XLim', [datenum(2023,4,12.5) datenum(2023,5,9.35)]); 
linkaxes(ha(1:2:end), 'x')
autodatetick(ha(1:2:end), 'x')
quietaxes(ha(1:2:end-2), 'x')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for ai=1:length(A)
    hi = ai; 
    if ai>Na; hi = ai-1; end
    
    VARx = A{2, ai}; VARy = A{1, ai}; lab = Alabs{ai}; col = cmap(ai,:);
    % 
    msk = ones(size(VARx)); msk(isnan(VARx + VARy)) = NaN;
    VARx = VARx.*msk; VARy = VARy.*msk;

    msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
    cc = arrayfun(@(i) ccfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    p = arrayfun(@(i) pfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    ccmsk = ones(size(X)); ccmsk(~isin(mag(X,Y), [30 50])) = NaN;

    n = sum(~isnan([VARx + VARy].*ccmsk), [1 2])./sum(~isnan(ccmsk), [1 2]); n = squeeze(n);
    cc(n<prcnanlim) = NaN;

    % tt = find(~isin(t, get(ha(end), 'XLim') + [-1 1]*0.4, 'inclusive')); 
    % cc(tt) = NaN;

    
    A_cc{ai} = cc'; A_n{ai} = n';
    % cc = abs(cc);
    cc_raw = cc;
    

    % if hi<=Na
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi*2); axes(thisha); hold on; 
    
        
        cc = A_cc{ai}; 
        

        % x = dth;
        % y = cc.^2; 
        % x = movmean(x,3, 'omitnan');
        % y = movmean(y,3, 'omitnan');
        % ii = find(~isnan(wsp + cc + dth));
        % ii = intersect(ii, find(wsp<=10));
        % plot(x(ii), y(ii), 'ko', ...
        %         'MarkerSize',5, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);
        % 
        % ii = find(~isnan(wsp + cc + dth));
        % ii = intersect(ii, find(wsp>10));
        % plot(x(ii), y(ii), 'k^', ...
        %     'MarkerSize',5, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));


        % wsplim = 8;
        % bins = [0:0.05:1];
        % II = find(~isnan(wsp + cc + dth));
        % histogram(y(II), bins, 'FaceColor', col, 'FaceAlpha',1);
        % xline(nanmedian(y(II)), 'Color', col, 'LineWidth',2)
        % ii = intersect(ii, find(abs(dth)<30));
        % histogram(y(ii), bins, 'FaceColor', col);
        % ii = intersect(ii, find(abs(dth)<45));

        % 
        % ii = intersect(II, find(abs(dth)<45));
        % histogram(y(ii), bins, 'FaceColor', col, 'FaceAlpha',1);
        % xline(nanmedian(y(ii)), 'Color', col, 'LineWidth',2)
        % 
        % ii = intersect(II, find(abs(dth)<45 & wsp<wsplim));
        % histogram(y(ii), bins, 'FaceColor', sqrt(col), 'FaceAlpha',1);
        % xline(nanmedian(y(ii)), 'Color', sqrt(col), 'LineWidth',2)
        % 
        % ii = intersect(II, find(abs(dth)>45));
        % histogram(y(ii), bins, 'FaceColor', sqrt(col).*0.4, 'FaceAlpha',1);
        % xline(nanmedian(y(ii)), 'Color', sqrt(col).*0.4, 'LineWidth',2)



        % ii = intersect(II, find(wsp<wsplim));
        % histogram(y(ii), bins, 'FaceColor', sqrt(col), 'FaceAlpha',1);
        % xline(nanmedian(y(ii)), 'Color', sqrt(col), 'LineWidth',2)
        % ylim([0 20])








        % ii = intersect(ii, find(abs(dth)<30));
        % plot(x(ii), y(ii), 'ko', ...
        %     'MarkerSize',5, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));

        x = dth; y = wsp; ylim([0 15]); xlim([0 180])

        ii = find(~isnan(wsp + cc + dth));
        ii = intersect(ii, find(isin(cc.^2, [0 0.25])));
        plot(x(ii), y(ii), 'kx', ...
            'MarkerSize',4, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));


        ii = intersect(ii, find(isin(cc.^2, [0.25 0.5])));
        plot(x(ii), y(ii), 'ko', ...
            'MarkerSize',4, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));
        
        ii = find(~isnan(wsp + cc + dth));
        ii = intersect(ii, find(isin(cc.^2, [0.5 1])));
        plot(x(ii), y(ii), 'k^', ...
            'MarkerSize',5, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);


        % 
        % scatter(dth(ii), wsp(ii), 30, cc(ii).^2, 'filled', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2); ylim([3 17])
        % caxis([0 1]);
        % colormap(flipud(cmocean('gray')))

    % end



    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(hi*2-1); axes(thisha); hold on; 

    x = t; y = cc; ylim([-0.2 1])
    x = t; y = cc.^2; ylim([0 1]); if ai>2; ylim([0 0.5]); end
    plot(x,y, 'kx', 'Color', col, 'LineWidth',1); 


    % y95
    % arrayfun(@())
    % y = y.^2; 
    ystd = movstd(y, 5, 'omitnan'); 
    y0 = movmean(y,5,'omitnan'); 
    y0(isnan(y)) = NaN;
    % lims =
    shaded_limits(x,y0,ystd.*0.5.*[-1 1]+y0, col, 'DisplayName', lab)
    
    if hi<3 | ai == 4
        leg = cleanLegend(gca, 'northeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
    end
end

% leg.Position(2) = leg.Position(2) + 0.02;
% leg.Position(1) = leg.Position(1) + 0.01;


%% [save] simple up down

% -------------------------------------------------------------------------

Na = length(A); Na = 3;
% -------------------------------------------------------------------------


wsp = squeeze(mag( nanmean(era5.u10,[1 2]), nanmean(era5.v10,[1 2]))); wsp = movmean(wsp, (4/24)./mode(diff(era5.time)));
dwspdt = gradient(wsp, mode(diff(era5.time))*24);
dwspdt = naninterp1(era5.time, dwspdt, t, 2/24);

dth = abs(squeeze(dirnanmean(TH_wnd_outer - TH_wv_outer, [1 2])));
wsp = squeeze(nanmean(WSP, [1 2]));
tp = squeeze(nanmean(TP, [1 2]));
swh = squeeze(nanmean(SWH, [1 2]));
clear pi; cp = [(9.81).*tp]./(2*pi);
wage = cp./wsp;
vwage = cp./wsp.*cosd(dth); % if phase speed of waves is larger than wind speed 
                            % they are old waves (long is fast and swelly)
                            % this means that young waves is equal to small
                            % wave age


% -------------------------------------------------------------------------
figure(5123); clf; 
setfigsize(gcf, [858   483])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
oldha = tight_subplot(2, Na, [0.01 0.015], [0.07 0.05], [0.06 0.03]);
[ha] = reorg_tightsubplot(oldha, [num2cell(1:Na) {[Na+1:Na*2]} ]);
dy = 0.2; 
for hi=1:Na
    ha(hi).Position([2 4]) = ha(hi).Position([2 4]) + [1 -1].*dy;
end
hi = length(ha); ha(hi).Position([2 4]) = ha(hi).Position([2 4])  + [0 1].*dy*0.5;


% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',10);
setaxes(ha(1:Na+1), 'YLabel', '$r^2$'); 
setaxes(ha(1:Na+1), 'YLim', [0 1]); 
setaxes(ha(1:Na), 'XLabel', '$\theta_{wnd} - \theta_{wv}$ [$^o$]'); 
setaxes(ha(1:Na), 'XLim', [0 180]); 
setaxes(ha(1:Na), 'XTick', [0 45 90 135 ]); 
quietaxes(ha(2:Na), 'y')
setaxes(ha(end:end), 'XLim', [datenum(2023,4,12.5) datenum(2023,5,9.35)]); 
autodatetick(ha(end:end), 'x')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for ai=1:length(A)
    hi = ai; 

    VARx = A{2, ai}; VARy = A{1, ai}; lab = Alabs{ai}; col = cmap(ai,:);
    % 
    msk = ones(size(VARx)); msk(isnan(VARx + VARy)) = NaN;
    VARx = VARx.*msk; VARy = VARy.*msk;

    msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
    cc = arrayfun(@(i) ccfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    p = arrayfun(@(i) pfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    ccmsk = ones(size(X)); ccmsk(~isin(mag(X,Y), [30 50])) = NaN;

    n = sum(~isnan([VARx + VARy].*ccmsk), [1 2])./sum(~isnan(ccmsk), [1 2]); n = squeeze(n);
    cc(n<0.8) = NaN;

    tt = find(~isin(t, get(ha(end), 'XLim') + [-1 1]*0.4, 'inclusive')); 
    cc(tt) = NaN;

    
    A_cc{ai} = cc; A_n{ai} = n';
    % cc = abs(cc);
    cc_raw = cc;
    

    if hi<=Na
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi); axes(thisha); hold on; 
    
        
        cc = A_cc{ai};
        
        
        % 
        plot(dth, cc.^2, 'kx', ...
            'MarkerSize',5, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);
    
        
    
        ii = find(wsp>10);
        plot(dth(ii), cc(ii).^2, 'k^','DisplayName', '$|U_{10}|>10$ m/s', ...
            'MarkerSize',3, 'LineWidth',2, 'Color', col*0.5, 'MarkerFaceColor', col);
    
        cleanLegend(gca, 'northeast', 'FontSize',8)
    end


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(end); axes(thisha); hold on; 
    plot(t, cc.^2, 'kx', 'Color', col, 'LineWidth',1); 

    x = t; y = cc.^2; 
    

    % y95
    % arrayfun(@())
    y = cc.^2; 
    ystd = movstd(y, 5, 'omitnan'); 
    y0 = movmedian(y,4,'omitnan'); 
    y0(isnan(y)) = NaN;
    
    shaded_limits(x,y0,ystd.*0.5.*[-1 1]'+y0, col, 'DisplayName', lab)
    
    

end

leg = cleanLegend(gca, 'northeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
leg.Position(2) = leg.Position(2) + 0.02;
leg.Position(1) = leg.Position(1) + 0.01;



savejpg(gcf, 'fig_timeseries_correlation', savepath, saveopt)


%% [save] even simpler LR
Rcc = 100;
prcnan = 0.25;
% -------------------------------------------------------------------------

Na = length(A); 
Lavg = 15; 
% -------------------------------------------------------------------------


wsp = squeeze(mag( nanmean(era5.u10,[1 2]), nanmean(era5.v10,[1 2]))); wsp = movmean(wsp, (4/24)./mode(diff(era5.time)));
dwspdt = gradient(wsp, mode(diff(era5.time))*24);
dwspdt = naninterp1(era5.time, dwspdt, t, 2/24);

dth = abs(squeeze(dirnanmean(TH_wnd_outer - TH_wv_outer, [1 2])));
wsp = squeeze(nanmean(WSP, [1 2]));
tp = squeeze(nanmean(TP, [1 2]));
swh = squeeze(nanmean(SWH, [1 2]));
clear pi; cp = [(9.81).*tp]./(2*pi);
wage = cp./wsp;
vwage = cp./wsp.*cosd(dth); % if phase speed of waves is larger than wind speed 
                            % they are old waves (long is fast and swelly)
                            % this means that young waves is equal to small
                            % wave age
                            



% -------------------------------------------------------------------------
figure(5123); clf; 
setfigsize(gcf, [947         400]*1.2)
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
oldha = tight_subplot(Na,2, [0.02 0.06], [0.09 0.025], [0.06 0.04]);
[ha] = reorg_tightsubplot(oldha, [num2cell(2:2:Na*2) {[1:2:Na*2]} ]);

dx = 0.3; 
for hi=1:Na
    ha(hi).Position([1 3]) = ha(hi).Position([1 3]) + [1 -1].*dx;
end
hi = length(ha); ha(hi).Position([1 3]) = ha(hi).Position([1 3])  + [0 1].*dx;


thisha = axes; thisha.Position = ha(end).Position; 
thisha.Position(4) = thisha.Position(4)*0.3;
% ha(end).Position([2 4]) = [sum(thisha.Position([2 4]))+0.05  ha(end).Position(4)*0.7];
ha(end).Position(2) = sum(thisha.Position([2 4]))+0.03;
ha(end).Position(4) = ha(end).Position(4) - (ha(end).Position(2) - thisha.Position(2));
ha(end+1) = thisha;


% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',9);
setaxes(ha(1:Na+1), 'YLabel', '$r^2$'); 
setaxes(ha(1:Na+1), 'YLim', [0 1]); 
setaxes(ha(1:Na), 'XLabel', '$\theta_{wnd} - \theta_{wv}$ [$^o$]'); 
setaxes(ha(1:Na), 'XLim', [0 180]); 
setaxes(ha(1:Na), 'YAxisLocation', 'right'); 
setaxes(ha(1:Na), 'XTick', [0 45 90 135 ]); 
quietaxes(ha(1:Na-1), 'x')
setaxes(ha(end-1:end), 'XLim', [datenum(2023,4,12.5) datenum(2023,5,9.35)]); 
autodatetick(ha(end-1:end), 'x')
quietaxes(ha(end-1), 'x')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for ai=1:length(A)
    hi = ai; 
    
    VARx = A{2, ai}; VARy = A{1, ai}; lab = Alabs{ai}; col = cmap(ai,:);
    % 

    % tt = find(~isin(t, get(ha(end), 'XLim'))); 
    % VARx(:,:,tt) = NaN; VARy(:,:,tt) = NaN;

    msk = ones(size(VARx)); msk(isnan(VARx + VARy)) = NaN;
    VARx = VARx.*msk; VARy = VARy.*msk;

    msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
    cc = arrayfun(@(i) ccfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    p = arrayfun(@(i) pfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    ccmsk = ones(size(X)); ccmsk(~isin(mag(X,Y), [30 50])) = NaN;

    n = sum(~isnan([VARx + VARy].*ccmsk), [1 2])./sum(~isnan(ccmsk), [1 2]); n = squeeze(n);
    cc(n<prcnan) = NaN;

    

    % cc(p>0.05) = NaN;
    
    A_cc{ai} = cc'; A_n{ai} = n';
    A_p{ai} = p';
    % cc = abs(cc);
    cc_raw = cc;
    

    if hi<=Na
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi); axes(thisha); hold on; 
    
        
        cc = A_cc{ai};
        tt = find(~isin(t, get(ha(end), 'XLim') + [-1 1]*0.4, 'inclusive')); 
        cc(tt) = NaN;
        
        
        % 
        plot(dth, cc.^2, 'kx', ...
            'MarkerSize',5, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);
        % 
        % 
        % 
        ii = find(wsp>10);
        plot(dth(ii), cc(ii).^2, 'k^','DisplayName', '$|U_{10}|>10$ m/s', ...
            'MarkerSize',3, 'LineWidth',2, 'Color', col*0.5, 'MarkerFaceColor', col);
        cleanLegend(gca, 'northeast', 'FontSize',8)



        % 
        % nmin = 3; 
        % bins = [0 45 90 135 180];
        % % bins = [0  90  180]; 
        % dbin = mode(diff(bins));
        % % boxplot_bybin(gca, cc.^2, dth, bins, 'nmin', 4, 'Color', col, 'FaceAlpha',0.4)
        % wsplim = 10;
        % ii = find(wsp>wsplim);
        % plot(dth(ii), cc(ii).^2, 'kx', ...
        %     'MarkerSize',5, 'LineWidth',1, 'Color', sqrt(col), 'MarkerFaceColor', sqrt(col));
        % boxplot_bybin(gca, cc(ii).^2, dth(ii)+dbin/5, bins+dbin/6, 'nmin', nmin, 'Color', sqrt(col), 'FaceAlpha',0.7, 'width',dbin/7, 'DisplayName', ['$|U_{10}|>$' num2str(wsplim) ' m/s'])
        % ii = find(wsp<=wsplim);
        % plot(dth(ii), cc(ii).^2, 'kx', ...
        %     'MarkerSize',5, 'LineWidth',1, 'Color', col.*0.7, 'MarkerFaceColor', col.*0.7);
        % boxplot_bybin(gca, cc(ii).^2, dth(ii)-dbin/5, bins-dbin/6, 'nmin', nmin, 'Color', col.*0.7, 'FaceAlpha',0.7, 'width',dbin/7, 'DisplayName', ['$|U_{10}|<$' num2str(wsplim) ' m/s'])
        % cleanLegend(gca, 'northeast', 'FontSize',8, 'NumColumns',2, 'box', 'off')

        % ii = find(~isnan(wsp + cc + dth));
        % scatter(dth(ii), wsp(ii), 30, cc(ii).^2, 'filled', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2); ylim([3 17])
        % caxis([0 1]);
        % colormap(cmocean('amp'))

    end



    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(end-1); axes(thisha); hold on; 
    cc = A_cc{ai};
    x = t; y = cc; ylim([-0.2 1])
    x = t; y = cc.^2; ylim([0 1])
    ii = find(p>0.05);
    plot(x(ii),y(ii), 'k.', 'Color', col.^(1/4), 'LineWidth',1, 'Markersize',8); 
    % plot(x(ii),y(ii), 'kx', 'Color', col, 'LineWidth',1); 
    ii = find(p<=0.05);
    plot(x(ii),y(ii), 'kx', 'Color', col, 'LineWidth',1); 
    % plot(x(ii),y(ii), 'k.', 'Color', col.*0.8, 'LineWidth',1, 'Markersize',10); 


    % y95
    % arrayfun(@())
    % y = y.^2; 
    y(p>0.05) = NaN;
    ystd = movstd(y, 4, 'omitnan'); 
    y0 = movmean(y,4,'omitnan'); 
    y0(isnan(y)) = NaN;
    % lims =
    shaded_limits(x,y0,ystd.*0.5.*[-1 1]+y0, col, 'DisplayName', lab)
    
    

end

leg = cleanLegend(gca, 'northeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
leg.Position(2) = leg.Position(2) + 0.02;
leg.Position(1) = leg.Position(1) + 0.01;




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(end); axes(thisha); hold on; 
yyaxis left; set(gca, 'YColor', 'k');
ylim([0 17]); ylabel('$|U_{10}|$ [m/s]'); yticks([0:5:20]);
x = era5.time; y = mag(era5.u10, era5.v10); y = squeeze(nanmean(y,[1 2])); 
plot(x,y,'k-', 'linewidth',1, 'DisplayName', 'ERA5')
% plot(t, wsp,'kx', 'linewidth',1)
% leg = cleanLegend(gca, 'southeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
col = [0.25 0.6 0.6];
yyaxis right; set(gca, 'YColor', col); ylabel('$\theta_{wnd} - \theta_{wv}$ [$^o$]', 'Interpreter', 'latex');
yticks([0:90:180]); ylim([0 180]);
% plot(t, abs(dth), 'rx', 'Color', col, 'linewidth',1)
x = era5.time; y_wnd = mod(dirnanmean(atan2d(era5.v10, era5.u10), [1 2]), 360); y_wv = mod(dirnanmean(era5.mwd, [1 2]), 360);
y_wv = met2oc(oc2pol(y_wv));
y =  y_wnd - y_wv; 
y = squeeze(y); y = mod(y, 360); y(y>180) = 360-y(y>180);
plot(x,y,'k:', 'linewidth',2, 'DisplayName', 'ERA5', 'Color', col)






savejpg(gcf, 'fig_timeseries_correlation_LR', savepath, saveopt)


%% [save] simple LR
Rcc = 100;
% -------------------------------------------------------------------------

Na = length(A); Na = 3; 
Lavg = 15; 
% -------------------------------------------------------------------------


wsp = squeeze(mag( nanmean(era5.u10,[1 2]), nanmean(era5.v10,[1 2]))); wsp = movmean(wsp, (4/24)./mode(diff(era5.time)));
dwspdt = gradient(wsp, mode(diff(era5.time))*24);
dwspdt = naninterp1(era5.time, dwspdt, t, 2/24);

dth = abs(squeeze(dirnanmean(TH_wnd_outer - TH_wv_outer, [1 2])));
wsp = squeeze(nanmean(WSP, [1 2]));
tp = squeeze(nanmean(TP, [1 2]));
swh = squeeze(nanmean(SWH, [1 2]));
clear pi; cp = [(9.81).*tp]./(2*pi);
wage = cp./wsp;
vwage = cp./wsp.*cosd(dth); % if phase speed of waves is larger than wind speed 
                            % they are old waves (long is fast and swelly)
                            % this means that young waves is equal to small
                            % wave age
                            



% -------------------------------------------------------------------------
figure(5123); clf; 
setfigsize(gcf, [947         400]*1.2)
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
oldha = tight_subplot(Na,2, [0.02 0.06], [0.09 0.025], [0.06 0.04]);
[ha] = reorg_tightsubplot(oldha, [num2cell(2:2:Na*2) {[1:2:Na*2]} ]);

dx = 0.2; 
for hi=1:Na
    ha(hi).Position([1 3]) = ha(hi).Position([1 3]) + [1 -1].*dx;
end
hi = length(ha); ha(hi).Position([1 3]) = ha(hi).Position([1 3])  + [0 1].*dx;


thisha = axes; thisha.Position = ha(end).Position; 
thisha.Position(4) = thisha.Position(4)*0.22;
% ha(end).Position([2 4]) = [sum(thisha.Position([2 4]))+0.05  ha(end).Position(4)*0.7];
ha(end).Position(2) = sum(thisha.Position([2 4]))+0.03;
ha(end).Position(4) = ha(end).Position(4) - (ha(end).Position(2) - thisha.Position(2));
ha(end+1) = thisha;


% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',9);
setaxes(ha(1:Na+1), 'YLabel', '$r^2$'); 
setaxes(ha(1:Na+1), 'YLim', [0 1]); 
setaxes(ha(1:Na), 'XLabel', '$\theta_{wnd} - \theta_{wv}$ [$^o$]'); 
setaxes(ha(1:Na), 'XLim', [0 180]); 
setaxes(ha(1:Na), 'YAxisLocation', 'right'); 
setaxes(ha(1:Na), 'XTick', [0 45 90 135 ]); 
quietaxes(ha(1:Na-1), 'x')
setaxes(ha(end-1:end), 'XLim', [datenum(2023,4,12.5) datenum(2023,5,9.35)]); 
autodatetick(ha(end-1:end), 'x')
quietaxes(ha(end-1), 'x')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for ai=1:length(A)
    hi = ai; 
    
    VARx = A{2, ai}; VARy = A{1, ai}; lab = Alabs{ai}; col = cmap(ai,:);
    % 
    msk = ones(size(VARx)); msk(isnan(VARx + VARy)) = NaN;
    VARx = VARx.*msk; VARy = VARy.*msk;

    msk = ones(size(X)); msk(mag(X,Y)>Rcc) = NaN;
    cc = arrayfun(@(i) ccfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    p = arrayfun(@(i) pfxn(VARx(:,:,i).*msk, VARy(:,:,i).*msk) , 1:size(VARx,3))    ;
    ccmsk = ones(size(X)); ccmsk(~isin(mag(X,Y), [30 50])) = NaN;

    n = sum(~isnan([VARx + VARy].*ccmsk), [1 2])./sum(~isnan(ccmsk), [1 2]); n = squeeze(n);
    cc(n<0.8) = NaN;

    tt = find(~isin(t, get(ha(end), 'XLim') + [-1 1]*0.4, 'inclusive')); 
    cc(tt) = NaN;

    
    A_cc{ai} = cc; A_n{ai} = n';
    % cc = abs(cc);
    cc_raw = cc;
    

    if hi<=Na
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        thisha = ha(hi); axes(thisha); hold on; 
    
        
        cc = A_cc{ai};
        
        
        % 
        plot(dth, cc.^2, 'kx', ...
            'MarkerSize',5, 'LineWidth',1, 'Color', col, 'MarkerFaceColor', col);
    
        
    
        ii = find(wsp>10);
        plot(dth(ii), cc(ii).^2, 'k^','DisplayName', '$|U_{10}|>10$ m/s', ...
            'MarkerSize',3, 'LineWidth',2, 'Color', col*0.5, 'MarkerFaceColor', col);
    
        cleanLegend(gca, 'northeast', 'FontSize',8)
    end



    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(end-1); axes(thisha); hold on; 
    plot(t, cc.^2, 'kx', 'Color', col, 'LineWidth',1); 

    x = t; y = cc.^2; 
    

    % y95
    % arrayfun(@())
    y = cc.^2; 
    ystd = movstd(y, 5, 'omitnan'); 
    y0 = movmedian(y,4,'omitnan'); 
    y0(isnan(y)) = NaN;
    
    shaded_limits(x,y0,ystd.*0.5.*[-1 1]'+y0, col, 'DisplayName', lab)
    
    

end

leg = cleanLegend(gca, 'northeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
leg.Position(2) = leg.Position(2) + 0.02;
leg.Position(1) = leg.Position(1) + 0.01;




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(end); axes(thisha); hold on; 
yyaxis left; set(gca, 'YColor', 'k');
ylim([0 17]); ylabel('$|U_{10}|$ [m/s]'); yticks([0:5:20]);
x = era5.time; y = mag(era5.u10, era5.v10); y = squeeze(nanmean(y,[1 2])); 
plot(x,y,'k-', 'linewidth',1, 'DisplayName', 'ERA5')
% plot(t, wsp,'kx', 'linewidth',1)
% leg = cleanLegend(gca, 'southeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
col = [0.25 0.6 0.6];
yyaxis right; set(gca, 'YColor', col); ylabel('$\theta_{wnd} - \theta_{wv}$ [$^o$]', 'Interpreter', 'latex');
yticks([0:90:180]); ylim([0 180]);
% plot(t, abs(dth), 'rx', 'Color', col, 'linewidth',1)
x = era5.time; y_wnd = mod(dirnanmean(atan2d(era5.v10, era5.u10), [1 2]), 360); y_wv = mod(dirnanmean(era5.mwd, [1 2]), 360);
y_wv = met2oc(oc2pol(y_wv));
y =  y_wnd - y_wv; 
y = squeeze(y); y = mod(y, 360); y(y>180) = 360-y(y>180);
plot(x,y,'k:', 'linewidth',2, 'DisplayName', 'ERA5', 'Color', col)






savejpg(gcf, 'fig_timeseries_correlation_LR', savepath, saveopt)


%%
% -------------------------------------------------------------------------


wsp = squeeze(mag( nanmean(era5.u10,[1 2]), nanmean(era5.v10,[1 2]))); wsp = movmean(wsp, (4/24)./mode(diff(era5.time)));
dwspdt = gradient(wsp, mode(diff(era5.time))*24);
dwspdt = naninterp1(era5.time, dwspdt, t, 2/24);

dth = abs(squeeze(dirnanmean(TH_wnd_outer - TH_wv_outer, [1 2])));
wsp = squeeze(nanmean(WSP, [1 2]));
tp = squeeze(nanmean(TP, [1 2]));
swh = squeeze(nanmean(SWH, [1 2]));
clear pi; cp = [(9.81).*tp]./(2*pi);
wage = cp./wsp;
vwage = cp./wsp.*cosd(dth); % if phase speed of waves is larger than wind speed 
                            % they are old waves (long is fast and swelly)
                            % this means that young waves is equal to small
                            % wave age
                            



% -------------------------------------------------------------------------
figure(5123); clf; 
% setfigsize(gcf, [858   453])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
oldha = tight_subplot(3, Na, [0.05 0.015], [0.07 0.05], [0.06 0.07]);
[ha] = reorg_tightsubplot(oldha, [num2cell(1:Na) {[Na+1:Na*2]} {[Na+1:Na*2]+Na}]);
dy = 0.05; 
for hi=1:Na
    ha(hi).Position([2 4]) = ha(hi).Position([2 4]) + [1 -1].*dy;
end
hi = length(ha); ha(hi).Position([2 4]) = ha(hi).Position([2 4])  + [0 1].*dy*0.5;
hi = length(ha)-1; ha(hi).Position([2 4]) = ha(hi).Position([2 4])  + [0 1].*dy*0.5;

dy = 0.1; 
hi = length(ha)-1; ha(hi).Position([2 4]) = ha(hi).Position([2 4])  + [-1 1].*dy;
hi = length(ha); ha(hi).Position([2 4]) = ha(hi).Position([2 4])  + [0 -1].*dy;


% ha = tight_subplot(length(A),2, [0.05 0.05], [0.05 0.05], [0.05 0.1]);
% dx = 0.4; 
% for hi=1:2:length(ha)
%     ha(hi).Position([1 3]) = ha(hi).Position([1 3]).*([1 1] + [0 1].*dx);
%     ha(hi+1).Position([1 3]) = ha(hi+1).Position([1 3]).*([1 1] + [1 -1].*dx);
% end

% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',10);
% setaxes(ha, 'DataAspectRatio', [1 1 1]); 
% setaxes(ha, 'XLim', lonlims); 
% setaxes(ha, 'YLim', latlims); 
setaxes(ha(1:length(A)+1), 'YLabel', '$r^2$'); 
setaxes(ha(1:length(A)+1), 'YLim', [0 1]); 
setaxes(ha(1:length(A)), 'XLabel', '$\theta_{wnd} - \theta_{wv}$'); 
% format_mapticks(ha(1:4), 'x')
setaxes(ha(1:length(A)), 'XLim', [0 180]); 
% setaxes(ha(1:4), 'XLim', [0 20]); 
setaxes(ha(1:length(A)), 'XTick', [0 45 90 135 ]); 
quietaxes(ha(2:length(A)), 'y')
setaxes(ha(end-1:end), 'XLim', [datenum(2023,4,12.5) datenum(2023,5,9.35)]); 
autodatetick(ha(end-1:end), 'x')
quietaxes(ha(end-1), 'x')
% linkaxes(ha)
% format_mapticks(ha, 'xy');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for ai=1:length(A)
    VARx = A{2, ai}; VARy = A{1, ai}; lab = Alabs{ai}; col = cmap(ai,:);
    % 
    msk = ones(size(VARx)); msk(isnan(VARx + VARy)) = NaN;
    VARx = VARx.*msk; VARy = VARy.*msk;

    
    cc = arrayfun(@(i) ccfxn(VARx(:,:,i).*cmsk, VARy(:,:,i).*cmsk) , 1:size(VARx,3))    ;
    p = arrayfun(@(i) pfxn(VARx(:,:,i).*cmsk, VARy(:,:,i).*cmsk) , 1:size(VARx,3))    ;
    ccmsk = ones(size(X)); ccmsk(~isin(mag(X,Y), [30 50])) = NaN;

    n = sum(~isnan([VARx + VARy].*ccmsk), [1 2])./sum(~isnan(ccmsk), [1 2]); n = squeeze(n);
    cc(n<0.8) = NaN;

    tt = find(~isin(t, get(ha(end), 'XLim') + [-1 1]*0.4, 'inclusive')); 
    cc(tt) = NaN;

    
    A_cc{ai} = cc; A_n{ai} = n';
    % cc = abs(cc);
    cc_raw = cc;
    

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(ai); axes(thisha); hold on; 

    
    cc = A_cc{ai};
    
    
    % 
    % plot(dth, cc.^2, 'k.', 'MarkerSize',10);

    

    ii = find(wsp>10);
    plot(dth(ii), cc(ii).^2, 'k^','DisplayName', '$|U_{10}|>10$ m/s', ...
        'MarkerSize',3, 'LineWidth',2, 'Color', col*0.5, 'MarkerFaceColor', col);

    ii = find(wsp<10);
    % ii = find(movmean(wsp,5)<10);
    plot(dth(ii), cc(ii).^2, 'k.','DisplayName', '$|U_{10}|<10$ m/s',...
        'MarkerSize',15, 'Color', col);
    % ii = find(movmean(wsp,5)>10);


    % ii = find(tp<6);
    % plot(dth(ii), cc(ii).^2, 'k*','DisplayName', '$|U_{10}|>10$ m/s', ...
    %     'MarkerSize',10, 'LineWidth',2, 'Color', col*0.1, 'MarkerFaceColor', col*0.1);


    cleanLegend(gca, 'northeast', 'FontSize',8)

    % bins = [0:30:180]
    % boxplot_bybin(gca, cc.^2, dth, bins, 'Color', col, 'FaceAlpha', 0.2, 'nmin',4, 'EdgeColor', col.*0.5, 'LineWidth',1.5, 'DisplayName', lab)
    % cleanLegend(gca, 'northeast', 'box', 'off')


    % bins = [0:0.1:1];
    % xlim(minmax(bins)); ylim([0 20])
    % histogram(cc.^2, bins, 'FaceColor', 'w'); 


    % ii = find(dth);
    % [N,edges] = histcounts(cc(ii).^2, bins);
    % plot(bins(1:end-1), movmean(N,2), 'k-', 'Linewidth',2, 'Color', col*0.7)
    % 
    % ii = find(dth<45);
    % [N,edges] = histcounts(cc(ii).^2, bins);
    % plot(bins(1:end-1), movmean(N,2), 'k-', 'Linewidth',2, 'Color', col*0.9)
    % 
    % ii = find(dth<45 & wsp<10 & wsp>4);
    % [N,edges] = histcounts(cc(ii).^2, bins);
    % plot(bins(1:end-1), movmean(N,2), 'k-', 'Linewidth',2, 'Color', col)

    
    % 
    % ii = find(dth<45);
    % histogram(cc(ii).^2, bins, 'FaceColor', col); 
    
    
    

    % scatter(dth, abs(dwspdt), 30, cc, 'filled'); ylim([0 2]); caxis([0 1]); colormap(nicejet)
    % scatter(dth,wsp, 30, cc, 'filled'); ylim([0 15]); caxis([0 1]); colormap(nicejet)

    
   
    % plot(dth, wsp, 'k.', 'MarkerSize',10, 'Color', [1 1 1]*0.6); ylim([0 15]); 
    % ii = find([cc.^2]>0.5);
    % plot(dth(ii), wsp(ii), 'k.', 'MarkerSize',15, 'Color', [1 1 1]*0.2);

    % scatter(dth, wage, 70, cc.^2, 'filled');
    % scatter(vwage, cc.^2, 70, cc.^2, 'filled');
    % ylim([0 3]); 
    % xlim([0 3]); xticks([0.5 1 ])
    % colormap(gca, nicejet)
    
    
    % plot(dth, cc.^2, 'k.', 'MarkerSize',10);
    % bins = [0:30:180]
    % boxplot_bybin(gca, cc.^2, dth, bins, 'Color', col, 'FaceAlpha', 0.2, 'nmin',4, 'EdgeColor', col.*0.5, 'LineWidth',1.5, 'DisplayName', lab)
    % cleanLegend(gca, 'northeast', 'box', 'off')




    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(end-1); axes(thisha); hold on; 
    plot(t, cc.^2, 'kx', 'Color', col, 'LineWidth',1); 

    % if contains(Alabs{ai}, 'WSPa') & contains(Alabs{ai}, 'hat')
    %     ii = find(wsp>10);
    %     plot(t(ii), cc(ii).^2, 'k^','DisplayName', '$|U_{10}|>10$ m/s', ...
    %         'MarkerSize',2, 'LineWidth',2, 'Color', col*0.5, 'MarkerFaceColor', col);
    % 
    %     ii = find(wsp<10);
    %     plot(t(ii), cc(ii).^2, 'k.','DisplayName', '$|U_{10}|<10$ m/s',...
    %         'MarkerSize',14, 'Color', col);
    % else
    %     plot(t, cc.^2, 'kx', 'Color', col, 'LineWidth',1); 
    % end



    % cc = movmedian(cc, 3, 'omitnan'); 
    % cc = movmean(cc, 3, 'omitnan'); 
    % cc(n<0.5) = NaN;
    % plot(t, cc.^2, '-', 'Linewidth',6,'Color', [col 0.2], 'DisplayName', lab)


    x = t; y = cc.^2; 
    
    % bins = datenum(2023,4,10):1.5:datenum(2023,5,15);
    % 
    % sts = boxplot_bybin(gca, y,x, bins, 'Color', col, 'FaceAlpha', 0.2, 'nmin',4, 'EdgeColor', col.*0.5, 'LineWidth',1.5, 'DisplayName', lab, 'Plot', 0);
    % shaded_limits(sts.x0,sts.q2,[sts.q1 sts.q3], col)
    % 


    

    % y95
    % arrayfun(@())
    y = cc.^2; 
    ystd = movstd(y, 5, 'omitnan'); 
    y0 = movmedian(y,4,'omitnan'); 
    y0(isnan(y)) = NaN;
    
    shaded_limits(x,y0,ystd.*0.5.*[-1 1]'+y0, col, 'DisplayName', lab)
    
    

end

leg = cleanLegend(gca, 'northeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
leg.Position(2) = leg.Position(2) + 0.02;
leg.Position(1) = leg.Position(1) + 0.01;



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(end); axes(thisha); hold on; 
yyaxis left; set(gca, 'YColor', 'k');
ylim([0 17]); ylabel('$|U_{10}|$ [m/s]'); yticks([0:5:20]);
x = era5.time; y = mag(era5.u10, era5.v10); y = squeeze(nanmean(y,[1 2])); 
plot(x,y,'k-', 'linewidth',1, 'DisplayName', 'ERA5')
% plot(t, wsp,'kx', 'linewidth',1)
% leg = cleanLegend(gca, 'southeast', 'NumColumns',4', 'Box', 'off', 'FontSize',8);
col = [0.25 0.6 0.6];
yyaxis right; set(gca, 'YColor', col); ylabel('$\theta_{wnd} - \theta_{wv}$', 'Interpreter', 'latex');
yticks([0:90:180]); ylim([0 180]);
% plot(t, abs(dth), 'rx', 'Color', col, 'linewidth',1)
x = era5.time; y_wnd = mod(dirnanmean(atan2d(era5.v10, era5.u10), [1 2]), 360); y_wv = mod(dirnanmean(era5.mwd, [1 2]), 360);
y_wv = met2oc(oc2pol(y_wv));
y =  y_wnd - y_wv; 
y = squeeze(y); y = mod(y, 360); y(y>180) = 360-y(y>180);
plot(x,y,'k:', 'linewidth',2, 'DisplayName', 'ERA5', 'Color', col)

% x = era5.time;  y = mod(dirnanmean(era5.dwps, [1 2]), 360); y = rad2deg(y);
% % y = met2oc(oc2pol(y)); 
% y = squeeze(y);
% plot(x,y,'k:', 'linewidth',2, 'DisplayName', 'ERA5', 'Color', col+0.3)



% set(gca, 'box', 'off')

% 
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% for ai=1:Na-1
%     % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     thisha = ha(ai); axes(thisha); hold on; 
% 
%     cc = A_cc{ai};
% 
%     ii = find([A_cc{3}.^2]'>0.5 & wsp<10);
%     plot(dth(ii), cc(ii).^2, 'k*','DisplayName', '$|U_{10}|>10$ m/s', ...
%         'MarkerSize',10, 'LineWidth',2, 'Color', 'k', 'MarkerFaceColor', 'k');
% 
% 
% end
% 

savejpg(gcf, 'fig_timeseries_correlation', savepath, saveopt)



%% CORRELATION PAIRS (WIND_CURRENT WAVE_CURRENt WIND_WAVE)

% -------------------------------------------------------------------------
t0 = datenum(2023,4,4.5); 
% daaterange = []
% t0 = datenum(2023,4,19.5); 
% % t0 = datenum(2023,4,13.5); 
% t0 = datenum(2023,4,30.5); 
% t0 = datenum(2023,5,8); 
% t0 = datenum(2023,5,13.5); 
[~, ti] = min(abs(t - t0)); t0 = t(ti);
% -------------------------------------------------------------------------
figure(5127); clf; 
% setfigsize(gcf, [500 600])
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right]) 
oldha = tight_subplot(3, Na+1, [0.01 0.015], [0.07 0.05], [0.06 0.03]);
[ha] = reorg_tightsubplot(oldha, [num2cell([(Na+1):Na+1:(Na+1)*3]) num2cell([((Na+1)*2+1):1:(length(oldha)-1)]) {[1:Na [(1:Na)+(Na+1)]]} ] );
% [ha] = reorg_tightsubplot(oldha, [num2cell(1:Na) {[Na+1:Na*2]}]);
dy = 0.15; 
hi = length(ha); 
% ha(hi).Position([2 4]) = ha(hi).Position([2 4]).*([1 1] + [1 -1].*dy);
ha(hi).Position([2 4]) = ha(hi).Position([2 4]) + [1 -1].*dy;
dy = 0.1 - 0.05;
for hi=Na+1:Na*2
    ha(hi).Position([2 4]) = ha(hi).Position([2 4]) + [0 1].*dy;
end
% format_fig4print(ha, 'fontsize',10);

% ha = tight_subplot(length(A),2, [0.05 0.05], [0.05 0.05], [0.05 0.1]);
% dx = 0.4; 
% for hi=1:2:length(ha)
%     ha(hi).Position([1 3]) = ha(hi).Position([1 3]).*([1 1] + [0 1].*dx);
%     ha(hi+1).Position([1 3]) = ha(hi+1).Position([1 3]).*([1 1] + [1 -1].*dx);
% end

% -------------------------------------------------------------------------
format_fig4print(ha, 'fontsize',10);
% setaxes(ha, 'DataAspectRatio', [1 1 1]); 
% setaxes(ha, 'XLim', lonlims); 
% setaxes(ha, 'YLim', latlims); 
setaxes(ha(1:Na), 'YLabel', '$r^2$'); 
setaxes(ha(1:Na), 'YLim', [-1 1]); 
setaxes(ha(1:Na), 'XLabel', '$\theta_{wnd} - \theta_{wv}$'); 
% format_mapticks(ha(1:4), 'x')
setaxes(ha(1:Na), 'XLim', [0 180]); 
% setaxes(ha(1:4), 'XLim', [0 20]); 
setaxes(ha(1:Na), 'XTick', [0 45 90 135 ]); 
% quietaxes(ha(2:length(A)), 'y')
setaxes(ha(end), 'XLim', [datenum(2023,4,3) datenum(2023,5,15)]); 
autodatetick(gca, 'x')
% linkaxes(ha)
% format_mapticks(ha, 'xy');
setaxes(ha(Na+1:end-1), 'XLim', [-1 1]*100); 
setaxes(ha(Na+1:end-1), 'YLim', [-1 1]*100); 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for ai=1:length(A)
    VARx = A{2, ai}; VARy = A{1, ai}; lab = Alabs{ai}; col = cmap(ai,:);
    % 
    msk = ones(size(VARx)); msk(isnan(VARx + VARy)) = NaN;
    VARx = VARx.*msk; VARy = VARy.*msk;

    
    cc = arrayfun(@(i) ccfxn(VARx(:,:,i).*cmsk, VARy(:,:,i).*cmsk) , 1:size(VARx,3))    ;
    p = arrayfun(@(i) pfxn(VARx(:,:,i).*cmsk, VARy(:,:,i).*cmsk) , 1:size(VARx,3))    ;
    
    ccmsk = ones(size(X)); ccmsk(~isin(mag(X,Y), [25 40])) = NaN;
    n = sum(~isnan([VARx + VARy].*ccmsk), [1 2])./sum(~isnan(ccmsk), [1 2]); n = squeeze(n);
    cc(n<0.80) = NaN;


    A_cc{ai} = cc; A_n{ai} = n';
    % cc = abs(cc);
    cc_raw = cc;
    

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(ai); axes(thisha); hold on; 

    
    cc = A_cc{ai};
    % plot(dth, cc.^2, 'k.', 'MarkerSize',10);
    plot(dth, cc, 'k.', 'MarkerSize',10);
    bins = [0:30:180]
    % boxplot_bybin(gca, cc.^2, dth, bins, 'Color', col, 'FaceAlpha', 0.2, 'nmin',4, 'EdgeColor', col.*0.5, 'LineWidth',1.5, 'DisplayName', lab)
    boxplot_bybin(gca, cc, dth, bins, 'Color', col, 'FaceAlpha', 0.2, 'nmin',4, 'EdgeColor', col.*0.5, 'LineWidth',1.5, 'DisplayName', lab)
    leg = cleanLegend(gca, 'northeast', 'box', 'off', 'FontSize',7)
    % leg.ItemTokenSize(2) = leg.ItemTokenSize(2)*0.1;

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(end); axes(thisha); hold on; 
    % plot(t, cc.^2, 'k^', 'Color', col, 'LineWidth',1, 'MarkerFaceColor', col, 'Markersize',4); 
    plot(t, cc, 'k^', 'Color', col, 'LineWidth',1, 'MarkerFaceColor', col, 'Markersize',4); 
    % plot(t, cc.^2, 'kx-', 'Color', col); 

    % cc = movmedian(cc, 3, 'omitnan'); 
    % cc = movmean(cc, 3, 'omitnan'); 
    % cc(n<0.5) = NaN;
    % plot(t, cc.^2, '-', 'Linewidth',6,'Color', [col 0.2], 'DisplayName', lab)
    % x = t; y = movmean(cc.^2, 4, 'omitnan'); ystd = movstd(cc.^2, 4, 'omitnan');
    x = t; y = movmean(cc, 4, 'omitnan'); ystd = movstd(cc.^2, 4, 'omitnan');
    shaded_limits(x,y,ystd.*[-1 1]'+y, col)

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thisha = ha(ai+Na); axes(thisha); hold on; 
    if ai==1; Z = A{2, ai};
    else; Z =  A{1, ai}; end
    
    % rotfldstr = '';
    % run load_vars_into_memory_universal
    % Z = U_wv; 
    % Z = Z - nanmedian(Z, [1 2]); Z = remove_planefit(X,Y,Z); 
    % tt = find(isin(t, t0 + [-1 1]*0.4));
    % Z = Z(:,:,tt); Z = nanmedian(Z,3);
    [~, ti] = min(abs(t - t0)); Z = Z(:,:,ti);

    if isnan(cc(ti)); title('no good'); end

    pcolor(X,Y,Z); shading flat;
    colormap(gca, niceredblue); caxis([-1 1]*1.5)


    
    
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thisha = ha(end); axes(thisha); hold on; 
xline(t0, 'k--')
    % 
    % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % thisha = ha(4); axes(thisha); hold on; 
    % % Z =  A{1, ai};
    % % rotfldstr = '';
    % % run load_vars_into_memory_universal
    % Z = U_wv; Z = Z - nanmedian(Z, [1 2]); Z = remove_planefit(X,Y,Z); 
    % % tt = find(isin(t, t0 + [-1 1]*0.4));
    % % Z = Z(:,:,tt); Z = nanmedian(Z,3);
    % [~, ti] = min(abs(t - t0)); Z = Z(:,:,ti);
    % 
    % pcolor(X,Y,Z); shading flat;
    % colormap(gca, flipud(cmocean('dense')))
    % 
    % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % thisha = ha(5); axes(thisha); hold on; 
    % % Z =  A{1, ai};
    % % rotfldstr = '_rot_wv';
    % % run load_vars_into_memory_universal
    % Z = SWHa; Z = Z - nanmedian(Z, [1 2]); Z = remove_planefit(X,Y,Z); 
    % % tt = find(isin(t, t0 + [-1 1]));
    % % % Z = Z(:,:,tt); Z = nanmedian(Z,3);
    % [~, ti] = min(abs(t - t0)); Z = Z(:,:,ti);
    % pcolor(X,Y,Z); shading flat;
    % colormap(gca, nicejet);
    % caxis([-1 1]*1)
    % 
    % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % thisha = ha(6); axes(thisha); hold on; 
    % % Z =  A{1, ai};
    % % rotfldstr = '_rot_wnd';
    % % run load_vars_into_memory_universal
    % Z = WSPa; Z = Z - nanmedian(Z, [1 2]); Z = remove_planefit(X,Y,Z); 
    % % tt = find(isin(t, t0 + [-1 1]));
    % % % Z = Z(:,:,tt); Z = nanmedian(Z,3);
    % [~, ti] = min(abs(t - t0)); Z = Z(:,:,ti);
    % pcolor(X,Y,Z); shading flat;
    % colormap(gca, flipud(cmocean('speed')))
    % caxis([-1 1]*1.5)
    % 




savejpg(gcf, ['fig_timeseries_correlation_w_snapshots_' datestr(t0, 'yyyymmddTHHMMSS')], savepath, saveopt)




%% Create analysis table

msk = ones(size(X)); msk(mag(X,Y)>75) = NaN;
cmsk = ones(size(X)); cmsk(mag(X,Y)>30) = NaN;
omsk = ones(size(X)); omsk(~isin( mag(X,Y), [60 120]  )) = NaN;
dth = abs(squeeze(dirnanmedian(TH_wnd_outer.*msk - TH_wv_outer.*msk, [1 2])));
wsp = squeeze(nanmedian(WSP.*omsk, [1 2]));
tp = squeeze(nanmedian(TP.*omsk, [1 2]));
swh = squeeze(nanmedian(SWH.*omsk, [1 2]));
umag = squeeze(nanmedian(movmean(mag(U,V),3,3,'omitnan').*msk, [1 2]));

ssta = squeeze(nanmean(SSTa_L4.*cmsk, [1 2])) - squeeze(nanmean(SSTa_L4.*omsk, [1 2]));



T = table();

% Coupling metrics
T.t = t';
T.cc_windwave    = A_cc{3}';
T.n_windwave    = A_n{3}';
T.cc_windcurrent = A_cc{2}';
T.n_windcurrent    = A_n{2}';
T.cc_currentwave = A_cc{1}';
T.n_currentwave    = A_n{1}';
T.cc_sstwind = A_cc{4}';
T.n_sstwind    = A_n{4}';

% Environmental / physical variables
T.swh = swh;
T.wsp = wsp;
T.tp  = tp;
T.dth = dth;
T.u_mag = umag;
T.ssta = ssta;

% save('/Users/alliho/Downloads/T.mat', 'T')
writetable(T, '/Users/alliho/Downloads/eddy_coupling_analysis.csv');



