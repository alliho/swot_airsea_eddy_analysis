
%% get rid of shitty ones

xgm = eddygrid.xgm; ygm = eddygrid.ygm; t = eddygrid.time; [X Y] = ndgrid(xgm,ygm); 
[~, fcor] = (calcinertialfreq(eddygrid.center_lat));
edgemsk = NaN(size(X)); edgemsk(isin(mag(X,Y), [75 125])) = 1;
cntrmsk = NaN(size(X)); cntrmsk(isin(mag(X,Y), [0 100])) = 1;
VAR = mag(eddygrid.dVdx,eddygrid.dUdy).*edgemsk; 
VAR = mag(eddygrid.U,eddygrid.V).*edgemsk; 

VAR = squeeze(nanmean(VAR, [1 2]));
bb = find(VAR>3);
% 
% figure(5); clf; hold on; 
% plot(eddygrid.time, VAR)
% plot(eddygrid.time(bb), VAR(bb), 'ro');
fldnms = fieldnames(eddygrid); 
flddims = cellfun(@(fldnm) size(eddygrid.(fldnm)), fldnms, 'Un', 0);
ff = find(cellfun(@length, flddims)==3 & ~contains(fldnms, 'wind') & ~contains(fldnms, 'wave') & ~contains(fldnms, 'swh') & ~contains(fldnms, 'sst')); 


for fi=ff'
    fldnm = fldnms{fi};
    VAR = eddygrid.(fldnm); 
    VAR(:,:,bb) = NaN; 
    eddygrid.(fldnm) = VAR;
end

%% assign variables

global rotfldstr

xgm = eddygrid.xgm; ygm = eddygrid.ygm; t = eddygrid.time; [X Y] = ndgrid(xgm,ygm);
[~, fcor] = calcinertialfreq(eddygrid.center_lat);
fcor = permute(fcor, [1 3 2]);


edgemsk = NaN(size(X)); edgemsk(isin(mag(X,Y), [75 155])) = 1;
cmsk = NaN(size(X)); cmsk(isin(mag(X,Y), [0 100], 'inclusive')) = 1;


LAT = eddygrid.(['latitude' rotfldstr]); LON = eddygrid.(['longitude' rotfldstr]);


U = eddygrid.(['U' rotfldstr]); V = eddygrid.(['V' rotfldstr]);
SSHa = eddygrid.(['ssha_karin_2_proc_patch' rotfldstr]); 
dUdx = eddygrid.(['dUdx' rotfldstr]); dVdx = eddygrid.(['dVdx' rotfldstr]);
dUdy = eddygrid.(['dUdy' rotfldstr]); dVdy = eddygrid.(['dVdy' rotfldstr]);
ZETA = eddygrid.(['zeta' rotfldstr]); dZETAdx = eddygrid.(['dzetadx' rotfldstr]);  dZETAdy = eddygrid.(['dzetady' rotfldstr]); 


%%% WAVES
SWH = eddygrid.(['swh_karin_proc_patch' rotfldstr]);
SWH_model = eddygrid.(['swh_model' rotfldstr]);
% TA_wv = eddygrid.(['mean_wave_period_t02' rotfldstr]); 
% dSWHdx = eddygrid.(['dSWHdx' rotfldstr]); dSWHdy = eddygrid.(['dSWHdy' rotfldstr]);
TP = eddygrid.(['mean_wave_period_t02' rotfldstr]); 
SWH(abs(SWH - nanmean(SWH.*edgemsk,[1 2]))./nanmean(SWH.*edgemsk,[1 2]) > 0.5) = NaN;


%%% WIND
WSP = eddygrid.(['wind_speed_karin_2_proc_patch' rotfldstr]);
WSP_U_model = eddygrid.(['wind_speed_model_u' rotfldstr]);
WSP_V_model = eddygrid.(['wind_speed_model_v' rotfldstr]);
dWSPdx = eddygrid.(['dwinddx' rotfldstr]); dWSPdy = eddygrid.(['dwinddy' rotfldstr]);
dWSP_Vdx = eddygrid.(['dwindvdx' rotfldstr]); dWSP_Vdy = eddygrid.(['dwindvdy' rotfldstr]);
dWSP_Udx = eddygrid.(['dwindudx' rotfldstr]); dWSP_Udy = eddygrid.(['dwindudy' rotfldstr]);

dWSP = mag(dWSPdx,dWSPdy);
WSP_model = mag(WSP_U_model, WSP_V_model);

TH_wnd = atan2d(WSP_V_model,WSP_U_model);
TH_wv = met2oc(oc2pol(eddygrid.(['mean_wave_direction' rotfldstr])));

if isfield(eddygrid, 'sst_ghrsst_l3'); 
    SST = eddygrid.(['sst_ghrsst_l3' rotfldstr]);
elseif isfield(eddygrid, 'sst_proc_patch');
    SST = eddygrid.(['sst_proc_patch' rotfldstr]);
else
    SST = NaN(size(SSHa));
end

if isfield(eddygrid, 'ghrsst_l4_sst_proc_patch'); 
    SST_L4 = eddygrid.(['ghrsst_l4_sst_proc_patch' rotfldstr]);
    dSSTdx_L4 = eddygrid.(['ghrsst_l4_dsstdx' rotfldstr]);
    dSSTdy_L4 = eddygrid.(['ghrsst_l4_dsstdy' rotfldstr]);
else
    SST_L4 = NaN(size(SSHa));
    dSSTdx_L4 = NaN(size(SSHa));
    dSSTdy_L4 = NaN(size(SSHa));

end

if isfield(eddygrid, 'era5_sst'); 
    SST_ERA5 = eddygrid.(['era5_sst' rotfldstr]);
else
    SST_ERA5 = NaN(size(SSHa));
end

dSSTdx = eddygrid.(['dsstdx' rotfldstr]);
dSSTdy = eddygrid.(['dsstdy' rotfldstr]);
dSST = mag(dSSTdx, dSSTdy);

if isfield(eddygrid, 'sig0_karin_2'); 
    SIG0 = eddygrid.(['sig0_karin_2' rotfldstr]);
else
    SIG0 = NaN(size(WSP));
end
if isfield(eddygrid, 'sea_state_bias_cor_2'); 
    SSB = eddygrid.(['sea_state_bias_cor_2' rotfldstr]);
else
    SSB = NaN(size(WSP));
end

if isfield(eddygrid, 'ascat_wind_u'); 
    WSP_U_ASCAT = eddygrid.(['ascat_wind_u' rotfldstr]);
    WSP_V_ASCAT = eddygrid.(['ascat_wind_v' rotfldstr]);
else
    WSP_U_ASCAT = NaN(size(X)); 
    WSP_V_ASCAT = NaN(size(X)); 
end

if isfield(eddygrid, 'duacs_adt'); 
    ADT_DUACS = eddygrid.(['duacs_adt' rotfldstr]);
    U_DUACS = eddygrid.(['duacs_ugos' rotfldstr]);
    V_DUACS = eddygrid.(['duacs_vgos' rotfldstr]);
    ZETA_DUACS = eddygrid.(['duacs_zeta' rotfldstr]);

else
    ADT_DUACS = NaN(size(X)); 
    U_DUACS = NaN(size(X)); 
    V_DUACS = NaN(size(X)); 
    ZETA_DUACS = NaN(size(X)); 
end

if isfield(eddygrid, 'neurost_adt'); 
    ADT_NEUROST = eddygrid.(['neurost_adt' rotfldstr]);
    U_NEUROST = eddygrid.(['neurost_ugos' rotfldstr]);
    V_NEUROST = eddygrid.(['neurost_vgos' rotfldstr]);
    ZETA_NEUROST = eddygrid.(['neurost_zeta' rotfldstr]);
else
    ADT_NEUROST = NaN(size(X)); 
    U_NEUROST = NaN(size(X)); 
    V_NEUROST = NaN(size(X)); 
    ZETA_NEUROST = NaN(size(X)); 
end

% A = {SSHa,U,V,dUdx, dUdy, dVdx,dVdy, SWH, WSP, SST, WSP_U_model, WSP_V_model, dWSP,dSST,dSSTdx,dSSTdy, dWSPdx,dWSPdy , dWSP_Udx, dWSP_Udy, dWSP_Vdx, dWSP_Vdy,SWH_model,CLOUD};
% for ai=1:length(A)
%     VAR = A{ai};
%     nmsk = ones(size(VAR)); nmsk(isnan(VAR)) = NaN;
%     VAR = cell2matn(arrayfun(@(i)    gaussfilt_2d(VAR(:,:,i),eddygrid.drkm, 6), 1:size(VAR,3), 'Un', 0),3);
%     % VAR = movmean(VAR,3,3,'omitnan'); 
%     VAR = VAR.*nmsk;
%     A{ai} = VAR;
% end
% [SSHa,U,V,dUdx, dUdy, dVdx,dVdy, SWH, WSP, SST, WSP_U_model, WSP_V_model, dWSP,dSST,dSSTdx,dSSTdy, dWSPdx,dWSPdy, dWSP_Udx, dWSP_Udy, dWSP_Vdx, dWSP_Vdy,SWH_model,CLOUD] = A{:};


%%% MORE CALCULATED VARIABLES
TH_wnd_outer = ones(size(TH_wnd)).*(dirnanmean(TH_wnd.*edgemsk,[1 2]));
TH_wv_outer = ones(size(TH_wv)).*(dirnanmean(TH_wv.*edgemsk,[1 2]));

[WSP_U, WSP_V] = pol2cart(deg2rad(TH_wnd_outer), WSP);
[U_wnd, V_wnd] = rotate_points(U,V,TH_wnd_outer);
[U_wv, V_wv] = rotate_points(U,V,TH_wv_outer);

[dWSPdx_wnd, dWSPdy_wnd] = rotate_points(dWSPdx,dWSPdy,TH_wnd_outer);
% [dSWHdx_wv, dSWHdy_wv] = rotate_points(dSWHdx,dSWHdy,TH_wv_outer);
[dUdx_wnd, dUdy_wnd] = rotate_points(dUdx,dUdy,TH_wnd_outer);
[dVdx_wnd, dVdy_wnd] = rotate_points(dVdx,dVdy,TH_wnd_outer);
[dSSTdx_wnd, dSSTdy_wnd] = rotate_points(dSSTdx,dSSTdy,TH_wnd_outer);
[dSSTdx_L4_wnd, dSSTdy_L4_wnd] = rotate_points(dSSTdx_L4,dSSTdy_L4,TH_wnd_outer);
[dZETAdx_wnd, dZETAdy_wnd] = rotate_points(dZETAdx,dZETAdy,TH_wnd_outer);


% [U_wnd, V_wnd] = rotate_points(U,V,TH_wnd_outer);
% [U_wv, V_wv] = rotate_points(U,V,TH_wv_outer);
% [dSSTdx_wnd, dSSTdy_wnd] = rotate_points(dSSTdx,dSSTdy,TH_wnd_outer);
% [dWSPdx_wnd, dWSPdy_wnd] = rotate_points(dWSPdx,dWSPdy,TH_wnd_outer);
% [dUdx_wnd, dUdy_wnd] = rotate_points(dUdx,dUdy,TH_wnd_outer);
% [dVdx_wnd, dVdy_wnd] = rotate_points(dVdx,dVdy,TH_wnd_outer);


SIG0a = SIG0 - nanmean(SIG0.*edgemsk,[1 2]);  
SSBa = SSB - nanmean(SSB.*edgemsk,[1 2]);  
SWHa = SWH - nanmean(SWH.*edgemsk,[1 2]);  
WSPa = WSP - nanmean(WSP.*edgemsk,[1 2]);  
WSPa_model = WSP_model - nanmean(WSP_model,[1 2]);  
SSTa = SST - nanmean(SST.*edgemsk,[1 2]);  
SSTa_L4 = SST_L4 - nanmean(SST_L4.*edgemsk,[1 2]);  
SSTa_ERA5 = SST_ERA5 - nanmean(SST_ERA5.*edgemsk,[1 2]);  

relvort = (ZETA)./fcor; 

%% load ERA5 


if ~exist('era5')
    load([datapath 'era5_point.mat'])
end
if ~exist('ERA5')
    load([datapath 'ERA5_grid.mat'])
end

%% load DUACS

if ~exist('DUACS')
    load([datapath 'DUACS.mat'])
end


%% load NEUROST

if ~exist('NEUROST')
    load([datapath 'NEUROST.mat'])
end

