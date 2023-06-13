function [time_percent, gain_value, coh_value, coh_map, F] = cp_wtc(x, y, time, th_wtc)

% Cardiopost analysis on signal x and y

fs = 1/(time(2)-time(1));
[Rsq, period, G, coi] = wtc([time x],[time y],'S0',2,'MaxScale',128);
outsidecoi=zeros(size(Rsq));
for s=1:length(period)
    outsidecoi(s,:)=(period(s)<=coi);
end

n = size(Rsq, 2);
WTC_TH = (th_wtc(:,2))*(ones(1,n));
WTC_TH = Rsq./WTC_TH;

coh_map = WTC_TH;
coh_map(coh_map >= 1) = 1;
coh_map(coh_map < 1) = 0;
coh_map = coh_map.*outsidecoi;

F = 1./period;

[mi, ind_hf_u] = min(abs(period-1/0.5));
[mi, ind_hf_l] = min(abs(period-1/0.15));
[mi, ind_lf_l] = min(abs(period-1/0.07));
[mi, ind_vlf_l] = min(abs(period-1/0.03));

%=========================== HF ====================================
len_hf = 0; wtc_hf = []; wtcsig_hf = []; gain_hf = []; gainsig_hf = []; ind_wtcsig_hf = [];

for i = ind_hf_u:ind_hf_l
    oscoi = outsidecoi(i,:);
    ind_oscoi = find(oscoi == 1);
    
    if ~isempty(ind_oscoi)
        wtc_all = Rsq(i, ind_oscoi);
        wgain = abs(G(i, ind_oscoi));
        wtc_hf = [wtc_hf wtc_all];
        gain_hf = [gain_hf wgain];

        wtc_th = WTC_TH(i, ind_oscoi);
        index = find(wtc_th >= 1);
        ind_wtcsig_hf = [ind_wtcsig_hf len_hf+index];
        wtcsig_hf = [wtcsig_hf wtc_all(index)];
        gainsig_hf = [gainsig_hf wgain(index)];
        
        len_hf = len_hf+length(ind_oscoi);
    end
end
mgain_hf = mean(gain_hf);
tp_hf = length(ind_wtcsig_hf)/len_hf;
if ~isempty(ind_wtcsig_hf)
    mgainsig_hf = mean(gainsig_hf);
    mcohsig_hf = mean(wtcsig_hf); 
else
    mgainsig_hf = NaN;
    mcohsig_hf = NaN;
end

%=========================== LF ====================================
len_lf = 0; wtc_lf = []; wtcsig_lf = []; gain_lf = []; gainsig_lf = []; ind_wtcsig_lf = [];

for i = ind_hf_l:ind_lf_l
    oscoi = outsidecoi(i,:);
    ind_oscoi = find(oscoi == 1);
    
    if ~isempty(ind_oscoi)
        wtc_all = Rsq(i, ind_oscoi);
        wgain = abs(G(i, ind_oscoi));
        wtc_lf = [wtc_lf wtc_all];
        gain_lf = [gain_lf wgain];

        wtc_th = WTC_TH(i, ind_oscoi);
        index = find(wtc_th >= 1);
        ind_wtcsig_lf = [ind_wtcsig_lf len_lf+index];
        wtcsig_lf = [wtcsig_lf wtc_all(index)];
        gainsig_lf = [gainsig_lf wgain(index)];
        
        len_lf = len_lf+length(ind_oscoi);
    end
end
mgain_lf = mean(gain_lf);
tp_lf = length(ind_wtcsig_lf)/len_lf;
if ~isempty(ind_wtcsig_lf)
    mgainsig_lf = mean(gainsig_lf);
    mcohsig_lf = mean(wtcsig_lf); 
else
    mgainsig_lf = NaN;
    mcohsig_lf = NaN; 
end

%=========================== VLF ====================================
len_vlf = 0; wtc_vlf = []; wtcsig_vlf = []; gain_vlf = []; gainsig_vlf = []; ind_wtcsig_vlf = [];

for i = ind_lf_l:ind_vlf_l
    oscoi = outsidecoi(i,:);
    ind_oscoi = find(oscoi == 1);
    
    if ~isempty(ind_oscoi)
        wtc_all = Rsq(i, ind_oscoi);
        wgain = abs(G(i, ind_oscoi));
        wtc_vlf = [wtc_vlf wtc_all];
        gain_vlf = [gain_vlf wgain];

        wtc_th = WTC_TH(i, ind_oscoi);
        index = find(wtc_th >= 1);
        ind_wtcsig_vlf = [ind_wtcsig_vlf len_vlf+index];
        wtcsig_vlf = [wtcsig_vlf wtc_all(index)];
        gainsig_vlf = [gainsig_vlf wgain(index)];
        
        len_vlf = len_vlf+length(ind_oscoi);
    end
end
mgain_vlf = mean(gain_vlf);
tp_vlf = length(ind_wtcsig_vlf)/len_vlf;
if ~isempty(ind_wtcsig_vlf)
    mgainsig_vlf = mean(gainsig_vlf);
    mcohsig_vlf = mean(wtcsig_vlf); 
else
    mgainsig_vlf = NaN;
    mcohsig_vlf = NaN; 
end

%================== Output ================================
gain_value = [mgain_hf mgain_lf mgain_vlf;...
          mgainsig_hf mgainsig_lf mgainsig_vlf];    % mean gain value
coh_value = [mcohsig_hf mcohsig_lf mcohsig_vlf];
time_percent = [tp_hf tp_lf tp_vlf];                % mean time percent of significant WTC
return

