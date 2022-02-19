function  TF_colormap(time_values,frequency_values,Cpower_values,color_lim,title_string)
%
%this function plot the colormap of power at each time and frequency
%obtained via CWT

delta=8;

min_freq=min(frequency_values); % oppure min(freq)
max_freq=max(frequency_values); % oppure max(freq)
mappa='parula';

pcolor(time_values,frequency_values,10*log10(Cpower_values)); %Cpower_values); %
shading interp
set(gca,'clim',color_lim,'xlim',[-200 750],'ylim',[min_freq max_freq],'yscale','log',...
    'ytick',logspace(log10(min_freq),log10(max_freq),delta),...
    'yticklabel',round(logspace(log10(min_freq),log10(max_freq),delta)*10)/10,...
    'xtick',[-200:100:700],'fontsize',8)
ylabel('Hz, log scale ')
xlabel('ms')
title([title_string,' power (dB)'],'fontsize',8)
colormap(mappa)
colorbar
%grid
set(gca,'layer','top')

end

