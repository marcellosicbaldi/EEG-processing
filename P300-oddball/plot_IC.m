function plot_IC(index_ICs,mix_matrix,Y_ICs,PSD_ICs,vector_f,chan_locs,Fs,t1,t2)
% plot_IC(index_ICs,mix_matrix,Y_ICs,PSD_ICs,vector_f,chan_locs,Fs,t1,t2):
% This function plots in the same figure 
% the time-pattern, the PSD and the topographic map of each indexed IC.
% 
% Input Arguments:
% index_ICs = index of the independent components to be plotted.
% mix_matrix = mixing matrix A.
% Y_ICs = the overall set of estimated independent components n x m (number
% of overall ICs x number of samples).
% PSD_ICs = the power spectral density of all independe components n_f x n
% (number of frequencies x overall number of independent components).
% vector_f = vector of frequencies for PSD n_f x 1.
% chan_locs = .locs file or structure containing channel locations. If file
% provide 'filename.locs' within quote marks.
% Fs = sampling frequency.
% t1 = initial time instant (in s) for representing time pattern.
% t2 = final time instant (in s) for representing time pattern.

A=mix_matrix;
good_locs=chan_locs;
f_IC=vector_f;
t=[0:1/Fs:(size(Y_ICs,2)-1)/Fs];
min_frame=find(t==t1);
max_frame=find(t==t2)

for c=1:length(index_ICs) %N_ICs

    figure
    subplot(2,2,1)
    topoplot(A(:,index_ICs(c)), good_locs, 'verbose', 'off', 'numcontour', 8);
    sgtitle(strcat('IC', num2str(index_ICs(c), '%03d')), 'fontsize', 24)
    
    subplot(2,2,2)
    plot(f_IC, PSD_ICs(:,index_ICs(c)),'k','linewidth',1);
    xlim([0.1 60])
    xlabel('Hz')
    ylabel('{\muV}^2/Hz')
    title('PSD')
    grid
    
    subplot(2,2,[3,4])
    plot(t(min_frame:max_frame), Y_ICs(index_ICs(c),min_frame:max_frame), 'k','linewidth',1);
    xlim([t(min_frame) t(max_frame)])
    ylim([-100, 100])
    xlabel('time (s)')
    ylabel('\muV')
    set(gca,'fontsize',11)
    title('time pattern')
    grid
    
%     savefig(strcat('IC_', num2str(index_ICs(c), '%03d'), '.fig'))
end

end

