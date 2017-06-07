function plot_wholebrain_roi_sim_results(subj_info, session_num, freq, snr, varargin)
% PLOT_WHOLEBRAIN_ROI_SIM_RESULTS  Plot free energy and t stastics for all
% simulations
%
% Use as
%   plot_wholebrain_roi_sim_results(subjects(1), 1, [10 30], -20)
% where the first argument is the subject info structure (from create_subjects),
% the second is the session number, the third is the frequency range, and
% the fourth is the SNR (db).
% 
%   plot_wholebrain_roi_sim_results(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or interger - moment of simulated
%    dipole
%    * surf_dir - directory containing subject surfaces
%    * sim_patch_size - 0 (default) or interger - simulated patch size
%    * reconstruct_patch_size - 0 (default) or interger - reconstruction patch size

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'surf_dir', 'd:\pred_coding\surf',...
    'sim_patch_size',0, 'reconstruct_patch_size',0);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
Nmeth=length(methodnames);
% Original and downsampled white matter surface
orig_white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_white.hires.deformed.surf.gii');

% Original and downsampled pial surface
orig_pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_pial.hires.deformed.surf.gii');

allmeshes=strvcat(white_mesh,pial_mesh);
Nmesh=size(allmeshes,1);

wholeBrainPialWhiteF=zeros(Nmeth,Nmesh,params.nsims);
% Load free energy results
if params.sim_patch_size==0 && params.reconstruct_patch_size==0
    fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',...
        freq(1),freq(2),snr,params.dipole_moment);
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
        num2str(session_num), fname);
else
    fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',...
        freq(1),freq(2),snr,params.dipole_moment, params.sim_patch_size,...
        params.reconstruct_patch_size);
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
        num2str(session_num), fname);
end
load(data_file);

for methind=1:Nmeth,       
    for simmeshind=1:Nmesh,           
        % F reconstructed on piwl - reconstructed on white
        % num simulations x number of folds
        pialF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind));
        whiteF=squeeze(allcrossF(simmeshind,1:params.nsims,1,methind));
        wholeBrainPialWhiteF(methind,simmeshind,:)=pialF-whiteF;
    end
end

roiWhitePialT=zeros(Nmeth,Nmesh*params.nsims);
data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id,...
    num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1),...
    freq(2), snr, params.dipole_moment));
for methind=1:Nmeth,    
    method=methodnames{methind};
    roiWhitePialT(methind,:)=get_wmpial_t(data_dir, method, params.nsims, allmeshes(2,:), ...
        allmeshes(1,:), orig_pial_mesh, orig_white_mesh, 'recompute_trials',false);
end

figure('Position',[1 1 1000 1800]);
mesh_idx=zeros(Nmeth,Nmesh);
mesh_idx(1,1)=3;
mesh_idx(1,2)=1;
mesh_idx(2,1)=7;
mesh_idx(2,2)=5;
mesh_idx(3,1)=11;
mesh_idx(3,2)=9;
mesh_idx(4,1)=15;
mesh_idx(4,2)=13;
for methind=1:Nmeth,       
    for simmeshind=1:Nmesh,    
        subplot(Nmeth*Nmesh,2,mesh_idx(methind,simmeshind));
        [path,file,ext]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        hold on
        for simind=1:params.nsims
            color='b';
            if wholeBrainPialWhiteF(methind,simmeshind,simind)<0
                color='r';
            end
            bar(simind,wholeBrainPialWhiteF(methind,simmeshind,simind),color,'EdgeColor','none');
        end
        plot([0 params.nsims+1],[3 3],'k--');
        plot([0 params.nsims+1],[-3 -3],'k--');
        xlim([0 params.nsims+1]);
        xlabel('Simulation')
        ylabel('Free energy diff (pial-white)');
        title(sprintf('Free energy, %s, %s',methodnames{methind},simmeshname));        
    end

end

mesh_idx=zeros(Nmeth,Nmesh);
mesh_idx(1,1)=4;
mesh_idx(1,2)=2;
mesh_idx(2,1)=8;
mesh_idx(2,2)=6;
mesh_idx(3,1)=12;
mesh_idx(3,2)=10;
mesh_idx(4,1)=16;
mesh_idx(4,2)=14;

dof=514;
alpha=1.0-(0.05/2);
t_thresh=tinv(alpha, dof);

for methind=1:Nmeth,    
    method=methodnames{methind};
    
    % For each simulated mesh
    for simmeshind=1:Nmesh,
        [path,file,ext]=fileparts(allmeshes(simmeshind,:));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        subplot(Nmeth*Nmesh,2,mesh_idx(methind,simmeshind));
        hold on
        for simind=1:params.nsims
            color='b';
            if roiWhitePialT(methind,(simmeshind-1)*params.nsims+simind)<0
                color='r';
            end
            bar(simind,roiWhitePialT(methind,(simmeshind-1)*params.nsims+simind),color,'EdgeColor','none');
        end
        hold on
        plot([0 params.nsims+1], [t_thresh t_thresh],'k--');
        plot([0 params.nsims+1], [-t_thresh -t_thresh],'k--');
        xlim([0 params.nsims+1]);
        %ylim([-40 40]);
        xlabel('Simulation')
        ylabel('Mean t (pial-white)');
        title(sprintf('t Stat, %s, %s',methodnames{methind},simmeshname)); 
    end
end
