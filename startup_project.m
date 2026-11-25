%% startup_project.m
% Interactive startup script for MD_ZT (Cone Micro-Motion Measurement & Identification).

addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
md_zt_paths = ensure_md_zt_paths(true);
project_root = md_zt_paths.projectRoot;
rv_isar_path_root = md_zt_paths.rvIsarPathRoot;
matlab_research_lib_root = md_zt_paths.matlabResearchLibRoot;
matlab_research_lib_dir = md_zt_paths.matlabResearchLibDir;
