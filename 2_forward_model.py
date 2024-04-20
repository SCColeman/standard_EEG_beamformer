# -*- coding: utf-8 -*-
"""
Create forward model with EEG BrainVision files.

@author: Sebastian C. Coleman, ppysc6@nottingham.ac.uk
"""

import os.path as op
import mne
import pandas as pd

### need the following line so 3d plotting works, for some reason
mne.viz.set_3d_options(depth_peeling=False, antialias=False)

#%% set up paths

root = r"R:\DRS-PSR\Seb\standard_EEG_beamformer\example_data"
data_path = op.join(root, "data")
deriv_path = op.join(root, "derivatives")
subjects_dir = op.join(root, "subjects_dir")

subject = "example"
fs_subject = "example"

data_fname = op.join(deriv_path, subject + "-raw.fif")
pos_fname = op.join(data_path, subject + ".pos")

#%% load pre-processed data

raw = mne.io.Raw(data_fname, preload=True)
print(raw.info)

#%% Get FS reconstruction for subject or use fsaverage for quick testing

#fs_dir = mne.datasets.fetch_fsaverage(verbose=True)  # for fsaverage
#subjects_dir = op.dirname(fs_dir)    # for fsaverage

plot_bem_kwargs = dict(
    subject=fs_subject,
    subjects_dir=subjects_dir,
    brain_surfaces="white",
    orientation="coronal",
    slices=[50, 100, 150, 200])

mne.viz.plot_bem(**plot_bem_kwargs)

#%% set proper montage using .pos file, CREATED FROM AN EINSCAN FILE, MODIFY ACCORDINGLY

elec_names = raw.ch_names

# load in pos file to pandas dataframe
df = pd.read_table(pos_fname, names=['point','x','y','z'], delim_whitespace=True)
pos = df.drop(df.index[0]).to_numpy()

# separate pos into fiducials, electrodes and headshape
# 3 fiducial points at the end
# 1 points for each channel (63 channels)
# the rest are headshape points

pos_fids = pos[-3:,1:] / 100  # change units to m for MNE
pos_elec = pos[-3-len(elec_names):-3,1:] / 100

pos_head = pos[0::100,1:] / 100   # downsample Einscan by 100 for speed

# zip electrode names with corresponding coordinates
elec_dict = dict(zip(elec_names, pos_elec))

nas = pos_fids[0,:].astype(float)
lpa = pos_fids[1,:].astype(float)
rpa = pos_fids[2,:].astype(float)
hsp = pos_head.astype(float)

# create head digitisation
digitisation = mne.channels.make_dig_montage(ch_pos=elec_dict, 
                         nasion=nas,
                         lpa=lpa,
                         rpa=rpa,
                         hsp=hsp)

raw.set_montage(digitisation, on_missing="ignore")
raw.save(data_fname, overwrite=True)  # overwrite with correct montage

#%% coregistration

plot_kwargs = dict(
    subject=fs_subject,
    subjects_dir=subjects_dir,
    surfaces="head-dense",
    dig=True,
    show_axes=True,
)

coreg = mne.coreg.Coregistration(raw.info, fs_subject, 
                            subjects_dir=subjects_dir)
mne.viz.plot_alignment(raw.info, trans=coreg.trans, **plot_kwargs)
coreg.fit_fiducials()
coreg.set_grow_hair(5)   # adjust as necessary for particular subject
coreg.fit_icp(20)
mne.viz.plot_alignment(raw.info, trans=coreg.trans, **plot_kwargs)

### UNCOMMENT FOLLOWING IF USING GUI RATHER THAN AUTOMATED ###
#coreg = mne.gui.coregistration(subject=fs_subject, subjects_dir=subjects_dir,
#                               scale_by_distance=False)
#coreg = coreg.coreg
###

trans_fname = subject + "-trans.fif"
coreg.trans.save(op.join(deriv_path, trans_fname), overwrite=True)

#%% compute source space

# can change oct6 to other surface source space
src = mne.setup_source_space(
    fs_subject, spacing="oct6", add_dist=False, subjects_dir=subjects_dir)
src.plot(subjects_dir=subjects_dir)

src_fname = subject + "-src.fif"
src.save(op.join(deriv_path, src_fname), overwrite=True)

#%% three layer bem conduction model

conductivity = (0.3, 0.006, 0.3)
model = mne.make_bem_model(
    subject=fs_subject, ico=4,
    conductivity=conductivity,
    subjects_dir=subjects_dir)
bem = mne.make_bem_solution(model)

bem_fname = subject + "-bem.fif"
mne.write_bem_solution(op.join(deriv_path, bem_fname), bem, overwrite=True)

#%% forward solution

fwd = mne.make_forward_solution(
    raw.info,
    trans=coreg.trans,
    src=src,
    bem=bem,
    meg=False,
    eeg=True
    )

fwd_fname = subject + "-fwd.fif"
mne.write_forward_solution(op.join(deriv_path, fwd_fname), fwd, overwrite=True)
