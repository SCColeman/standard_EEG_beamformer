# -*- coding: utf-8 -*-
"""
LCMV source reconstruction for finger abductions.

@author: Sebastian C. Coleman, ppysc6@nottingham.ac.uk
"""

import os.path as op
import mne
import numpy as np
from matplotlib import pyplot as plt

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
fwd_fname = op.join(deriv_path, subject + "-fwd.fif")
events_fname = op.join(deriv_path, subject + "-eve.fif")

#%% load pre-processed data and forward solution

raw = mne.io.Raw(data_fname, preload=True)
fwd = mne.read_forward_solution(fwd_fname)
events = mne.read_events(events_fname)
print(raw.info)

#%% epoch

event_id = 96  # trigger of interest
tmin, tmax = 0.4, 2.4
epochs = mne.Epochs(
    raw,
    events,
    event_id,
    tmin,
    tmax,
    baseline=(0.4, 0.6),
    preload=True,
    reject_by_annotation=False,
    reject=None)

#%% reject bad epochs by eye (good practice for EEG), LEFT CLICK BAD EPOCHS

epochs.plot("all")

#%% compute covariance from unfiltered data

cov = mne.compute_covariance(epochs)
cov.plot(epochs.info)

#%% compute LCMV weights

filters = mne.beamformer.make_lcmv(
    epochs.info,
    fwd,
    cov,
    reg=0.05, # regularisation, 5% is fairly standard but consider adjusting
    noise_cov=None,
    pick_ori="max-power",
    weight_norm="unit-noise-gain",
    rank=None)

#%% get pseudo-T

# filter epochs
fband = [13, 30]   # beta
epochs_filt = epochs.copy().filter(fband[0], fband[1])

# compute active and control covariance of filtered data        
act_min, act_max = 1.8, 2.2
con_min, con_max = 0.6, 1

active_cov = mne.compute_covariance(epochs_filt, tmin=act_min, tmax=act_max)
control_cov= mne.compute_covariance(epochs_filt, tmin=con_min, tmax=con_max)

stc_active = mne.beamformer.apply_lcmv_cov(active_cov, filters)
stc_control = mne.beamformer.apply_lcmv_cov(control_cov, filters)
pseudoT = (stc_active - stc_control) / (stc_active + stc_control)

# default plot, interactive
pseudoT.plot(src=fwd['src'], subject=fs_subject,
            subjects_dir=subjects_dir, hemi="both")

#%% morph pseudoT to fsaverage for group averaging and better plotting

fname_fsaverage_src = op.join(subjects_dir, "fsaverage", "bem", 
                              "fsaverage-ico-5-src.fif")
src_to = mne.read_source_spaces(fname_fsaverage_src)

morph = mne.compute_source_morph(
    pseudoT,
    subject_from=fs_subject,
    subject_to="fsaverage",
    src_to=src_to,
    subjects_dir=subjects_dir,
)
pseudoT_morphed = morph.apply(pseudoT)

# more publication-ready plot
pseudoT_morphed.plot(src=src_to, subject="fsaverage",
            subjects_dir=subjects_dir,
            surface="pial",
            cortex=[0.54, 0.52, 0.505],  # a pleasing neutral
            background='white',
            views=["dorsal", "med", "lat"],
            view_layout="horizontal",
            size=[1000, 320],
            hemi="both",
            smoothing_steps=10,
            time_viewer=False,
            show_traces=False,
            colorbar=False,
            )

#%% peak timecourse in label

# create LCMV generator
stc_epochs = mne.beamformer.apply_lcmv_epochs(epochs, filters,
                                              return_generator=True)
# get labels from parcellation
parc = "aparc"
labels = mne.read_labels_from_annot(fs_subject, parc=parc, subjects_dir=subjects_dir)
names = [i.name for i in labels]

# make atlas label from combination of other labels (i.e., make a bigger ROI)
label_list = [32, 44, 48]
hemi = "lh"
label_name = hemi + " motor"   # CHANGE THIS TO MATCH LABEL_LIST!!!!!!

# combine vertices and pos from labels
vertices = []
pos = []
for l in label_list:
    vertices.append(labels[l].vertices)
    pos.append(labels[l].pos)
vertices = np.concatenate(vertices, axis=0)   
pos = np.concatenate(pos, axis=0)

# sort vertices and pos
vert_order = np.argsort(vertices)
vertices_ordered = vertices[vert_order]
pos_ordered = pos[vert_order,:]

new_label = mne.Label(vertices_ordered, pos_ordered, hemi=hemi, 
                      name=label_name, subject=fs_subject)

# extract peak from pseudoT in your new label
stc_inlabel = pseudoT.in_label(new_label)
label_peak = stc_inlabel.get_peak(mode="abs", vert_as_index=True)[0]

# extract timecourse of peak
n_epochs = len(epochs_filt)
epoch_len = np.shape(epochs[0])[2]
epoch_peak_data = np.zeros((n_epochs,1,epoch_len))
for s,stc_epoch in enumerate(stc_epochs):
    stc_epoch_label = mne.extract_label_time_course(stc_epoch, new_label, 
                                                    fwd['src'], mode=None)
    epoch_peak_data[s,0,:] = stc_epoch_label[0][label_peak,:]
    
# make source epoch object
ch_names = ["peak"]
ch_types = ["misc"]
source_info = mne.create_info(ch_names=ch_names, sfreq=epochs.info["sfreq"],
                              ch_types=ch_types)
source_epochs = mne.EpochsArray(epoch_peak_data, source_info,
                                tmin=epochs.tmin)

# TFR
baseline = (0.4, 0.6)
freqs = np.arange(1,40)
n_cycles = freqs/2
power = mne.time_frequency.tfr_morlet(source_epochs, freqs=freqs, n_cycles=n_cycles,
                                           use_fft=True, picks="all"
                                           )
power[0].plot(picks="all", baseline=baseline)

# timecourse
source_epochs_filt = source_epochs.copy().filter(fband[0], fband[1], picks="all")
source_epochs_hilb = source_epochs_filt.copy().apply_hilbert(envelope=True, picks="all")
peak_timecourse = source_epochs_hilb.average(picks="all").apply_baseline(baseline)

plt.figure()
plt.plot(peak_timecourse.times, peak_timecourse.get_data()[0], color="black")
plt.ylabel("Oscillatory Power (A.U)")
plt.xlabel("Time (s)")
plt.title(new_label.name)
