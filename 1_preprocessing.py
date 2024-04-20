# -*- coding: utf-8 -*-
"""
Pre-process BrainVision EEG files

@author: Sebastian C. Coleman, ppysc6@nottingham.ac.uk
"""

import os.path as op
import numpy as np
import mne

### need the following line so 3d plotting works, for some reason
mne.viz.set_3d_options(depth_peeling=False, antialias=False)

#%% set up paths

root = r"R:\DRS-PSR\Seb\standard_EEG_beamformer\example_data"
data_path = op.join(root, "data")
deriv_path = op.join(root, "derivatives")

subject = "example"
data_fname = op.join(data_path, subject + "_faces_circles")

#%% load data

raw = mne.io.read_raw_brainvision(data_fname + '.vhdr', preload=True)
print(raw.info)

#%% make appropriate montage, this is replaced in 2_forward_model by .pos file

montage = mne.channels.make_standard_montage("easycap-M1")
montage.plot()
raw.set_montage(montage, on_missing="ignore")

#%% downsample

sfreq = 500
raw.resample(sfreq=sfreq)

#%% get events from annotations

events, IDs = mne.events_from_annotations(raw)

#%% drop ECG

raw.drop_channels("ECG") # consider keeping in a separate data object if using 
                         # for automated ICA
                         
#%% broadband filter

raw.filter(l_freq=1, h_freq=45)

#%% drop bad channels, LEFT CLICK IN PLOT

raw.plot()

#%% interpolate bads

raw.interpolate_bads()

#%% set average reference

raw.set_eeg_reference('average', projection=True)

#%% annotate out muscle artefacts for cleaner ICA - we aren't actually removing
# the segments at this stage, we reject bad data later when we epoch. We can 
# therefore be overly conservative here.

threshold_muscle = 6  # z-score
annot_muscle, scores_muscle = mne.preprocessing.annotate_muscle_zscore(
    raw,
    ch_type="eeg",
    threshold=threshold_muscle,
    min_length_good=0.2,
    filter_freq=[35, 45],
)

raw.set_annotations(annot_muscle)

#%% plot psd (this should look pretty clean at this point)

raw.copy().plot_psd(fmax=45, picks='eeg').show()

#%% ICA (LEFT CLICK ON CARDIAC AND BLINKING TIMECOURSES)

ica = mne.preprocessing.ICA(n_components=20)
ica.fit(raw, reject_by_annotation=True)
ica.plot_components()
ica.plot_sources(raw)

#%% remove bad components (THE ONES YOU CLICKED ON IN PREVIOUS PLOT)

ica.apply(raw)

#%% remove annotations, bad segments are removed later when we epoch

raw.annotations.delete(np.arange(len(raw.annotations)))

#%% save out

preproc_fname = subject + "-raw.fif"
events_fname = subject + "-eve.fif"
raw.save(op.join(deriv_path, preproc_fname), overwrite=True)
mne.write_events(op.join(deriv_path, events_fname), events, overwrite=True)
