# Praat_voice_analysis_sliding_window
Perform various voice analyses on a set of sound files, using a sliding window of fixed duration.
Typical use: voice analysis of sustained vowels.

This script extracts selected features among f0, intensity and/or voice perturbation parameters (jitter, shimmer, HNR, ZCR, CPP and/or CPPS) and/or energy in frequency bands 0-1kHz, 1-4kHz, 0-4kHz, 4-8kHz, 0-5kHz, 5-8kHz on overlapping frames (duration and overlap in milliseconds set by the user) for each sound file in a given folder.
Note that the computation of CPP and CPPS requires much more computation power compared to other features, include measures of cepstral peak prominence only if needed.
Parameters used in the analysis are defined in an external text file.

Note that some analyses (typically the extraction of intensity and to a lesser extent f0) may require longer frames, pay attention to error messages issued by Praat if any and adjust parameters if needed.

Optionnally, .TextGrid files with matching names may be used to label output frames for further filtering of comparison, according to the labels in the tier set as reference. If not specified, .TextGrid files are loaded from the same folder as sound files.
