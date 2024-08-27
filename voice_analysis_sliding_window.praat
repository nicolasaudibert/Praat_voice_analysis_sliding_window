# voice_analysis_sliding_window_paired_sound_textgrid_files.praat
#
# Requires Praat 6.2.04 or newer.
#
# Extract f0, intensity and/or voice perturbation parameters (jitter, shimmer, HNR, ZCR, CPP and CPPS)
# and/or energy in frequency bands 0-1kHz, 1-4kHz, 0-4kHz, 4-8kHz, 0-5kHz, 5-8kHz
# on overlapping frames for each sound file in a given folder.
# Parameters used in the analysis are defined in an external text file.
# Note that some analyses (typically the extraction of intensity and to a lesser extent f0) may require longer frames,
# pay attention to error messages issued by Praat if any and adjust parameters if needed.
# Optionnally, .TextGrid files with matching names may be used to label output frames according to the labels found in the tier set as reference.
# Typical use: voice analysis of sustained vowels.
#
# Author: Nicolas Audibert, LPP UMR7018 CNRS & Sorbonne Nouvelle, February 2022 - last modified August 2024
form voice_analysis_sliding_window
	folder sndFilesFolder .
	word sndFilesExtension .wav
	folder tgFilesFolder_leave_empty_if_irrelevant
	natural targetTierIndex 1
	sentence sndFilesSuffix 
	sentence resultsFileBasename sliding_window_analysis_f0_jitter_shimmer_HNR
	infile parameters_file voice_analysis_sliding_window_default_parameters.txt
	natural windowLengthMs 30
	natural windowOverlapMs 5
	boolean extract_F0 1
	boolean extract_intensity 0
	boolean extract_jitter 1
	boolean extract_shimmer 1
	boolean extract_HNR 1
	boolean extract_ZCR 0
	boolean extract_CPP 0
	boolean extract_CPPS 0
	boolean extract_energy_bands 0
	optionmenu windowPositionStrategy: 1
		button start
		button mid
		button end
	optionmenu windowShape: 1
		button rectangular
		button parabolic
		button Hanning
		button Hamming
		button Gaussian1
		button Gaussian2
		button Gaussian3
		button Gaussian4
		button Gaussian5
		button Kaiser1
		button Kaiser2
endform

# Check input parameters
if windowPositionStrategy$<>"start" and windowPositionStrategy$<>"mid" and windowPositionStrategy$<>"end"
	writeInfoLine: "Unknown windowPositionStrategy parameter value: ", windowPositionStrategy$
else
	# Clean out the info window
	clearinfo

	# Read the external parameters file and create variables listed in the file
	parametersTable = Read Table from tab-separated file: parameters_file$
	nParameters = Get number of rows
	for iParam from 1 to nParameters
		currentParamName$ = Get value: iParam, "Parameter"
		currentParamType$ = Get value: iParam, "Type"
		if currentParamType$ = "Txt"
			currentParamName$ = currentParamName$ + "$"
			'currentParamName$' = Get value: iParam, "Value"
		elsif currentParamType$ = "Num"
			'currentParamName$' = Get value: iParam, "Value"
		else
			appendInfoLine: "File ",parameters_file$ ,", line ", iParam+1," - unknown type for parameter ", currentParamName$
		endif
	endfor
	removeObject: parametersTable

	# Get the list of sound files in the specified folder that match the regular expression
	sndFolderRegex$ = sndFilesFolder$+ "/*" + sndFilesSuffix$ + sndFilesExtension$
	filesList = Create Strings as file list: "fileslist", sndFolderRegex$
	nFiles = Get number of strings

	# Get the results file path from the basename and window duration parameters
	resultsFilePath$ = resultsFileBasename$ + "_" + fixed$(windowLengthMs,0) + "ms_overlap" + fixed$(windowOverlapMs,0) + "ms.txt"

	# Write the header of the extracts reference file
	writeFile: resultsFilePath$, "soundFile", tab$, "windowStartTime", tab$, "windowEndTime", tab$, "intervalLabel"
	if extract_F0
		appendFile: resultsFilePath$, tab$, "F0"
	endif
	if extract_intensity
		appendFile: resultsFilePath$, tab$, "intensity"
	endif
	if extract_jitter
		appendFile: resultsFilePath$, tab$, "jitter_local", tab$, "jitter_local_abs", tab$, "jitter_rap", tab$, "jitter_ppq5", tab$, "jitter_ddp"
	endif
	if extract_shimmer
		appendFile: resultsFilePath$, tab$, "shimmer_local", tab$, "shimmer_local_dB", tab$, "shimmer_apq3", tab$, "shimmer_apq5", tab$, "shimmer_apq11", tab$, "shimmer_dda"
	endif
	if extract_HNR
		appendFile: resultsFilePath$, tab$, "HNR"
	endif
	if extract_ZCR
		appendFile: resultsFilePath$, tab$, "ZCRwindow", tab$, "ZCRperiod"
	endif
	if extract_CPP
		appendFile: resultsFilePath$, tab$, "CPP"
	endif
	if extract_CPPS
		appendFile: resultsFilePath$, tab$, "CPPS"
	endif
	if extract_energy_bands
		appendFile: resultsFilePath$, tab$, "energy0_1kHz", tab$, "energy1_4kHz", tab$, "energy0_4kHz", tab$, "energy4_8kHz", tab$, "energy0_5kHz", tab$, "energy5_8kHz"
	endif
	appendFile: resultsFilePath$, newline$

	# Convert times to seconds
	windowLengthSec = windowLengthMs/1000
	windowOverlapSec = windowOverlapMs/1000

	# Sound pressure-related constant needed to convert energy values in frequency bands (in Pa^2.s-1) to intensity values in dB
	p0 = 2e-5
	p0squared = p0*p0

	# Loop over every sound files
	for iFile from 1 to nFiles
		# Read the sound
		selectObject: filesList
		currentSnd$ = Get string: iFile
		appendInfoLine: "Processing file ", currentSnd$
		currentSndPath$ = sndFilesFolder$+ "/" + currentSnd$

		currentTG$ = currentSnd$ - ".wav" + ".TextGrid"
		currentTGpath$ = tgFilesFolder_leave_empty_if_irrelevant$ + "/" + currentTG$

		# Nowarn to avoid warning display in case mp3 files are used (use wav files instead if possible)
		snd = nowarn Read from file: currentSndPath$

		tg = Read from file: currentTGpath$

		# Get sampling frequency and number of channels (needed for zero padding) and create silence to be added before and after each frame in intensity analysis
		selectObject: snd
		currentSndSamplingFrequency = Get sampling frequency
   		currentSndNchannels = Get number of channels
   		silence = Create Sound from formula: "silence", currentSndNchannels, 0, zeroPaddingDurationIntensityAnalysisMs/1000, currentSndSamplingFrequency, "0"

		# Define times according to selected strategy
		selectObject: snd
		sndDuration = Get total duration
		sndDurationMs = sndDuration * 1000
		nFrames = 1 + floor((sndDurationMs-windowLengthMs)/windowOverlapMs)
		nFramesProcessedBetweenDotDisplay = floor(nFrames/5)
		appendInfo: tab$, sndDuration, " sec, ", nFrames, " frames "
		analyzedSignalTimeSec = windowLengthSec + (nFrames-1)*windowOverlapSec
		if windowPositionStrategy$="start"
			currentFrameStartTime = 0
		elsif windowPositionStrategy$="mid"
			currentFrameStartTime = (sndDuration-analyzedSignalTimeSec)/2
		else
			# windowPositionStrategy$="end"
			currentFrameStartTime = analyzedSignalTimeSec-sndDuration
		endif

		for iFrame from 1 to nFrames
			currentFrameEndTime = currentFrameStartTime + windowLengthSec
			selectObject: snd
			currentFrameSignal = noprogress Extract part: currentFrameStartTime, currentFrameEndTime, windowShape$, 1, "no"

			currentFrameMidTime = currentFrameStartTime + windowLengthSec/2
			selectObject: tg
			currentIntervalIndex = Get interval at time: targetTierIndex, currentFrameMidTime
			currentIntervLabel$ = Get label of interval: targetTierIndex, currentIntervalIndex

			###################
			# Perform acoustic analyses and extract values from computed objects

			if extract_F0 or extract_jitter or extract_shimmer
				# F0
				selectObject: currentFrameSignal
				currentFramePitchCC = noprogress To Pitch (cc): timeStepF0detection, minF0, f0DetectionMaxNumberCandidates, f0DetectionVeryAccurateMode$, f0DetectionSilenceThreshold, f0DetectionVoicingThreshold, f0DetectionOctaveCost, f0DetectionOctaveJumpCost, f0DetectionVoicedUnvoicedCost, maxF0
				selectObject: currentFramePitchCC
				currentFrameF0 = Get mean: 0, 0, "Hertz"

				selectObject: currentFrameSignal
				plusObject: currentFramePitchCC
				currentFramePP = noprogress To PointProcess (cc)
			endif

			if extract_intensity
				# Intensity
				# first, add silence (zero padding) before and after the target frame to make sure it's long enough for analysis
				selectObject: silence
				silenceCopy = Copy: "currentFrameSignalCopy"
				selectObject: currentFrameSignal
				currentFrameSignalCopy = Copy: "currentFrameSignalCopy"
				selectObject: silence
				silenceCopy2 = Copy: "currentFrameSignalCopy2"
				selectObject: silenceCopy
				plusObject: currentFrameSignalCopy
				plusObject: silenceCopy2
				currentFrameSignalZeroPadded = Concatenate
				# extract intensity in zero-padded frame
				selectObject: currentFrameSignalZeroPadded
				currentFrameIntensityObject = To Intensity: minF0, timeStepIntensityComputation, "no"
				selectObject: currentFrameIntensityObject
				currentFrameIntensity = Get mean: 0, 0, averagingMethodIntensity$
				removeObject: silenceCopy, currentFrameSignalCopy, silenceCopy2, currentFrameSignalZeroPadded, currentFrameIntensityObject
			endif

			if extract_shimmer
				# Jitter and shimmer
				selectObject: currentFramePP
				currentFrameJitterLocal = Get jitter (local): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor
				currentFrameJitterLocalAbs = Get jitter (local, absolute): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor
				currentFrameJitterRap = Get jitter (rap): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor
				currentFrameJitterPpq5 = Get jitter (ppq5): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor
				currentFrameJitterDdp = Get jitter (ddp): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor
			endif
			
			if extract_shimmer
				selectObject: snd
				plusObject: currentFramePP
				currentFrameShimmerLocal = Get shimmer (local): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor, shimmerMaximumAmplitudeFactor
				currentFrameShimmerLocal_dB = Get shimmer (local_dB): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor, shimmerMaximumAmplitudeFactor
				currentFrameShimmerApq3 = Get shimmer (apq3): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor, shimmerMaximumAmplitudeFactor
				currentFrameShimmerApq5 = Get shimmer (apq5): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor, shimmerMaximumAmplitudeFactor
				currentFrameShimmerApq11 = Get shimmer (apq11): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor, shimmerMaximumAmplitudeFactor
				currentFrameShimmerDda = Get shimmer (dda): 0, 0, jitterShimmerShortestPeriod, jitterShimmerLongestPeriod, jitterShimmerMaximumPeriodFactor, shimmerMaximumAmplitudeFactor
			endif
			
			if extract_HNR
				# HNR
				selectObject: currentFrameSignal
				currentFrameHarmonicity = noprogress To Harmonicity (cc): timeStepHNRextraction, minF0, silenceTresholdHNRextraction, periodsPerWindowHNRextraction

				selectObject: currentFrameHarmonicity
				currentFrameHNR = Get mean: 0, 0

				removeObject: currentFrameHarmonicity
			endif
			
			if extract_ZCR
				# ZCR
				selectObject: currentFrameSignal
				currentFrameZCR_PP = noprogress To PointProcess (zeroes): 1, includeRaisingPartsInZeroCrossingRateComputation$, includeFallingPartsInZeroCrossingRateComputation$
				currentFrameZCR_TT = Up to TextTier: ""
				currentFrameZCR_TG = Into TextGrid
				removeObject: currentFrameZCR_PP, currentFrameZCR_TT

				selectObject: currentFrameZCR_TG
				nCrossingPoints = Get number of points: 1
				currentFrameZCR = nCrossingPoints / windowLengthSec

				# Get F0 on a longer window
				currentF0WinStartTime = currentFrameMidTime - windowLengthF0computationMilliseconds
				currentF0WinEndTime = currentFrameMidTime + windowLengthF0computationMilliseconds
				if currentF0WinStartTime<0
					currentF0WinStartTime = 0
				endif
				if currentF0WinEndTime>sndDuration
					currentF0WinEndTime = sndDuration
				endif
				selectObject: snd
				currentF0WinSignal = noprogress Extract part: currentF0WinStartTime, currentF0WinEndTime, windowShape$, 1, "yes"
				pitchF0Win  = noprogress To Pitch: timeStepF0detection, minF0, maxF0
				f0ValInF0WinCenter = Get value at time: currentFrameMidTime, "Hertz", "linear"

				nPeriodsInWindow = windowLengthSec * f0ValInF0WinCenter
				currentFramePeriodwiseZCR = nCrossingPoints / nPeriodsInWindow

				# removeObject: currentF0WinSignal, pitchF0Win
				removeObject: currentFrameZCR_TG, currentF0WinSignal, pitchF0Win
			endif
			
			if extract_CPP
				# CPP
				selectObject: currentFrameSignal
				currentFrameSpectrum = noprogress To Spectrum: "yes"
				currentFramePowerCepstrum = noprogress To PowerCepstrum

				selectObject: currentFramePowerCepstrum
				currentFrameCPP = Get peak prominence: minF0forPeakProminenceComputation, maxF0forPeakProminenceComputation, interpolationMethodPeakProminenceComputation$, trendLineQuefrencyMinValuePeakProminenceComputation, trendLineQuefrencyMaxValuePeakProminenceComputation, trendTypePeakProminenceComputation$, fitMethodPeakProminenceComputation$
				removeObject: currentFrameSpectrum, currentFramePowerCepstrum
			endif
			
			if extract_CPPS
				# CPPS
				selectObject: currentFrameSignal
				currentFramePowerCepstrogram = noprogress To PowerCepstrogram: minF0, timeStepPowerCepstrogram, maxFrequencyPowerCepstrogram, preEmphasisStartFrequencyPowerCepstrogram

				selectObject: currentFramePowerCepstrogram
				currentFrameCPPS = Get CPPS: subtractTrendBeforeSmoothingCPPS$, timeAveragingWindowCPPS, quefrencyAveragingWindowCPPS, minF0forPeakProminenceComputation, maxF0forPeakProminenceComputation, peakSearchToleranceFactorCPPS, interpolationMethodPeakProminenceComputation$, trendLineQuefrencyMinValuePeakProminenceComputation, trendLineQuefrencyMaxValuePeakProminenceComputation, trendTypePeakProminenceComputation$, fitMethodPeakProminenceComputation$
				removeObject: currentFramePowerCepstrogram
			endif
			
			if extract_energy_bands
				# Energy in frequency bands
				selectObject: currentFrameSignal
				currentFrameSpectrum = noprogress To Spectrum: "yes"
				selectObject: currentFrameSpectrum
				current_energy_value_Pa2s_0_1kHz = Get band energy: 0, 1000
				currentFrameEnergy0_1kHz = 10*log10(current_energy_value_Pa2s_0_1kHz/(windowLengthSec*p0squared))
				current_energy_value_Pa2s_1_4kHz = Get band energy: 1000, 4000
				currentFrameEnergy1_4kHz = 10*log10(current_energy_value_Pa2s_1_4kHz/(windowLengthSec*p0squared))
				current_energy_value_Pa2s_0_4kHz = Get band energy: 0, 4000
				currentFrameEnergy0_4kHz = 10*log10(current_energy_value_Pa2s_0_4kHz/(windowLengthSec*p0squared))
				current_energy_value_Pa2s_4_8kHz = Get band energy: 4000, 8000
				currentFrameEnergy4_8kHz = 10*log10(current_energy_value_Pa2s_4_8kHz/(windowLengthSec*p0squared))
				current_energy_value_Pa2s_0_5kHz = Get band energy: 0, 5000
				currentFrameEnergy0_5kHz = 10*log10(current_energy_value_Pa2s_0_5kHz/(windowLengthSec*p0squared))
				current_energy_value_Pa2s_5_8kHz = Get band energy: 5000, 8000
				currentFrameEnergy5_8kHz = 10*log10(current_energy_value_Pa2s_5_8kHz/(windowLengthSec*p0squared))
				removeObject: currentFrameSpectrum
			endif

			###################
			# Write extracted values to the results file
			appendFile: resultsFilePath$, currentSnd$, tab$, currentFrameStartTime, tab$, currentFrameEndTime, tab$, currentIntervLabel$
			if extract_F0
				appendFile: resultsFilePath$, tab$, currentFrameF0
			endif
			if extract_intensity
				appendFile: resultsFilePath$, tab$, currentFrameIntensity
			endif
			if extract_jitter
				appendFile: resultsFilePath$, tab$, currentFrameJitterLocal, tab$, currentFrameJitterLocalAbs, tab$, currentFrameJitterRap, tab$, currentFrameJitterPpq5, tab$, currentFrameJitterDdp
			endif
			if extract_shimmer
				appendFile: resultsFilePath$, tab$, currentFrameShimmerLocal, tab$, currentFrameShimmerLocal_dB, tab$, currentFrameShimmerApq3, tab$, currentFrameShimmerApq5, tab$, currentFrameShimmerApq11, tab$, currentFrameShimmerDda
			endif
			if extract_HNR
				appendFile: resultsFilePath$, tab$, currentFrameHNR
			endif
			if extract_ZCR
				appendFile: resultsFilePath$, tab$, currentFrameZCR, tab$, currentFramePeriodwiseZCR
			endif
			if extract_CPP
				appendFile: resultsFilePath$, tab$, currentFrameCPP
			endif
			if extract_CPPS
				appendFile: resultsFilePath$, tab$, currentFrameCPPS
			endif
			if extract_energy_bands
				appendFileLine: resultsFilePath$, tab$, currentFrameEnergy0_1kHz, tab$, currentFrameEnergy1_4kHz, tab$, currentFrameEnergy0_4kHz, tab$, currentFrameEnergy4_8kHz, tab$, currentFrameEnergy0_5kHz, tab$, currentFrameEnergy5_8kHz
			endif
			appendFile: resultsFilePath$, newline$

			# Clean-up: remove temporary objects
			removeObject: currentFrameSignal
			if extract_F0 or extract_jitter or extract_shimmer
				removeObject: currentFramePitchCC, currentFramePP
			endif
			
			# Update frame start time to process next frame
			currentFrameStartTime = currentFrameStartTime + windowOverlapSec
			if iFrame mod nFramesProcessedBetweenDotDisplay = 0
				appendInfo: "."
			endif
		endfor
		appendInfo: newline$

		removeObject: snd, silence, tg
	endfor

	removeObject: filesList
