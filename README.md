# preprocessing
#do single station preprocessing: (as it says)
#              executable file:pre_processing_2017_02_20 in directory cc-love/preparation/src
#              source file    :pre_processing_2017_02_20.f90 in directory cc-love/preparation/src. compile with a single "make", see the makefile for subroutines been called.
#                              This code does single station preprocceing for N,E component together (Lin et al., 2008):
#                               a, time domain normalisation with time window length of 128s
#                               b, frequency domain whitening with window length of 0.01Hz
#                              The output file would be with append of "im" and "re" for imaginary and real part of the waveform respectively.
#              how to run     :pre_processing_2017_02_20 year_begin day_begin year_end day_end
#                              The parameter file "for_pre" is needed, with the form:
#                              "station.list
#                               length_of_data_file_in_seconds multiplication_in_precent component
#                               corner_frequency_fa fb f1 f2
#                               SAC_file_directory"
#                              The station.list is the same as previous one
#                              See the example: for_pre
#              attention      :corner frequency fa and fb here is performed to the raw data, and result is taken as the denominater in the time domain normalization,
#                              in which frequency band the earthquake signal is usually strong ( e.g. 0.02 0.0667)
#                              f1 and f2 is the frequency band for bandpass filter applying to the raw data (as numerater in time domain normalisation)
#                              try for your self.

