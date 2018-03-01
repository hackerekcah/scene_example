The included C++ program, “segment.cpp” is a software implementation of the auditory segmentation algorithm described in the following paper:

	G. Hu and D.L. Wang, "Auditory segmentation based on onset and offset analysis," IEEE Transactions on Audio, Speech, and Language Processing, in press, 2006.

To use this program, first build the C++ program into an executable file, say "segment.exe". Then type
	
	segment x y

	x is the input data.
	y is the output file describing the resulting segments. In this file, every line corresponds to a time frame, from the first to the last. Every line has 128 integers corresponding to 128 time-frequency units from low frequency to high frequency. The label of each unit indicates the corresponding segment, i.e., the units belonging to one segment have the same label and different segments are indicated by different labels. All the units with label 0 form the background.

We have also included some sample files:

	"input.dat": the mixture of a male utterance and crowd noise in a playground. The sampling frequency is 20 kHz.
	"output.dat": the resulting segments.

We distribute this program freely, but please cite the paper if you have made any use of this program.


By Guoning Hu and DeLiang Wang
September 2006