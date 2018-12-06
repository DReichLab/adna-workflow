all: imports.zip imports2.zip

imports.zip: demultiplex.wdl
	zip imports.zip $^
	
imports2.zip: demultiplex.wdl analysis.wdl analysis_clipping.wdl release_and_pulldown.wdl
	zip imports2.zip $^

