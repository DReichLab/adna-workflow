all: imports.zip imports2.zip

graphs: analysis.pdf demultiplex.pdf sample_merge.pdf

imports.zip: demultiplex.wdl
	zip imports.zip $^
	
imports2.zip: demultiplex.wdl analysis.wdl analysis_clipping.wdl release_and_pulldown.wdl
	zip imports2.zip $^

analysis.pdf: analysis.wdl
	java -jar womtool-36.jar graph $^ |dot -Tpdf -o $@ 

demultiplex.pdf: demultiplex.wdl
	java -jar womtool-36.jar graph $^ |dot -Tpdf -o $@

sample_merge.pdf: sample_merge.wdl
	java -jar womtool-36.jar graph $^ |dot -Tpdf -o $@
