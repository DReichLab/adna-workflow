imports.zip: demultiplex.wdl
	zip imports.zip $<

all: imports
