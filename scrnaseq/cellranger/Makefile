all: merge
	./make_ref_genome

merge: tdt
	./append_tdt_to_mm10

tdt: srcs
	python3 get_tdtomato_sequence.py

srcs:
	./get_ref_sources
