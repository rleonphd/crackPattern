all: 
	+$(MAKE) -C 0-D-STIT
	+$(MAKE) -C 1-D-GAUSS
	+$(MAKE) -C 2-D-STIT-ASA
	+$(MAKE) -C 3-D-GAUSS-ASA
	+$(MAKE) -C 4-D-RDMIN-ASA
	+$(MAKE) -C 5-D-RDSSQ-ASA

clean:
	rm 0-D-STIT/cdt
	rm 1-D-GAUSS/cdt
	rm 2-D-STIT-ASA/cdt
	rm 3-D-GAUSS-ASA/cdt
	rm 4-D-RDMIN-ASA/cdt
	rm 5-D-RDSSQ-ASA/cdt
