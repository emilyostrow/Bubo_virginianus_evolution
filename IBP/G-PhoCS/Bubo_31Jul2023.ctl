GENERAL-INFO-START

	seq-file            Bubo_seqfile.txt
	trace-file          Bubo_31July2023_mcmc.log				
	num-loci            5000
	locus-mut-rate      VAR 1.0
	
	burn-in	          50000
	mcmc-iterations	  500000
	iterations-per-log  10
	logs-per-line       10
#	start-mig           10000

	find-finetunes		FALSE
	finetune-coal-time	0.01		
	finetune-mig-time	0.3		
	finetune-theta		0.04
	finetune-mig-rate	0.02
	finetune-tau		0.0000008
	finetune-mixing		0.003
	finetune-locus-rate 0.3
	
	tau-theta-print		1.0
	tau-theta-alpha		3.0			# for STD/mean ratio of 100%
	tau-theta-beta		10000.0		# for mean of 1e-4

	mig-rate-print		1.0
	mig-rate-alpha		1
	mig-rate-beta		10

GENERAL-INFO-END

CURRENT-POPS-START	

	POP-START
		name		Pop1
		samples		ENO324_Allele-0 h ENO324_Allele-1 h 48729_Allele-0 h 48729_Allele-1 h ENO318_Allele-0 h ENO318_Allele-1 h ENO319_Allele-0 h ENO319_Allele-1 h ENO305_Allele-0 h ENO305_Allele-1 h ENO303_Allele-0 h ENO303_Allele-1 h LHD1893_Allele-0 h LHD1893_Allele-1 h ENO306_Allele-0 h ENO306_Allele-1 h ENO308_Allele-0 h ENO308_Allele-1 h OCGR12152_Allele-0 h OCGR12152_Allele-1 h 
	POP-END

	POP-START
		name		Pop2
		samples		179928_Allele-0 h 179928_Allele-1 h 186079_Allele-0 h 186079_Allele-1 h ENO328_Allele-0 h ENO328_Allele-1 h ENO331_Allele-0 h ENO331_Allele-1 h ENO326_Allele-0 h ENO326_Allele-1 h 23256_Allele-0 h 23256_Allele-1 h 23272_Allele-0 h 23272_Allele-1 h 74120_Allele-0 h 74120_Allele-1 h 86276_Allele-0 h 86276_Allele-1 h 42797_Allele-0 h 42797_Allele-1 h 
			POP-END

	POP-START
		name		Pop3
		samples		236217_Allele-0 h 236217_Allele-1 h 18275_Allele-0 h 18275_Allele-1 h 18278_Allele-0 h 18278_Allele-1 h 23170_Allele-0 h 23170_Allele-1 h 23829_Allele-0 h 23829_Allele-1 h 19811_Allele-0 h 19811_Allele-1 h 22174_Allele-0 h 22174_Allele-1 h 24309_Allele-0 h 24309_Allele-1 h 22173_Allele-0 h 22173_Allele-1 h 24379_Allele-0 h 24379_Allele-1 h 
	POP-END

CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			Pop3_2
		children		Pop3		Pop2
		tau-initial	0.001
#		tau-alpha		3.0	
#		tau-beta		10000.0	
		finetune-tau			0.0000008
#		theta-alpha	5.0
#		theta-beta		1000.0	
	POP-END


	POP-START
		name			root
		children		Pop1	Pop3_2
		tau-initial	0.003
#		tau-alpha	3.0
#		tau-beta		10000.0	
#		theta-alpha	5.0
#		theta-beta		1000.0	
		finetune-tau			0.00000286
	POP-END

ANCESTRAL-POPS-END

MIG-BANDS-START	
	BAND-START		
       source  Pop1
       target  Pop2
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Pop2
       target  Pop1
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Pop2
       target  Pop3
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Pop3
       target  Pop2
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Pop1
       target  Pop3
       mig-rate-print 0.1
	BAND-END

	BAND-START		
       source  Pop3
       target  Pop1
       mig-rate-print 0.1
	BAND-END
	
MIG-BANDS-END

