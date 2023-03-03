GENERAL-INFO-START

	seq-file            rice_7_indv_Omeri.txt
	trace-file		rice_7_indv_Omeri_mig_3.log
	locus-mut-rate          CONST

	mcmc-iterations	  200000
	iterations-per-log  10
	logs-per-line       10


	find-finetunes		TRUE
	#find-finetunes-samples-per-step		10000
	#finetune-coal-time	0.01		
	#finetune-mig-time	0.3		
	#finetune-theta		0.04
	#finetune-mig-rate	0.02
	#finetune-tau		0.0000008
	#finetune-mixing		0.003
#   finetune-locus-rate 0.3
	
	tau-theta-print		10000.0
	tau-theta-alpha		1.0			# for STD/mean ratio of 100%
	tau-theta-beta		10000.0		# for mean of 1e-4

	mig-rate-print		0.001
	mig-rate-alpha		0.002
	mig-rate-beta		0.00001

GENERAL-INFO-END

CURRENT-POPS-START	

	POP-START
		name		A
		samples		Osjap h
	POP-END

	POP-START
		name		B
		samples		Orufi h
	POP-END

	POP-START
		name		C
		samples		Osind h
	POP-END

	POP-START
		name		D
		samples		Oniva h
	POP-END

	POP-START
		name		E
		samples		Oglab h
	POP-END

	POP-START
		name		F
		samples		Obart h
	POP-END

	POP-START
		name		G
		samples		Omeri h
	POP-END


CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			AB
		children		A		B
		#tau-initial	0.000005
		#tau-beta		20000.0	
		#finetune-tau			0.0000008
	POP-END


	POP-START
		name			CD
		children		C		D
		#tau-initial	0.000005
		#tau-beta		20000.0	
		#finetune-tau			0.0000008
	POP-END

	POP-START
		name			EF
		children		E		F
		#tau-initial	0.00005
		#tau-beta		20000.0	
		#finetune-tau			0.00000286
	POP-END

	POP-START
		name			ABCD
		children		AB		CD
		#tau-initial	0.00005
		#tau-beta		20000.0	
		#finetune-tau			0.00000286
	POP-END

	POP-START
		name			ABCDEF
		children		ABCD		EF
		#tau-initial	0.00005
		#tau-beta		20000.0	
		#finetune-tau			0.00000286
	POP-END
	POP-START
		name			root
		children		ABCDEF		G
		#tau-initial	0.00005
		#tau-beta		20000.0	
		#finetune-tau			0.00000286
	POP-END



ANCESTRAL-POPS-END

	
MIG-BANDS-START
	BAND-START
	source	A
	target	C
	BAND-END

	BAND-START
	source	B
	target	C
	BAND-END

	BAND-START
	source	C
	target	A
	BAND-END



MIG-BANDS-END


