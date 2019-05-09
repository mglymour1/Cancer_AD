/***************************************************************************************************************************************************
	Title:	Rate of Memory decline Before and After Cancer diagnosis 
		
	Project overview: 	We evaluate whether the pre-cancer slope and post-cancer slopes are different to the change in memory with each additional
						year of age in the cancer free group.
						We evaluate whether the slope of memory decline changes after cancer diagnosis.
	
	Programmer:	Monica Ospina-Romero
	
	Section 1: Merge data sets, rename some variables for flip of data from wide to long 
	Section 2: Clean relevant variables, apply eligibility and exclusion criteria to obtain the sample, censor observations after missing 
			   interviews(reshaping wide to long)
	Section 3: Clean the exposure variable (canceriwave, cancerever, and date of diagnosis)
	Section 4: Results of the study (follow-up information and Table 1)
	Section 5: Reshape data from wide to long, cleaning some variables for the analysis (time variables related to cancer, canceryet)
	Section 6: Analysis using mixed models. Table 2. and eTable 1. 
	Section 7: Figure 1.
	Section 8: Sensitivity analysis restricting the sample to those who were in more than 3,4 or 5 wave
	
    Note: includes interaction terms with baseline covariates	4/25/2019
	
*****************************************************************************************************************************************************/ 
	global data "/Users/monic/Dropbox"  // you can change this global to your dropbox
	
	cd "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory"
	
	capture log close
	quietly log using cancer_memory_survival.smcl, replace

	version 15
	clear all
	set linesize 80

/**************************************************************
	SECTION 1: Merge data sets, rename some variables.
**************************************************************/

	
	/* import and append data */ 
		* import tracker
		use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/trk2014tr_r.dta"
		sort HHID PN
		
		*Rename interview wave for each wave
			*rename AAGE age1992
			*rename BAGE age1993
			*rename CAGE age1994
			*rename DAGE age1995
			*rename EAGE age1996
			rename FAGE age1998
			rename GAGE age2000
			rename HAGE age2002
			rename JAGE age2004
			rename KAGE age2006
			rename LAGE age2008
			rename MAGE age2010
			rename NAGE age2012
			rename OAGE age2014
		
		*Rename interview month in for each wave 
			rename AIWMONTH iwmo1992
			rename BIWMONTH iwmo1993
			rename CIWMONTH iwmo1994
			rename DIWMONTH iwmo1995
			rename EIWMONTH iwmo1996
			rename FIWMONTH iwmo1998
			rename GIWMONTH iwmo2000
			rename HIWMONTH iwmo2002
			rename JIWMONTH iwmo2004
			rename KIWMONTH iwmo2006
			rename LIWMONTH iwmo2008
			rename MIWMONTH iwmo2010
			rename NIWMONTH iwmo2012
			rename OIWMONTH iwmo2014
		*Rename interview year for each wave 
			rename AIWYEAR iwyr1992
			rename BIWYEAR iwyr1993
			rename CIWYEAR iwyr1994
			rename DIWYEAR iwyr1995
			rename EIWYEAR iwyr1996
			rename FIWYEAR iwyr1998
			rename GIWYEAR iwyr2000
			rename HIWYEAR iwyr2002
			rename JIWYEAR iwyr2004
			rename KIWYEAR iwyr2006
			rename LIWYEAR iwyr2008
			rename MIWYEAR iwyr2010
			rename NIWYEAR iwyr2012
			rename OIWYEAR iwyr2014
		
		* clean alive - rename for reshape of data to long format 
			*rename aalive alive1992
			*rename balive alive1993
			*rename calive alive1994
			*rename dalive alive1995
			*rename ealive alive1996
			rename FALIVE alive1998
			rename GALIVE alive2000
			rename HALIVE alive2002
			rename JALIVE alive2004
			rename KALIVE alive2006
			rename LALIVE alive2008
			rename MALIVE alive2010
			rename NALIVE alive2012
			rename OALIVE alive2014		
			
		keep HHID PN BIRTHYR FIWTYPE FIRSTIW iw* alive* age* RACE GENDER HISPANIC DEGREE OVRESULT OVHHID OVPN
		save trk, replace 

		* import H14C_R to get date of cancer diagnosis reported in 2014
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H14C_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2014	
			rename OC028 yrcancer2014
			replace yrcancer2014=. if yrcancer2014==9998 | yrcancer2014==9999

		* rename month of cancer diagnosis reported in 2014
			rename OC029 mocancer2014
			replace mocancer2014=. if mocancer2014==98 | mocancer2014==99
			
			keep HHID PN yrcancer2014 mocancer2014
			save HRS14, replace	
			
		* import H12C_R to get date of cancer diagnosis reported in 2012
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H12C_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2012	
			rename NC028 yrcancer2012
			replace yrcancer2012=. if yrcancer2012==9998 | yrcancer2012==9999

		* rename month of cancer diagnosis reported in 2012
			rename NC029 mocancer2012
			replace mocancer2012=. if mocancer2012==98 | mocancer2012==99
			
			keep HHID PN yrcancer2012 mocancer2012
			save HRS12, replace	
		
		* import H10C_R to get date of cancer diagnosis reported in 2010
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H10C_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2010	
			rename MC028 yrcancer2010
			replace yrcancer2010=. if yrcancer2010==9998 | yrcancer2010==9999

		* rename month of cancer diagnosis reported in 2010
			rename MC029 mocancer2010
			replace mocancer2010=. if mocancer2010==98 | mocancer2010==99
			
			keep HHID PN yrcancer2010 mocancer2010
			save HRS10, replace	

		* import H08C_R to get date of cancer diagnosis reported in 2008
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H08C_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2008	
			rename LC028 yrcancer2008
			replace yrcancer2008=. if yrcancer2008==9998 | yrcancer2008==9999

		* rename month of cancer diagnosis reported in 2008
			rename LC029 mocancer2008
			replace mocancer2008=. if mocancer2008==98 | mocancer2008==99
			
			keep HHID PN yrcancer2008 mocancer2008
			save HRS08, replace	
		
		* import H06C_R to get date of cancer diagnosis reported in 2006
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H06C_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2006	
			rename KC028 yrcancer2006
			replace yrcancer2006=. if yrcancer2006==9998 | yrcancer2006==9999

		* rename month of cancer diagnosis reported in 2006
			rename KC029 mocancer2006
			replace mocancer2006=. if mocancer2006==98 | mocancer2006==99
			
			keep HHID PN yrcancer2006 mocancer2006
			save HRS06, replace	
		
		* import H04C_R to get date of cancer diagnosis reported in 2004
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H04C_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2004	
			rename JC028 yrcancer2004
			replace yrcancer2004=. if yrcancer2004==9998 | yrcancer2004==9999

		* rename month of cancer diagnosis reported in 2004
			rename JC029 mocancer2004
			replace mocancer2004=. if mocancer2004==98 | mocancer2004==99
			
			keep HHID PN yrcancer2004 mocancer2004
			save HRS04, replace	
		
		* import H02C_R to get date of cancer diagnosis reported in 2002
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H02C_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2002	
			rename HC028 yrcancer2002
			replace yrcancer2002=. if yrcancer2002==9998 | yrcancer2002==9999

		* rename month of cancer diagnosis reported in 2002
			rename HC029 mocancer2002
			replace mocancer2002=. if mocancer2002==98 | mocancer2002==99
			
			keep HHID PN yrcancer2002 mocancer2002
			save HRS02, replace			
		
		* import H00B_R to get date of cancer diagnosis reported in 2000
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H00B_R"
			sort HHID PN
		
		* rename year of cancer diagnosis reported in 2000	
			rename G1274 yrcancer2000
			replace yrcancer2000=. if yrcancer2000==9998 | yrcancer2000==9999
			
		* rename month of cancer diagnosis reported in 2000
			rename G1275 mocancer2000
			replace mocancer2000=. if mocancer2000==98 | mocancer2000==99
			
			keep HHID PN yrcancer2000 mocancer2000
			save HRS00, replace
			
		*****************
			
		* import memory and dementia imputations__Data set provided by Maria
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/dpmemimp_apr2017.dta"
			rename hhid HHID
			rename pn PN
			sort HHID PN
				
		* rename memimp (for reshaping later)
		
					*rename memimp92 memimp1992
					*rename memimp93 memimp1993
					*rename memimp94 memimp1994
					rename memimp95 memimp1995
					rename memimp96 memimp1996
					rename memimp98 memimp1998
					rename memimp00 memimp2000
					rename memimp02 memimp2002
					rename memimp04 memimp2004
					rename memimp06 memimp2006
					rename memimp08 memimp2008
					rename memimp10 memimp2010
					rename memimp12 memimp2012
					rename memimp14 memimp2014
					
		* rename dementpimp 
		
					*rename dementpimp92 dementpimp1992
					*rename dementpimp93 dementpimp1993
					*rename dementpimp94 dementpimp1994
					rename dementpimp95 dementpimp1995
					rename dementpimp96 dementpimp1996
					rename dementpimp98 dementpimp1998
					rename dementpimp00 dementpimp2000
					rename dementpimp02 dementpimp2002
					rename dementpimp04 dementpimp2004
					rename dementpimp06 dementpimp2006
					rename dementpimp08 dementpimp2008
					rename dementpimp10 dementpimp2010
					rename dementpimp12 dementpimp2012
					rename dementpimp14 dementpimp2014
				
				keep HHID PN memimp* dementpimp*
				save memdemimpu, replace
		
		* Import H98A to get Childhood self-rated health (covariate)
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/H98A"
			sort HHID PN
			rename F992 childhe
						
			keep HHID PN childhe
			save HRS98A, replace
			
		* Import cses_measures to get cses_index (covariate)
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/cses_measures"
			rename hhid HHID
			rename pn PN
			sort HHID PN
			keep HHID PN cses_index
			save cses_cam,replace	
				
		* Import RAND: to get cancer diagnosis categorical variable and covariates
			use "$data/CancerDementia/Data_HRS_PublicUse/cancer_memory/rndhrs_p.dta"
			rename hhid HHID
			rename pn PN
			sort HHID PN
								
		* rename cancer diagnosis
			rename r4cancre  canceryet1998
			rename r5cancre  canceryet2000
			rename r6cancre  canceryet2002
			rename r7cancre  canceryet2004
			rename r8cancre  canceryet2006
			rename r9cancre  canceryet2008
			rename r10cancre canceryet2010
			rename r11cancre canceryet2012
			rename r12cancre canceryet2014
			
			keep HHID PN hhidpn canceryet* r4smoken r4drinkn r4vigact r4hibpe ///
				 r4diabe r4lunge r4hearte r4stroke r4arthre r4bmi h4itot h4atotb rabplace raedyrs
				 
			save rand, replace
		
		* append data sets 
			use trk, clear
			
				merge 1:1 HHID PN using HRS00,		gen(m6)
				merge 1:1 HHID PN using HRS02,		gen(m7)
				merge 1:1 HHID PN using HRS04,		gen(m8)
				merge 1:1 HHID PN using HRS06, 		gen(m9) 
				merge 1:1 HHID PN using HRS08, 		gen(m10) 
				merge 1:1 HHID PN using HRS10, 		gen(m11) 
				merge 1:1 HHID PN using HRS12, 		gen(m12) 
				merge 1:1 HHID PN using HRS14, 		gen(m13) 
				merge 1:1 HHID PN using memdemimpu, gen(m14) 
				merge 1:1 HHID PN using HRS98A, 	gen(m15)
				merge 1:1 HHID PN using cses_cam,	gen(m16)
				merge 1:1 HHID PN using rand, 		gen(m17) 
				
			replace BIRTHYR = . if BIRTHYR == 0 
			cou
			save camem_merged_sur, replace

/************************************************************************************************************************
	Section 2: clean variables, apply eligibility and exclusion criteria, censor observations after missing 2 follow-ups
***************************************************************************************************************************/ 	
	use camem_merged_sur, clear

	* Clean the interview date: 
	
		gen iwdate1998 = mdy(iwmo1998, 15 , iwyr1998)
		gen iwdate2000 = mdy(iwmo2000, 15 , iwyr2000)
		gen iwdate2002 = mdy(iwmo2002, 15 , iwyr2002)
		gen iwdate2004 = mdy(iwmo2004, 15 , iwyr2004)
		gen iwdate2006 = mdy(iwmo2006, 15 , iwyr2006)
		gen iwdate2008 = mdy(iwmo2008, 15 , iwyr2008)
		gen iwdate2010 = mdy(iwmo2010, 15 , iwyr2010)
		gen iwdate2012 = mdy(iwmo2012, 15 , iwyr2012)
		gen iwdate2014 = mdy(iwmo2014, 15 , iwyr2014)
		
		
	/***self-reported cancer diagnosis in 1998. This variable is cancer ever vs never in 1998 and was taken from rand dataset.
		I used hiscancer98 to exclude participants with prevalent cancer */
	
			*Recode missing in hiscancer98
			gen hiscancer98 = canceryet1998
			replace hiscancer98 = . if hiscancer98 == .d 
			replace hiscancer98 = . if hiscancer98 == .m
			replace hiscancer98 = . if hiscancer98 == .r
		
		* Identify people with no data in hiscancer98
					gen no_can = 0
					replace no_can = 1 if hiscancer98 == .
					tab no_can
	
		
	*** outcomes: memory score.  I got these imputed values from Maria, based on this paper: doi: 10.1097/WAD.0b013e31826cfe90
				
				* identify people with no outcome data to exclude from analysis 
				* memory 
					gen no_mem = 0
					replace no_mem = 1 if memimp1998 == . 
									
					tab no_mem 	
		
	/* Clean Covariates: Recoded for this study. */

		* Gender is coded 1= male 2=female in Tracker
			rename GENDER male
			recode male 2=0
						
		* Race is coded 0= missing 1=white 2= black 7=other in Tracker
			gen nonwhite = RACE
			replace nonwhite = . if RACE == 0 
			replace nonwhite = 0 if RACE == 1
			replace nonwhite = 1 if RACE == 2
			replace nonwhite = 1 if RACE == 7
			
		* Hispanic is coded 0. Not obtained 1. Hispanic, mexican 2. Hispanic, other 3. Hispanic, type unkwon 5. Non-hispanic
			
			gen hisp = .
			replace hisp = 1 if HISPANIC == 1
			replace hisp = 1 if HISPANIC == 2
			replace hisp = 1 if HISPANIC == 3
			replace hisp = 0 if HISPANIC == 5
			
							
		*Age at baseline is age at wave F (1998)
			gen baseage = age1998 
			replace baseage = . if baseage == 999
			sum baseage , d
			
		/* Education: DEGREE: HRS is coded 0.  No degree 1.  GED 2.  High school diploma 3.  Two year college degree 4.  Four year college degree 
		5.  Master degree  6.  Professional degree (Ph.D., M.D., J.D.)  9.  Degree unknown/Some College
		Used Jennifer Karas in Demography. 2012 Feb; 49(1): 315–336. model 13:lths + edu + edu*lths */
			
			recode DEGREE 5/6=3 9=2 3/4=2 2=1   
			label define edulabel 3 "Master/Professional" 2 "Some/completed college" 1 "High school" 0 "less than high school"
			label values DEGREE edulabel
			
			gen hsdip = 0
			replace hsdip = 1 if DEGREE>= 1
			
			gen coldip = 0
			replace coldip = 1 if DEGREE>= 2
			
			gen masdip = 0 
			replace masdip = 1 if DEGREE== 3
			
		* Years of education: continuous variable but colapse at 17 years +
			rename raedyrs edu  
			replace edu = . if edu== .m
			gen c_edu = edu - 12 //centered variable to 12
			gen lths = 0 // variable less than hs
			replace lths = 1 if DEGREE == 0 
			gen int_edu = c_edu*lths //interaction education and lths 
									
		* Birth place: respondent birth place. 
		/**	    1. Northeast Region: New England Division (ME, NH, VT, MA, RI, CT)
				2. Northeast Region: Middle Atlantic Division (NY, NJ, PA)
				3. Midwest Region: East North Central Division (OH, IN, IL, MI, WI)                         
				4. Midwest Region: West North Central Division (MN, IA, MO, ND,SD, NE, KS)                          
				5. South Region: South Atlantic Division (DE, MD, DC, VA, West VA,NC, SC, GA, FL)                         
				6. South Region: East South Central Division (KY, TN, AL, MS)
				7. South Region: West South Central Division (AR, LA, OK, TX)
				8. West Region: Mountain Division (MT, ID, WY, CO, NM, AZ, UT, NV)
				9. West Region: Pacific Division (WA, OR, CA, AK, HI)
				10. U.S., NA state (without census division info)
				11. Foreign Country: Not in a Census Division (includes U.S. territories)
		**/	
			rename rabplace birthp
			recode birthp 1/2=1 3/4=2 5/7=3 8/9=4 10=. 10= .m 11=5 
			label define birthplabel 1 "Northeast" 2 "Midwest" 3 "South" 4 "West" 5 "Foreign country"
			label values birthp birthplabel
			
			gen SUS = . 
			replace SUS = 1 if birthp == 3
			replace SUS = 0 if birthp == 1
			replace SUS = 0 if birthp == 2
			replace SUS = 0 if birthp == 4
			replace SUS = 0 if birthp == 5		
	
			
		* Childhood self-rate health is coded 1. excellent 2. very good 3. good 4. fair 5. poor 8. DK 9. RF
			
			recode childhe 1/2=0 3=1 4/5=2 8/9=.
			
			label define childhelabel 2 "low" 1 "good" 0 "high" 
			label values childhe childhelabel
			
			
		* Tobaco use current
			rename r4smoken csmoker
			replace csmoker = .  if csmoker == .d
			replace csmoker = .  if csmoker == .r
			tab csmoker, m
			
		* Alcohol use : Classification using the definition from the National Institute on Alcohol Abuse and Alcoholism
			rename r4drinkn etohw

			gen 	etoh_cat = .
			replace etoh_cat = 2 if etohw >=4 & etohw !=. & male ==0
			replace etoh_cat = 2 if etohw >=5 & etohw !=. & male ==1
			replace etoh_cat = 1 if etohw >=1 & etohw <=3 & male ==0
			replace etoh_cat = 1 if etohw >=1 & etohw <=4 & male ==1
			replace etoh_cat = 0 if etohw ==0
			
			label define etohlabel 0 "none" 1 "low-risk" 2 "Binge"
			label values etoh_cat etohlabel
			tab etoh_cat, missing
			
		* Physical activity: R4VIGACT:W4 R Wtr vigorus phys act 3+/wk
			rename r4vigact pacti
			replace pacti = . if pacti == .d
			replace pacti = . if pacti == .r
			tab pacti, m
		
		* comorbidities
		*** Hypertension: R4HIBPE:W4 R ever had high blood pressure
			rename r4hibpe HTN
			replace HTN = . if HTN == .d
			replace HTN = . if HTN == .m
			tab HTN, m
		*** Diabetes: R4DIABE:W4 R ever had diabetes
			rename r4diabe DM2
			replace DM2 = . if DM2 == .d
			replace DM2 = . if DM2 == .m
			tab DM2, m
		*** Lung disease: R4LUNGE:W4 R ever had lung disease
			rename r4lunge lungd
			replace lungd = . if lungd == .d
			replace lungd = . if lungd == .m
			tab lungd, m
		*** Heart disease: R4HEARTE:W4 R ever had heart problems
			rename r4hearte CVD
			replace CVD = . if CVD == .d
			replace CVD = . if CVD == .m
			tab CVD, m
		*** Stroke: R4STROKE:W4 R ever had stroke
			rename r4stroke stroke
			replace stroke = . if stroke == .d
			replace stroke = . if stroke == .m
			tab stroke, m			
		*** Athritis: R4ARTHRE:W4 R ever had arthritis
			rename r4arthre arths
			replace arths = . if arths == .d
			replace arths = . if arths == .m
			replace arths = . if arths == .r
			tab arths, m
		*** BMI
			rename r4bmi BMI
			replace BMI = .  if BMI == .d
			replace BMI = .  if BMI == .m
											
		* Wealth: H4ATOTB
			rename h4atotb wealth 
			gen wea_t = wealth/10000
			
				* Identify participants with missing data (I use this variables to exclude people with missing data)
					gen missing = 0
					replace missing = 1 if male == .
					replace missing = 1 if nonwhite == .
					replace missing = 1 if edu == .
					replace missing = 1 if cses_index == .
					replace missing = 1 if SUS == .
					replace missing = 1 if wea_t == .
					replace missing = 1 if childhe == .
					replace missing = 1 if BMI == .
					replace missing = 1 if pacti == .
					replace missing = 1 if csmoker == .
					replace missing = 1 if etoh_cat == .
					replace missing = 1 if HTN == .
					replace missing = 1 if DM2 == .
					replace missing = 1 if CVD == .
					replace missing = 1 if stroke == .
					replace missing = 1 if lungd == .
					replace missing = 1 if arths == .
		
					tab missing	
				
				
	/** Drop ineligible
	Inclusion criteria: Born before 1948 ( age > 50)
						Completed core interview in 1998
						No history of cancer in 1998
						At least one follow-up interview
	***/				
			
	drop if BIRTHYR == .				// missing data on birthyr
	drop if BIRTHYR > 1948				// drop ages < 50 yrs at 1998
	drop if FIWTYPE != 1				// drop participants with no core interview
	drop if hiscancer98==1
	
	cou // Here there are 18,493 participants 
	
	*** Exclusion criteria: Missing baseline data on memory score, cancer diagnosis in 1998 & covariates
	
	tab1 male RACE hisp DEGREE  if no_mem==0
	tab1 male RACE hisp DEGREE  if no_mem==1  //75% Latinos were missing memory function
	
	*** To build eTable 3. Comparison complete data with missing
	gen anymiss = 0
	replace anymiss = 1 if no_mem == 1
	replace anymiss = 1 if no_can == 1
	replace anymiss = 1 if missing == 1  
	
	tab anymiss // this number should be 2270
	
	tab1 male RACE hisp DEGREE  if anymiss==0
	sum baseage if anymiss==0
	tab1 male RACE hisp DEGREE  if anymiss==1
	sum baseage if anymiss==1
	
		
	*** drop eligible, but with missing data

	drop if no_mem == 1				// no data on memory
	drop if no_can == 1				// no data on cancer in 1998
	drop if missing == 1			// no data on baseline covarites 
	
	cou // 16,223
	
	save camem_cleaning_wide, replace
		
	/***** 	Reshape to drop observations with missing in 2+ waves(censoring at the last follow-up before missing 2 or more interviews)
			and select participants with at least one follow-up ****/
	
	*Reshape data from wide to long
	
	use camem_cleaning_wide, clear
	
	reshape long memimp alive iwdate age canceryet, i(hhidpn) j(year)
	
	order hhidpn year alive memimp age iwdate canceryet 
	sort hhidpn  year
	
	* determine sample size for longitudinal analysis 
		* drop observations of memory score and dementimp before 1998
			drop if year == 1995
			drop if year == 1996
			
		* drop observations for when people are dead
			drop if alive == 5
			drop if alive == 6
		
		* drop observations with no outcome data
			drop if memimp == . 
							
		* drop observations with no follow-up data in canceryet ( I should do a sensitivity analysis for observations with missing data in canceryet)
			drop if canceryet == .d
			drop if canceryet == .r
			drop if canceryet == .m
			drop if canceryet == .t
					
		* indicator for first observation 
			gen obs_num = 1 if hhidpn != hhidpn[_n-1]
			tab obs_num 		// number of people n =   16,223, as expected 
			replace obs_num = obs_num[_n-1] + 1 if hhidpn == hhidpn[_n-1]
			tab obs_num 		// number of people is obs_num = 1 
					
		* Identify nonconsecutive observations. This variable tells the number of years between observations
			gen gap_fup =  0 if obs_num==1
			replace gap_fup = year - year[_n-1] if hhidpn == hhidpn[_n-1]  
			
		* I want to exclude observation after a gap of 6 years of more, then I need to identify observations that come after this gap.
			gen cens_obs = .
			replace cens_obs = 1 if gap_fup>=6
			replace cens_obs = 1 if hhidpn == hhidpn[_n-1] & cens_obs[_n-1]==1

		* drop observations with a gap >= 6 years
			
			drop if cens_obs==1
			
		* determine number of observations per person 
		
		cou 						// number of observations after censoring
		replace obs_num = obs_num[_n-1]+1 if obs_num == . 	// number of obs per person 
		bysort hhidpn: egen obs_tot = max(obs_num) 			// maximum number of observations per person 
		tab obs_tot 
		
			* I used obs_tot to exclude participants with only baseline observation (obs_tot==1)
			
		drop if obs_tot==1
		
		cou		// 96,585 observations
		cou if obs_num==1 //  14,626 number of eligible participants (before excluding those with cancer diagnosis before 1998 reported during follow-up)
	
					
		save camem_cleaning_long_censored, replace
					
	*** Reshape data back to wide format
		use camem_cleaning_long_censored, clear
	
		reshape wide memimp alive iwdate age canceryet obs_num gap_fup, i(hhidpn) j(year)
			
			cou // sample 14,626 (before excluding people who reported cancer diagnosis before 1998 in the waves 2000 to 2014)
			
			save camem_cleaning_wide_censored, replace
			

/************************************************************************************************************************
	Section 3: clean exposure (cancerever, canceryet, date of cancer), exclude participants with cancer before 1998 
***************************************************************************************************************************/ 
		use camem_cleaning_wide_censored, clear
	
	*********Clean exposure: 
	
	*Incident cancer (from 2000 to 2014) for survival analysis
		* Run the loop across each wave from 2000 to 2014
				
		gen canceriwave=. 		// wave in which it was first reported cancer. Used the clean variables from RAND (canceryet*)
		gen canceriyear=.		// self-reported year of diagnosis. Used variables from core interview (yrcancer*)
		gen cancerimonth=.		// self-reported month of diagnosis. Used variables from core interview (mocancer*)
				
			
		forvalues year = 2000(2)2014 { 
		replace canceriwave =`year' if canceryet`year'==1 & canceriwave==. 
		replace canceriyear = yrcancer`year' if canceriwave==`year' & canceriyear==. 
		replace cancerimonth= mocancer`year' if canceriwave==`year'& cancerimonth==.
		}
		
		cou if canceriwave!=. & canceryet1998==0 // 2,293 total cancer cases first reported from 2000 to 2014
		
		tab canceriyear if canceriyear<1998  // 43 cancer cases reported during 2000 to 2014 were diagnosed before 1998
		
		******** Drop participants with cancer before 1998
		
		drop if canceriyear<1998 // exclude cancer diagnosed before 1998
		
		cou  // 14,583 participants included in the sample 
		
		* Flag variable for people with no self-reported cancer year
		
		gen myricancer = .
		replace myricancer = 1 if canceriyear==. & canceriwave!=.  // number of cases with missing year of diagnosis
		tab myricancer if canceriwave!=. , m
		
		
		* Variable that shows the difference in canceriyear and canceriwave
		
		gen diffwayr = canceriwave-canceriyear
		tab diffwayr 
		
		/* flag variable for individuals with missing canceriyear who has at least one reported year of cancer diagnosis in core dataset which
		is not in the same wave as canceryet*  */
		
		gen missyearca = .
				
		forvalues year = 2000(2)2014 {
		replace missyearca = 1 if canceriwave!=. & canceriyear==. & yrcancer`year'!=. & missyearca==. 
		}
		tab missyearca
		
		*** Imputation for those with missing canceriyear 
		
		/*** asumptions to assign the date of diagnosis: 
		1. Canceriwave which is a variable created using RAND data is the variable with more accurate information, 
		I trust the information in yrcancer20* (year of diagnosis from the core dataset) when is from the same wave as the canceriwave 
		(even if canceriwave-canceriyear>2, flag var: diffwayr). Those with missing in yrcancer20* (core data) in the wave of canceriwave 
		will have missing in the canceriyear (myricancer==1)
		2. Those with missing in canceriyear will be replace with midpoint between the date the first reported cancer and the last cancer-free iw.
		3. There should never be a difference of canceriwave-canceriyear<-1
		***/
		
				* First I created a variable that tells me the date when they first reported cancer
		
				gen date_repca = .
		
				forvalues year = 1998(2)2014 {
				replace date_repca = iwdate`year' if canceriwave==`year'
				}
				* Then I created a variable that tells me the interview date of the last wave (when the participant was cancer free)
		
				gen date_lastcafree = .
		
				forvalues year = 1998(2)2014 {
				replace date_lastcafree = iwdate`year' if canceriwave==`year'+2
				}
		
				forvalues year = 1998(2)2014 {
				replace date_lastcafree = iwdate`year' if canceriwave==`year'+4 & date_lastcafree==.
				}
			
					
				****Determine the midpoint between these two dates
					
					**this date is used as the date of diagnosis for those with missing year of diagnosis
		
				gen midpoint = (date_repca+date_lastcafree)/2
				
		* Missing month of diagnosis: 
		
		gen no_moicancer = .
		replace no_moicancer = 1 if canceriwave!=. & cancerimonth==.
		
		tab no_moicancer
		
		* Assigned month 1 (january) to those with missing month
		
		gen impucaimonth = cancerimonth
		replace impucaimonth = 1 if cancerimonth==. & canceriwave!=.
		
		* Date of incident cancer diagnosis (for those with year of diagnosis)
		
		gen dateicancer = mdy(impucaimonth, 1 ,canceriyear)
		
		* Date of incident cancer diagnosis (includes participants with missing values in the date of diagnosis)
		
		gen impudateica = dateicancer
		replace impudateica = midpoint if myricancer==1
		
		*** Age at cancer diagnosis ****
			* Calculate the years since baseline to diagnosis
				gen yrbase = (impudateica - iwdate1998)/365
				gen canceriage= age1998 + yrbase		// age at cancer diagnosis			
		
	* Indicator variable for incident cancer during the entire follow-up. Classifies participants in expose and unexposed.
	gen 	cancerever = 0 
	replace cancerever = 1 if canceriwave!=.
	tab cancerever
	
	* To standardize memory variable (used these variables later (after fliping the data) to standarized all memory assessments)
	egen meanmem98 = mean(memimp1998)
	egen sdmem98 = sd(memimp1998)
	
	* Centering variables for the analysis
		*** BMI
			gen cBMI = BMI-25
			
		* Wealth: H4ATOTB
			egen mean_wea =mean(wea_t)
			gen c_wea = wea_t - mean_wea
			
		*Childhood SES
			egen mean_cses = mean(cses)
			gen c_cses = cses - mean_cses
			
	**** Follow-up time: this variable is to report the average follow-up time.
	
	gen fut1998 = (iwdate1998 - iwdate1998)/365
	gen fut2000 = (iwdate2000 - iwdate1998)/365
	gen fut2002 = (iwdate2002 - iwdate1998)/365
	gen fut2004 = (iwdate2004 - iwdate1998)/365
	gen fut2006 = (iwdate2006 - iwdate1998)/365
	gen fut2008 = (iwdate2008 - iwdate1998)/365
	gen fut2010 = (iwdate2010 - iwdate1998)/365
	gen fut2012 = (iwdate2012 - iwdate1998)/365
	gen fut2014 = (iwdate2014 - iwdate1998)/365
	
	egen maxfut = rmax(fut2000 fut2002 fut2004 fut2006 fut2008 fut2010 fut2012 fut2014)
	
			
	***Time at risk of cancer
	
	gen carisktime=maxfut if cancerever==0
	replace carisktime=yrbase if cancerever==1
	
	******* Crude incidence rate of cancer
	*** Declare dataset as survival data
		
		stset carisktime, failure(cancerever) id(hhidpn)
		
		stptime	
	

	save camem_analysis_wide, replace
	

/************************************************************************************************************************
	Section 4: Results of the study (follow-up information and Table 1)  
***************************************************************************************************************************/ 
	use camem_analysis_wide, clear

	
	***** Years of follow-up 
	sum maxfut, d
	tab cancerever, sum(maxfut)
	
	***** Number of interviews
	sum obs_tot, d
	tab cancerever, sum(obs_tot)

	***** Average age at diagnosis of cancer
	sum canceriage, d
	
	
	****Table 1:
	* Total sample
	
	sum baseage	   , d
	tab male	 
	tab nonwhite  
	sum edu 	   , d
	sum cses_index  , d
	tab SUS    	
	sum wea_t	   , d
	* General health status:
	   
	sum BMI		   , d
	tab pacti	   
	tab csmoker	   
	tab etoh_cat
	tab childhe
	tab HTN 	  
	tab DM2 	  
	tab CVD 	  
	tab stroke 	  
	tab lungd 	  
	tab arths	  
	
	
	*** Participants with no incident cancer during follow-up
	
	* Demographics
	sum baseage	  if cancerever == 0 , d
	tab male	  if cancerever == 0
	tab nonwhite  if cancerever == 0
	sum edu 	  if cancerever == 0 , d
	sum cses_index if cancerever == 0 , d
	tab SUS    if cancerever == 0	
	sum wea_t	  if cancerever == 0 , d
	* General health status:
	sum BMI		  if cancerever == 0 , d
	tab pacti	  if cancerever == 0 
	tab csmoker	  if cancerever == 0 
	tab etoh_cat  if cancerever == 0 
	tab childhe   if cancerever == 0
	tab HTN 	  if cancerever == 0
	tab DM2 	  if cancerever == 0
	tab CVD 	  if cancerever == 0
	tab stroke 	  if cancerever == 0
	tab lungd 	  if cancerever == 0
	tab arths	  if cancerever == 0

	
	*** Participants with incident cancer during follow-up
	
	* Demographics:
	sum baseage	  if cancerever == 1 , d
	tab male	  if cancerever == 1
	tab nonwhite	  if cancerever == 1
	sum edu 	  if cancerever == 1 , d
	sum cses_index if cancerever == 1 , d
	tab SUS		   if cancerever == 1
	sum wea_t	  if cancerever == 1 , d
	* General health status:	
	sum BMI		  if cancerever == 1 , d
	tab pacti	  if cancerever == 1 
	tab csmoker	  if cancerever == 1 
	tab etoh_cat   if cancerever == 1
	tab childhe   if cancerever == 1
	tab HTN 	  if cancerever == 1
	tab DM2 	  if cancerever == 1
	tab CVD 	  if cancerever == 1
	tab stroke 	  if cancerever == 1
	tab lungd 	  if cancerever == 1
	tab arths	  if cancerever == 1

	*Check for overlap with continuous variables
	
	*Baseline age
	twoway ///
 	(kdensity baseage if cancerever==1, bw(0.8) lpattern(solid)) ///
 	(kdensity baseage if cancerever==0, bw(0.8) lpattern(longdash)), ///
 	ytitle("Density") xtitle("Baseline Age (years)") ///
 	legend(order(1 "Cancer" 2 "Comparison")) ///
 	name(bage_overlap, replace)
	
	* Years of education
	twoway ///
 	(kdensity edu if cancerever==1, bw(0.8) lpattern(solid)) ///
 	(kdensity edu if cancerever==0, bw(0.8) lpattern(longdash)), ///
 	ytitle("Density") xtitle("Education (years)") ///
 	legend(order(1 "Cancer" 2 "Comparison")) ///
 	name(bage_overlap, replace)
	
	* Childhood SES score
	
	twoway ///
 	(kdensity cses_index if cancerever==1, bw(0.8) lpattern(solid)) ///
 	(kdensity cses_index if cancerever==0, bw(0.8) lpattern(longdash)), ///
 	ytitle("Density") xtitle("Childhood SES (units)") ///
 	legend(order(1 "Cancer" 2 "Comparison")) ///
 	name(cses_overlap, replace)
	
	* Total household weath
	
	twoway ///
 	(kdensity wea_t if cancerever==1, bw(0.8) lpattern(solid)) ///
 	(kdensity wea_t if cancerever==0, bw(0.8) lpattern(longdash)), ///
 	ytitle("Density") xtitle("Total Household Weath (10,000 dollars)") ///
 	legend(order(1 "Cancer" 2 "Comparison")) ///
 	name(wealth_overlap, replace)
	
	*BMI
	twoway ///
 	(kdensity BMI if cancerever==1, bw(0.8) lpattern(solid)) ///
 	(kdensity BMI if cancerever==0, bw(0.8) lpattern(longdash)), ///
 	ytitle("Density") xtitle("BMI (kg/m2)") ///
 	legend(order(1 "Cancer" 2 "Comparison")) ///
 	name(BMI_overlap, replace)
	
	* To estimate p-values of table 1
	
	foreach x in male nonwhite SUS pacti csmoker etoh_cat childhe HTN DM2 CVD stroke lungd arths {
		tab `x' cancerever, col chi
		}
		
	foreach x in baseage cses_index BMI {
		ttest `x', by(cancerever)
		}
		
	foreach x in edu wea_t {
		ranksum `x', by(cancerever)
		}
	

/**************************************************************************************************
	Section 5: Reshaped data from wide to long for the analysis, creates time variables related to cancer
***************************************************************************************************/ 
		
* Reshape data from wide to long
	use camem_analysis_wide, clear
	
	reshape long memimp alive iwdate age canceryet gap_fup obs_num, i(hhidpn) j(year)
	
	order hhidpn year alive memimp age iwdate canceryet 
	sort hhidpn  year
	
		* drop observations with no exposure or outcome data
				* drop observations with no outcome data
				drop if memimp == . 
							
				* drop observations with no follow-up data in canceryet 
				drop if canceryet == .d
				drop if canceryet == .r
				drop if canceryet == .m
				drop if canceryet == .t
 

	**** Predictor: time with respect to cancer diagnosis. 
	
	gen 	timecancer = (iwdate - impudateica)/365 	
	replace timecancer =  0 if cancerever==0  // cancer-free group is set as 0
		
		** time since cancer diagnosis. This variable is the same as canceryet*timecancer
	
		gen tsca = 0
		replace tsca = timecancer if timecancer>0 
	
	/*** Clean the variable canceryet: Replace canceryet 0 to 1 when timecancer>0. Canceryet is built from the time-varying variable 
	self-reported cancer each wave. Participants might have not reported cancer in the wave they were diagnosed so this code fix that. 
	Example: someone who reported cancer in 2008 and the year of diagnosis was 1999, canceryet should be recode to 1 from 2000.
	**/
	
	replace canceryet = 1 if tsca>0	& tsca!=.
	
		* determine number of observation per person after cancer diagnosis
		
		gen after_obs = 0 if canceryet == 0
		replace after_obs = 1 if canceryet == 1 & canceryet[_n-1] == 0 & hhidpn == hhidpn[_n-1]
		replace after_obs = 1 if canceryet == 1 & year == 1998 & hhidpn != hhidpn[_n-1]
		replace after_obs = after_obs[_n-1] + 1 if canceryet == 1 & hhidpn == hhidpn[_n-1]
		bysort hhidpn: egen after_tot = max(after_obs) if cancerever==1			// maximum number of observations per person 
		
		tab after_tot if obs_num==1
		
	**** Time variable: Age at each interview. Centered at 75
	
		gen cageiw = age-75
		
	**** Age at time of cancer diagnosis
	
	gen 	cagedx = canceriage - 75 if cancerever==1 
	replace cagedx = 0 if cancerever==0  // cancer-free participants is set at 0
	
		
	**** Quadratic age for curvilinear trends
		
	gen cageiw2 = cageiw*cageiw	
	
	**** Dummy variables for childhood self-rated health and etoh categories
	
	gen ch_go = 0 if childhe!=.
	replace ch_go = 1 if childhe == 1
	
	gen ch_lo = 0 if childhe!=.
	replace ch_lo = 1 if childhe == 2	
	
	gen etoh_lo = 0  if etoh_cat!=.
	replace etoh_lo = 1 if etoh_cat==1
	
	gen etoh_bin = 0 if etoh_cat!=.
	replace etoh_bin = 1 if etoh_cat==2
	
	**** Interaction terms baseline covariates and age
	** early life
	gen edu_t = c_edu*cageiw  
	gen nonw_t = nonwhite*cageiw
	gen m_t = male*cageiw
	gen sus_t = SUS*cageiw
	gen ses_t = c_cses*cageiw
	gen c_wea_t = c_wea*cageiw
	gen ch_go_t = ch_go*cageiw
	gen ch_lo_t = ch_lo*cageiw
	
	** baseline covariates	
	gen smok_t = csmoker*cageiw
	gen etohl_t = etoh_lo*cageiw
	gen etohb_t = etoh_bin*cageiw
	gen pact_t = pacti*cageiw
	gen bmi_t = cBMI*cageiw
	gen htn_t = HTN*cageiw
	gen dm_t = DM2*cageiw
	gen lu_t = lungd*cageiw
	gen cvd_t = CVD*cageiw
	gen st_t = stroke*cageiw
	gen art_t = arths*cageiw
	
	**** To standardize memory scores
	gen zmemimp = (memimp - meanmem98)/sdmem98 
	
	**** Lowess memory score by age
	
	twoway ///
	(lowess zmemimp cageiw if cancerever==0 , bw(.2) lcolor(blue))  ///
	(lowess zmemimp cageiw if cancerever==1 & canceryet==0, bw(.2) lcolor(red)) ///
	,   ytitle("Observed Memory Score, SD Units")  xtitle("Age, centered at 75 years")  ///
		legend(order( 1 "No cancer" 2 "Cancer group before diagnosis") row(2) position(1) ring(0))     ///
		plotregion(lcolor(white)) graphregion(margin(large) fcolor(white)) scheme(s2mono)

	twoway ///
	(lowess zmemimp cageiw if cancerever==0 , bw(.2) lcolor(blue))  ///
	(lowess zmemimp cageiw if cancerever==1 & canceryet==1, bw(.2) lcolor(red)) ///
	,   ytitle("Observed Memory Score, SD Units")  xtitle("Age, centered at 75 years")  ///
		legend(order( 1 "No cancer" 2 "Cancer group after diagnosis") row(2) position(1) ring(0))     ///
		plotregion(lcolor(white)) graphregion(margin(large) fcolor(white)) scheme(s2mono)

	save camem_analysis_long, replace
	
/***************************************************
	Section 6: Survival analysis with incident cases
****************************************************/ 
	
	use camem_analysis_long, clear

	* Macros for different confounder combinations
	
	* Early life confounders 
	global conf1 c_edu nonwhite  male SUS c_cses  c_wea ch_go ch_lo //earlylife confounders

	* Baseline confounders related to behavior and health conditions
	global conf3 csmoker etoh_lo etoh_bin pacti cBMI HTN DM2 lungd CVD stroke arths   // General health in 1998
	
	* interactions 1	
	global int1 edu_t nonw_t m_t sus_t ses_t c_wea_t ch_go_t ch_lo_t  // interaction early life confounders
	
	* interactions 2
	global int2 smok_t etohl_t etohb_t pact_t bmi_t htn_t dm_t lu_t cvd_t st_t art_t // interactions with baseline confounders

	*** Modified on April 19/19
	
	* model 0: 
	mixed zmemimp cancerever canceryet timecancer tsca cageiw cagedx cageiw2 || hhidpn: cageiw , reml covariance(unstructured) residuals(indep)
			
			*** for table 2: Model 0. 
			* Memory slope with linear age 
				nlcom  _b[cageiw]*(10)
			* Memory slope with quadratic age
				nlcom _b[cageiw2]*100
			*Difference:
			* Precancer memory slope (linear)
				nlcom _b[timecancer]*10
			* Postcancer memory slope slope (linear)
				nlcom _b[timecancer]*10 + _b[tsca]*10
		
	* model 1: simplest model only sex, race,  and SUS adjusted
	mixed zmemimp cancerever canceryet timecancer tsca cageiw cagedx male nonwhite SUS cageiw2 || hhidpn: cageiw , reml covariance(unstructured) residuals(indep)
	
			*** for table 2: Model 1. 
			* Memory slope with linear age 
				nlcom  _b[cageiw]*(10)
			* Memory slope with quadratic age
				nlcom _b[cageiw2]*100
			*Difference:
			* Precancer memory slope (linear)
				nlcom _b[timecancer]*10
			* Postcancer memory slope slope (linear)
				nlcom _b[timecancer]*10 + _b[tsca]*10
	
	* model 2 : all confounders with education with quadratic time and age
	mixed zmemimp cancerever canceryet timecancer tsca cageiw cagedx cageiw2 $conf1 $conf3 || hhidpn: cageiw , reml covariance(unstructured) residuals(indep)
	
			*** for table 2: Model 2. 
			* Memory slope with linear age 
				nlcom  _b[cageiw]*(10)
			* Memory slope with quadratic age
				nlcom _b[cageiw2]*100
			*Difference:
			* Precancer memory slope (linear)
				nlcom _b[timecancer]*10
			* Postcancer memory slope slope (linear)
				nlcom _b[timecancer]*10 + _b[tsca]*10
								

	*** In the text:
	
		*** Memory function at age 75 for the cancer group
			nlcom _b[_cons] + _b[cancerever]		
		
		*** Memory decline from 65 to 75
			nlcom -1*(_b[cageiw]*(-10) + _b[cageiw2]*100)  //cancer-free person
					// difference in decline with the cancer group when diagnosis was at age 75
			nlcom -1*(_b[timecancer]*(-10))
						
	
		*** Memory decline from 75 to 85
			nlcom _b[cageiw]*(10) + _b[cageiw2]*100  //cancer-free person
					// difference in decline after cancer when cancer dx was at 75
			nlcom _b[timecancer]*(10) + _b[tsca]*(10)
				
		*** Percent difference in slope 
			nlcom  100*((_b[timecancer]*(-10))/(_b[cageiw]*(-10) + _b[cageiw2]*100)) // percent difference in slopes yrs before cancer diagnosis
			nlcom  100*((_b[timecancer]*(10) + _b[tsca]*(10))/(_b[cageiw]*(10) + _b[cageiw2]*100))  // percent difference in slopes after cancer diagnosis
	
		*** For Risk Ratio for cognitive impairment 2 and 10 years after cancer diagnosis
		
		*** Memory function (level) 2 years after cancer diagnosis at age 75
			nlcom _b[_cons] + _b[cageiw]*(2) + _b[cageiw2]*4   // Memory at age 77 for cancer free person
			nlcom _b[_cons] + _b[cancerever] +_b[canceryet] + _b[cageiw]*(2) + _b[cageiw2]*4 + _b[timecancer]*2 + _b[tsca]*2 
			*memory for cancer person at age 77 with cancer diagnosis age 75
			
			**Difference in memory function between the cancer and comparison groups 2 years after diagnosis at age 75
			nlcom _b[cancerever] +_b[canceryet] + _b[timecancer]*2 + _b[tsca]*2
			
		*** Memory function (level) 10 years after cancer diagnosis at age 75
			nlcom _b[_cons] + _b[cageiw]*(10)+ _b[cageiw2]*100  // Memory at age 85 for cancer free person
			nlcom _b[_cons] + _b[cancerever] +_b[canceryet]+ _b[cageiw]*(10) + _b[cageiw2]*100 + _b[timecancer]*10 + _b[tsca]*10
			*memory for cancer person at age 85 with cancer diagnosis age 75
			
			**Difference in memory function between the cancer and comparison groups 10 years after diagnosis at age 75
			nlcom _b[cancerever] +_b[canceryet] + _b[timecancer]*2 + _b[tsca]*2
			
	* model 3 : all confounders and interactions with age
	mixed zmemimp cancerever canceryet timecancer tsca cageiw cagedx cageiw2 $conf1 $conf3 $int1 $int2 || hhidpn: cageiw , reml covariance(unstructured) residuals(indep)
			
				*Difference:
			* Precancer memory slope (linear)
				nlcom _b[timecancer]
			* Postcancer memory slope slope (linear)
				nlcom _b[timecancer] + _b[tsca]
				
/****************************************************************
	Section 7: Figure 1. 
****************************************************************/ 
*** used model 2 for the paper		 		
	*Never cancer prediction
	forval i=1/5 {
	predictnl NC_yu`i'=_b[_cons] + _b[cageiw]*(-`i') + _b[cageiw2]*(`i')*(`i')  , ci(NC_yu_LCL`i' NC_yu_UCL`i')
	predictnl NC_ys`i'=_b[_cons] + _b[cageiw]*`i'    + _b[cageiw2]*(`i'*`i')	, ci(NC_ys_LCL`i' NC_ys_UCL`i')
	}
	predictnl NC_ya0 =_b[_cons] 	 , ci(NC_ya_LCL0 NC_ya_UCL0)

	
	*Prediction before cancer
	forval i=1/5 {
	predictnl CG1_yu`i'=_b[_cons] + _b[cancerever]*1 + _b[cageiw]*(-`i') + _b[timecancer]*(-`i') + _b[cageiw2]*(`i'*`i') , ci(CG1_yu_LCL`i' CG1_yu_UCL`i')
	gen CG1_ys`i' = .
	gen CG1_ys_LCL`i' = .
	gen CG1_ys_UCL`i' = .
	}
	predictnl CG1_ya0 = _b[_cons] + _b[cancerever]*1   , ci(CG1_ya_LCL0 CG1_ya_UCL0)  
	
	*Prediction after cancer 
	forval i=1/5 {
	gen CG2_yu`i'= .
	gen CG2_yu_LCL`i' = .
	gen CG2_yu_UCL`i' = .
	predictnl CG2_ys`i'=_b[_cons] + _b[cancerever]*1 + _b[canceryet]*1 + _b[cageiw]*(`i') + _b[tsca]*(`i') + _b[timecancer]*(`i') + _b[cageiw2]*(`i'*`i')	, ci(CG2_ys_LCL`i' CG2_ys_UCL`i')
	}
	predictnl CG2_ya0 = _b[_cons] + _b[cancerever]*1 + _b[canceryet]*1  , ci(CG2_ya_LCL0 CG2_ya_UCL0)
	
	* The code above should show missing values for some variables that were set as missing
	
	
	* Rename variable to as years
	****No cancer group********************
	rename NC_yu5 NC_70
	rename NC_yu4 NC_71
	rename NC_yu3 NC_72
	rename NC_yu2 NC_73
	rename NC_yu1 NC_74
	rename NC_ya0 NC_75
	rename NC_ys1 NC_76
	rename NC_ys2 NC_77
	rename NC_ys3 NC_78
	rename NC_ys4 NC_79
	rename NC_ys5 NC_80
	*NC_Lower limit
	rename NC_yu_LCL5 NC_LCL70
	rename NC_yu_LCL4 NC_LCL71
	rename NC_yu_LCL3 NC_LCL72
	rename NC_yu_LCL2 NC_LCL73
	rename NC_yu_LCL1 NC_LCL74
	rename NC_ya_LCL0 NC_LCL75
	rename NC_ys_LCL1 NC_LCL76
	rename NC_ys_LCL2 NC_LCL77
	rename NC_ys_LCL3 NC_LCL78
	rename NC_ys_LCL4 NC_LCL79
	rename NC_ys_LCL5 NC_LCL80
	
	*NC_upper limit
	rename NC_yu_UCL5 NC_UCL70
	rename NC_yu_UCL4 NC_UCL71
	rename NC_yu_UCL3 NC_UCL72
	rename NC_yu_UCL2 NC_UCL73
	rename NC_yu_UCL1 NC_UCL74
	rename NC_ya_UCL0 NC_UCL75
	rename NC_ys_UCL1 NC_UCL76
	rename NC_ys_UCL2 NC_UCL77
	rename NC_ys_UCL3 NC_UCL78
	rename NC_ys_UCL4 NC_UCL79
	rename NC_ys_UCL5 NC_UCL80
	****Before cancer*******************
	rename CG1_yu5 CG1_70
	rename CG1_yu4 CG1_71
	rename CG1_yu3 CG1_72
	rename CG1_yu2 CG1_73
	rename CG1_yu1 CG1_74
	rename CG1_ya0 CG1_75
	rename CG1_ys1 CG1_76
	rename CG1_ys2 CG1_77
	rename CG1_ys3 CG1_78
	rename CG1_ys4 CG1_79
	rename CG1_ys5 CG1_80
	*CG1_Lower limit
	rename CG1_yu_LCL5 CG1_LCL70
	rename CG1_yu_LCL4 CG1_LCL71
	rename CG1_yu_LCL3 CG1_LCL72
	rename CG1_yu_LCL2 CG1_LCL73
	rename CG1_yu_LCL1 CG1_LCL74
	rename CG1_ya_LCL0 CG1_LCL75
	rename CG1_ys_LCL1 CG1_LCL76
	rename CG1_ys_LCL2 CG1_LCL77
	rename CG1_ys_LCL3 CG1_LCL78
	rename CG1_ys_LCL4 CG1_LCL79
	rename CG1_ys_LCL5 CG1_LCL80
	*CG1_upper limit
	rename CG1_yu_UCL5 CG1_UCL70
	rename CG1_yu_UCL4 CG1_UCL71
	rename CG1_yu_UCL3 CG1_UCL72
	rename CG1_yu_UCL2 CG1_UCL73
	rename CG1_yu_UCL1 CG1_UCL74
	rename CG1_ya_UCL0 CG1_UCL75
	rename CG1_ys_UCL1 CG1_UCL76
	rename CG1_ys_UCL2 CG1_UCL77
	rename CG1_ys_UCL3 CG1_UCL78
	rename CG1_ys_UCL4 CG1_UCL79
	rename CG1_ys_UCL5 CG1_UCL80
	****After cancer*************************
	rename CG2_yu5 CG2_70
	rename CG2_yu4 CG2_71
	rename CG2_yu3 CG2_72
	rename CG2_yu2 CG2_73
	rename CG2_yu1 CG2_74
	rename CG2_ya0 CG2_75
	rename CG2_ys1 CG2_76
	rename CG2_ys2 CG2_77
	rename CG2_ys3 CG2_78
	rename CG2_ys4 CG2_79
	rename CG2_ys5 CG2_80
	*CG2_Lower limit
	rename CG2_yu_LCL5 CG2_LCL70
	rename CG2_yu_LCL4 CG2_LCL71
	rename CG2_yu_LCL3 CG2_LCL72
	rename CG2_yu_LCL2 CG2_LCL73
	rename CG2_yu_LCL1 CG2_LCL74
	rename CG2_ya_LCL0 CG2_LCL75
	rename CG2_ys_LCL1 CG2_LCL76
	rename CG2_ys_LCL2 CG2_LCL77
	rename CG2_ys_LCL3 CG2_LCL78
	rename CG2_ys_LCL4 CG2_LCL79
	rename CG2_ys_LCL5 CG2_LCL80
	*CG2_upper limit
	rename CG2_yu_UCL5 CG2_UCL70
	rename CG2_yu_UCL4 CG2_UCL71
	rename CG2_yu_UCL3 CG2_UCL72
	rename CG2_yu_UCL2 CG2_UCL73
	rename CG2_yu_UCL1 CG2_UCL74
	rename CG2_ya_UCL0 CG2_UCL75
	rename CG2_ys_UCL1 CG2_UCL76
	rename CG2_ys_UCL2 CG2_UCL77
	rename CG2_ys_UCL3 CG2_UCL78
	rename CG2_ys_UCL4 CG2_UCL79
	rename CG2_ys_UCL5 CG2_UCL80
	
	*** Here it creates a data set for the predicted values for ages 70 to 80 for a group with no cancer (NC) and cancer group (CG) with diagnosis at age 75
		
	collapse (lastnm) NC_* CG1_* CG2_* //This takes the last row/observation (could be any row/observation -- the values are fixed over individuals)
	gen predcancer=1
	save predmem_incidentca.dta, replace
	
	*I use reshape to have one row for each age 

	use predmem_incidentca.dta,clear
	
	reshape long NC_ NC_LCL NC_UCL CG1_ CG1_LCL CG1_UCL CG2_ CG2_LCL CG2_UCL , i(predcancer) j(year) //If I want the graph with age as x-axis, I need to rename the variables
	save predmem_incidentca_figure, replace
	
	use predmem_incidentca_figure, clear
	
	*Figure
	 		 twoway  ///	 
	 (rarea CG1_LCL CG1_UCL year, fcolor(blue*0.2) lcolor(blue*0.2) lwidth(vvvthin))   ///
	 (rarea NC_LCL NC_UCL year, fcolor(orange*0.2) lcolor(orange*0.2) lwidth(vvvthin)) ///
	 (rarea CG2_LCL CG2_UCL year, fcolor(blue*0.2) lcolor(blue*0.2) lwidth(vvvthin))   ///
	 (connected NC_ year, msize(zero) lcolor(orange) lpattern(solid)) (connected CG1_ year, msize(zero) lcolor(blue) lpattern(solid)) ///
	 (connected CG2_ year, msize(zero) lcolor(blue) lpattern(solid)) ///
	  , xline(75, lcolor(black)) ytitle("Estimated Memory Score, SD Units")  xtitle(Age)  ///
		legend(order(1 "Cancer" 2 "No cancer") row(2) position(1) ring(0))     ///
		xlabel(70(1)80) ylabel(-.8(0.2).4) ///
		plotregion(lcolor(white)) graphregion(margin(large) fcolor(white)) scheme(s2mono)
	  
	  graph export "C:\Users\monic\Dropbox\cancer_memory\Figure 2..eps", as(eps) preview(on) replace


/**************************************************************************************************************	  
	  Section 8: Sensitivity analysis restricting the sample to those who were in more than 3,4 or 5 wave
*****************************************************************************************************************/
	use camem_analysis_long, clear
	
	* Macros for different confounder combinations
	
	* Early life confounders 
	global conf1 c_edu nonwhite  male SUS c_cses  c_wea i.childhe //earlylife confounders

	* Baseline confounders related to behavior and health conditions
	global conf3 csmoker i.etoh_cat pacti cBMI HTN DM2 lungd CVD stroke arths   // General health in 1998

	* Model 3: people with at least 3 waves (interviews)
	mixed zmemimp cancerever canceryet timecancer tsca cageiw cagedx cageiw2 $conf1 $conf3 if obs_tot>=3 || hhidpn: cageiw , reml covariance(unstructured) residuals(indep)
		
			*** for etable 3: Model 3. Postcancer rate of memory decline (in years)
							
			* Memory slope with linear age 
				nlcom  _b[cageiw]*(10)
			* Memory slope with quadratic age
				nlcom _b[cageiw2]*100
			*Difference:
			* Precancer memory slope (linear)
				nlcom _b[timecancer]*10
			* Postcancer memory slope slope (linear)
				nlcom _b[timecancer]*10 + _b[tsca]*10
				
	* Model 4: people with at least 4 waves (interviews)
	mixed zmemimp cancerever canceryet timecancer tsca cageiw cagedx cageiw2 $conf1 $conf3 if obs_tot>=4 || hhidpn: cageiw , reml covariance(unstructured) residuals(indep)
		
			*** for etable 3: Model 4. Postcancer rate of memory decline (in years)
			
			* Memory slope with linear age 
				nlcom  _b[cageiw]*(10)
			* Memory slope with quadratic age
				nlcom _b[cageiw2]*100
			*Difference:
			* Precancer memory slope (linear)
				nlcom _b[timecancer]*10
			* Postcancer memory slope slope (linear)
				nlcom _b[timecancer]*10 + _b[tsca]*10

	* Model 5: people with at least 5 waves (interviews)
	mixed zmemimp cancerever canceryet timecancer tsca cageiw cagedx cageiw2 $conf1 $conf3 if obs_tot>=5 || hhidpn: cageiw , reml covariance(unstructured) residuals(indep)
		
			*** for etable 3: Model 5. Postcancer rate of memory decline (in years)
			
			* Memory slope with linear age 
				nlcom  _b[cageiw]*(10)
			* Memory slope with quadratic age
				nlcom _b[cageiw2]*100
			*Difference:
			* Precancer memory slope (linear)
				nlcom _b[timecancer]*10
			* Postcancer memory slope slope (linear)
				nlcom _b[timecancer]*10 + _b[tsca]*10

	
**************************************************************************			
			log close
			exit
		
		
		
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		
		
		
		
		
		
		
		
		
		
		
		
		
		