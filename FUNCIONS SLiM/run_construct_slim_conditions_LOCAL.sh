#### RUN 15 DIFFERENT CONDITIONS. OBTAIN SFS IN SYN AND NONSYN. CALCULATE MKTasymptotic. ####
#### GENERAK CONDITIONS:
#### Two populations: ONE target plus ONE outgroup.
#### 1.The initial population run for 5Ne generations to achieve some equilibrium mutation-selection-drift
#### 2.Split target and outgroup for 10Ne generations
#### 3.Possible change in number of individuals in target population after 5Ne+10Ne

for((iter=1; iter<=100; i++))
do
	#fixed paraneters
	Ne=400; L=500000; ngenes=100;
	mut_rate=1e-6;
	ind_sample_size=25; out_sample_size=1;

	# CONDITION 0:
	#Neutral. No change Ne.
	FILEOUT="'./00_slim_SFS_SNM.txt'"
	Neb=400;
	rec_rate=1e-4;
	rate_ben=0; s_backg_ben=0;
	rate_del=0; s_backg_del=0;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\> slim_time_memory.txt \& > ./run_slim_conditions.sh

	# CONDITION 1:
	# BACKGROUND SELECTION. No beneficial selection. No change Ne. 
	FILEOUT="'./01_slim_SFS_BGS.txt'"
	Neb=400;
	rec_rate=1e-4
	rate_ben=0; s_backg_ben=0;
	rate_del=8; s_backg_del=-0.05;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\>\>  slim_time_memory.txt \& >> ./run_slim_conditions.sh

	# CONDITION 2:
	#No background selection. BENEFICIAL SELECTION. No change Ne. 
	FILEOUT="'./02_slim_SFS_PSL.txt'"
	Neb=400;
	rec_rate=1e-4
	rate_ben=0.05; s_backg_ben=0.02;
	rate_del=0; s_backg_del=0;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\>\>  slim_time_memory.txt \& >> ./run_slim_conditions.sh

	# CONDITION 3:
	# BACKGROUND SELECTION. BENEFICIAL SELECTION. POPULATION REDUCTION. 
	FILEOUT="'./03_slim_SFS_BGS_PSL_RED.txt'"
	Neb=100;
	rec_rate=1e-4
	rate_ben=0.5; s_backg_ben=0.02;
	rate_del=8; s_backg_del=-0.05;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\>\>  slim_time_memory.txt \& >> ./run_slim_conditions.sh

	# CONDITION 4:
	# BACKGROUND SELECTION. BENEFICIAL SELECTION. POPULATION EXPANSION. 
	FILEOUT="'./04_slim_SFS_BGS_PSL_EXP.txt'"
	Neb=1000;
	rec_rate=1e-4
	rate_ben=0.5; s_backg_ben=0.02;
	rate_del=8; s_backg_del=-0.05;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\>\>  slim_time_memory.txt \& >> ./run_slim_conditions.sh

	# CONDITION 5:
	# BACKGROUND SELECTION. SMALL PROPORTION BENEFICIAL SELECTION. No change Ne. 
	FILEOUT="'./05_slim_SFS_BGS_PSL.txt'"
	Neb=400;
	rec_rate=1e-4
	rate_ben=0.05; s_backg_ben=0.02;
	rate_del=8; s_backg_del=-0.05;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\>\>  slim_time_memory.txt \& >> ./run_slim_conditions.sh

	# CONDITION 6:
	# BACKGROUND SELECTION. MIDDLE PROPORTION BENEFICIAL SELECTION. No change Ne. 
	FILEOUT="'./06_slim_SFS_BGS_PSM.txt'"
	Neb=400;
	rec_rate=1e-4
	rate_ben=0.5; s_backg_ben=0.02;
	rate_del=8; s_backg_del=-0.05;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\>\>  slim_time_memory.txt \& >> ./run_slim_conditions.sh

	# CONDITION 7:
	# BACKGROUND SELECTION. HIGH PROPORTION BENEFICIAL SELECTION. No change Ne. 
	FILEOUT="'./07_slim_SFS_BGS_PSH.txt'"
	Neb=400;
	rec_rate=1e-4
	rate_ben=2; s_backg_ben=0.02;
	rate_del=8; s_backg_del=-0.05;

	echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT_$iter\" ./slim_template.slim 2\>\>  slim_time_memory.txt >> ./run_slim_conditions.sh
done
