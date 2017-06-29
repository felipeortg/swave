CONFIG="config_files/config_g_1.txt"

COUNTER=0
TOTAL=5

echo ------- BEGIN SPECTRUM CALCULATION --------

while [  $COUNTER -lt $TOTAL ]; do

	
	echo -------------------------------
	echo --- Modifying configuration file ---

	python changeL.py $CONFIG $COUNTER $TOTAL

	let COUNTER=COUNTER+1

	echo -------------------------------
	echo ----- Loop $COUNTER out of $TOTAL -------- 
	echo -------------------------------
	echo -------------------------------
	echo --- Evaluating 'zdlm.py' ---
	python ../felipe/zdlm.py $CONFIG

	sleep 1

	echo -------------------------------
	echo --- Evaluating 'quicklook.py' ---
	python quicklook.py $CONFIG

	sleep 1

	echo -------------------------------
	echo --- Evaluating 'fzdlm_g.py' ---
	python fzdlm_g.py $CONFIG

	sleep 1

	echo -------------------------------
	echo --- Evaluating 'findspec.py' ---
	python findspec.py $CONFIG

done

python final_plot.py $CONFIG