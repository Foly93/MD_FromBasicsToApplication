#!/bin/bash
echo "total/s per.CPU/s Acc" | awk '{printf "%12s %11s %6s\n", $1,$2,$3}' > timings.txt

for i in $(ls _LOG/*)
do 
	grep -A1 'Core t' $i | tail -1 | awk '{printf "%12.3f %11.3f %6.1f\n", $2,$3,$4}' >> timings.txt
done

sum_total=$(awk 'NR > 1{print $1}' timings.txt | paste -sd+ | bc)
sum_pCPU=$(awk 'NR > 1{print $2}' timings.txt | paste -sd+ | bc)

echo "$sum_total $sum_pCPU" | awk '{printf "%12.3f %11.3f",$1,$2}' >> timings.txt

