cd 000mapping; perl bat.pl | bash ; cd ..
cd 010besthit; bash bat.sh; cd ..
cd 020getmutation; perl bat1.pl | bash; bash bat2.sh; cd ..
cd 020getmutation; grep uniq 78G10_S8_L001_result.txt | grep SNP > final_result.txt; cd ..
