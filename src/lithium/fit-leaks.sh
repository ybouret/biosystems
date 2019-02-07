rm -f data_fuite*.fit.dat
rm -f fit-leaks.dat
../../bin/fit-leak data/fuite*.txt
cp data_fuite*.fit.dat supmat/leaks
cp fit-leaks.dat       supmat/leaks
