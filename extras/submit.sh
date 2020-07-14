# This scripts allows for submission of the Bayes-factor computation code #
# in multiple batches. It is a very simple tools to quickly compute the #
# Bayes-factor with multiple trials. For example if you need to run 100 trials #
# you can simply set N = 10 and n = 10. This will result in running of 10 parallel #
# processes, each running 10 trials, leading to a total of 100 trial runs. #

# Sample command:
# ./submit.sh <target eos> <reference eos> <path to posterior_samples.dat file>
#----------------------------------------------------------------------------------#

N=3  # Number of trials
n=7 # Number of processes

if test -f "$3"; then
    echo "$3 exists."
else
    echo "\n$3 does not exists."
    echo "Use correct option for submitting the run"
    echo "./submit.sh <target-eos> <reference-eos> <path-to-posterior-samples>"
    echo
    exit
fi

now=$(date +"%T")
echo "Run started at "$now
echo "Total number of processes running = "$n
echo "Each job running "$N "trials"
results="results"_$1
mkdir -p outfiles
rm -rf $results # Just in case previous instance of run exits.
mkdir -p $results
# Loop over all the processes
for ((i = 1; i <= n; i++));
do
  ./driver -i $3 -T $1 -r $2 -N $N -o $results/testing$i.json > outfiles/stdout$i.out &
done 

./aggregator -d $results -p $n


