# this script will run FastSMC on a simulated data as described in the paper (in FILES/FASTSMC_EXAMPLE/)
# parameters can be changed if desired

# This example will run multiple jobs in different threads on the same machine. If you are running FastSMC on a cluster
# then it may be more appropriate to instead use the job scheduler such as `qsub`.

# The total number of jobs you want to run in parallel (this should be a square number).
# Note that the standard output will be messy as information will be printed from every job simultaneously.
total_num_jobs=4

cd ../FASTSMC_BUILD_DIR/ || exit

run_single_job() {
  local job=$1
  ./FastSMC_exe --inFileRoot ../FILES/FASTSMC_EXAMPLE/example \
    --outFileRoot ../cpp_example/FastSMC_output_example \
    --decodingQuantFile ../FILES/FASTSMC_EXAMPLE/example.decodingQuantities.gz \
    --mode array \
    --time 50 \
    --min_m 1.5 \
    --segmentLength \
    --hashing \
    --jobs "${total_num_jobs}" \
    --jobInd "${job}" \
    --perPairPosteriorMeans \
    --perPairMAP \
    --noConditionalAgeEstimates \
    --bin
  echo Finished job "$1"
}

# Run the jobs in parallel
for ((i = 1; i <= total_num_jobs; i++)); do
  # substitute the command below with the appropriate command if running this on a cluster, e.g. "qsub run_single_job $i"
  run_single_job "$i" &
done
wait

# Binary output file can be converted with the following command line
echo 'Showing first lines of the binary output...'
./convertBinary_exe ../cpp_example/FastSMC_output_example.1.4.FastSMC.bibd.gz | head

# Note that there will be the same number of output files as jobs, which will need to be concatenated.
