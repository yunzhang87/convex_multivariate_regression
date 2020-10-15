cd GGM_code
make clean
make

wait

cd ../cGGM_code
make clean
make
wait

cd ../simulation
mkdir result
mkdir out

k=1
R CMD BATCH "--args $k" sim_example.R sim_example$k.out &
wait

mv *.out out

