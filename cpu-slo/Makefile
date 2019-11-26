main: main.cpp
	g++ -std=c++11 -O3 -fopenmp main.cpp -o main	

qsub: main
	qsub -q mamba -l walltime=08:00:00 -l nodes=1:ppn=16 -d `pwd` run.sh
