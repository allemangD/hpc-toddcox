main : main.cu util.h
	nvcc -o main -std=c++11 -O3 main.cu

clean :
	rm main

queue : main
	qsub -q mamba -l walltime=24:00:00 -l nodes=1:ppn=1:gpus=1 -d `pwd` run_gpu_opt.sh

