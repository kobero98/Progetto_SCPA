#module load mpi;mpicc MPI-prod-kobero.c -o MPI-prod-kobero;
#mpicc MPI-prod-fanfa.c -o MPI-prod-fanfa;
all:
	echo "Bash version ${BASH_VERSION}...";
	module load mpi;
	mkdir target
	mkdir result
	mpicc MPI-prod-fanfa.c -o target/MPI-prod-fanfa;
	mpicc MPI-prod-kobero.c -o target/MPI-prod-kobero;
	bash start.sh
clean:
	rm -d -r target;
	rm -d -r result;
