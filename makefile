#Makefile
main: fdtd3d.o memory_allocate3d.o memory_allocate4d.o sigma_calc.o D_update.o D_update_pml.o E_update.o H_update.o H_update_pml.o pml_class.o
	g++ -o main fdtd3d.o memory_allocate3d.o memory_allocate4d.o sigma_calc.o D_update.o D_update_pml.o E_update.o H_update.o H_update_pml.o pml_class.o -fopenmp

fdtd3d.o: fdtd3d.cpp fdtd3d.h
	g++ -c fdtd3d.cpp

memory_allocate3d.o: memory_allocate3d.cpp fdtd3d.h
	g++ -c memory_allocate3d.cpp

memory_allocate4d.o: memory_allocate4d.cpp fdtd3d.h
	g++ -c memory_allocate4d.cpp

sigma_calc.o: sigma_calc.cpp fdtd3d.h
	g++ -c sigma_calc.cpp

D_update.o: D_update.cpp fdtd3d.h
	g++ -c D_update.cpp

D_update_pml.o: D_update_pml.cpp fdtd3d.h
	g++ -c D_update_pml.cpp

E_update.o: E_update.cpp fdtd3d.h
	g++ -c E_update.cpp

H_update.o: H_update.cpp fdtd3d.h
	g++ -c H_update.cpp

H_update_pml.o: H_update_pml.cpp fdtd3d.h
	g++ -c H_update_pml.cpp

pml_class: pml_class.cpp fdtd3d.h
	g++ -c pml_class.cpp

main.o: main.cpp fdtd3d.h
	g++ -c main.cpp









