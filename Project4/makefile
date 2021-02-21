all: compile execute

compile:
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_c.x CodeBase/main_c.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_d.x CodeBase/main_d.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_e.x CodeBase/main_e.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -fopenmp -Wall -o CodeBase/main_f.x CodeBase/main_f.cpp CodeBase/solver.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/unit.x CodeBase/Unittest.cpp CodeBase/solver.cpp -larmadillo -llapack -lblas
	# Mac computers needs this instead: -Xpreprocessor -fopenmp -lomp
execute:
	#CodeBase/./main_c.x
	#CodeBase/./main_d.x
	#CodeBase/./main_e.x
	#CodeBase/./main_f.x
	#CodeBase/./unit.x
	#python3 CodeBase/EMC_visual.py
	#python3 CodeBase/E_histogram.py
	#python3 Codebase/T_critical.py