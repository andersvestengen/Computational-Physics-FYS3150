all: compile execute

compile:
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_a.x CodeBase/main_a.cpp CodeBase/Finance.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_c.x CodeBase/main_c.cpp CodeBase/Finance.cpp   -larmadillo -llapack -lblas
	#g++ -std=c++17 -O3 -flto -Wall -o CodeBase/main_d.x CodeBase/main_d.cpp CodeBase/Finance.cpp   -larmadillo -llapack -lblas
	g++ -std=c++17  -O2 -pipe -g0 -fstack-protector -flto  -Wall -o CodeBase/main_e.x CodeBase/main_e.cpp CodeBase/Finance.cpp   -larmadillo -llapack -lblas
	#For running Python in c++ : -framework Python -I/opt/anaconda3/include/python3.7m -I/opt/anaconda3/include/python3.7m -Wno-unused-result -Wsign-compare -Wunreachable-code -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -I/opt/anaconda3/include -arch x86_64 -I/opt/anaconda3/include -arch x86_64
	# For parallallized code : -Xpreprocessor -fopenmp -lomp
execute:
	#CodeBase/./main_a.x
	#CodeBase/./main_c.x
	#CodeBase/./main_d.x
	CodeBase/./main_e.x
	#python3 Codebase/Variance.py
	#python3 Codebase/histo.py
	python3 Codebase/hist_e.py
	#python3 Codebase/Omega.py
	#python3 Codebase/diff_savings.py
	python3 Codebase/Alpha_vary.py
