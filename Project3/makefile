all: compile execute

compile:
  #g++ -std=c++17 -O2 -flto -Wall -o CodeBase/main_a.x CodeBase/main_a.cpp CodeBase/object.cpp CodeBase/solving.cpp  -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/main_c.x CodeBase/main_c.cpp CodeBase/object.cpp CodeBase/solving.cpp  -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/main_f.x CodeBase/main_f.cpp Codebase/object.cpp CodeBase/solving.cpp  -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/main_g.x CodeBase/main_g.cpp CodeBase/object.cpp CodeBase/solving.cpp  -larmadillo -llapack -lblas
	g++ -std=c++17 -O2 -flto -Wall -o CodeBase/main_h.x CodeBase/main_h.cpp CodeBase/object.cpp CodeBase/solving.cpp  -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/main_i.x CodeBase/main_i.cpp CodeBase/object.cpp CodeBase/solving.cpp  -larmadillo -llapack -lblas
	#g++ -std=c++17 -O2 -flto -Wall -o CodeBase/unit.x CodeBase/unittest.cpp CodeBase/solving.cpp CodeBase/object.cpp -larmadillo -llapack -lblas
execute:
	#CodeBase/./main_a.x
	#python3 CodeBase/Plotte_code.py
	#CodeBase/./main_c.x
	#CodeBase/./main_f.x
	#CodeBase/./main_g.x
	CodeBase/./main_h.x
	#CodeBase/./main_i.x
	#CodeBase/./unit.x
	python3 CodeBase/Plotte_code_MB.py
