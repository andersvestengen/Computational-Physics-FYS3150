all: compile execute

compile:
	#g++ -std=c++11 -O3 -flto -Wall -o main_B.x main_B.cpp classtuff.cpp  -larmadillo -llapack -lblas
	#g++ -std=c++11 -O3 -flto -Wall -o unit.x unittest.cpp classtuff.cpp -larmadillo -llapack -lblas
	#g++ -std=c++11 -O3 -flto -Wall -o main_D.x main_D.cpp classtuff.cpp  -larmadillo -llapack -lblas
	g++ -std=c++11 -O3 -flto -Wall -o main_E.x main_E.cpp classtuff.cpp  -larmadillo -llapack -lblas
execute:
	#./main_B.x 1
	#./unit.x
	#./main_D.x
	./main_E.x

	python3 Plotte_code.py
