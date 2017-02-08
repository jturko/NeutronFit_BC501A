all: fit fit_foreward fit_backward fit_middle

fit: 
	g++ -Wall fit.cc -o fit `root-config --cflags --libs --glibs`

fit_foreward: 
	g++ -Wall fit_foreward.cc -o fit_foreward `root-config --cflags --libs --glibs`

fit_backward: 
	g++ -Wall fit_backward.cc -o fit_backward `root-config --cflags --libs --glibs`

fit_middle: 
	g++ -Wall fit_middle.cc -o fit_middle `root-config --cflags --libs --glibs`

clean:
	rm -f fit
	rm -f fit_foreward
	rm -f fit_backward
	rm -f fit_middle



