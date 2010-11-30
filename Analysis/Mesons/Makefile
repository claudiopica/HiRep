DEFAULT: bs_mesons

ranlxs.o: ranlxs.cpp ranlxs.h Makefile
	g++ -Wall -O2 -c ranlxs.cpp

datasample.o: datasample.cpp datasample.h Makefile
	g++ -Wall -O2 -c datasample.cpp

effective_mass.o: effective_mass.cpp effective_mass.h Makefile
	g++ -Wall -O2 -c effective_mass.cpp

bootstrap.o: bootstrap.cpp bootstrap.h Makefile
	g++ -Wall -O2 -c bootstrap.cpp

progbar.o: progbar.cpp progbar.h Makefile
	g++ -Wall -O2 -c progbar.cpp

bs_reader.o: bs_reader.cpp bs_reader.h bs_type.h Makefile
	g++ -Wall -O2 -c bs_reader.cpp

bs_ctrl.o: bs_ctrl.cpp bs_ctrl.h Makefile
	g++ -Wall -O2 -c bs_ctrl.cpp

bs_observables.o: bs_observables.cpp bs_observables.h Makefile
	g++ -Wall -O2 -c bs_observables.cpp

bs_mesons: bs_mesons.cpp bs_type.h bs_ctrl.o bs_reader.o progbar.o datasample.o ranlxs.o bootstrap.o effective_mass.o bs_observables.o Makefile
	g++ -Wall -O2 bs_mesons.cpp bs_ctrl.o bs_reader.o progbar.o datasample.o ranlxs.o bootstrap.o effective_mass.o bs_observables.o -o bs_mesons

clean:
	rm -f *.o bs_mesons

.PHONY: clean
