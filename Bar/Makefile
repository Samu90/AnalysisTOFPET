CC = g++
CFLAGS = -Wall -c -g -fpermissive 


INCLUDES = -I. -I.. -I$(ROOTSYS)/include -I$(ROOFIT_INCLUDE)/ 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG  -lPhysics -lz -lTreePlayer -lRooFitCore -lRooFit


bar: BarAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o BarAnalysis.exe BarAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

temp: AnalysisTemp.cpp
	$(CC) -Wall $(INCLUDES) -o AnalysisTemp.exe AnalysisTemp.cpp $(ROOTFLAG) $(EXTRALIBS)

clean:
	rm BarAnalysis.exe
run: 
	./macro
