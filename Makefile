
CXX = g++
CXXFLAGS = -Wall -fPIC
ROOTLIBS = -lMathMore -lProof
ROOTCONFIG=`root-config --cflags --glibs`
HEADERS = NeutronFit_BC501A.hh ProtonFitter.hh vec.hh
SOURCES = NeutronFit_BC501A.cc ProtonFitter.cc
FITTERLIB = -L. -lProtonFitter

fit_proton: fit_proton.cc libProtonFitter.so
	$(CXX) $(CXXFLAGS) -o fit_proton $(FITTERLIB) $(ROOTLIBS) $(ROOTCONFIG) $^

fit_protonNM: fit_protonNM.cc libProtonFitter.so
	$(CXX) $(CXXFLAGS) -o fit_protonNM $(FITTERLIB) $(ROOTLIBS) $(ROOTCONFIG) $^

draw: draw.cc libProtonFitter.so
	$(CXX) $(CXXFLAGS) -o draw $(FITTERLIB) $(ROOTLIBS) $(ROOTCONFIG) $^

libProtonFitter.so: ProtonFitterDict.cxx $(SOURCES)
	g++ -shared -o $@ $(ROOTCONFIG) $(CXXFLAGS) $(ROOTLIBS) -I$(ROOTSYS)/include $^

ProtonFitterDict.cxx: $(HEADERS) LinkDef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

NeutronFit_BC501A.o: NeutronFit_BC501A.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) NeutronFit_BC501A.cc $(ROOTCONFIG)                

ProtonFitter.o: ProtonFitter.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) ProtonFitter.cc $(ROOTCONFIG)                

clean:
	rm -f *.o *.pcm ProtonFitterDict.cxx libProtonFitter.so fit_proton fit_protonNM draw 

all: fit_proton draw
