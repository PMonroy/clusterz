CPP=g++ 
CPPFLAGS=-Wall -O2 -fopenmp

LDFLAGS= -lnetcdf_c++ -fopenmp
RM=rm -rf

DESTDIR=$(HOME)/bin

common_src=  rparameters.cpp gridconstruction.cpp lagrangian_engine.cpp constants.cpp vflow.cpp eqdate.cpp vectorXYZ.cpp VTK.cpp mt19937ar.cpp
common_obj=$(common_src:.cpp=.o) 
common_dep=$(common_obj:.o=.d)  # one dependency file for each source

#CLUSTERz
clusterz_src= clusterz.cpp
clusterz_obj=$(clusterz_src:.cpp=.o) 
clusterz_dep=$(clusterz_obj:.o=.d)  # one dependency file for each source

#SINKTRAJ
sinktraj_src= sinktraj.cpp
sinktraj_obj=$(sinktraj_src:.cpp=.o) 
sinktraj_dep=$(sinktraj_obj:.o=.d)  # one dependency file for each source

#AGGREGATE
aggregate_src= aggregate.cpp
aggregate_obj=$(aggregate_src:.cpp=.o) 
aggregate_dep=$(aggregate_obj:.o=.d)  # one dependency file for each source

.PHONY: all clusterz sinktraj aggregate

all:  clusterz sinktraj aggregate

clusterz: $(common_obj) $(clusterz_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

sinktraj: $(common_obj) $(sinktraj_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

aggregate: $(common_obj) $(aggregate_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

-include $(common_dep) # include all dep files in makefile
-include $(clusterz_dep)
-include $(sinktraj_dep) 
-include $(aggregate_dep) 

%.d: %.cpp 	# rule to generate a dep file by using the g++ prepocesor
	$(CPP) $(CPPFLAGS) -MM -MT $(@:.d=.o) $< -MF $@

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -o $@ -c $<

.PHONY: debug clusterz sinktraj aggregate
debug: CPPFLAGS+= -DDEBUG -ggdb # debug with gdb
debug: clusterz sinktraj aggregate

.PHONY: clean
clean:
	$(RM) $(common_obj) $(clusterz_obj) $(aggregate_obj)  $(sinktraj_obj) *.d *~ *# clusterz sinktraj aggregate

.PHONY: install
install: clusterz sinktraj aggregate
	install $^ $(DESTDIR)

# Trick to run comfortably in emacs (emulating compile command)
#ARGS="--default"
#run:	# For include other args use -> make run ARGS="-arg1 -arg2=foo" 
#	rtime $(ARGS)
