Sources=b2ddsj_evt.cc b2ddsj_mcevt.cc

INC=-I../include -I/home/belle/vitaly/belle/include

b2ddsj_evt_dict.cpp: $(Sources)
	rootcint -f $@ -c $(INC) -p $^

