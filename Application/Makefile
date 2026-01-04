
.PHONY: all clean

all: script_application

#all: bos_application poh_application sn_application rs_application tc_application script_application bos_opt1

BOS_OBJS = BOS_Application.o gcm.o

POH_OBJS = POH_Application.o gcm.o

SN_OBJS = SN_Application.o gcm.o

RS_OBJS = RS_Application.o gcm.o

TC_OBJS = TC_Test_Application.o gcm.o

SA_OBJS = Script_Application.o gcm.o

BOSO1_OBJS = BOS_Optimizer_1.o

%.o: %.cpp
	g++ -std=c++11 -O0 -ggdb -c $<

bos_application: $(BOS_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -lOSort -lcrypto -Wl,--rpath=$(CURDIR) -o $@

poh_application: $(POH_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -lOSort -lcrypto -Wl,--rpath=$(CURDIR) -o $@

sn_application: $(SN_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -lOSort -lcrypto -Wl,--rpath=$(CURDIR) -o $@

rs_application: $(RS_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -lOSort -lcrypto -Wl,--rpath=$(CURDIR) -o $@

tc_application: $(TC_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -lOSort -lcrypto -Wl,--rpath=$(CURDIR) -o $@

script_application: $(SA_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -lOSort -lcrypto -Wl,--rpath=$(CURDIR) -o $@

bos_opt1: $(BOSO1_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -lOSort -lcrypto -Wl,--rpath=$(CURDIR) -o $@

borp_fp_tester_app: $(BORP_OBJS)
	g++ -std=c++11 $^ -L=$(CURDIR) -Wl,--rpath=$(CURDIR) -o $@

clean:
	rm -f $(OBJS) $(BOS_OBJS) $(POH_OBJS) $(SN_OBJS) $(RS_OBJS) $(TC_OBJS) $(SA_OBJS) $(BOSO1_OBJS) *_application *.so
