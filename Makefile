#TARGET=mu_pi_ratio
#TARGET=connect_track
TARGET=draw_mu_tracks

FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion -lEDA -lShower -lScan -lMLP -lSpectrum

$TARGET: $(TARGET).cpp
	g++ $(TARGET).cpp -w -Iinclude src/*.cpp `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs` -o $(TARGET)

clean:
	$(RM) $(TARGET)
