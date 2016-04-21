#!/usr/bin/make -f

PREFIX ?= /usr
LIBDIR ?= lib
LV2DIR ?= $(PREFIX)/$(LIBDIR)/lv2

#Basic Flags

OPTIMIZATIONS ?= -msse -msse2 -mfpmath=sse -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only

LDFLAGS ?= -Wl,--as-needed -shared -Wl,-Bstatic -Wl,-Bdynamic
CXXFLAGS ?= $(OPTIMIZATIONS) -Wall -fPIC -DPIC -lfftw3 -lfftw3f -lm

BUNDLE = nrepel.lv2
LIB_EXT=.so
###############################################################################

#library detection
ifeq ($(shell pkg-config --exists lv2 lv2core lv2-plugin || echo no), no)
  $(error "LV2 SDK was not found")
else
  LV2FLAGS=`pkg-config --cflags --libs lv2 lv2core lv2-plugin`
endif

#directory creation
$(BUNDLE): manifest.ttl nrepel.ttl nrepel$(LIB_EXT)
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp manifest.ttl nrepel.ttl $(BUNDLE)
	mv nrepel$(LIB_EXT) $(BUNDLE)

#file compiling
nrepel$(LIB_EXT): nrepel.c denoise.c nestim.c
	$(CXX) -o nrepel$(LIB_EXT) \
		$(CXXFLAGS) \
		nrepel.c \
		denoise.c \
		nestim.c \
		$(LV2FLAGS) $(LDFLAGS)

#ttl files control
nrepel.peg: nrepel.ttl
	lv2peg nrepel.ttl nrepel.peg

#make recipes
install: $(BUNDLE)
	install -d $(DESTDIR)$(LV2DIR)/$(BUNDLE)
	install -t $(DESTDIR)$(LV2DIR)/$(BUNDLE) $(BUNDLE)/*

uninstall:
	rm -rf $(DESTDIR)$(LV2DIR)/$(BUNDLE)

clean:
	rm -rf $(BUNDLE) nrepel$(LIB_EXT) nrepel.peg

.PHONY: clean install uninstall
