#!/usr/bin/make -f

PREFIX ?= /usr
LIBDIR ?= lib
LV2DIR ?= $(PREFIX)/$(LIBDIR)/lv2

OPTIMIZATIONS ?= -msse -msse2 -mfpmath=sse -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only

LDFLAGS ?= -Wl,--as-needed
CXXFLAGS ?= $(OPTIMIZATIONS) -Wall
CFLAGS ?= $(OPTIMIZATIONS) -Wall

###############################################################################
BUNDLE = nrepel.lv2

CXXFLAGS += -fPIC -DPIC
CFLAGS += -fPIC -DPIC

UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
  LIB_EXT=.dylib
  LDFLAGS += -dynamiclib
else
  LDFLAGS += -shared -Wl,-Bstatic -Wl,-Bdynamic
  LIB_EXT=.so
endif


ifeq ($(shell pkg-config --exists lv2 lv2core lv2-plugin || echo no), no)
  $(error "LV2 SDK was not found")
else
  LV2FLAGS=`pkg-config --cflags --libs lv2 lv2core lv2-plugin`
endif

ifeq ($(shell pkg-config --exists lv2-gui || echo no), no)
  $(error "LV2-GUI is required ")
else
  LV2GUIFLAGS=`pkg-config --cflags --libs lv2-gui lv2 lv2core lv2-plugin`
endif


$(BUNDLE): manifest.ttl nrepel.ttl nrepel$(LIB_EXT)
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp manifest.ttl nrepel.ttl nrepel$(LIB_EXT) $(BUNDLE)

nrepel$(LIB_EXT): nrepel.c
	$(CXX) -o nrepel$(LIB_EXT) \
		$(CXXFLAGS) \
		nrepel.c \
		$(LV2FLAGS) $(LDFLAGS)

zameq2.peg: nrepel.ttl
	lv2peg nrepel.ttl nrepel.peg

install: $(BUNDLE)
	install -d $(DESTDIR)$(LV2DIR)/$(BUNDLE)
	install -t $(DESTDIR)$(LV2DIR)/$(BUNDLE) $(BUNDLE)/*

uninstall:
	rm -rf $(DESTDIR)$(LV2DIR)/$(BUNDLE)

clean:
	rm -rf $(BUNDLE) nrepel$(LIB_EXT) nrepel.peg

.PHONY: clean install uninstall
