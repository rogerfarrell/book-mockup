# Makefile for book-mockup
#
# SPDX-License-Identifier: GPL-3.0-or-later

PREFIX     ?= /usr/local
BINDIR     ?= $(PREFIX)/bin
CC         ?= cc
PKG_CONFIG ?= pkg-config
CFLAGS     ?= -O2 -Wall

GEGL_CFLAGS := $(shell $(PKG_CONFIG) --cflags gegl-0.4)
GEGL_LIBS   := $(shell $(PKG_CONFIG) --libs gegl-0.4)

PROGRAM = book-mockup

.PHONY: all install uninstall clean

all: $(PROGRAM)

$(PROGRAM): book-mockup.c
	$(CC) -std=c11 $(CFLAGS) $(GEGL_CFLAGS) -o $@ book-mockup.c \
	  $(LDFLAGS) $(GEGL_LIBS) -lm

install: $(PROGRAM)
	install -d $(DESTDIR)$(BINDIR)
	install -m 755 $(PROGRAM) $(DESTDIR)$(BINDIR)/$(PROGRAM)

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(PROGRAM)

clean:
	rm -f $(PROGRAM)
