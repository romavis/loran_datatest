TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CFLAGS += -std=gnu99 -O2
#QMAKE_CFLAGS_DEBUG += -pg
#QMAKE_LFLAGS_DEBUG += -pg

SOURCES += main.c \
    loran.c \
    loran_reference.c \
    loran_corr.c

HEADERS += \
    loran.h \
    loran_reference.h

