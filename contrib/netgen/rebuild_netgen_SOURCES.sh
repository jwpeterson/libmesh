#!/bin/sh

# Hopefully this gets everything libMesh uses
sources=$(find netgen/libsrc/ netgen/nglib/ -name '*.h' -o -name '*.hpp' -o -name '*.c' -o -name '*.cpp' -type f | LC_COLLATE=POSIX sort)

echo "# Do not edit - automatically generated from $0" > netgen_SOURCES
printf '%s' "netgen_SOURCE_FILES = " >> netgen_SOURCES
for source_with_path in $sources ; do

    echo " \\" >> netgen_SOURCES
    printf '%s' "        "$source_with_path >> netgen_SOURCES

done

echo " " >> netgen_SOURCES

