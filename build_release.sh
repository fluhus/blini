# Compiles blini executables for all platforms.

set -e

VERSION=v1.1.0
VERSION2=v2.0.0

EXE=blini
OUTDIR=../release
TAGS=
SUFFIX=

if [ "$1" == "2" ]; then
  VERSION=$VERSION2
  TAGS="-tags=blini2"
  SUFFIX=2
fi

FLAGS="-ldflags=-s -X main.version=$VERSION"

rm -fr $OUTDIR
mkdir $OUTDIR

# Linux
go build $TAGS "$FLAGS" -o $OUTDIR ./$EXE
zip -j $OUTDIR/${EXE}${SUFFIX}_linux.zip $OUTDIR/$EXE
rm $OUTDIR/$EXE

# Mac
GOOS=darwin go build $TAGS "$FLAGS" -o $OUTDIR ./$EXE
zip -j $OUTDIR/${EXE}${SUFFIX}_mac.zip $OUTDIR/$EXE
rm $OUTDIR/$EXE

# Windows
GOOS=windows go build $TAGS "$FLAGS" -o $OUTDIR ./$EXE
zip -j $OUTDIR/${EXE}${SUFFIX}_win.zip $OUTDIR/$EXE.exe
rm $OUTDIR/$EXE.exe

ls -lh $OUTDIR
