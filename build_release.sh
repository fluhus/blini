# Compiles blini executables for all platforms.

set -e

VERSION=v2.0.0
VERSION1=v1.1.0

EXE=blini
OUTDIR=../release
TAGS=
SUFFIX=2

if [ "$1" == "1" ]; then
  VERSION=$VERSION1
  TAGS="-tags=blini1"
  SUFFIX=
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
