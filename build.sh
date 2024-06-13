# Compiles blini executables for all platforms.
set -e

rm -fr build
mkdir build

# Linux
go build -o build ./blini
zip -j build/blini_linux_amd64.zip build/blini

# Mac
GOOS=darwin GOARCH=arm64 go build -o build ./blini
zip -j build/blini_macos_arm64.zip build/blini

# Windows
GOOS=windows go build -o build ./blini
zip -j build/blini_win_amd64.zip build/blini.exe

rm build/blini build/blini.exe
