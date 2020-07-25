if dpkg --compare-versions $(dpkg -s libitpp-dev | awk '/^Version:/ { print $2}') lt 4.3.1-9;then
    cd /usr/include
    wget -O - -q - https://sources.debian.org/data/main/libi/libitpp/4.3.1-9/debian/patches/sse2-immediate-fix.diff | patch -p1
fi
