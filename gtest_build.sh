cd /usr/src/gtest
cmake CMakeLists.txt
make
cp *.a /usr/lib
ln -s /usr/lib/libgtest.a /usr/local/lib/gtest/libgtest.a
ln -s /usr/lib/libgtest_main.a /usr/local/lib/gtest/libgtest_main.a
