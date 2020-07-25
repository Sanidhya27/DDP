cd /usr/src/gtest
sudo -i cmake CMakeLists.txt
sudo -i make
sudo cp *.a /usr/lib
sudo ln -s /usr/lib/libgtest.a /usr/local/lib/gtest/libgtest.a
sudo -i ln -s /usr/lib/libgtest_main.a /usr/local/lib/gtest/libgtest_main.a
