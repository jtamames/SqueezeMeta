# fix libs for samtools and mothur
cd $PREFIX/lib
ln -s libncurses.so.6 libncurses.so.5
ln -s libtinfo.so.6 libtinfo.so.5
ln -s libreadline.so.8 libreadline.so.6
cp libboost_iostreams.so libboost_iostreams.so.1.85.0
cp libboost_system.so libboost_system.so.1.85.0
cp libboost_filesystem.so libboost_filesystem.so.1.85.0
