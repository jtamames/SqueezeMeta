# fix libs for samtools and mothur
cd $PREFIX/lib
ln -s libncurses.so.6 libncurses.so.5
#ln -s libtinfo.so.6 libtinfo.so.5 # this makes the more and less commands fail in some systems, remove and see if it's still needed
ln -s libreadline.so.8 libreadline.so.6
cp libboost_iostreams.so libboost_iostreams.so.1.85.0
cp libboost_system.so libboost_system.so.1.85.0
cp libboost_filesystem.so libboost_filesystem.so.1.85.0
# Install CONCOCT locally, ensures that we have the same math libraries that we'll have at runtime
cd $PREFIX/SqueezeMeta/bin/CONCOCT-1.1.0/
$PYTHON setup.py install
