mkdir -f fa;
cd fa;
for f in `cat download.list`; do wget $f; done;
cd ..
