set -xe

make -j ${CPU_COUNT}
mkdir -p $PREFIX/bin
cp bin/unitig_distance $PREFIX/bin/unitig_distance
cp bin/gfa1_parser $PREFIX/bin/gfa1_parser
cp scripts/gwes_plot.r $PREFIX/bin/gwes_plot.r
