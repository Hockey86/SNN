cd $HOME
rm -r JointSTDPModule_build
rm -r JointSTDPModule
rm -rfv $HOME/opt/nest/lib/nest && mkdir $HOME/opt/nest/lib/nest
cp -r $HOME/nest-2.6.0/examples/JointSTDPModule $HOME/JointSTDPModule
cd JointSTDPModule
./bootstrap.sh
cd ..
mkdir JointSTDPModule_build
cd JointSTDPModule_build
../JointSTDPModule/configure --with-nest=$HOME/opt/nest/bin/nest-config
make
make install
cd $HOME/nest-2.6.0-build
../nest-2.6.0/configure --prefix=$HOME/opt/nest --with-modules="jointstdpmodule"
make
make install
