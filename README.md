## Introduction is unavailable temporarily until something good happens.

## Installation

### Dependencies

* Ubuntu 14.04;
* Python 2.7;
* NEST 2.6.0.

### Steps

1. Make sure that NEST is correctly installed in $HOME
    * $HOME/opt/nest exists;
    * $HOME/nest-2.6.0 exists;
2. Extract and copy the folder JointSTDPModule to $HOME/nest-2.6.0/examples;
3. Move install\_jointstdpmodule\_nest.sh and jointstdpmodule\_example.py to $HOME;
4. cd $HOME;
5. chmod +x install\_jointstdpmodule\_nest.sh;
6. sudo ./install\_jointstdpmodule\_nest.sh.

## Usage

In your python script, use the following code to load the module.

```
import nest
nest.Install('jointstdpmodule')
```

