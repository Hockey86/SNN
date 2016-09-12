# Joint Weight-Delay Spike-Timing Dependent Plasticity: A NEST Module

This is a [NEural Simulation Tool (NEST)](http://www.nest-simulator.org/) module for joint weight-delay spike-timing dependent plasticity (joint STDP). The joint STDP allows to update weights and delays jointly, so that

* synapses carrying causal spikes are strengthened, while others are weakened.
* the axonal delays become complementary to the spike timing, which maximizes the spike synchrony

With the joint STDP, we can use the weights and delays to recover the spatiotemporal pattern of cell assemblies (or polychronous neuronal groups) at millisecond scale.

![The joint STDP windows](figuers/fig_joint_stdp.png "The joint STDP windows")

## Neuron and Synapse Models

* Synapse model
  * `joint_stdp_synapse`: Synapse with joint plasticity;
  * `weight_stdp_synapse`: Synapse with the conventional STDP which only updates weights;
  * `immediate_inhibition_synapse`: Synapse which can inhibit postsynaptic neurons with zero delay (useful for lateral inhibition).

* Neuron model
  * `spike_response_neuron`: Same with `iaf_chs_2007` in NEST, but additionally stores presynaptic delays and spike history;
  * `parrot_neuron2`: Same with `parrot_neuron` in NEST, but additionally stores presynaptic delays and spike history.

These synapse models can only be used to connect neuron models in this module.

## Important Note

### Dendritic and axonal delays

NEST handles spike event at presynaptic site and uses dendritic delay. However, joint STDP is based on axonal delay, which handles spike event at postsynaptic site. To implement handling spike event at postsynaptic site, when a spike event happens at presynaptic site,
we call `event.source.update_upstream_weights/delays`.

This is the reason we need to define neuron models which store presynaptic information.

We are also considering whether it is equivalent to use dendritic delay with sign-reversed delay STDP, which allows an easier implementation in NEST.

### NEST versions

The code used in the paper was written with NEST 2.6.0, while latest NEST version is 2.10.0 (as of Sep 12, 2016). This is the reason that we provide two versions here. But they should be equivalent. We will try to update this module when new versions of NEST is released.

## Installation

### Dependencies

* Ubuntu 14.04 or later;
* NEST 2.10.0 or NEST 2.6.0.

### Steps

1. Make sure that NEST is correctly installed. Assume
    * the installation destination is $HOME/opt/nest;
    * the source folder is $HOME/nest-X.Y.Z;
    * the build folder is $HOME/nest-X.Y.Z-build.
2. Copy the folder JointSTDPModule to $HOME/JointSTDPModule.
3. Execute the following code in the terminal.
```bash
cd $HOME/JointSTDPModule
chmod +x bootstrap.sh && ./bootstrap.sh
mkdir $HOME/JointSTDPModule-build
cd $HOME/JointSTDPModule-build
../JointSTDPModule/configure --with-nest=$HOME/opt/nest/bin/nest-config
make
make install
cd $HOME/nest-X.Y.Z-build
../nest-X.Y.Z/configure --prefix=$HOME/opt/nest --with-modules="jointstdpmodule"
make
make install
```

## Example

In your Python script, use the following code to load the module.

```python
import nest
nest.Install('jointstdpmodule')
```

Detailed examples can be found in
* joint_stdp_example.py: An example to demonstrate the joint STDP under input patterns without noise;
* validate_STDP.py: Validation of `weight_stdp_synapse` and `joint_stdp_synapse` in this module.

## How to cite

### MLA format

> Sun, Haoqi, Olga Sourina, and Guang-Bin Huang. "Learning Polychronous Neuronal Groups Using Joint Weight-Delay Spike-Timing-Dependent Plasticity." *Neural Computation* (2016).

### BibTex
>  @article{sun2016learning,  
>  title={Learning Polychronous Neuronal Groups Using Joint Weight-Delay Spike-Timing-Dependent Plasticity},  
>  author={Sun, Haoqi and Sourina, Olga and Huang, Guang-Bin},  
>  journal={Neural Computation},  
>  year={2016},  
>  publisher={MIT Press}  
>  }

## Contact

hockeysun86@gmail.com

## License

This module is provided under [NEST GNU General Public License 2](http://www.nest-simulator.org/license/) or later.

