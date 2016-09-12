#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

validate_STDP.py

Validation of the weight/joint STDP in JointSTDP Module.
Modify variable "synapse_type" to switch between weight/joint STDP.

Haoqi Sun, Sep 12, 2016
"""


import numpy as np
import nest
nest.Install('jointstdpmodule')

## define constants

WEIGHT_STDP_TYPE = 0
JOINT_STDP_TYPE = 1
NO_CONNECTION = np.inf
    
def get_weights(sources, targets=None):
    """helper function to get weights as numpy array, where no connection is denoted by NO_CONNECTION"""
    if targets is None:
        targets = sources
    weights = np.zeros((len(sources),len(targets)))+NO_CONNECTION
    conn = nest.GetConnections(source=sources,target=targets)
    nest_weights = nest.GetStatus(conn,'weight')
    for i in range(len(conn)):
        weights[sources.index(conn[i][0]),targets.index(conn[i][1])] = nest_weights[i]
    return weights


def get_delays(sources, targets=None):
    """helper function to get delays as numpy array, where no connection is denoted by NO_CONNECTION"""
    if targets is None:
        targets = sources
    delays = np.zeros((len(sources),len(targets)))+d_min
    conn = nest.GetConnections(source=sources,target=targets)
    nest_delays = nest.GetStatus(conn,'delay')
    for i in range(len(conn)):
        delays[sources.index(conn[i][0]),targets.index(conn[i][1])] = nest_delays[i]
    return delays


def reset_all():
    """reset the simulator kernel and define the models"""
    nest.ResetKernel()
    nest.set_verbosity(verbose)
    np.random.seed(random_seed)
    pyrngs = [np.random.RandomState(i) for i in range(random_seed, random_seed+vp_num)]
    nest.SetKernelStatus({'resolution':time_step,'print_time':print_time,"local_num_threads":vp_num,'grng_seed':random_seed+vp_num,'rng_seeds':range(random_seed+vp_num+1, random_seed+2*vp_num+1)})
    
    ## neuron models

    nest.SetDefaults('spike_response_neuron',{
            'tau_epsp':tau_epsp,
            'tau_ipsp':tau_ipsp,
            'tau_reset':tau_reset,
            'V_reset':u_reset,
            't_ref':t_abs_refrac,
            #'tau_minus':tau_n,  # time constant in the negative side of the weight stdp learning window
        })

    nest.CopyModel('parrot_neuron2','input_neuron')
    nest.CopyModel('spike_response_neuron','readout_neuron', {'flush_PSP_after_firing':True})
    
    # auxiliary_neuron:
    # use it to make sure every neuron in the network has at least
    # 1 downstream neuron connected with 'weight/joint_stdp_synapse'.
    
    # This is because jointstdpmodule uses a different mechanism to update weights and delays.
    # Note that NEST handles spike event at presynaptic site and uses dendritic delay.
    # However, joint stdp is based on axonal delay, which handles spike event at postsynaptic site.
    # To implement handling spike event at postsynaptic site, when a spike event happens at presynaptic site,
    # we call event.source.update_upstream_weights/delays.
    # As a result, if a neuron has no downstream neuron, the call will not be executed, hence the weights/delays
    # will not be updated.
    nest.CopyModel('spike_response_neuron','auxiliary_neuron', { 'V_epsp':0.01 })  # V_epsp is small enough, never fires

    ## synapse models

    nest.SetDefaults('static_synapse',{
            'max_delay':d_max,
            'min_delay':d_min,
        })
        
    nest.SetDefaults('stdp_synapse',{ # additive STDP
        'Wmax':w_max,
        'tau_plus':tau_p,
        'lambda':A_n,
        'alpha':A_p/A_n,
        'mu_plus':0.0,
        'mu_minus':0.0
    })

    # synapse with weight STDP (the ``normal'' STDP), but with axonal delay
    nest.SetDefaults('weight_stdp_synapse',
        {
            #'Wmax':w_max,  # Wmax is hard-coded to 1.0
            'tau_plus':tau_p,
            'tau_minus':tau_n,
            'A_plus':A_p,
            'A_minus':A_n,
            'w_in':w_in,
            'w_out':w_out,
            'STDP_weight_effective_len':STDP_effective_len,  # the effective length of weight STDP window is STDP_weight_effective_len*tau
        })

    # synapse with joint STDP, with axonal delay
    nest.SetDefaults('joint_stdp_synapse',
        {
            'tau_plus':tau_p,
            'tau_minus':tau_n,
            'A_plus':A_p,
            'A_minus':A_n,
            #'Wmax':w_max,  # Wmax is hard-coded to 1.0
            'w_in':w_in,
            'w_out':w_out,
            'STDP_weight_effective_len':STDP_effective_len,  # the effective length of weight STDP window is STDP_weight_effective_len*tau
            'delay_min':d_min,'min_delay':d_min,###### yes they are duplicate, but pls do like this
            'max_delay':d_max,'delay_max':d_max,
            'sigma_plus':sigma_p,
            'sigma_minus':sigma_n,
            'B_plus':B_p,
            'B_minus':B_n,
            'delta_plus':delta_p,
            'delta_minus':delta_n,
            'nearest_approximation':True,
        })


if __name__=='__main__':
    
    ## define parameters
    
    #synapse_type = 'weight_stdp_synapse'
    synapse_type = 'joint_stdp_synapse'
    
    # weight STDP parameters
    tau_p = 16.8  # [ms], weight STDP time constant on the positive side
    tau_n = 33.7  # [ms], weight STDP time constant on the negative side
    A_p = 0.85*2**-6  # weight STDP amplitude on the positive side
    A_n = 2**-6  # weight STDP amplitude on the negative side
    w_in = 2**-10
    w_out = -w_in
    w_max = 1.0  # weight is dimensionless scaling factor of PSP
    #w_min = 0.0  # minimum weight is always 0
    STDP_effective_len = 5.
    
    # delay STDP parameters
    sigma_p = 16.8  # [ms], delay STDP time constant on the positive side
    sigma_n = 16.8  # [ms], delay STDP time constant on the negative side
    
    B_p = 1.  # [ms], delay STDP amplitude on the positive side
    B_n = B_p  # [ms], delay STDP amplitude on the negative side
    
    d_max = 35.0  # [ms], maximum axonal delay
    d_min = 1.0  # [ms], minimum axonal delay
    delta_p = -0.01
    delta_n = 0.01
    
    w_init = w_max/2.  # initial weight
    d_init = (d_min+d_max)/2.  # initial delay
    
    # neuron parameters
    tau_epsp = 3.  # [ms], EPSP time constant
    tau_ipsp = 3.  # [ms], IPSP time constant
    tau_reset = 12.0  # [ms], reset time constant
    
    thres = 3.0*w_max  # [mV], firing threshold
    inh = 3.  # lateral inhibition weight
    
    u_reset = thres  # reset potential is -thres (resting potential is 0)
    t_abs_refrac = 5.0  # [ms], period of absolute refractory
        
    # simulation parameters
    print_time = False
    verbose = 'M_ERROR'  # NEST verbose level
    vp_num = 2  # number of virtual processes to run the simulation
    time_step = 1.  # [ms]
    random_seed = 1
    sim_time = 100.0  # [ms]

    ## apply the above parameters
    
    reset_all()
    np.random.seed(random_seed)
        
    ## create two neurons, spike sources and spike detectors

    # the structure is:
    # spike_source           spike_source
    #     |                       |
    #  (static)                (static)
    #     |                       |
    #     v                       v
    # neuron1 --(weight STDP)--> neuron2 --(static)--> auxiliary_neuron
    #(parrot_neuron)      (spike_response_neuron)
    spike_sources = nest.Create('spike_generator',2)
    neuron1 = nest.Create('input_neuron',1)
    neuron2 = nest.Create('readout_neuron',1)
    nest.SetStatus(neuron2,{'V_epsp':0.1, 'V_ipsp':0.1})
    auxiliary_neuron = nest.Create('auxiliary_neuron',1)  # why auxiliary neuron? see comments in reset_all()
    spike_detectors = nest.Create('spike_detector', 2, params={'withgid':True,'withtime':True})
    
    # connect spike sources to neurons
    nest.Connect(spike_sources, [neuron1[0],neuron2[0]],
            conn_spec={'rule':'one_to_one','autapses':False,'multapses':False},
            syn_spec={'model':'static_synapse','weight':w_max,'delay':d_min})

    # connect neuron to neuron
    nest.Connect(neuron1, neuron2,
            conn_spec={'rule':'one_to_one','autapses':False,'multapses':False},
            syn_spec={'model':synapse_type,'weight':w_init,'delay':d_init})

    # this is a compulsory step!!
    # make sure every neuron in the network has at least
    # 1 downstream neuron connected with 'weight/joint_stdp_synapse'.
    
    # This is because jointstdpmodule uses a different mechanism to update weights and delays.
    # Note that NEST handles spike event at presynaptic site and uses dendritic delay.
    # However, joint stdp is based on axonal delay, which handles spike event at postsynaptic site.
    # To implement handling spike event at postsynaptic site, when a spike event happens at presynaptic site,
    # we call event.source.update_upstream_weights/delays.
    # As a result, if a neuron has no downstream neuron, the call will not be executed, hence the weights/delays
    # will not be updated.
    nest.Connect(neuron2, auxiliary_neuron,
            conn_spec={'rule':'all_to_all','autapses':False,'multapses':False},
            syn_spec={'model':synapse_type,'weight':w_max,'delay':d_min})
    
    # this is a compulsory step!!
    # tell the upstream neuron ids and connection types
    # to all neurons with incoming synpases of weight_stdp_synapse or joint_stdp_synapse
    
    # because NEST performs STDP in presynaptic site.
    # In order to implement postsynaptic update and axonal delays,
    # we have to store a copy of its upstream neuron ids.
    if synapse_type=='joint_stdp_synapse' or synapse_type=='weight_stdp_synapse':
        if synapse_type=='joint_stdp_synapse':
            st = JOINT_STDP_TYPE
        else:
            st = WEIGHT_STDP_TYPE
        nest.SetStatus(neuron2,{
            'upstream_neuron_ids':neuron1, 'upstream_connection_types':[st]
        })
        nest.SetStatus(auxiliary_neuron,{
            'upstream_neuron_ids':neuron2, 'upstream_connection_types':[st]
        })
    
    # connect spike detectors to neurons
    nest.Connect( [neuron1[0],neuron2[0]], spike_detectors, conn_spec={'rule':'one_to_one'})

    ## set input spikes
    
    input_spikes = [[1.],
                    [20.]]
    nest.SetStatus([spike_sources[0]],{'spike_times':input_spikes[0],'spike_weights':[w_max]*len(input_spikes[0])})
    nest.SetStatus([spike_sources[1]],{'spike_times':input_spikes[1],'spike_weights':[100*w_max]*len(input_spikes[1])})
    
    ## run it!!
    
    weights_history = []
    delays_history = []
    # record weights and delays at the beginning
    weights_history.append(get_weights(neuron1, neuron2))
    delays_history.append(get_delays(neuron1, neuron2))

    nest.Simulate(sim_time)
    weights_history.append(get_weights(neuron1, neuron2))
    delays_history.append(get_delays(neuron1, neuron2))
            
    # record spikes
    ss = nest.GetStatus(spike_detectors, 'events')
    spike_history = [ss[i]['times'] for i in range(2)]
    
    ## print weights and delays at the end of simulation
    print('\n*** validation of \'%s\' ***'%synapse_type)
    
    print('\nneuron1 fired at %gms, neuron2 fired at %gms'%(spike_history[0][0],spike_history[1][0]))
    print('the spike from neuron1 arrived at neuron2 %gms later at %gms'%(delays_history[0][0,0],spike_history[0][0]+delays_history[0][0,0]))
    dt = spike_history[0][0]+delays_history[0][0,0]-spike_history[1][0]
    print('the difference (t_1+delay-t_2) is %gms'%dt)
    
    print('\nweight at the beginning of simulation:')
    print(weights_history[0][0,0])
    
    print('\nweight at the end of simulation should be:')
    if dt<0:
        ww = weights_history[0][0,0]+w_in+w_out+A_p*np.exp(dt/tau_p)
        print('w_init + w_in + w_out + A_p*exp(dt/tau_p) = %g'%ww)
    else:
        ww = weights_history[0][0,0]+w_in+w_out-A_n*np.exp(-dt/tau_n)
        print('w_init + w_in + w_out - A_n*exp(-dt/tau_n) = %g'%ww)
    print('weight from JointSTDP Module is: %g'%weights_history[1][0,0])
    print('diff < 1e-10?: %s'%(np.abs(ww-weights_history[1][0,0])<1e-10,))
    
    if synapse_type=='joint_stdp_synapse':
        print('\ndelay at the beginning of simulation:')
        print(delays_history[0][0,0])
        print('\ndelay at the end of simulation should be:')
        if dt<0:
            dd = delays_history[0][0,0]-(w_max-weights_history[0][0,0]+delta_n)/w_max*B_n*(sigma_n/np.e+dt*np.exp(dt/sigma_n))
            print('d_init - (w_max - w + delta_n) / w_max * B_n * (sigma_n / e + dt * exp(dt / sigma_n)) = %g'%dd)
        else:
            dd = delays_history[0][0,0]+(weights_history[0][0,0]-delta_p)/w_max*B_p*(sigma_p/np.e-dt*np.exp(-dt/sigma_p))
            print('d_init + (w - delta_p) / w_max * B_p * (sigma_p / e - dt * exp( -dt / sigma_p)) = %g'%dd)
        dd_r = np.round(dd/time_step)*time_step
        print('after rounding to simulation time step (%gms): %g'%(time_step, dd_r))
        print('delay from JointSTDP Module is: %g'%delays_history[1][0,0])
        print('equal?: %s'%(dd_r==delays_history[1][0,0]))
    else:
        print('\ndelay at the beginning of simulation:')
        print(delays_history[0][0,0])
        print('delay at the end of simulation:')
        print(delays_history[1][0,0])
    
    # expected output for synapse_type=='weight_stdp_type':
    """
    *** validation of 'weight_stdp_synapse' ***

    neuron1 fired at 2ms, neuron2 fired at 22ms
    the spike from neuron1 arrived at neuron2 18ms later at 20ms
    the difference (t_1+delay-t_2) is -2ms

    weight at the beginning of simulation:
    0.5

    weight at the end of simulation should be:
    w_init + w_in + w_out + A_p*exp(dt/tau_p) = 0.511791
    weight from JointSTDP Module is: 0.511791
    diff < 1e-10?: True

    delay at the beginning of simulation:
    18.0
    delay at the end of simulation:
    18.0
    """
    
    # expected output for synapse_type=='joint_stdp_type':
    """
    *** validation of 'joint_stdp_synapse' ***

    neuron1 fired at 2ms, neuron2 fired at 22ms
    the spike from neuron1 arrived at neuron2 18ms later at 20ms
    the difference (t_1+delay-t_2) is -2ms

    weight at the beginning of simulation:
    0.5

    weight at the end of simulation should be:
    w_init + w_in + w_out + A_p*exp(dt/tau_p) = 0.511791
    weight from JointSTDP Module is: 0.511791
    diff < 1e-10?: True

    delay at the beginning of simulation:
    18.0

    delay at the end of simulation should be:
    d_init - (w_max - w + delta_n) / w_max * B_n * (sigma_n / e + dt * exp(dt / sigma_n)) = 15.7535
    after rounding to simulation time step (1ms): 16
    delay from JointSTDP Module is: 16
    equal?: True
    """

    # that's all folks.

