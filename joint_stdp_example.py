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

joint_stdp_example.py

An example to demonstrate the joint STDP under input patterns without noise.
The network structure consists of input and readout layers, which are fully connected using joint STDP,
and there are static lateral inhibition in the readout layer (see Figure 2 in Section 3).

Haoqi Sun, Sep 10, 2016
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
    
    ##### the relativity between B_p and A_p needs careful tuning !!!
    B_p = 0.1  # [ms], delay STDP amplitude on the positive side
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
    
    ##### thres and inh need careful tuning !!!
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
    trial_time = 100.0  # [ms]

    ## apply the above parameters
    
    reset_all()
    np.random.seed(random_seed)
    
    ## define input pattern(s)

    # [
    #  [neuron id, firing time (ms)],
    #  [neuron id, firing time (ms)],
    #  ...]
    # note: this is timing of spike source, the firing of input neurons is spike source + d_min
    input_pattern1 = np.array([
            [0,1.0],
            [1,3.0],
            [2,5.0],
            [3,7.0],
            [4,9.0],
            [5,11.0],
            [6,13.0],
            [7,15.0],
            [8,17.0],
            [9,19.0],
        ])
    input_pattern2 = np.array([
            [0,9.0],
            [1,5.0],
            [2,1.0],
            [3,5.0],
            [4,9.0],
            [5,5.0],
            [6,1.0],
            [7,5.0],
            [8,9.0],
            [9,5.0],
        ])
    input_patterns_mapping = [input_pattern1,input_pattern2]  # for more patterns, this should be [input_pattern1, input_pattern2, ...]
    pattern_num = len(input_patterns_mapping)
    trial_num = 400
    
    # generate the random order of patterns, in case of single pattern, it has no effect
    pattern_order = []
    for i in range(int(round(1.0*trial_num/pattern_num))):
        pattern_order.extend(range(pattern_num))
    pattern_order = np.array(pattern_order)
    np.random.shuffle(pattern_order)
    patterns = [np.array(input_patterns_mapping[pattern_order[i]],copy=True) for i in range(trial_num)]
        
    ## create network

    # the structure is:
    #                                    (static lateral inhibition)
    #                                             +-----+
    #             (static)         (joint STDP)   v     v    (static)
    # spike_sources ----> input_neurons ----> readout_neurons ----> auxiliary_neurons
    
    # why not spike_sources --> readout_neurons?
    # because you cannot use plastic synapse on spike_sources, which is a device in NEST
    input_num = 10
    readout_num = 2
    
    spike_sources = nest.Create('spike_generator',input_num)
    input_neurons = nest.Create('input_neuron',input_num)
    readout_neurons = nest.Create('readout_neuron',readout_num)
    auxiliary_neuron = nest.Create('auxiliary_neuron',1)  # why auxiliary neuron? see comments in reset_all()

    # set readout threshold
    # SRM in NEST has normalized threshold of 1.0, so we set V_epsp instead
    nest.SetStatus(readout_neurons,
        {'V_epsp':1./thres, 'V_ipsp':1./thres})

    # connect spike sources to input layer
    nest.Connect(spike_sources, input_neurons,
            conn_spec={'rule':'one_to_one','autapses':False,'multapses':False},
            syn_spec={'model':'static_synapse','weight':w_max,'delay':d_min})

    # connect readout layer to input layer
    nest.Connect(input_neurons, readout_neurons,
            conn_spec={'rule':'all_to_all','autapses':False,'multapses':False},
            syn_spec={'model':synapse_type,'weight':w_init,'delay':d_init})

    # lateral inhibition connection within readout layer
    # the synapses for lateral inhibition must have zero delay,
    # otherwise it cannot inhibit other neurons fired at the same time
    # however, every synapse in NEST must have a positive minimum delay,
    # therefore we created 'immediate_inhibition_synapse'.
    # In order to enable this synapse,
    # you must specify 'downstream_neuron_ids'
    # See the example below.
    if readout_num>1:
        nest.Connect(readout_neurons, readout_neurons,
                conn_spec={'rule':'all_to_all','autapses':False,'multapses':False},
               syn_spec={'model':'immediate_inhibition_synapse','weight':-abs(inh)})
    
        # this is a compulsory step!!
        # 'downstream_neuron_ids' indicate which neurons to immediately inhibit
        for i in range(readout_num):
            nest.SetStatus([readout_neurons[i]],{
                'downstream_neuron_ids':list(readout_neurons[:i]+readout_neurons[i+1:])
            })

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
    nest.Connect(readout_neurons, auxiliary_neuron,
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
        for i in range(readout_num):
            nest.SetStatus([readout_neurons[i]],{
                'upstream_neuron_ids':list(input_neurons), 'upstream_connection_types':[st]*input_num
            })
        nest.SetStatus(auxiliary_neuron,{
            'upstream_neuron_ids':readout_neurons, 'upstream_connection_types':[st]*readout_num
        })
    
    # connect spike detectors to input neurons and readout neurons
    spike_detectors = nest.Create('spike_detector', input_num+readout_num, params={'withgid':True,'withtime':True})
    nest.Connect(input_neurons+readout_neurons, spike_detectors, conn_spec={'rule':'one_to_one'})

    ## generate input spikes according to patterns
    
    # it looks like
    # | p0    p1    p0
    # | *     *     *
    # |  *     *     *      ...
    # |   *   *       *
    # +------------------------> time
    # |trial1|trial2|trial3|...
    input_spikes = [[] for i in range(input_num)]
    for i in range(trial_num):
        pp = patterns[i]
        for j in range(pp.shape[0]):
            if int(pp[j,0])<input_num:
                input_spikes[int(pp[j,0])].append(pp[j,1]+i*trial_time)
    for i in range(input_num):
        nest.SetStatus([spike_sources[i]],{'spike_times':input_spikes[i],'spike_weights':[w_max]*len(input_spikes[i])})
    
    ## run it!!
    
    weights_history = []
    delays_history = []
    # record weights and delays at the beginning
    weights_history.append(get_weights(input_neurons, readout_neurons))
    delays_history.append(get_delays(input_neurons, readout_neurons))

    for k in range(trial_num):
        nest.Simulate(trial_time)  # run the network for a trial
        weights_history.append(get_weights(input_neurons, readout_neurons))
        delays_history.append(get_delays(input_neurons, readout_neurons))
        if k%100==0:
            print('trial %d/%d'%(k+1,trial_num))
            
    # record spikes
    ss = nest.GetStatus(spike_detectors, 'events')
    spike_history = [ss[i]['times'] for i in range(input_num+readout_num)]
    
    ## print weights and delays at the end of simulation
    
    np.set_printoptions(suppress=True,precision=4)
    
    print('\nweights at the end of simulation:')
    print(weights_history[-1].T)  # each row is the weights of a readout neuron
    print('\ndelays at the end of simulation:')
    print(delays_history[-1].T)  # each row is the delays of a readout neuron
    
    print('\npresynaptic spike arrival time (delay + pattern1) for synapses with weight>0.9,\nwhich row (readout neuron) has similar value?')
    presynaptic_arrival_time1 = (delays_history[-1].T+[input_pattern1[input_pattern1[:,0].tolist().index(i),1] if i in input_pattern1[:,0] else np.nan for i in range(input_num)])
    presynaptic_arrival_time1[weights_history[-1].T<=0.9]=np.nan
    print(presynaptic_arrival_time1)
    
    print('\npresynaptic spike arrival time (delay + pattern2) for synapses with weight>0.9\nwhich row (readout neuron) has similar value?')
    presynaptic_arrival_time2 = (delays_history[-1].T+[input_pattern2[input_pattern2[:,0].tolist().index(i),1] if i in input_pattern2[:,0] else np.nan for i in range(input_num)])
    presynaptic_arrival_time2[weights_history[-1].T<=0.9]=np.nan
    print(presynaptic_arrival_time2)
        
    ## hit rate and false alarm frequency in the period of the last 100 trials
    
    num_last_trial = 100
    T = num_last_trial*trial_time
    spike_history_T = [spike_history[k][spike_history[k]>=trial_num*trial_time-T] for k in range(len(spike_history))]
    
    # hit rate is the fraction of patterns with readout spike within its valid response interval
    # over all patterns during the period T.
    # hit_rate: the i-th row and the j-th column: the hit rate of pattern j for readout neuron i
    hit_rate = np.zeros((readout_num,pattern_num))
    for i in range(pattern_num):
        pattern_ids = np.where(np.logical_and(pattern_order==i,[False]*(trial_num-num_last_trial)+[True]*num_last_trial))[0]
        for k in range(readout_num):
            # here, the valid response interval of the i-th appearing pattern is [i*trial_time, (i+1)*trial_time]
            # so we can compare i*trial_time and round(spike_history[]/trial_time)*trial_time
            hit_ids = np.where(np.in1d(pattern_ids*trial_time, np.round(spike_history_T[input_num+k]/trial_time)*trial_time))[0]
            hit_rate[k,i] = hit_ids.shape[0]*1.0/pattern_ids.shape[0]
    print('\nhit rate of the last 100 trials:')
    print(hit_rate)
    
    # false alarm frequency is the number of responses out of valid response interval divided by a period T.
    # false_alarm_frequency: the i-th row and the j-th column: the false alarm frequency of pattern j for readout neuron i
    false_alarm_frequency = np.zeros((readout_num,pattern_num))
    for i in range(pattern_num):
        pattern_ids = np.where(np.logical_and(pattern_order==i,[False]*(trial_num-num_last_trial)+[True]*num_last_trial))[0]
        for k in range(readout_num):
            false_alarm_ids = np.where(np.in1d(np.round(spike_history_T[input_num+k]/trial_time)*trial_time, pattern_ids*trial_time, invert=True))[0]
            false_alarm_frequency[k,i] = false_alarm_ids.shape[0]*1.0/T
    print('\nfalse alarm frequency (Hz) of the last 100 trials:')
    print(false_alarm_frequency)

    # expected output:
    """
    trial 1/400
    trial 101/400
    trial 201/400
    trial 301/400

    weights at the end of simulation:
    [[ 0.  1.  1.  1.  0.      1.      1.      1.      0.      0.    ]
     [ 1.  1.  1.  1.  0.0013  0.0012  0.0011  0.0012  0.0013  0.0011]]

    delays at the end of simulation:
    [[ 18.  12.  14.  12.  13.  12.  14.  12.  20.  21.]
     [ 16.  15.  16.  15.  22.  22.  20.  19.  21.  16.]]

    presynaptic spike arrival time (delay + pattern1) for synapses with weight>0.9,
    which row (readout neuron) has similar value?
    [[ nan  15.  19.  19.  nan  23.  27.  27.  nan  nan]
     [ 17.  18.  21.  22.  nan  nan  nan  nan  nan  nan]] <---

    presynaptic spike arrival time (delay + pattern2) for synapses with weight>0.9
    which row (readout neuron) has similar value?
    [[ nan  17.  15.  17.  nan  17.  15.  17.  nan  nan] <---
     [ 25.  20.  17.  20.  nan  nan  nan  nan  nan  nan]]

    hit rate of the last 100 trials:
    [[ 0.  1.]
     [ 1.  0.]]

    false alarm frequency (Hz) of the last 100 trials:
    [[ 0.0054  0.    ]
     [ 0.      0.0046]]
    """

    # The expected output indicates that the first readout neuron learned the second pattern,
    # and the seoncd readout neuron learned the first pattern.
    # The presynaptic spike arrival time has similar values means that the weights and delays achieved maximum spike arrival synchrony, from which we can recover the input.
    # The fact that the presynaptic arrival times are not exactly same is due to the long tail of PSPs, which leads to some tolerance to the spike timing.

    # that's all folks.
