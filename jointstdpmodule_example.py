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
"""

import numpy as np
import nest
# the jointstdpmodule does not support nest.topology module yet, you have to arrange the topology of neurons by your own
#import nest.topology as tp
nest.Install('jointstdpmodule')


def compute_PSP_norm(taumem, cmem, tausyn):
    """ compute the \varepsilon_max, which is
    the maximum of the postsynaptic response induced by a single maximally strengthened spike) """
    a = taumem / tausyn
    b = 1.0 / tausyn - 1.0 / taumem
    # time of maximum
    nest.sli_push(-np.exp(-1.0/a)/a)
    nest.sli_run('LambertWm1')
    t_max = 1.0/b * ( -nest.sli_pop() - 1.0/a )
    # maximum of PSP for current of unit amplitude
    return np.e/(tausyn*cmem*b) * ((np.exp(-t_max/taumem) - np.exp(-t_max/tausyn)) / b - t_max*np.exp(-t_max/tausyn))
    
    
def get_weights(sources, targets=None):
    """helper function to get weights as numpy array, where no connection is denoted by NO_CONNECTION"""
    if targets is None:
        targets = sources
    weights = np.zeros((len(sources),len(targets)))+NO_CONNECTION
    conn = nest.GetConnections(source=sources,target=targets)
    nest_weights = nest.GetStatus(conn,'weight')
    for i in range(len(conn)):
        weights[conn[i][0]-1,conn[i][1]-1] = nest_weights[i]
    return weights


def get_delays(sources, targets=None):
    """helper function to get delays as numpy array, where no connection is denoted by NO_CONNECTION"""
    if targets is None:
        targets = sources
    delays = np.zeros((len(sources),len(targets)))+d_min
    conn = nest.GetConnections(source=sources,target=targets)
    nest_delays = nest.GetStatus(conn,'delay')
    for i in range(len(conn)):
        delays[conn[i][0]-1,conn[i][1]-1] = nest_delays[i]
    return delays


def reset_all():
    """reset the simulator kernel and define the models"""
    nest.set_verbosity(verbose)
    nest.ResetKernel()
    global print_time
    global vp_num
    if debug:
        print_time = False
        vp_num = 1
    np.random.seed(random_seed)
    pyrngs = [np.random.RandomState(i) for i in range(random_seed, random_seed+vp_num)]
    nest.SetKernelStatus({'resolution':time_step,'print_time':print_time,"local_num_threads": vp_num,'grng_seed':random_seed+vp_num,'rng_seeds':range(random_seed+vp_num+1, random_seed+2*vp_num+1)})

    nest.SetDefaults('static_synapse',{
            'max_delay':d_max,  # [ms]
            'min_delay':time_step,  # [ms]
        })
    nest.SetDefaults('iaf_psc_alpha2',{ # this is same with iaf_psc_alpha, but stores presynaptic delays and spike history
            'C_m':C_m,  # membrane capacity [pF]
            'tau_m':tau_m,  # membrane time constant [ms]
            'tau_syn_ex':tau_epsp,  # rise time of the excitatory synaptic alpha function [ms]
            'tau_syn_in':tau_ipsp,  # rise time of the inhibitory synaptic alpha function [ms]
            'E_L':0.0,  # resting potential [mV], for simplicity we use 0
            'V_m':0.0,  # potential [mV]
            'tau_minus':tau_n,  # time constant in the negative side of the weight stdp learning window
        })

    # if we directly connect spike sources to readout neurons, the synapse cannot be plastic,
    # therefore the spike sources are first connected to input neurons (which has very low threshold),
    # and then the input neurons are connected to readout neurons, where the synapses are joint_stdp_synapse
    nest.CopyModel('iaf_psc_alpha2','input_neuron',
        {
            't_ref':1.0,  # efractory period [ms]
            'V_reset':0.0,  # reset potential [mV]
            'V_th':0.9*w_max*epsilon_max  # threshold [mV]
        })

    nest.CopyModel('iaf_psc_alpha2','readout_neuron', # readout neurons are also normal iaf neurons
        {
            't_ref':readout_abs_refrac,
            'V_reset':readout_reset_potential,
            'V_th':15.0
        })

    # because NEST only process synaptic changes in presynaptic site, if a neuron's axon is not connected to any other neuron, its synapses will not be updated,
    # therefore we append an auxiliary_neuron to neurons without downstream neurons, to make sure every neuron in the network has downstream neurons.
    nest.CopyModel('iaf_psc_alpha2','auxiliary_neuron',
        { 'V_th':15.0 })

    # synapse with weight (the classic) STDP. Different from NEST, the delays are all axonal.
    nest.CopyModel('weight_stdp_synapse','readout_synapse_without_delay_plasticity',
        {
            'Wmax':w_max,  # [a.u.]
            'tau_plus':tau_p,  # time constant in the postive side of the weight learning window
            'tau_minus':tau_n,  # time constant in the negative side of the weight learning window
            'A_plus':readout_A_p,  # amplitude in the postive side of the weight learning window
            'A_minus':readout_A_n,  # amplitude in the negative side of the weight learning window
            'a_pre':readout_w_in,  # first order Hebbian term in weight for presynaptic spike
            'a_post':readout_w_out,  # first order Hebbian term in weight for postsynaptic spike
            'STDP_weight_window_plus':readout_STDP_weight_window_plus,  # the effective length of the positive side of the weight learning window is STDP_weight_window_plus*tau_plus
            'STDP_weight_window_minus':readout_STDP_weight_window_minus,  # the effective length of the negative side of the weight learning window is STDP_weight_window_minus*tau_minus
            'debug':debug  # debug mode, but the debug output is commented out in model source file, you can recover them and rebuild the model
        })

    # synapse with joint plasticity. Different from NEST, the delays are all axonal.
    nest.CopyModel('joint_stdp_synapse','readout_synapse',
        {
            'tau_plus':tau_p,
            'tau_minus':tau_n,
            'A_plus':readout_A_p,
            'A_minus':readout_A_n,
            'Wmax':w_max,
            'a_pre':readout_w_in,
            'a_post':readout_w_out,
            'STDP_weight_window_plus':readout_STDP_weight_window_plus,
            'STDP_weight_window_minus':readout_STDP_weight_window_minus,
            'delay_min':d_min,######### these two seem duplicate, but pls do like this
            'min_delay':d_min,######### 
            'max_delay':d_max,
            'sigma_plus':sigma_p,  # time constant in the postive side of the delay learning window
            'sigma_minus':sigma_n,  # time constant in the negative side of the delay learning window
            'B_plus':readout_B_p,  # amplitude in the postive side of the delay learning window
            'B_minus':readout_B_n,  # amplitude in the negative side of the delay learning window
            'STDP_delay_window_plus':readout_STDP_delay_window_plus,  # the effective length of the positive side of the delay learning window is STDP_delay_window_plus*sigma_plus
            'STDP_delay_window_minus':readout_STDP_delay_window_minus,  # the effective length of the negative side of the delay learning window is STDP_delay_window_minus*sigma_minus
            'debug':debug
        })
        
def run_snn(weights=None,delays=None,synapse_type='readout_synapse',record_history=False,record_spike=None,record_potential=None):
    """build the network and run"""
    global debug, verbose, print_time
    if record_history:
        print_time = False
    if record_spike is None:
        record_spike = record_history
    if record_potential is None:
        record_potential = record_history
    verbose = 'M_ERROR'
    
    reset_all()

    # build input and readout neurons
    input_sources = nest.Create('spike_generator',input_num)  # spike sources
    input_neurons = nest.Create('input_neuron',input_num)  # input neurons, see comments in reset_all()
    readout_neurons = nest.Create('readout_neuron',readout_num)
    auxiliary_neurons = nest.Create('auxiliary_neuron',readout_num)  # auxiliary neurons, see comments in reset_all()

    # set readout threshold
    for i in range(readout_num):
        nest.SetStatus([readout_neurons[i]], {'V_th':readout_thresholds[i]})

    # connect sources to input layer
    nest.Connect(input_sources, input_neurons,
            conn_spec={'rule':'one_to_one','autapses':False,'multapses':False},
            syn_spec={'model':'static_synapse','weight':w_max,'delay':d_min})

    # connect readout layer to input layer
    if weights is None or delays is None:
        # if weights or delays is not specified,
        # the weights take w_init=w_max
        # the delays take uniform random value from [d_min,d_max/2]
        for i in range(input_num):
            for j in range(readout_num):
                nest.Connect([input_neurons[i]], [readout_neurons[j]],
                        conn_spec={'rule':'all_to_all','autapses':False,'multapses':False},
                        syn_spec={'model':synapse_type,'weight':w_init,'delay':{'distribution':'uniform','low':d_min,'high':d_max/2}})
    else:
        for i in range(input_num):
            for j in range(readout_num):
                nest.Connect([input_neurons[i]], [readout_neurons[j]],
                        conn_spec={'autapses':False,'multapses':False},
                        syn_spec={'model':synapse_type,'weight':weights[i,j],'delay':delays[i,j]})

    # lateral inhibition connection within readout layer
    if inhibitory_strength > 0.0:
        nest.Connect(readout_neurons, readout_neurons,
                conn_spec={'rule':'all_to_all','autapses':False,'multapses':False},
               syn_spec={'model':'static_synapse','weight':-inhibitory_strength*w_max,'delay':time_step})

    # connect readout layer to auxiliary neurons
    nest.Connect(readout_neurons, auxiliary_neurons,
            conn_spec={'rule':'one_to_one','autapses':False,'multapses':False},
            syn_spec={'model':synapse_type,'weight':w_max,'delay':d_min})
            
    weights_readout = get_weights(input_sources+input_neurons+readout_neurons)
    delays_readout = get_delays(input_sources+input_neurons+readout_neurons)
    
    # tell readout neurons its upstream neuron ids and delays
    # this is a very important step!!
    # because NEST updates the weights in presynaptic sites, and the delays are all dendritic,
    # to implement postsynaptic update and axonal delays, we have to use iaf_psc_neuron2 to store a copy of its presynaptic neuron id and delays,
    # so that they can store the presynaptic spike history
    # during the update made by NEST in presynaptic sites, it updates its presynaptic neurons according to the stored values
    # (I am the postsynaptic neuron of my presynaptic neurons)
    for i in range(readout_num):
        nest.SetStatus([readout_neurons[i]],{
            'upstream_neuron_ids':np.where(weights_readout[:2*input_num,2*input_num+i]!=NO_CONNECTION)[0],
            'upstream_delays':delays_readout[weights_readout[:2*input_num,2*input_num+i]!=NO_CONNECTION,2*input_num+i]
        })
        nest.SetStatus([auxiliary_neurons[i]],{
            'upstream_neuron_ids':[2*input_num+i], 'upstream_delays':[d_min]
        })
        
    if record_spike:
        sds = nest.Create('spike_detector',input_readout_num,params={'withgid': True, 'withtime': True})
        nest.Connect(input_neurons+readout_neurons, sds,conn_spec={'rule':'one_to_one'})
    if record_potential:
        mm = nest.Create('multimeter',1,params={'withtime':True,'interval':time_step,'record_from':['V_m']})[0]
        nest.Connect(mm,readout_neurons[0])

    # generate the input spikes according to training_patterns
    input_spikes = [[] for i in range(input_num)]
    for i in range(sample_num):
        pp = training_patterns[i]
        for j in range(len(pp)):
            input_spikes[int(pp[j,0])].append(pp[j,1]+i*trial_time)
    for i in range(input_num):
        nest.SetStatus([input_sources[i]],{'spike_times':input_spikes[i],'spike_weights':[w_max]*len(input_spikes[i])})
    
    spike_history = []
    potential_history = []
    weights_history = []
    delays_history = []
    weights_history.append(get_weights(input_sources+input_neurons+readout_neurons)[input_num:2*input_num,2*input_num:])
    delays_history.append(get_delays(input_sources+input_neurons+readout_neurons)[input_num:2*input_num,2*input_num:])
    
    if record_history:  # record the weights and delays after each trial
        #pdb.set_trace()
        #print(0,'before')
        #print('\nweights')
        #print(weights_history[-1].T)
        #print('\ndelays'  )
        #print(delays_history[-1].T)
        for k in range(sample_num):
            nest.Simulate(trial_time)  # run the network for a trial
            """
            print(k,'after')
            ss = nest.GetStatus(sds,'events')
            for i in range(input_num,len(ss)):
                print('%d %s'%(i+1, np.round([j for j in ss[i]['times'] if trial_time*k<=j<trial_time*(k+1)],1)))
            """
            weights_history.append(get_weights(input_sources+input_neurons+readout_neurons)[input_num:2*input_num,2*input_num:])
            delays_history.append(get_delays(input_sources+input_neurons+readout_neurons)[input_num:2*input_num,2*input_num:])
            
            """
            print('\nweights')
            print(weights_history[-1].T)
            print('\ndelays'  )
            print(delays_history[-1].T)
            print('\narrival time, pattern #%d: %s'%(pattern_order[k],np.round([ss[i]['times'][k]-k*trial_time for i in range(input_num)],1)))
            print(delays_history[-1].T+[ss[i]['times'][k]-k*trial_time for i in range(input_num)])
            """
            if k%50==0:
                print('trial %d/%d'%(k+1,sample_num))
                
    else:  # if not record history, weights/delays_history has 2 elements, the first is before simulation, the second is after simulation
        nest.Simulate(trial_time*sample_num)
        weights_history.append(get_weights(input_sources+input_neurons+readout_neurons)[input_num:2*input_num,2*input_num:])
        delays_history.append(get_delays(input_sources+input_neurons+readout_neurons)[input_num:2*input_num,2*input_num:])

    if record_potential:
        dmm = nest.GetStatus([mm],keys='events')[0]
        potential_history = (dmm['times'],dmm['V_m'])
    if record_spike:
        ss = nest.GetStatus(sds,'events')
        spike_history = [ss[i]['times'] for i in range(input_readout_num)]
        
    return weights_history,delays_history,spike_history,potential_history
    

print_time = False
debug = False  # but the debug output is commented out in model source file, you can recover them and rebuild the model
verbose = 'M_ERROR'
vp_num = 4
time_step = 0.1  # [ms]

random_seed = 0
np.random.seed(random_seed)
pyrngs = [np.random.RandomState(i) for i in range(random_seed, random_seed+vp_num)]

NO_CONNECTION = np.inf
d_min = 1.0
d_max = 30.0
w_max = 1.0
w_init = w_max

tau_m = 10.0
C_m = 250.0
tau_epsp = 1.0
tau_ipsp = 5.0
tau_p = 10.0
tau_n = 10.0
sigma_p = 10.0
sigma_n = 10.0
epsilon_max = compute_PSP_norm(tau_m,C_m,tau_epsp)

readout_w_in = 0.0
readout_w_out = 0.0
readout_A_p = 0.2
readout_A_n = 0.1
readout_B_p = 6.0
readout_B_n = 0.6
readout_STDP_weight_window_minus = readout_STDP_weight_window_plus = 5.0
readout_STDP_delay_window_minus = readout_STDP_delay_window_plus = 5.0

trial_time = 200.0  # [ms]

readout_abs_refrac = trial_time/2.0  # refractory period is set long enough to allow only one readout spike in each trial

training_patterns = []


########### now we are ready, the example is to show the convergence of weights and delays under single input pattern

# this is the input pattern
# [[neuron id #1, firing time #1 (ms)],
#  [neuron id #2, firing time #2 (ms)],
#  ...]
# note: this is timing of spike source, the firing of input neurons is spike source + rise time of input neuron (~4ms)
input_pattern1 = np.array([
        [0,21.0],
        [1,16.0],
        [2,11.0],
        [3,6.0],
        [4,1.0],
    ])
class_num = 1  # here we show single input pattern
input_patterns_mapping = [input_pattern1]  # for more patterns, this should be [input_pattern1, input_pattern2, ...]
    
input_num = 5  # number of input neurons
# number of readout neurons, but why put it in []?
# because you can specify readout neurons with different thresholds,
# for example, readout_each_threshold_num = [5] and readout_threshold_values = [readout_threshold1] means there are 5 readout neurons with readout_threshold1
# readout_each_threshold_num = [5,3] and readout_threshold_values = [readout_threshold1,readout_threshold2] means there are 5 readout neurons with readout_threshold1, and 3 readout neurons with readout_threshold2
readout_each_threshold_num = [5]
readout_num = sum(readout_each_threshold_num)
input_readout_num = input_num + readout_num

readout_threshold = 2.9*w_max*epsilon_max
readout_threshold_values = [readout_threshold]
readout_thresholds = []
for i in range(len(readout_threshold_values)):
    readout_thresholds.extend([readout_threshold_values[i]]*readout_each_threshold_num[i])
readout_reset_potential = -readout_threshold  # rest potential is 0, reset potential is -threshold

sample_num = 100  # 100 trials

inhibitory_strength = 0.1  # lateral inhibition strength, the actual inhibitory weight = -inhibitory_strength*w_max*epsilon_max

# generate the random order of patterns, in case of single pattern, it has no effect
pattern_order = []
for i in range(int(round(1.0*sample_num/class_num))):
    pattern_order.extend(range(class_num))
np.random.shuffle(pattern_order)

training_patterns = [np.array(input_patterns_mapping[pattern_order[i]],copy=True) for i in range(sample_num)]

# run it!!
weights_history,delays_history,_,_ = run_snn(record_history=True)

# now we have weights and delays after each trial
print('')
print('weights:')
print(weights_history[-1].T)  # the final weights, each row is the weights of a readout neuron
print('')
print('delays:')
print(delays_history[-1].T)  # the final delays, each row is the delays of a readout neuron

# expected output:
"""
trial 1/100
trial 51/100

weights:
[[ 1.  1.  1.  1.  1.]  # weights of readout neuron #1
 [ 0.  1.  1.  1.  1.]  # weights of readout neuron #2
 [ 1.  1.  1.  1.  1.]  # weights of readout neuron #3
 [ 1.  1.  1.  1.  1.]  # weights of readout neuron #4
 [ 1.  1.  1.  1.  1.]] # weights of readout neuron #5

delays:
[[  8.1  14.6   9.6  14.2  11.7]  # delays of readout neuron #1
 [ 30.    1.    6.   11.   16. ]  # delays of readout neuron #2
 [  6.8   1.4  13.8   8.6   4.8]  # delays of readout neuron #3
 [  7.6   8.    7.9   4.6  10.1]  # delays of readout neuron #4
 [  1.    6.   11.   16.   21. ]] # delays of readout neuron #5
 
# Therefore readout neuron #2 and #5 developed weights and delays with maximum spike arrival synchrony, from which we can recover the input.
# Other readout neurons are inhibited by lateral inhibition.
"""

# that's all folks.
    
