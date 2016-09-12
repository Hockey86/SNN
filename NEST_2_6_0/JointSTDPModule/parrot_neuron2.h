/*
 *  parrot_neuron2.h
 *
 *  This file is part of JointSTDPModule.
 *
 */


  /* BeginDocumentation
Name: parrot_neuron2 - Neuron that repeats incoming spikes, inherited from my_base_neuron.

Description:

  The parrot neuron simply emits one spike for every incoming spike.
  One possible application for this is to create different channels
  for the output of generator devices such as the poisson_generator
  or the mip_generator.
  
  Inherited from my_base_neuron, it stores its upstream connections
  and downstream neurons with immediate inhibition connection.

Remarks:

  Network-wise the parrot neuron behaves like other neuron models
  regarding connections and communication. While the number of
  outgoing spikes equals that of incoming ones, the weigth of the 
  outgoing spikes solely depends on the weigth of outgoing connections.

  A Poisson generator that would send multiple spikes during a single
  time step due to a high rate will send single spikes with 
  multiple synaptic strength instead, for effiacy reasons.
  This can be realized because of the way devices are implemented
  in the threaded environment. A parrot neuron on the other 
  hand always emits single spikes. Hence, in a situation where for 
  example a poisson generator with a high firing rate is connected
  to a parrot neuron, the communication cost associated with outgoing
  spikes is much bigger for the latter.

Receives: SpikeEvent
  
Sends: SpikeEvent
  
Parameters: 

  No parameters to be set in the status dictionary.

References:
  No references

Author:  May 2006, Reichert, Morrison

Modified: Nov 2015, Sun, Haoqi
*/



/**
 * The parrot neuron emits one spike for every incoming spike.
 * It is a (strongly) simplified version of the iaf_neuron class,
 * stripped of the dynamics and unneeded features.
 * Instead of the accumulated weigths of the incoming spikes the
 * number of the spikes is stored within a ring buffer.
 * \author David Reichert
 * \date may 2006 
 */

#ifndef PARROT_NEURON2_H
#define PARROT_NEURON2_H

#include "nest.h"
#include "event.h"
#include "my_base_neuron.h"
#include "ring_buffer.h"
#include "connection.h"

//#include "nestmodule.h"
//#include "target_identifier.h"

namespace mynest
{
  //class nest::Network;

  class parrot_neuron2: public my_base_neuron
  {
    
  public:        
    
    parrot_neuron2();
    parrot_neuron2(const parrot_neuron2&);

    /**
     * Import sets of overloaded virtual functions.
     * @see Technical Issues / Virtual Functions: Overriding,
     * Overloading, and Hiding
     */
    using nest::Node::handle;
    using nest::Node::handles_test_event;
    bool has_proxies()    const { return false; }
    bool local_receiver() const { return true;  }

    nest::port send_test_event(nest::Node&, nest::rport, nest::synindex, bool);
    
    void handle(nest::SpikeEvent &);
    nest::port handles_test_event(nest::SpikeEvent &, nest::rport);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
    void immediate_inhibit(nest::double_t) {} // parrot_neuron cannot be immediately inhibited, because it's super strong.
    void immediate_excite(nest::double_t) {} // it cannot be immediately excited either

  private:
      
    void init_state_(const nest::Node&){}  // no state
    void init_buffers_();
    void calibrate(){}  // no variables
    
    void update(nest::Time const &, const nest::long_t, const nest::long_t);
    
    /**
       Buffers and accumulates the number of incoming spikes per time step; 
       RingBuffer stores doubles; for now the numbers are casted.
    */
    struct Buffers_ {
      nest::RingBuffer n_spikes_;
    };
    
    Buffers_ B_;
  };

  inline
  nest::port parrot_neuron2::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
  {
  nest::SpikeEvent e;
  e.set_sender(*this);
  
  return target.handles_test_event(e, receptor_type);
  }

  inline
  nest::port parrot_neuron2::handles_test_event(nest::SpikeEvent&, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
  
} // namespace

#endif //PARROT_NEURON2_H
