/*
 *  parrot_neuron2.h
 *
 *  This file is part of JointSTDPModule in NEST.
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Sun Haoqi, Jun 29, 2016, NTU, Singapore
 */


/* BeginDocumentation
Name: parrot_neuron2 - Neuron that repeats incoming spikes, inherited
from base_neuron in JointSTDPModule.

Description:

Same with parrot_neuron in NEST. The only difference is that it is
inherited from base_neuron, so it stores its upstream connections, and
downstream neurons for immediate inhibition connection.

The purpose of this model is to function as input neurons in a network
connected by weight_stdp_synapse and joint_stdp_synapse, because these
synapses can only connect neurons inherited from base_neuron in JointSTDPModule.

Receives: SpikeEvent

Sends: SpikeEvent

Parameters:
No parameters to be set in the status dictionary.

Author: David Reichert, Abigail Morrison, Alexander Seeholzer, Hans Ekkehard
Plesser
FirstVersion: May 2006

Modified:  Nov 2015, Sun, Haoqi
*/


/**
 * The parrot neuron emits one spike for every incoming spike,
 * but may use multiplicity to indicate number of spikes in a single
 * time step.
 * Instead of the accumulated weigths of the incoming spikes, the
 * number of the spikes is stored within a ring buffer.
 */

#ifndef PARROT_NEURON2_H
#define PARROT_NEURON2_H

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"

#include "base_neuron.h"

namespace mynest
{

class nest::Network;

//TODO should better use multiple inheritance or mixin interface?
//class parrot_neuron2 : public nest::Archiving_Node, public base_neuron
class parrot_neuron2 : public base_neuron
{

public:
  parrot_neuron2();
  //parrot_neuron2(const parrot_neuron2&);

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding,
   * Overloading, and Hiding
   */
  using nest::Node::handle;
  using nest::Node::handles_test_event;
  using nest::Node::sends_signal;
  using nest::Node::receives_signal;

  nest::port send_test_event( nest::Node&, nest::rport, nest::synindex, bool );
  nest::SignalType sends_signal() const;
  nest::SignalType receives_signal() const;

  void handle( nest::SpikeEvent& );
  nest::port handles_test_event( nest::SpikeEvent&, nest::rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );
  
  //bool has_proxies()    const { return false; }
  //bool local_receiver() const { return true;  }

  void immediate_inhibit(nest::double_t) {} // parrot_neuron cannot be immediately inhibited, because it's super strong.
  void immediate_excite(nest::double_t) {} // it cannot be immediately excited either

private:
  void
  init_state_( const nest::Node& )
  {
  } // no state
  void init_buffers_();
  void
  calibrate()
  {
  } // no variables

  void update( nest::Time const&, const nest::long_t, const nest::long_t );

  /**
     Buffers and accumulates the number of incoming spikes per time step;
     RingBuffer stores doubles; for now the numbers are casted.
  */
  struct Buffers_
  {
    nest::RingBuffer n_spikes_;
  };

  Buffers_ B_;
};

inline nest::port
parrot_neuron2::send_test_event( nest::Node& target,
  nest::rport receptor_type,
  nest::synindex,
  bool )
{
  nest::SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline nest::port
parrot_neuron2::handles_test_event( nest::SpikeEvent&, nest::rport receptor_type )
{
  // Allow connections to port 0 (spikes to be repeated)
  // and port 1 (spikes to be ignored).
  if ( receptor_type == 0 or receptor_type == 1 )
  {
    return receptor_type;
  }
  else
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }
}

inline nest::SignalType
parrot_neuron2::sends_signal() const
{
  return nest::ALL;
}

inline nest::SignalType
parrot_neuron2::receives_signal() const
{
  return nest::ALL;
}

} // namespace

#endif // PARROT_NEURON2_H
