/*
 *  parrot_neuron2.cpp
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

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"

#include <limits>

#include "parrot_neuron2.h"


namespace mynest
{

parrot_neuron2::parrot_neuron2()
  : base_neuron()
{
}

/*parrot_neuron2::parrot_neuron2( const parrot_neuron2& n )
  : base_neuron( n )
{
}*/

void
parrot_neuron2::init_buffers_()
{
  B_.n_spikes_.clear(); // includes resize
  base_neuron::clear_history();
  base_neuron::clear_upstream_spike_history();
}

void
parrot_neuron2::update( nest::Time const& origin, const nest::long_t from, const nest::long_t to )
{
  assert(
    to >= 0 && ( nest::delay ) from < nest::Scheduler::get_min_delay() );
  assert( from < to );

  for ( nest::long_t lag = from; lag < to; ++lag )
  {
    const nest::ulong_t current_spikes_n =
      static_cast< nest::ulong_t >( B_.n_spikes_.get_value( lag ) );
    if ( current_spikes_n > 0 )
    {
      // create a new SpikeEvent, set its multiplicity and send it
      nest::SpikeEvent se;
      se.set_multiplicity( current_spikes_n );
      network()->send( *this, se, lag );

      // set the spike times, respecting the multiplicity
      for ( nest::ulong_t i = 0; i < current_spikes_n; i++ )
      {
        set_spiketime( nest::Time::step( origin.get_steps() + lag + 1 ) );
      }
      
      // if there are downstream neurons for immediate inhibition, inhibit them
      for ( nest::ulong_t i = 0; i < baseparam_.downstream_immediate_inhibition_neuron_num_; ++i )
      {
        if ( baseparam_.downstream_immediate_inhibition_neurons_[i] != NULL )
        {
          //nest::SpikeEvent* myse = new nest::SpikeEvent(se);
          //myse->set_weight(baseparam_.downstream_weights_[it-baseparam_.downstream_immediate_inhibition_neurons_.begin()]);
          //myse->set_delay(0);
          //(*it)->handle(*myse);
          //delete myse;            
          baseparam_.downstream_immediate_inhibition_neurons_[i]->immediate_inhibit(baseparam_.downstream_weights_[i]);
        }
      }
    }
  }
}


inline void
parrot_neuron2::get_status( DictionaryDatum& d ) const
{
  // get our own parameter and state data
  def<double>(d, nest::names::t_spike, get_spiketime_ms());
  //baseparam_.get( d );

  // get information managed by parent class
  base_neuron::get_status( d ); ////// base_neuron::get_status not implemented, roll back to nest::Archiving_Node::get_status( d )
}

inline void
parrot_neuron2::set_status( const DictionaryDatum& d )
{
    //BaseParam_ bptmp = baseparam_;  // temporary copy in case of errors
    //bptmp.set( d );                 // throws if BadProperty
    
    base_neuron::set_status( d );   ////// base_neuron::sget_status not implemented, roll back to nest::Archiving_Node::set_status( d )
  
    //baseparam_ = bptmp;
}

void
parrot_neuron2::handle( nest::SpikeEvent& e )
{
  // Repeat only spikes incoming on port 0, port 1 will be ignored
  if ( 0 == e.get_rport() )
  {
    B_.n_spikes_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
      static_cast< nest::double_t >( e.get_multiplicity() ) );
  }
}

} // namespace
