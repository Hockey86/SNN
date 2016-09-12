/*
 *  base_neuron.cpp
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
 *  Sun Haoqi, Jun 28, 2016, NTU, Singapore
 */

#include "base_neuron.h"

// Includes from nestkernel:
#include "exceptions.h"

// Includes from sli:
#include "dict.h"

using namespace nest;


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

mynest::base_neuron::base_neuron()
  : Archiving_Node()
  , baseparam_()
{
}

mynest::base_neuron::base_neuron( const mynest::base_neuron& n )
  : Archiving_Node( n )
  , baseparam_( n.baseparam_ )
{
}


/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

mynest::base_neuron::BaseParam_::BaseParam_()
{
}


/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */
 
void mynest::base_neuron::BaseParam_::get(DictionaryDatum &d) const
{
  //( *d )[ "upstream_delays" ] = upstream_delays_;
  ( *d )[ "upstream_neuron_ids" ] = upstream_neuron_ids_;
  ( *d )[ "upstream_connection_types" ] = upstream_connection_types_;
  ( *d )[ "downstream_neuron_ids" ] = downstream_neuron_ids_;
}

void mynest::base_neuron::BaseParam_::set(const DictionaryDatum& d)
{
  //updateValue< std::vector< double_t > >( d, "upstream_delays", upstream_delays_ );
  updateValue< std::vector< long_t > >( d, "upstream_neuron_ids", upstream_neuron_ids_ );
  updateValue< std::vector< long_t > >( d, "upstream_connection_types", upstream_connection_types_ );
  updateValue< std::vector< long_t > >( d, "downstream_neuron_ids", downstream_neuron_ids_ );

  //if ( upstream_delays_.size() != upstream_neuron_ids_.size() )
  //  throw BadProperty( "Inconsistent upstream neuron number and delay number." );
      
  if ( upstream_neuron_ids_.size() > 0 && *std::min_element(upstream_neuron_ids_.begin(), upstream_neuron_ids_.end()) < 0 )
    throw BadProperty( "Negative upstream neurons id???" );
    
  upstream_neuron_num_ = upstream_neuron_ids_.size();
  
  if ( upstream_neuron_num_ != upstream_connection_types_.size() )
    throw BadProperty( "Inconsistent number of upstream neurons and connection types???" );
  
  for ( size_t i = 0; i < upstream_neuron_num_; ++i )
    if ( upstream_connection_types_[i] != 0 && upstream_connection_types_[i] != 1 )
      throw BadProperty( "Undefined upstream connection type! Use 0 for weight_stdp_synapse, 1 for joint_stdp_synapse." );
      
  upstream_spike_history_ = std::vector< std::list< double_t > >( upstream_neuron_num_, std::list< double_t >( 0 ) );
  upstream_connections_ = std::vector< Connection< TargetIdentifierPtrRport >* >( upstream_neuron_num_, NULL );
  downstream_immediate_inhibition_neuron_num_ = downstream_neuron_ids_.size();
  downstream_immediate_inhibition_neurons_ = std::vector< base_neuron* >( downstream_immediate_inhibition_neuron_num_, NULL );
  downstream_weights_ = std::vector< double_t >( downstream_immediate_inhibition_neuron_num_, 0.0 );
  downstream_delay_steps_ = std::vector< long_t >( downstream_immediate_inhibition_neuron_num_, 0 );
  updated_spikes_ = std::list< delay >( 0 );
}
