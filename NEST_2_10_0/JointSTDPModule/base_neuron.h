/*
 *  base_neuron.h
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

#ifndef BASE_NEURON_H
#define BASE_NEURON_H

// Includes from libnestutil:
#include "compose.hpp"

// Includes from nestkernel:
#include "nest.h"
#include "archiving_node.h"
#include "connection.h"
#include "nestmodule.h"
#include "target_identifier.h"

namespace mynest
{
/* BeginDocumentation
Name: base_neuron - base class for parrot_neuron2 and spike_response_neuron in JointSTDPModule,
  inherited from nest::Archiving_Node.

Description:
This is an abstract class. Do not instantiate it.
This class of neuron stores its upstream connections,
upstream neuron ids, upstream delays and upstream spike history,
as well as downstream neurons for immediate inhibition connection.

Parameters:
upstream_neuron_ids        std::vector<nest::long_t> - upstream neuron indices
upstream_connection_types  std::vector<nest::long_t> - 0: weight_stdp_synapse, 1: joint_stdp_synapse
downstream_neuron_ids      std::vector<nest::long_t> - downstream neuron indices

Author:  Nov 2015, Sun, Haoqi

SeeAlso: parrot_neuron2, spike_response_neuron
*/

class base_neuron: public nest::Archiving_Node
{
public:        
  /**
   * The constructor is only used to create the model prototype in the model
   * manager.
   */
  base_neuron();

  /**
   * The copy constructor is used to create model copies and instances of the
   * model.
   * @node The copy constructor needs to initialize the parameters and the
   * state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c calibrate().
   */
  base_neuron( const base_neuron& );
    
  bool has_proxies()    const { return false; }
  bool local_receiver() const { return true;  }
  
  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );
    
  std::size_t
  get_downstream_immediate_inhibition_neuron_num()
  {
    return baseparam_.downstream_immediate_inhibition_neuron_num_;
  }

  std::vector< base_neuron* >
  get_downstream_immediate_inhibition_neurons() const
  {
    return baseparam_.downstream_immediate_inhibition_neurons_;
  }

  base_neuron*
  get_downstream_immediate_inhibition_neuron( nest::size_t id ) const
  {
    return baseparam_.downstream_immediate_inhibition_neurons_[id];
  }

  std::vector< nest::Connection<nest::TargetIdentifierPtrRport>* >
  get_upstream_connections()
  {
    return baseparam_.upstream_connections_;
  }

  nest::Connection< nest::TargetIdentifierPtrRport >*
  get_upstream_connection( nest::size_t id )
  {
    return baseparam_.upstream_connections_[id];
  }

  nest::long_t
  get_upstream_connection_type( nest::size_t id ) const
  {
    return baseparam_.upstream_connection_types_[id];
  }

  nest::long_t
  get_upstream_neuron_id( nest::size_t id ) const
  {
    return baseparam_.upstream_neuron_ids_[id];
  }

  /*std::vector<nest::double_t>
  get_upstream_delays() const
  {
    return baseparam_.upstream_delays_;
  }*/

  std::list< nest::double_t >
  get_upstream_spike_history( std::size_t id ) const
  {
    return baseparam_.upstream_spike_history_[id];
  }

  /*nest::double_t
  get_upstream_delay( std::size_t id ) const
  {
    return baseparam_.upstream_delays_[id];
  }*/

  std::size_t
  convert_upstream_id( nest::size_t id )
  {
    std::vector< nest::long_t >::iterator loc = std::find( baseparam_.upstream_neuron_ids_.begin(), baseparam_.upstream_neuron_ids_.end(), id);
    if ( loc != baseparam_.upstream_neuron_ids_.end() )
      return ( std::size_t )( loc - baseparam_.upstream_neuron_ids_.begin() );
    else
      throw nest::BadProperty( "Upstream id not found." );
  }

  std::size_t
  convert_downstream_id( nest::size_t id )
  {
    std::vector< nest::long_t >::iterator loc = std::find( baseparam_.downstream_neuron_ids_.begin(), baseparam_.downstream_neuron_ids_.end(), id);
    if ( loc != baseparam_.downstream_neuron_ids_.end() )
      return ( std::size_t )( loc - baseparam_.downstream_neuron_ids_.begin() );
    else
      throw nest::BadProperty( "Downstream id not found." );
  }
  
  std::size_t
  get_upstream_num()
  {
    return baseparam_.upstream_neuron_num_;
  }

  void
  set_downstream_immediate_inhibition_neuron( nest::size_t id, base_neuron* nn )
  {
    baseparam_.downstream_immediate_inhibition_neurons_[id] = nn;
  }
  
  void
  set_downstream_weight_and_delay_step( nest::size_t id, nest::double_t w, nest::long_t ds )
  {
    baseparam_.downstream_weights_[id] = w; baseparam_.downstream_delay_steps_[id] = ds;
  }

  void
  set_upstream_connection( nest::size_t id, nest::Connection<nest::TargetIdentifierPtrRport>* conn )
  {
    baseparam_.upstream_connections_[id] = conn;
  }

  void
  set_upstream_spike_history( nest::size_t id, std::list<nest::double_t> ush )
  {
    baseparam_.upstream_spike_history_[id] = ush;
  }

  void
  append_upstream_spike_history( nest::size_t id, nest::double_t t, nest::double_t lower_bound )
  {      
    // remove old spikes to keep the history span being [lower_bound,t]
    // where lower_bound = t - (stdp_effective_len+max_delay)
    std::list< nest::double_t >::iterator it;
    std::list< nest::double_t >::iterator itb = baseparam_.upstream_spike_history_[id].begin();
    for ( it = itb; ( it != baseparam_.upstream_spike_history_[id].end() && *it < lower_bound ); ++it );
    baseparam_.upstream_spike_history_[id].erase( itb, it );
    
    // append t to upstream spike history
    baseparam_.upstream_spike_history_[id].push_back( t );
  }
  
  bool
  find_updated_spike( nest::delay step )
  {
    std::list< nest::delay >::iterator it = std::find( baseparam_.updated_spikes_.begin(), baseparam_.updated_spikes_.end(), step );
    return ( it != baseparam_.updated_spikes_.end() );
  }

  void
  mark_updated_spike( nest::delay step, nest::delay lower_bound_step )
  {      
    // remove old spikes to keep the history span being [lower_bound_step,step]
    // where lower_bound_step = step - STEP((stdp_effective_len+max_delay))
    std::list< nest::delay >::iterator it;
    for ( it = baseparam_.updated_spikes_.begin(); ( it != baseparam_.updated_spikes_.end() && *it < lower_bound_step ); ++it );
    baseparam_.updated_spikes_.erase( baseparam_.updated_spikes_.begin(), it );
    
    // append t to upstream spike history
    baseparam_.updated_spikes_.push_back( step );
  }

  /*void
  set_upstream_delay( nest::size_t id, nest::double_t dv )
  {
    baseparam_.upstream_delays_[id] = dv;
  }*/
  
  void
  clear_upstream_spike_history()
  {
    for ( nest::size_t i = 0; i < baseparam_.upstream_neuron_num_; ++i )
      baseparam_.upstream_spike_history_[i].clear();
  }
  
  virtual void immediate_inhibit( nest::double_t ) = 0;
  virtual void immediate_excite( nest::double_t ) = 0;

protected:
  
  struct BaseParam_
  {
    nest::size_t downstream_immediate_inhibition_neuron_num_;
    std::vector< base_neuron*> downstream_immediate_inhibition_neurons_;
    std::vector<nest::double_t> downstream_weights_;
    std::vector<nest::long_t> downstream_delay_steps_;
    std::vector<nest::long_t> downstream_neuron_ids_;
    nest::size_t upstream_neuron_num_;
    std::vector<nest::long_t> upstream_connection_types_;
    std::vector<nest::Connection<nest::TargetIdentifierPtrRport>*> upstream_connections_;
    //std::vector<nest::double_t> upstream_delays_;
    std::vector<std::list<nest::double_t> > upstream_spike_history_;
    std::vector<nest::long_t> upstream_neuron_ids_;
    std::list<nest::delay> updated_spikes_;
        
    //! Initialize parameters to their default values.
    BaseParam_();

    //! Store parameter values in dictionary.
    void get( DictionaryDatum& ) const;

    //! Set parameter values from dictionary.
    void set( const DictionaryDatum& );
  };
  
  BaseParam_ baseparam_;
};

inline void
base_neuron::get_status( DictionaryDatum& d ) const
{
  // get our own parameter and state data
  baseparam_.get( d );

  // get information managed by parent class
  nest::Archiving_Node::get_status( d );
}

inline void
base_neuron::set_status( const DictionaryDatum& d )
{
  BaseParam_ bptmp = baseparam_;  // temporary copy in case of errors
  bptmp.set(d);                   // throws if BadProperty
  
  nest::Archiving_Node::set_status(d);

  baseparam_ = bptmp;
}

} // namespace

#endif /* #ifndef BASE_NEURON_H */
