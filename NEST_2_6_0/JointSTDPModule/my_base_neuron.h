/*
 *  my_base_neuron.h
 *
 *  This file is part of JointSTDPModule.
 *
 */


  /* BeginDocumentation
Name: my_base_neuron.

Description:

  This is the base class for all neurons in JointSTDPModule,
  inherited from nest::Archiving_Node.
  
  This neuron stores its upstream connections, including
  upstream connection instances, upstream neuron ids, upstream delays
  and upstream spike history, as well as downstream neurons with
  immediate inhibition connection.

Author:  Nov 2015, Sun, Haoqi
*/

#ifndef MY_BASE_NEURON_H
#define MY_BASE_NEURON_H

#include "nest.h"
//#include "event.h"
#include "archiving_node.h"
//#include "ring_buffer.h"
#include "connection.h"

#include "nestmodule.h"
#include "target_identifier.h"

namespace mynest
{
  //class nest::Network;

  class my_base_neuron: public nest::Archiving_Node
  {
    
  public:        
    
    my_base_neuron();
    
    bool has_proxies()    const { return false; }
    bool local_receiver() const { return true;  }
    virtual void immediate_inhibit(nest::double_t) = 0;
    virtual void immediate_excite(nest::double_t) = 0;
    
    std::size_t get_downstream_immediate_inhibition_neuron_num()
    { return baseparam_.downstream_immediate_inhibition_neuron_num_; }

    std::vector<my_base_neuron*> get_downstream_immediate_inhibition_neurons() const
    { return baseparam_.downstream_immediate_inhibition_neurons_; }

    my_base_neuron* get_downstream_immediate_inhibition_neuron(nest::size_t id) const
    { return baseparam_.downstream_immediate_inhibition_neurons_[id]; }

    std::vector<nest::Connection<nest::TargetIdentifierPtrRport>* > get_upstream_connections()
    { return baseparam_.upstream_connections_; }

    nest::Connection<nest::TargetIdentifierPtrRport>* get_upstream_connection(nest::size_t id)
    { return baseparam_.upstream_connections_[id]; }

    nest::long_t get_upstream_connection_type(nest::size_t id) const
    { return baseparam_.upstream_connection_types_[id]; }

    nest::long_t get_upstream_neuron_id(nest::size_t id) const
    { return baseparam_.upstream_neuron_ids_[id]; }

    //std::vector<nest::double_t> get_upstream_delays() const
    //{ return baseparam_.upstream_delays_; }

    std::list<nest::double_t> get_upstream_spike_history(std::size_t id) const
    { return baseparam_.upstream_spike_history_[id]; }

    //nest::double_t get_upstream_delay(std::size_t id) const
    //{ return baseparam_.upstream_delays_[id]; }

    std::size_t convert_upstream_id(nest::size_t id)
    {
      std::vector<nest::long_t>::iterator loc = std::find( baseparam_.upstream_neuron_ids_.begin(), baseparam_.upstream_neuron_ids_.end(), id);
      if ( loc != baseparam_.upstream_neuron_ids_.end() )
        return (std::size_t)(loc-baseparam_.upstream_neuron_ids_.begin());
      else
        throw nest::BadProperty("Upstream id not found.");
    }

    std::size_t convert_downstream_id(nest::size_t id)
    {
      std::vector<nest::long_t>::iterator loc = std::find( baseparam_.downstream_neuron_ids_.begin(), baseparam_.downstream_neuron_ids_.end(), id);
      if ( loc != baseparam_.downstream_neuron_ids_.end() )
        return (std::size_t)(loc-baseparam_.downstream_neuron_ids_.begin());
      else
        throw nest::BadProperty("Downstream id not found.");
    }
    
    std::size_t get_upstream_num() { return baseparam_.upstream_neuron_num_; }

    void set_downstream_immediate_inhibition_neuron(nest::size_t id, my_base_neuron* nn)
    { baseparam_.downstream_immediate_inhibition_neurons_[id] = nn; }
    
    void set_downstream_weight_and_delay_step(nest::size_t id, nest::double_t w, nest::long_t ds)
    { baseparam_.downstream_weights_[id] = w; baseparam_.downstream_delay_steps_[id] = ds; }

    void set_upstream_connection(nest::size_t id, nest::Connection<nest::TargetIdentifierPtrRport>* conn)
    { baseparam_.upstream_connections_[id] = conn; }

    void set_upstream_spike_history(nest::size_t id, std::list<nest::double_t> ush)
    { baseparam_.upstream_spike_history_[id] = ush; }

    void append_upstream_spike_history(nest::size_t id, nest::double_t t, nest::double_t lower_bound)
    {      
      // remove old spikes to keep the history span being [lower_bound,t]
      // where lower_bound = t - (stdp_effective_len+max_delay)
      std::list<nest::double_t>::iterator it;
      std::list<nest::double_t>::iterator itb = baseparam_.upstream_spike_history_[id].begin();
      for (it=itb; (it!=baseparam_.upstream_spike_history_[id].end() && *it < lower_bound) ; ++it);
      baseparam_.upstream_spike_history_[id].erase(itb,it);
      
      // append t to upstream spike history
      baseparam_.upstream_spike_history_[id].push_back(t);
    }
    
    bool find_updated_spike(nest::delay step)
    {
      std::list<nest::delay>::iterator it = std::find(baseparam_.updated_spikes_.begin(), baseparam_.updated_spikes_.end(), step);
      return (it!=baseparam_.updated_spikes_.end());
    }

    void mark_updated_spike(nest::delay step, nest::delay lower_bound_step)
    {      
      // remove old spikes to keep the history span being [lower_bound_step,step]
      // where lower_bound_step = step - STEP((stdp_effective_len+max_delay))
      std::list<nest::delay>::iterator it;
      for (it=baseparam_.updated_spikes_.begin(); (it!=baseparam_.updated_spikes_.end() && *it < lower_bound_step) ; ++it);
      baseparam_.updated_spikes_.erase(baseparam_.updated_spikes_.begin(),it);
      
      // append t to upstream spike history
      baseparam_.updated_spikes_.push_back(step);
    }

    //void set_upstream_delay(nest::size_t id, nest::double_t dv)
    //{ baseparam_.upstream_delays_[id] = dv; }
    
    void clear_upstream_spike_history()
    {
      for (nest::size_t i=0; i<baseparam_.upstream_neuron_num_; ++i)
        baseparam_.upstream_spike_history_[i].clear();
    }

  protected:
    
    struct BaseParam_ {
      
      /**downstream neuron number */
      nest::size_t downstream_immediate_inhibition_neuron_num_;
    
      /** Downstream neurons with immediate inhibition */
      std::vector< my_base_neuron*> downstream_immediate_inhibition_neurons_;
      
      /** Downstream weights of the immediate inhibition neurons */
      std::vector<nest::double_t> downstream_weights_;
      
      /** Downstream delays of the immediate inhibition neurons */
      std::vector<nest::long_t> downstream_delay_steps_;

      /** downstream immediate inhibition neuron ids */
      std::vector<nest::long_t> downstream_neuron_ids_;
      
      /**upstream neuron number */
      nest::size_t upstream_neuron_num_;
      
      /**upstream connection types */
      std::vector<nest::long_t> upstream_connection_types_; // 0 for weight_stdp_synapse, 1 for joint_stdp_synapse

      /** upstream connections */
      std::vector<nest::Connection<nest::TargetIdentifierPtrRport>*> upstream_connections_;
      
      /** upstream delays */
      //std::vector<nest::double_t> upstream_delays_;

      /** upstream spike history */
      std::vector<std::list<nest::double_t> > upstream_spike_history_;

      /** upstream neuron ids */
      std::vector<nest::long_t> upstream_neuron_ids_;
      
      /** spikes already used for STDP update by upstream neurons */
      std::list<nest::delay> updated_spikes_;
      
      BaseParam_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary

      /** Set values from dictionary */
      void set(const DictionaryDatum&);
      
    };
    
    BaseParam_ baseparam_;
  };

} // namespace

#endif //MY_BASE_NEURON_H
