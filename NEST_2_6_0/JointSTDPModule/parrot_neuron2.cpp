/*
 *  parrot_neuron2.h
 *
 *  This file is part of JointSTDPModule.
 *
 */

#include "exceptions.h"
#include "parrot_neuron2.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"

#include <limits>
namespace mynest
{

parrot_neuron2::parrot_neuron2()
  : my_base_neuron()
{}

parrot_neuron2::parrot_neuron2(const parrot_neuron2& n)
: my_base_neuron(n)
{}
  
void parrot_neuron2::init_buffers_()
{
  B_.n_spikes_.clear();  // includes resize
  my_base_neuron::clear_history();
  my_base_neuron::clear_upstream_spike_history();
}

void parrot_neuron2::update(nest::Time const & origin, 
			   const nest::long_t from, const nest::long_t to)
{
  assert(to >= 0 && (nest::delay) from < nest::Scheduler::get_min_delay());
  assert(from < to);

  nest::SpikeEvent se;

  for ( nest::long_t lag = from ; lag < to ; ++lag )
  {
    const nest::ulong_t current_spikes_n
      = static_cast<nest::ulong_t>(B_.n_spikes_.get_value(lag));      

    if (current_spikes_n > 0)
    {
      for (nest::ulong_t i_spike=0; i_spike < current_spikes_n; i_spike++)
        network()->send(*this, se, lag);
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      
      for (int i=0; i<baseparam_.downstream_immediate_inhibition_neuron_num_; ++i)
      {
        if (baseparam_.downstream_immediate_inhibition_neurons_[i]!=NULL)
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

void parrot_neuron2::get_status(DictionaryDatum &d) const
{
  def<double>(d, nest::names::t_spike, get_spiketime_ms());
  baseparam_.get(d);
  
}

void parrot_neuron2::set_status(const DictionaryDatum& d)
{

    BaseParam_ bptmp = baseparam_;            // temporary copy in case of errors
    bptmp.set(d);         // throws if BadProperty
    
    my_base_neuron::set_status(d);
  
    baseparam_ = bptmp;
}

void parrot_neuron2::handle(nest::SpikeEvent & e)
{
  B_.n_spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
                         static_cast<nest::double_t>(e.get_multiplicity()));
}
 
} // namespace
