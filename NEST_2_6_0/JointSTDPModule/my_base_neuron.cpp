/*
 *  my_base_neuron.cpp
 *
 *  This file is part of JointSTDPModule.
 *
 */

#include "exceptions.h"
#include "my_base_neuron.h"
//#include "network.h"
#include "dict.h"
//#include "integerdatum.h"
//#include "doubledatum.h"
//#include "dictutils.h"
//#include "numerics.h"

//#include <limits>
namespace mynest
{

my_base_neuron::my_base_neuron()
  : nest::Archiving_Node(),
    baseparam_()
{}

my_base_neuron::BaseParam_::BaseParam_(){}

void my_base_neuron::BaseParam_::get(DictionaryDatum &d) const
{
  //def<std::vector<nest::double_t> >(d, "upstream_delays", upstream_delays_);
  def<std::vector<nest::long_t> >(d, "upstream_neuron_ids", upstream_neuron_ids_);
  def<std::vector<nest::long_t> >(d, "upstream_connection_types", upstream_connection_types_);
  def<std::vector<nest::long_t> >(d, "downstream_neuron_ids", downstream_neuron_ids_);
}

void my_base_neuron::BaseParam_::set(const DictionaryDatum& d)
{
  //updateValue<std::vector<nest::double_t> >(d, "upstream_delays", upstream_delays_);
  updateValue<std::vector<nest::long_t> >(d, "upstream_neuron_ids", upstream_neuron_ids_);
  updateValue<std::vector<nest::long_t> >(d, "upstream_connection_types", upstream_connection_types_);
  updateValue<std::vector<nest::long_t> >(d, "downstream_neuron_ids", downstream_neuron_ids_);

  //if ( upstream_delays_.size() != upstream_neuron_ids_.size() )
  //  throw nest::BadProperty("Inconsistent upstream neuron number and delay number.");
      
  if ( upstream_neuron_ids_.size() > 0 && *std::min_element(upstream_neuron_ids_.begin(), upstream_neuron_ids_.end()) < 0 )
    throw nest::BadProperty("Negative upstream neurons id???");
    
  upstream_neuron_num_ = upstream_neuron_ids_.size();
  
  if (upstream_neuron_num_ != upstream_connection_types_.size())
    throw nest::BadProperty("Inconsistent number of upstream neurons and connection types???");
  
  for (int i=0;i<upstream_neuron_num_;++i)
    if (upstream_connection_types_[i]!=0 && upstream_connection_types_[i]!=1)
      throw nest::BadProperty("Undefined upstream connection type! Use 0 for weight_stdp_synapse, 1 for joint_stdp_synapse.");
      
  upstream_spike_history_ = std::vector<std::list<nest::double_t> >(upstream_neuron_num_,std::list<nest::double_t>(0));
  upstream_connections_ = std::vector<nest::Connection<nest::TargetIdentifierPtrRport>* >(upstream_neuron_num_, NULL);
  downstream_immediate_inhibition_neuron_num_ = downstream_neuron_ids_.size();
  downstream_immediate_inhibition_neurons_ = std::vector< my_base_neuron*>(downstream_immediate_inhibition_neuron_num_, NULL);
  downstream_weights_ = std::vector< nest::double_t>(downstream_immediate_inhibition_neuron_num_, 0.0);
  downstream_delay_steps_ = std::vector< nest::long_t>(downstream_immediate_inhibition_neuron_num_, 0);
  updated_spikes_ = std::list<nest::delay>(0);
}
 
} // namespace
