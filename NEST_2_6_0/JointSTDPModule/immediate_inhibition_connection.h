/*
 *  immediate_inhibition_connection.h
 *
 *  This file is part of JointSTDPModule.
 *
 */

#ifndef IMMEDIATE_INHIBITION_CONNECTION_H
#define IMMEDIATE_INHIBITION_CONNECTION_H


/* BeginDocumentation
  Name: immediate_inhibition_synapse - Static inhibitory synapse with NO delays.

  Description:
   It induces an immediate IPSP in its target neuron. However the IPSP kernel
   is different compared with other synapses: this synapse induces exponential
   IPSP (immediate rising), while for other synapses induces alpha-function
   shaped IPSP (finite rising time).
   
   Note: to use this synapse, you must set the '' property in its presynaptic neuron.
   
   This kind of synapse is especially useful to implement lateral inhibition for
   simulators requiring min delay > 0, so that the inhibitory spike can be transmitted
   to neighbor neurons to prevent simultaneous firing.

  Author:  Nov 2015, Sun, Haoqi
*/

#include "connection.h"
#include "my_base_neuron.h"

namespace mynest
{

/**
 * Class representing a static connection. A static connection has the properties weight, delay and receiver port.
 * A suitable Connector containing these connections can be obtained from the template GenericConnector.
 */


template<typename targetidentifierT>
class ImmediateInhibitionConnection : public nest::Connection<targetidentifierT>
{
  nest::double_t weight_;

 public:

  // this line determines which common properties to use
  typedef nest::CommonSynapseProperties CommonPropertiesType;
  typedef nest::Connection<targetidentifierT> ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  ImmediateInhibitionConnection() : ConnectionBase(), weight_(-1.0)
  { }

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  ImmediateInhibitionConnection(const ImmediateInhibitionConnection& rhs) : ConnectionBase(rhs), weight_(rhs.weight_)
  { }

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;


  class ConnTestDummyNode: public nest::ConnTestDummyNodeBase
  {
  public:
	// Ensure proper overriding of overloaded virtual functions.
	// Return values from functions are ignored.
	using nest::ConnTestDummyNodeBase::handles_test_event;
    nest::port handles_test_event(nest::SpikeEvent&, nest::rport) { return nest::invalid_port_; }
    nest::port handles_test_event(nest::RateEvent&, nest::rport) { return nest::invalid_port_; }
    nest::port handles_test_event(nest::DataLoggingRequest&, nest::rport) { return nest::invalid_port_; }
    nest::port handles_test_event(nest::CurrentEvent&, nest::rport) { return nest::invalid_port_; }
    nest::port handles_test_event(nest::ConductanceEvent&, nest::rport) { return nest::invalid_port_; }
    nest::port handles_test_event(nest::DoubleDataEvent&, nest::rport) { return nest::invalid_port_; }
    nest::port handles_test_event(nest::DSSpikeEvent&, nest::rport) { return nest::invalid_port_; }
    nest::port handles_test_event(nest::DSCurrentEvent&, nest::rport) { return nest::invalid_port_; }
  };

  void check_connection(nest::Node & s, nest::Node & t, nest::rport receptor_type, nest::double_t, const CommonPropertiesType &)
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_(dummy_target, s, t, receptor_type);
  }

  void send(nest::Event& e, nest::thread t, nest::double_t, const nest::CommonSynapseProperties &)
  {
    // get source and target
    my_base_neuron* mytarget = static_cast<my_base_neuron*>(get_target(t));
    my_base_neuron* mysource = static_cast<my_base_neuron*>(&(e.get_sender()));
    if (mytarget==NULL || mysource==NULL)
    {
      std::cout<<"Connection 'immediate_inhibition_synapse' can only be used together with neurons inherited from my_base_neuron!"<<std::endl;
      throw nest::BadProperty("Connection 'immediate_inhibition_synapse' can only be used together with neurons inherited from my_base_neuron!");
    }
    
    // if you want to use inhibitory plasticity,
    // you can update weight and delay here.
    
    // set downstream immediate neuron
    std::size_t converted_target_id = mysource->convert_downstream_id(mytarget->get_gid());
    if (mysource->get_downstream_immediate_inhibition_neuron(converted_target_id) == NULL)
      mysource->set_downstream_immediate_inhibition_neuron(converted_target_id, mytarget);
    mysource->set_downstream_weight_and_delay_step(converted_target_id, weight_, get_delay_steps());
    
    // Here is the trick (pretty cool):
    // immediate_inhibition_connection actually does not send anything,
    // it just tell the source neuron its downstream neurons, and their weights and delays.
    // The actual update is done in the source neuron using these downstream information,
    // for example in spike_response_neuron::update.
    /*e.set_weight(weight_);
    // you can also achieve immediate inhibition by explicitly setting delays to 0 here,
    // but this calls spike_response_neuron::handle and induces a alpha-shaped IPSP,
    // where the inhibition gradually rise with finite time, not immediate as we want.
    e.set_delay(0);
    e.set_receiver(*get_target(t));
    e.set_rport(get_rport());
    e();*/
  }

  void get_status(DictionaryDatum & d) const;

  void set_status(const DictionaryDatum & d, nest::ConnectorModel& cm);

  void set_weight (double_t w) { weight_ = w; }
};

template<typename targetidentifierT>
void ImmediateInhibitionConnection<targetidentifierT>::get_status(DictionaryDatum & d) const
{

  ConnectionBase::get_status(d);
  def<nest::double_t>(d, nest::names::weight, weight_);
  def<nest::long_t>(d, nest::names::size_of, sizeof(*this));
}

template<typename targetidentifierT>
void ImmediateInhibitionConnection<targetidentifierT>::set_status(const DictionaryDatum & d, nest::ConnectorModel& cm)
{
  ConnectionBase::set_status(d, cm);
  updateValue<nest::double_t>(d, nest::names::weight, weight_);
  if (weight_>=0.0)
    throw nest::BadProperty("Weight must be negative.");
}

} // namespace

#endif /* #ifndef IMMEDIATE_INHIBITION_CONNECTION_H */
