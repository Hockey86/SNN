/*
 *  my_base_connection.h
 *
 *  This file is part of JointSTDPModule.
 *
 */

#ifndef MY_BASE_CONNECTION_H
#define MY_BASE_CONNECTION_H

/* BeginDocumentation

  This is the base class for all connections in JointSTDPModule,
  inherited from nest::Connection.
  
  Parameters:
   nearest_approximation  bool - whether only consider nearest spike pairs in STDP
   //Wmax        double - Maximum allowed weight, always 1.0
   debug         bool - whether display debug information.
                        For efficiency, debugging is commented out. To recover it,
                        uncomment all blocks starting with if (debug_)...

  Transmits: SpikeEvent

  Author: Nov 2015, Sun, Haoqi

  SeeAlso: weight_stdp_synapse, joint_stdp_synapse
  
*/

#include "connection.h"

#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "event.h"
#include "my_base_neuron.h"


namespace mynest
{

  // connections are templates of target identifier type (used for pointer / target index addressing)
  // derived from generic connection template
  template<typename targetidentifierT>
  class MyBaseConnection : public nest::Connection<targetidentifierT>
  {

  public:

    typedef nest::CommonSynapseProperties CommonPropertiesType;
    typedef nest::Connection<targetidentifierT> ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
    MyBaseConnection();

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
    MyBaseConnection(const MyBaseConnection &);

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, nest::ConnectorModel &cm);

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  virtual void send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const CommonPropertiesType &cp) = 0;


  class ConnTestDummyNode: public nest::ConnTestDummyNodeBase
  {
  public:
	// Ensure proper overriding of overloaded virtual functions.
	// Return values from functions are ignored.
	using nest::ConnTestDummyNodeBase::handles_test_event;
    nest::port handles_test_event(nest::SpikeEvent&, nest::rport) { return nest::invalid_port_; }
  };

  void check_connection(nest::Node & s, nest::Node & t, nest::rport receptor_type, nest::double_t t_lastspike, const CommonPropertiesType &)
  {
    ConnTestDummyNode dummy_target;    
    ConnectionBase::check_connection_(dummy_target, s, t, receptor_type);
    t.register_stdp_connection(t_lastspike - get_delay());
  }

  void set_weight(nest::double_t w) { weight_ = w; }

  nest::double_t get_weight() { return weight_; }
  
  virtual void update_source_upstream(my_base_neuron*, nest::double_t, std::size_t) = 0;
  
  void update_weight(nest::double_t dw) { weight_ += dw;
  if (weight_ >= Wmax_) weight_ = Wmax_;
  else if (weight_ <= 0.0) weight_ = 0.0;
  }
  
  virtual void update_delay(nest::double_t dd) = 0;
 
 protected:
  
  static const nest::double_t DT_PRECISION = 100000.;
  
  // data members of each connection
  nest::double_t weight_;
  nest::double_t Wmax_;  
  bool nearest_approximation_;
  bool debug_;
};

  template<typename targetidentifierT>
  MyBaseConnection<targetidentifierT>::MyBaseConnection() :
    ConnectionBase(),
    weight_(1.0),
    Wmax_(1.0),
    nearest_approximation_(true),
    debug_(false)
  { }

  template<typename targetidentifierT>
  MyBaseConnection<targetidentifierT>::MyBaseConnection(const MyBaseConnection<targetidentifierT> &rhs) :
    ConnectionBase(rhs),
    weight_(rhs.weight_),
    Wmax_(rhs.Wmax_),
    nearest_approximation_(rhs.nearest_approximation_),
    debug_(rhs.debug_)
  {  }

  template<typename targetidentifierT>
  void MyBaseConnection<targetidentifierT>::get_status(DictionaryDatum & d) const
  {
    ConnectionBase::get_status(d);
    def<nest::double_t>(d, nest::names::weight, weight_);
    def<nest::double_t>(d, "Wmax", Wmax_);
    def<bool>(d, "nearest_approximation", nearest_approximation_);
    def<bool>(d, "debug", debug_);
    def<nest::long_t>(d, nest::names::size_of, sizeof(*this));
  }

  template<typename targetidentifierT>
  void MyBaseConnection<targetidentifierT>::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    ConnectionBase::set_status(d, cm);
    updateValue<nest::double_t>(d, nest::names::weight, weight_);
    //updateValue<nest::double_t>(d, "Wmax", Wmax_);
    Wmax_ = 1.0;
    updateValue<bool>(d, "nearest_approximation", nearest_approximation_);
    updateValue<bool>(d, "debug", debug_);
  }

} // of namespace mynest

#endif // my_base_connection.h
