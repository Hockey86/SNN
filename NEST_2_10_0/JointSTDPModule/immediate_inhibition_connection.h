/*
 *  immediate_inhibition_connection.h
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

#ifndef IMMEDIATE_INHIBITION_CONNECTION_H
#define IMMEDIATE_INHIBITION_CONNECTION_H

// Includes from nestkernel:
#include "connection.h"

#include "base_neuron.h"


/* BeginDocumentation
  Name: immediate_inhibition_synapse - Static inhibitory synapse with NO delays.

  Description:
   The minimum delay implemented in NEST is 1 time step (usually 1ms). However,
   in computer simulations, for lateral inhibition, 0 delay is needed so that the
   inhibition can happen immediately to prevent simultaneous firing.
   
   The IPSP kernel is an exponential function, not alpha function.
   
   To use this synapse, you must set the 'downstream_neuron_ids' property
   in its presynaptic neuron.

  Transmits: SpikeEvent   

  Author:  Nov 2015, Sun, Haoqi
*/

namespace mynest
{

template< typename targetidentifierT >
class ImmediateInhibitionConnection : public nest::Connection< targetidentifierT >
{
private:
  nest::double_t weight_;

public:
  //! Type to use for representing common synapse properties
  typedef nest::CommonSynapseProperties CommonPropertiesType;

  //! Shortcut for base class
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  ImmediateInhibitionConnection()
    : ConnectionBase()
    , weight_( -1.0 )
  {
  }

  //! Default Destructor.
  ~ImmediateInhibitionConnection()
  {
  }
  
  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  ImmediateInhibitionConnection( const ImmediateInhibitionConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  {
  }
  
  /**
   * Helper class defining which types of events can be transmitted.
   *
   * These methods are only used to test whether a certain type of connection
   * can be created.
   *
   * `handles_test_event()` should be added for all event types that the
   * synapse can transmit. The methods shall return `invalid_port_`; the
   * return value will be ignored.
   *
   * Since this is a synapse model dropping spikes, it is only for spikes,
   * therefore we only implement `handles_test_event()` only for spike
   * events.
   *
   * See Kunkel et al (2014), Sec 3.3.1, for background information.
   */
  class ConnTestDummyNode : public nest::ConnTestDummyNodeBase
  {
  public:
    using nest::ConnTestDummyNodeBase::handles_test_event;
    nest::port
    handles_test_event( nest::SpikeEvent&, nest::rport )
    {
      return nest::invalid_port_;
    }
  };

  /**
   * Check that requested connection can be created.
   *
   * This function is a boilerplate function that should be included unchanged
   * in all synapse models. It is called before a connection is added to check
   * that the connection is legal. It is a wrapper that allows us to call
   * the "real" `check_connection_()` method with the `ConnTestDummyNode
   * dummy_target;` class for this connection type. This avoids a virtual
   * function call for better performance.
   *
   * @param s  Source node for connection
   * @param t  Target node for connection
   * @param receptor_type  Receptor type for connection
   * @param lastspike Time of most recent spike of presynaptic (sender) neuron,
   *                  not used here
   */
  void
  check_connection( nest::Node& s,
    nest::Node& t,
    nest::rport receptor_type,
    nest::double_t,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
  }
  
  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_target;
  //using ConnectionBase::get_rport;
  
  /**
   * Send an event to the receiver of this connection.
   * @param e The event to send
   * @param t Thread
   * @param t_lastspike Point in time of last spike sent.
   * @param cp Common properties to all synapses.
   */
  void send( nest::Event& e,
    nest::thread t,
    nest::double_t t_lastspike,
    const CommonPropertiesType& cp );
    
  // The following methods contain mostly fixed code to forward the
  // corresponding tasks to corresponding methods in the base class and the w_
  // data member holding the weight.

  //! Store connection status information in dictionary
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set connection status.
   *
   * @param d Dictionary with new parameter values
   * @param cm ConnectorModel is passed along to validate new delay values
   */
  void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );

  //! Allows efficient initialization on contstruction
  void
  set_weight( nest::double_t w )
  {
    weight_ = w;
  }
};

template < typename targetidentifierT >
inline void
ImmediateInhibitionConnection< targetidentifierT >::send( nest::Event& e,
  nest::thread t,
  nest::double_t,
  const CommonPropertiesType& )
{
  // get source and target
  base_neuron* mytarget = static_cast< base_neuron* >( get_target( t ) );
  base_neuron* mysource = static_cast< base_neuron* >( &( e.get_sender() ) );
  if ( mytarget==NULL || mysource==NULL )
  {
    std::cout << "Connection 'immediate_inhibition_synapse' can only connect neurons inherited from mynest::base_neuron!" << std::endl;
    throw nest::BadProperty( "Connection 'immediate_inhibition_synapse' can only connect neurons inherited from mynest::base_neuron!" );
  }
  
  // TODO inhibitory plasticity here
  
  // set downstream neuron for immediate inhibition
  std::size_t converted_target_id = mysource->convert_downstream_id( mytarget->get_gid() );
  if ( mysource->get_downstream_immediate_inhibition_neuron( converted_target_id ) == NULL )
    mysource->set_downstream_immediate_inhibition_neuron( converted_target_id, mytarget );
  mysource->set_downstream_weight_and_delay_step( converted_target_id, weight_, get_delay_steps() );
  
  // Here is the trick:
  // immediate_inhibition_connection actually does not send anything,
  // it just tell the source neuron who are the downstream neurons for immediate inhibition, and their weights and delays.
  // The actual update is done in the source neuron using these downstream information,
  // for example, see spike_response_neuron::update.
  
  // you can also implement immediate inhibition by explicitly setting delays to 0 here,
  // but this calls spike_response_neuron::handle and induces a alpha-shaped IPSP,
  // where the inhibition gradually rise with finite time, not immediate as we want.
  /*
  e.set_weight( weight_ );
  e.set_delay( 0 );
  e.set_receiver( *ConnectionBase::get_target( t ) );
  e.set_rport( ConnectionBase::get_rport() );
  e(); // this sends the event
  */
}
  
template < typename targetidentifierT >
void
ImmediateInhibitionConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< nest::double_t >( d, nest::names::weight, weight_ );
  def< nest::long_t >( d, nest::names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
ImmediateInhibitionConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  nest::ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< nest::double_t >( d, nest::names::weight, weight_ );
  if (weight_>=0.0)
    throw nest::BadProperty("Weight must be negative.");
}

} // namespace

#endif /* #ifndef IMMEDIATE_INHIBITION_CONNECTION_H */
