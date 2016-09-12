/*
 *  base_connection.h
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
 
#ifndef BASE_CONNECTION_H
#define BASE_CONNECTION_H

// Includes from nestkernel:
#include "connection.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

#include "base_neuron.h"



/* BeginDocumentation
Name: base - base class for JointSTDPConnection and WeightSTDPConnection in JointSTDPModule,
inherited from nest::Connection.

Description:
This is an abstract class. Do not instantiate it.
Descriptions are available in its inherited classes.

Parameters:
nearest_approximation  bool - whether only consider nearest spike pairs in STDP
// Wmax        double - Maximum allowed weight, always 1.0
debug          bool - whether display debug information.
                    For efficiency, debugging is commented out. To recover it,
                    uncomment all blocks starting with "if (debug_)"

Author: Nov 2015, Sun, Haoqi

SeeAlso: weight_stdp_synapse, joint_stdp_synapse
*/

namespace mynest
{

template < typename targetidentifierT >
class BaseConnection : public nest::Connection< targetidentifierT >
{
protected:
  // used to round delay values,
  // for example, DT_PRECISION = 100000. means the delays are rounded to 5 digits after decimal point
  // keep it constant, 100000. is enough
  static const nest::double_t DT_PRECISION = 100000.;
  
  nest::double_t weight_; //!< Synaptic weight
  nest::double_t Wmax_;  
  bool nearest_approximation_;
  bool debug_;

public:
  //! Type to use for representing common synapse properties
  typedef nest::CommonSynapseProperties CommonPropertiesType;

  //! Shortcut for base class
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  BaseConnection()
    : ConnectionBase()
    , weight_( 1.0 )
    , Wmax_( 1.0 )
    , nearest_approximation_( true )
    , debug_( false )
  {
  }

  //! Default Destructor.
  ~BaseConnection()
  {
  }
  
  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  BaseConnection( const BaseConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , Wmax_( rhs.Wmax_ )
  , nearest_approximation_( rhs.nearest_approximation_ )
  , debug_( rhs.debug_ )
  {
  }

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  //using ConnectionBase::get_delay;

  /**
   * Send an event to the receiver of this connection.
   * @param e The event to send
   * @param t Thread
   * @param t_lastspike Point in time of last spike sent.
   * @param cp Common properties to all synapses.
   */
  virtual void send( nest::Event& e,
    nest::thread t,
    nest::double_t t_lastspike,
    const CommonPropertiesType& cp ) = 0;
  
  virtual void update_source_upstream( base_neuron*,
    nest::double_t,
    std::size_t ) = 0;
    
  virtual void update_delay( nest::double_t dd ) = 0;

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
  
  nest::double_t
  get_weight()
  {
    return weight_;
  }
  
  void
  update_weight( nest::double_t dw )
  {
    weight_ += dw;
    if ( weight_ >= Wmax_ ) weight_ = Wmax_;
    else if ( weight_ <= 0.0 ) weight_ = 0.0;
  } 
};

template < typename targetidentifierT >
void
BaseConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< nest::double_t >( d, nest::names::weight, weight_ );
  def< nest::double_t >( d, "Wmax", Wmax_ );
  def< bool >( d, "nearest_approximation", nearest_approximation_ );
  def< bool >( d, "debug", debug_ );
  def< nest::long_t >( d, nest::names::size_of, sizeof( *this ) );
}
  
template < typename targetidentifierT >
void
BaseConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  nest::ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< nest::double_t >( d, nest::names::weight, weight_ );
  //updateValue<nest::double_t>(d, "Wmax", Wmax_);
  Wmax_ = 1.0;
  updateValue<bool>(d, "nearest_approximation", nearest_approximation_);
  updateValue<bool>(d, "debug", debug_);
}

} // namespace

#endif // base_connection.h
