/*
 *  joint_stdp_connection.h
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
 
#ifndef JOINT_STDP_CONNECTION_H
#define JOINT_STDP_CONNECTION_H

// C++ includes:
#include <cmath>
#include <limits>

// Includes from nestkernel:
//#include "connection.h" // included in base_connection.h
#include "nest_time.h"

#include "base_connection.h"


/* BeginDocumentation
  Name: joint_stdp_synapse - Synapse with spike-timing dependent plasticity
  which updates weights and delays jointly.
  
  Description:
    This synpase is similar with stdp_synapse in NEST. The differences are
    1) updates weights and delays jointly
    2) axonal delays (dendritic delay in stdp_synapse);
    3) weight and delay updates happen at postsynaptic site (presynaptic in stdp_synapse);
    4) homeostatic nonhebbian terms;
    5) addictive only.

  Parameters:   
   tau_plus      double - Time constant of weight STDP window (positive amp. side), in ms, > 0
   tau_minus     double - Time constant of weight STDP window (negative amp. side), in ms, > 0
   A_plus        double - Amplitude of of weight STDP window (positive amp. side), dimensionless, >= 0
   A_minus       double - Amplitude of of weight STDP window (negative amp. side), dimensionless, >= 0
   w_in          double - Amount of weight change when receiving a presynaptic spike, dimensionless, >=0
   w_out         double - Amount of weight change when firing a spike, dimensionless, <=0
   STDP_weight_effective_len  double - Effective weight STDP window length L_plus = STDP_weight_effective_len*tau_plus, > 0
   
   sigma_minus   double - Time constant of delay STDP window (positive amp. side), in ms, > 0
   sigma_plus    double - Time constant of delay STDP window (negative amp. side), in ms, > 0
   B_minus       double - Amplitude of of delay STDP window (positive amp. side), in ms, >= 0
   B_plus        double - Amplitude of of delay STDP window (negative amp. side), in ms, >= 0
   delta_minus   double - weight-related gain offset (negative amp. side), dimensionless, >= 0
   delta_plus    double - weight-related gain offset (positive amp. side), dimensionless, <= 0
   delay_min     double - minimum delay, in ms, > 0
   delay_max     double - maximum delay, in ms, >= delay_min

  Transmits: SpikeEvent

  Author: Moritz Helias, Abigail Morrison
  
  Modified: Nov 2015, Sun, Haoqi

  SeeAlso: weight_stdp_synapse, base_connection
*/

namespace mynest
{
template<typename targetidentifierT>
class WeightSTDPConnection;

// connections are templates of target identifier type (used for pointer / target index addressing)
// derived from generic connection template
template<typename targetidentifierT>
class JointSTDPConnection : public BaseConnection<targetidentifierT>
{

public:
  //! Type to use for representing common synapse properties
  typedef nest::CommonSynapseProperties CommonPropertiesType;

  //! Shortcut for base class
  typedef BaseConnection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  JointSTDPConnection()
  : ConnectionBase()
  , internal_delay_( -1.0 )
  , tau_plus_( 20.0 )
  , tau_minus_( 20.0 )
  , A_plus_( 0.01 )
  , A_minus_( 0.02 )
  //, mu_plus_( 1.0 )
  //, mu_minus_( 1.0 )
  //, Kplus_( 0.0 )
  , w_in_( 0.0 )
  , w_out_( 0.0 )
  , STDP_weight_effective_len_( 5.0 )
  //, eta_d_( 0.001 )
  , delay_max_( 40.0 )
  , delay_min_( 1.0 )
  , sigma_minus_( 20.0 )
  , sigma_plus_( 20.0 )
  , B_minus_( 0.01 )
  , B_plus_( 0.01 )
  , delta_minus_( 0.0 )
  , delta_plus_( 0.0 )
  //, STDP_delay_effective_len_( 5.0 )
  {
  }

  //! Default Destructor.
  ~JointSTDPConnection()
  {
  }
  
  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  JointSTDPConnection( const JointSTDPConnection< targetidentifierT >& rhs )
    : ConnectionBase( rhs )
    , internal_delay_( rhs.internal_delay_ )
    , tau_plus_( rhs.tau_plus_ )
    , tau_minus_( rhs.tau_minus_ )
    , A_plus_( rhs.A_plus_ )
    , A_minus_( rhs.A_minus_ )
    //, mu_plus_( rhs.mu_plus_ )
    //, mu_minus_( rhs.mu_minus_ )
    //, Kplus_( rhs.Kplus_ )
    , w_in_( rhs.w_in_ )
    , w_out_( rhs.w_out_ )
    , STDP_weight_effective_len_( rhs.STDP_weight_effective_len_ )
    //, eta_d_( rhs.eta_d_ )
    , delay_max_( rhs.delay_max_ )
    , delay_min_( rhs.delay_min_ )
    , sigma_minus_( rhs.sigma_minus_ )
    , sigma_plus_( rhs.sigma_plus_ )
    , B_minus_( rhs.B_minus_ )
    , B_plus_( rhs.B_plus_ )
    , delta_minus_( rhs.delta_minus_ )
    , delta_plus_( rhs.delta_plus_ )
    //, STDP_delay_effective_len_( rhs.STDP_delay_effective_len_ )
  {
  }

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;
  using ConnectionBase::set_delay_steps;
  using ConnectionBase::update_weight;
  using ConnectionBase::DT_PRECISION;
  using ConnectionBase::weight_;
  using ConnectionBase::Wmax_;  
  using ConnectionBase::nearest_approximation_;
  using ConnectionBase::debug_;

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
    nest::double_t t_lastspike,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
    t.register_stdp_connection( t_lastspike - get_delay() );
  }

  //! Store connection status information in dictionary
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set connection status.
   *
   * @param d Dictionary with new parameter values
   * @param cm ConnectorModel is passed along to validate new delay values
   */
  void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );

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

  void update_source_upstream( base_neuron*,
    nest::double_t,
    std::size_t );

  void update_delay( nest::double_t dd)
  {
    if ( internal_delay_ <= 0.0 ) {
      //std::cout << "internal_delay_ not set, set to " << get_delay() << std::endl;
      internal_delay_ = get_delay();
    }
    internal_delay_ += dd;
    if ( internal_delay_ > delay_max_ ) internal_delay_ = delay_max_;
    else if ( internal_delay_ < delay_min_ ) internal_delay_ = delay_min_;
    set_delay_steps( nest::Time::delay_ms_to_steps( internal_delay_ ) );
  }

protected:
  // There is variable 'delay' in ConnectionBase, but it is rounded to time_step
  // we need to store the precise delay to avoid rounding error
  nest::double_t internal_delay_;

private:
  /*
  double_t
  facilitate_( double_t w, double_t kplus )
  {
    double_t norm_w = ( w / Wmax_ )
      + ( lambda_ * std::pow( 1.0 - ( w / Wmax_ ), mu_plus_ ) * kplus );
    return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
  }

  double_t
  depress_( double_t w, double_t kminus )
  {
    double_t norm_w = ( w / Wmax_ )
      - ( alpha_ * lambda_ * std::pow( w / Wmax_, mu_minus_ ) * kminus );
    return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
  }*/

  // data members of each connection
  nest::double_t tau_plus_;
  nest::double_t tau_minus_;
  nest::double_t A_plus_;
  nest::double_t A_minus_;
  //nest::double_t mu_plus_;
  //nest::double_t mu_minus_;  
  //nest::double_t Kplus_;
  nest::double_t w_in_;
  nest::double_t w_out_;
  nest::double_t STDP_weight_effective_len_;

  //nest::double_t eta_d_;
  nest::double_t delay_max_;
  nest::double_t delay_min_;
  nest::double_t sigma_minus_;
  nest::double_t sigma_plus_;
  nest::double_t B_minus_;
  nest::double_t B_plus_;
  nest::double_t delta_minus_;
  nest::double_t delta_plus_;
  //nest::double_t STDP_delay_effective_len_;
};


template < typename targetidentifierT >
inline void
JointSTDPConnection< targetidentifierT >::send( nest::Event& e,
  nest::thread t,
  nest::double_t t_lastspike,
  const CommonPropertiesType&)
{
  /* downstream stuff */
  
  // get source and target
  base_neuron* mytarget = static_cast< base_neuron* >( get_target( t ) );
  base_neuron* mysource = static_cast< base_neuron* >( &( e.get_sender() ) );
  if ( mytarget==NULL || mysource==NULL )
  {
    std::cout << std::endl << std::endl << "Connection 'joint_stdp_synapse' can only connect neurons inherited from base_neuron!" << std::endl << std::endl;
    throw nest::BadProperty( "Connection 'joint_stdp_synapse' can only connect neurons inherited from base_neuron!" );
  }
  
  std::deque< nest::histentry >::iterator start;
  std::deque< nest::histentry >::iterator finish;
  nest::double_t t_spike = e.get_stamp().get_ms(); // t_lastspike_ = 0 initially
  nest::double_t dt, tmp1, tmp2;
  if ( internal_delay_ <= 0.0 )
    internal_delay_ = get_delay();
  //nest::double_t delay_max_ = nest::Time::delay_steps_to_ms( e.get_max_delay() );
  //nest::double_t original_dendritic_delay = dendritic_delay;

  std::size_t converted_sender_id = mytarget->convert_upstream_id( mysource->get_gid() ), cc;
  /*if ( debug_ ) {
    std::cout << "downstream" << std::endl;
    std::cout << "from neuron #" << mysource->get_gid() << " to #" << mytarget->get_gid() << " original weight: " << weight_ << " original delay: " << internal_delay_ << std::endl;
  }*/
  
  // tell downstream neurons my spike history (I'm the upstream of downstream neurons)
  if ( mytarget->get_upstream_connection( converted_sender_id ) == NULL ) {
    //if ( debug_ ) std::cout << "target #" << mytarget->get_gid() << "does not have upstream connection, create." << std::endl;
    mytarget->set_upstream_connection( converted_sender_id, this );
  }
  mytarget->append_upstream_spike_history( converted_sender_id, t_spike, t_spike-tau_plus_*STDP_weight_effective_len_-delay_max_ );
  /*if ( debug_ ) {
    std::list< nest::double_t > ttt = mytarget->get_upstream_spike_history( converted_sender_id );
    for ( std::list<nest::double_t>::iterator it = ttt.begin(); it != ttt.end(); ++it )
      std::cout << *it << std::endl;
  }*/
	        
  // update downstream weights
  // use mytarget->get_history to get downstream neuron spike history
  tmp1 = B_plus_*(weight_-delta_plus_)/Wmax_;
  tmp2 = tmp1*sigma_plus_/std::exp(1.0);
  bool delay_not_updated = true; // delay plasticity always uses nearest spike scheme, i.e. update once
  mytarget->get_history( std::min( t_spike+internal_delay_-STDP_weight_effective_len_*tau_minus_, t_spike+internal_delay_-sigma_plus_ ), t_spike, &start, &finish );
  while ( start != finish )
  {
    --finish;
    dt = t_spike - finish->t_ + internal_delay_;
    dt = dt < 0 ? std::ceil( dt*DT_PRECISION-0.5 )/DT_PRECISION : std::floor( dt*DT_PRECISION+0.5 )/DT_PRECISION;
    //if ( debug_ ) std::cout << t_spike << "-" << finish->t_ << "+" << internal_delay_ << "=" << dt << std::endl;
    
    if ( dt <= STDP_weight_effective_len_*tau_minus_ )  // dt must be within (delay, STDP_weight_effective_len_*tau_minus_]
      weight_ -= ( A_minus_ * std::exp( -dt/tau_minus_ ) );
    //else  // dt cannot be within [-STDP_weight_effective_len_*tau_plus_, 0], so comment out
    //  weight_ += ( A_plus_ * std::exp( dt / tau_plus_ ) );
    
    if ( delay_not_updated && dt <= sigma_plus_ ) // dt must be within (delay, STDP_weight_effective_len_*sigma_plus_]
    {
      internal_delay_ -= ( tmp1 * dt * std::exp( -dt / sigma_plus_ ) - tmp2 );
      delay_not_updated = false; // delay plasticity always uses nearest spike scheme, i.e. update once
    }
    //else  // dt cannot be within [-STDP_weight_effective_len_*sigma_minus_, 0], so comment out
    //  internal_delay_ += ??;
    
    if ( nearest_approximation_ ) break;
  }
  
  //if ( debug_ ) std::cout << "after downstream dW: " << weight_ << std::endl;
  // first order non-Hebbian term on weight, w_in
  weight_ += w_in_;
  //if ( debug_ ) std::cout << "after downstream w_in: " << weight_ << std::endl;
  // cut the weight
  if ( weight_ > Wmax_ ) weight_ = Wmax_;
  else if ( weight_ < 0.0 ) weight_ = 0.0;
  //if ( debug_ ) std::cout << "after cut: " << weight_ << std::endl;
  
  //if ( debug_ ) std::cout << "after downstream dD: " << internal_delay_ << std::endl;
  // cut the delay
  if ( internal_delay_ > delay_max_ ) internal_delay_ = delay_max_;
  else if ( internal_delay_ < delay_min_ ) internal_delay_ = delay_min_;
  //if (debug_) std::cout << "after cut: " << internal_delay_ << std::endl;
  
  
  /* upstream stuff */
  
  std::size_t source_upstream_num = mysource->get_upstream_num();
  if ( source_upstream_num > 0 )
  {
    //if ( debug_ ) std::cout << "upstream" << std::endl;
    
    // mark this spike is updated, so that different downstream connctions will not update its upstream multiple times
    nest::delay spike_step = e.get_stamp().get_steps();
    if ( !mysource->find_updated_spike( spike_step ) )  // not updated yet
    {        
      mysource->mark_updated_spike( spike_step, spike_step-3 );  // mark it as updated, keep a short memory of 3 steps
      
      ConnectionBase* myconn;
      
      /**** the following part was commented out because
      we decided to use local delay plasticity rule ****/
      
      /*tmp8 = 0.; // max(last spike, current spike-affective_len)
      mysource->get_history(t_spike-STDP_delay_effective_len_*sigma_minus_, t_spike-0.01, &start, &finish);
      if (start != finish)
        tmp8 = (--finish)->t_;
      tmp8 = std::max(tmp8,t_spike-sigma_minus_*STDP_delay_effective_len_);
        
      // first, downcast upstream connection vector to type std::vector<JointSTDPConnection<targetidentifierT>*>
      // and get all presynaptic spike arriving times between max(last spike, current spike-affective_len) and current spike
      std::vector<JointSTDPConnection<targetidentifierT>*> myconns;
      std::vector<nest::double_t> presynaptic_arriving_times(source_upstream_num,0.0);//zero means no spike between..., "meaningless"
      std::vector<nest::double_t> dts(source_upstream_num,0.0);
      //tmp2 = 0;
      if (debug_) std::cout<<"upstream spike arrival between "<<tmp8<<" and "<<t_spike<<std::endl;
      for (cc = 0; cc < source_upstream_num; ++cc)
      {
        myconn = static_cast<JointSTDPConnection<targetidentifierT>*>(mysource->get_upstream_connection(cc));
        myconns.push_back(myconn);
        if (myconn!=NULL)
        {
          upstream_spike_history = mysource->get_upstream_spike_history(cc);
          for (std::list<nest::double_t>::reverse_iterator rit = upstream_spike_history.rbegin() ; rit != upstream_spike_history.rend(); ++rit)
          {
            tmp1 = *rit + myconn->internal_delay_;
            dt = tmp1 - t_spike; // presynaptic spike arriving time, t_k^f+d_ik
            dt = dt < 0 ? std::ceil(dt*DT_PRECISION-0.5)/DT_PRECISION : std::floor(dt*DT_PRECISION+0.5)/DT_PRECISION;
            //if (tmp1 < t_spike && *rit >= tmp8) // if both presynaptic spike and its arriving time are between last and current postsynaptic spike
            if (dt < 0 && dt >= -sigma_minus_)
            {
              presynaptic_arriving_times[cc] = tmp1;
              dts[cc] = dt;
              //tmp2 += 1; // tmp2 is the number of "meaningful" values
              if (debug_) std::cout<<"k="<<cc<<": t_k^f+d_ik="<<*rit<<"+"<<myconn->internal_delay_<<"="<<tmp1<<std::endl;
              break;
            }
            else if (tmp1 < tmp8)
            {
              if (debug_) std::cout << "neuron #" << mysource->get_gid() << " has no upstream spike history." << std::endl;
              break;  // correct only if spike history is sorted ascendly
            }
          }          
        }
        else
          if (debug_)
            std::cout << "neuron #" << mysource->get_gid() << " upstream neuron #" << cc << " connection unknown." << std::endl;
      }
      
      // max of all "meaningful" values
      tmp2 = *std::max_element(presynaptic_arriving_times.begin(),presynaptic_arriving_times.end()); // max presynaptic arriving time
      //tmp2 = std::accumulate(presynaptic_arriving_times.begin(), presynaptic_arriving_times.end(), 0.0)/tmp2;
      if (debug_) std::cout<<"max upstream spike arrival: "<<tmp2<<std::endl;*/
      
      // make upstream changes on delay and weight
      for ( cc = 0; cc < source_upstream_num; ++cc )
      {
        if ( mysource->get_upstream_connection_type( cc ) == 0 )  // weight_stdp_synapse
          myconn = static_cast< WeightSTDPConnection< targetidentifierT >* >( mysource->get_upstream_connection( cc ) );
        else  // joint_stdp_synapse
          myconn = static_cast< JointSTDPConnection< targetidentifierT >* >( mysource->get_upstream_connection( cc ) );
        if ( myconn!=NULL )
          myconn->update_source_upstream( mysource, t_spike, cc );
      }
    }
  }

  // send event
  e.set_receiver( *mytarget );
  e.set_weight( weight_ );
  e.set_delay( nest::Time::delay_ms_to_steps( internal_delay_ ) );
  e.set_rport( get_rport() );
  e(); // this sends the event
  //if (debug_) std::cout << std::endl;
}

template < typename targetidentifierT >
inline void
JointSTDPConnection< targetidentifierT >::update_source_upstream( base_neuron* source,
  nest::double_t t_spike,
  std::size_t cc)
{
  std::list< nest::double_t > upstream_spike_history = source->get_upstream_spike_history( cc );
  nest::double_t dt = std::exp( 1.0 ), tmp1 = w_out_;  // weight change
  nest::double_t tmp2 = 0.0;  // delay change
  nest::double_t tmp3 = sigma_minus_/dt, tmp4 = sigma_plus_/dt;
  nest::double_t tmp5 = (Wmax_-weight_+delta_minus_)/Wmax_, tmp6=(weight_-delta_plus_)/Wmax_;
  nest::double_t tmp7 = tau_minus_*STDP_weight_effective_len_, tmp8 = -tau_plus_*STDP_weight_effective_len_;
  bool delay_not_updated = true; // delay always update once for each postsynaptic spike, regardless of STDP all-to-all relation
  
  for ( std::list< nest::double_t >::reverse_iterator rit = upstream_spike_history.rbegin(); rit != upstream_spike_history.rend(); ++rit )
  {
    dt = *rit + internal_delay_ - t_spike;
    dt = dt < 0 ? std::ceil( dt*DT_PRECISION-0.5 )/DT_PRECISION : std::floor( dt*DT_PRECISION+0.5 )/DT_PRECISION;
    
    // if outside STDP effective length, break
    // correct only if spike history is sorted ascendly
    //if ( dt < tmp8 ) break;
    // spikes outside STDP effective length is already removed in mytarget->append_upstream_spike_history
    
    // update weight
    if ( dt >= 0.0 && dt <= tmp7 )
      tmp1 -= ( A_minus_ * std::exp( -dt / tau_minus_ ) );
    else if ( dt < 0.0 && dt >= tmp8 )
      tmp1 += ( A_plus_ * std::exp( dt / tau_plus_ ) );
      
    // update delay
    if ( delay_not_updated )
    {         
      if ( dt >= 0.0 && dt <= sigma_plus_ )
      {
        tmp2 -= (B_plus_*tmp6*(dt*std::exp(-dt/sigma_plus_)-tmp4));
        delay_not_updated = false; // delay plasticity always uses nearest spike scheme, i.e. update once
      }
      else if ( dt < 0.0 && dt >= -sigma_minus_ )
      {
        tmp2 -= (B_minus_*tmp5*(dt*std::exp(dt/sigma_minus_)+tmp3));
        delay_not_updated = false; // delay plasticity always uses nearest spike scheme, i.e. update once
      }
    }
    //if ( debug_ )
    //  std::cout << cc << ": dt=" << *rit << "-" << t_spike << "+" << internal_delay_ << "=" << dt << ", dW=" << tmp1 << ", dd=" << tmp2 << std::endl;
    
    if ( nearest_approximation_ ) break;
  }
    
  /**** the following part was commented out because
  we decided to use local delay plasticity rule ****/
  
  /*// update delay for dt<0
  // delay plasticity always uses nearest spike scheme, i.e. update once
  if (tmp6>0.0 && presynaptic_arriving_times[cc]>0.0)// && delay_once_flag (dts[cc] < 0.0 always true)
  {
    //tmp2 += (tmp6-presynaptic_arriving_times[cc])*B_minus_;
    //if (debug_) std::cout << cc << ": dd=(" << tmp6 << "-" << presynaptic_arriving_times[cc] << ")*" << B_minus_ << "=" << (tmp6-presynaptic_arriving_times[cc])*B_minus_ <<std::endl;
    tmp2 += (tmp6-presynaptic_arriving_times[cc])*B_minus_*std::exp(dts[cc]/sigma_minus_);
    if (debug_) std::cout << cc << ": dd=(" << tmp6 << "-" << presynaptic_arriving_times[cc] << ")*" << B_minus_ << "*(exp("<<dts[cc]<<"/"<<sigma_minus_<<")="<<(tmp6-presynaptic_arriving_times[cc])*B_minus_*std::exp(dts[cc]/sigma_minus_)<<std::endl;
  }*/
  update_weight( tmp1 );
  update_delay( tmp2 );
}

template < typename targetidentifierT >
void
JointSTDPConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< nest::double_t >( d, "tau_plus", tau_plus_ );
  def< nest::double_t >( d, "tau_minus", tau_minus_ );
  def< nest::double_t >( d, "A_plus", A_plus_ );
  def< nest::double_t >( d, "A_minus", A_minus_ );
  //def< nest::double_t >( d, "mu_plus", mu_plus_ );
  //def< nest::double_t >( d, "mu_minus", mu_minus_ );
  def< nest::double_t >( d, "w_in", w_in_ );
  def< nest::double_t >( d, "w_out", w_out_ );
  def< nest::double_t >( d, "STDP_weight_effective_len", STDP_weight_effective_len_ );
  //def< nest::double_t >( d, "eta_d", eta_d_ );
  def< nest::double_t >( d, "delay_max", delay_max_ );
  def< nest::double_t >( d, "delay_min", delay_min_ );
  def< nest::double_t >( d, "sigma_minus", sigma_minus_ );
  def< nest::double_t >( d, "sigma_plus", sigma_plus_ );
  def< nest::double_t >( d, "B_minus", B_minus_ );
  def< nest::double_t >( d, "B_plus", B_plus_ );
  def< nest::double_t >( d, "delta_minus", delta_minus_ );
  def< nest::double_t >( d, "delta_plus", delta_plus_ );
  //def< nest::double_t >( d, "STDP_delay_effective_len", STDP_delay_effective_len_ );
  def< nest::long_t >( d, nest::names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
JointSTDPConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  nest::ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< nest::double_t >( d, "tau_plus", tau_plus_ );
  updateValue< nest::double_t >( d, "tau_minus", tau_minus_ );
  updateValue< nest::double_t >( d, "A_plus", A_plus_ );
  updateValue< nest::double_t >( d, "A_minus", A_minus_ );
  //updateValue< nest::double_t >( d, "mu_plus", mu_plus_ );
  //updateValue< nest::double_t >( d, "mu_minus", mu_minus_ );
  updateValue< nest::double_t >( d, "w_in", w_in_ );
  updateValue< nest::double_t >( d, "w_out", w_out_ );
  updateValue< nest::double_t >( d, "STDP_weight_effective_len", STDP_weight_effective_len_ );
  //updateValue< nest::double_t >( d, "eta_d", eta_d_ );
  updateValue< nest::double_t >( d, "delay_max", delay_max_ );
  updateValue< nest::double_t >( d, "delay_min", delay_min_ );
  updateValue< nest::double_t >( d, "sigma_minus", sigma_minus_ );
  updateValue< nest::double_t >( d, "sigma_plus", sigma_plus_ );
  updateValue< nest::double_t >( d, "B_minus", B_minus_ );
  updateValue< nest::double_t >( d, "B_plus", B_plus_ );
  updateValue< nest::double_t >( d, "delta_minus", delta_minus_ );
  updateValue< nest::double_t >( d, "delta_plus", delta_plus_ );
  //updateValue< nest::double_t >( d, "STDP_delay_effective_len", STDP_delay_effective_len_ );
  
  if ( tau_plus_ <= 0.0 ) throw nest::BadProperty( "tau_plus should be > 0." );
  if ( tau_minus_ <= 0.0 ) throw nest::BadProperty( "tau_minus should be > 0." );
  if ( A_plus_ < 0.0 ) throw nest::BadProperty( "A_plus should be >= 0." );
  if ( A_minus_ < 0.0 ) throw nest::BadProperty( "A_minus should be >= 0." );
  if ( w_out_ > 0.0 ) throw nest::BadProperty( "w_out should be <= 0." );
  if ( w_in_ < 0.0 ) throw nest::BadProperty( "w_in should be >= 0." );
  if ( STDP_weight_effective_len_ <= 0.0 ) throw nest::BadProperty( "STDP_weight_effective_len should be > 0." );
  if ( delay_min_ <= 0.0 ) throw nest::BadProperty( "delay_min should be > 0." );
  if ( delay_max_ < delay_min_ ) throw nest::BadProperty( "delay_max should be >= delay_min." );
  if ( sigma_plus_ <= 0.0 ) throw nest::BadProperty( "sigma_plus should be > 0." );
  if ( sigma_minus_ <= 0.0 ) throw nest::BadProperty( "sigma_minus should be > 0." );
  if ( B_plus_ < 0.0 ) throw nest::BadProperty( "B_plus should be >= 0." );
  if ( B_minus_ < 0.0 ) throw nest::BadProperty( "B_minus should be >= 0." );
  if ( delta_plus_ > 0.0 ) throw nest::BadProperty( "delta_plus should be <= 0." );
  if ( delta_minus_ < 0.0 ) throw nest::BadProperty( "delta_minus should be >= 0." );
}

} // namespace

#endif // joint_stdp_connection.h
